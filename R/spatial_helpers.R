#' Make a marked point pattern for cells in a single field
#'
#' Makes a marked point pattern ([`spatstat::ppp`]) for the cells
#' of the specified phenotype in the
#' specified tissue categories and field.
#'
#' @param csd Cell seg data, may contain multiple fields
#' @param export_path Path to a directory containing composite and component
#'   image files from inForm
#' @param pheno Phenotype definition. Either a (possibly named) character
#' vector to be parsed by [parse_phenotypes()] or a named list as
#' used by [select_rows()].
#' @param field_name Sample Name or Annotation ID for the field of interest.
#' May be omitted if `csd` contains data for only one field.
#' @param tissue_categories Tissue categories of interest. If supplied, the
#' returned `ppp` will contain only cells in these categories and the
#' associated window will be restricted to the extent of these categories.
#' @return Returns a marked point pattern ([`spatstat::ppp`] object)
#' with a single mark value.
#' @examples
#' # ppp for CD8+ cells in the sample data
#' suppressPackageStartupMessages(library(spatstat))
#' pp <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
#'   "CD8+", tissue_categories="Tumor")
#' plot(pp, show.window=FALSE, main='')
#'
#' # To include multiple phenotypes in a single point pattern,
#' # create them separately and join them with [spatstat::superimpose()].
#'
#' pp2 <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
#'   "CK+", tissue_categories="Tumor")
#' plot(superimpose(pp2, pp), cols=c('red', 'blue'), show.window=FALSE, main='')
#' @export
#' @md
make_ppp = function(csd, export_path, pheno,
                    field_name=NULL, tissue_categories=NULL) {
  stopifnot(!is.null(csd),
            unit_is_microns(csd),
            dir.exists(export_path),
            !is.null(pheno),
            length(pheno) == 1)

  # Get data for the correct field
  field_col = field_column(csd)
  if (is.null(field_name)) {
    # If field_name is not provided, csd must be for a single field
    stopifnot(dplyr::n_distinct(csd[[field_col]]) == 1)
    field_name = csd[[field_col]][[1]]
    field_data = csd
  } else {
    field_data = csd[csd[[field_col]]==field_name, ]
  }

  # Subset to tissue category if provided
  if (!is.null(tissue_categories)) {
    stopifnot('Tissue Category' %in% names(field_data))
    field_data = field_data %>%
      dplyr::filter(`Tissue Category` %in% tissue_categories)
  }

  # Process the phenotype
  if (is.character(pheno))
    pheno = parse_phenotypes(as.list(pheno))
  if (!rlang::is_named(pheno))
    pheno = rlang::set_names(pheno)
  field_data = field_data[select_rows(field_data, pheno), ]

  if (nrow(field_data) == 0)
    warning(stringr::str_glue(
      'Field "{field_name}" contains no {names(pheno)[[1]]} cells.'))

  field_info = read_field_info(field_name, export_path)
  stopifnot(!is.null(field_info))

  field_data = correct_for_pixel_data(field_data, field_info)

  # Get the spatstat window
  xrange=c(field_info$location[1], field_info$location[1]+field_info$field_size[1])
  yrange=c(field_info$location[2], field_info$location[2]+field_info$field_size[2])
  if (is.null(tissue_categories)) {
    # Window is the full field
    wind = spatstat::owin(xrange=xrange, yrange=yrange, unitname='micron')
  } else {
    # Window based on tissue categories
    map_path = get_map_path(field_name, export_path)
    maps = read_maps(map_path)
    stopifnot('Tissue' %in% names(maps))

    tissue = maps[['Tissue']]
    layers = parse_tissue_description(tissue)
    layer_nums = layers[tissue_categories]

    mask = tissue %in% layer_nums
    dim(mask) = dim(tissue) # Convert back to matrix
    wind = spatstat::owin(mask=mask, xrange=xrange, yrange=yrange,
                          unitname='micron')
  }

  # Finally create the actual point pattern
  pp = spatstat::ppp(field_data$`Cell X Position`, field_data$`Cell Y Position`,
                window=wind,
                marks=factor(rep(names(pheno)[[1]], nrow(field_data))))
  pp
}


# Operations with ROIs
# These require {sf}

#' Read polygons and tags from a Phenochart annotation.
#'
#' Reads the polygons and included rectangles for all ROI annotations
#' in a single file.
#' @param xml_path Path to an annotations file.
#' @return An `sf::st_sf` object with columns `tags`, `geometry` and `rects`.
#' Multiple tags for a single ROI are separated by spaces in a single string.
#' @export
#' @importFrom magrittr %>%
read_phenochart_polygons = function(xml_path) {
  xml = xml2::read_xml(xml_path)
  rois = xml2::xml_find_all(xml, './/Annotations-i[@subtype="ROIAnnotation"]')
  result = purrr::map_dfr(rois, function(roi) {
    tags = xml2::xml_find_all(roi, './/Tags-i') %>%
      xml2::xml_text() %>%
      paste(collapse=' ')
    locs = xml2::xml_find_all(roi, './/Perimeter-i')
    x = xml2::xml_find_all(locs, './/X') %>% xml2::xml_double()
    y = xml2::xml_find_all(locs, './/Y') %>% xml2::xml_double()

    # A valid polygon has at least four points, making a triangle
    if (length(x)<4) return(NULL)

    poly = sf::st_polygon(list(matrix(c(x, y), ncol=2)))

    # At this point poly may not contain a valid geometry. There are two
    # ways to fix this: st_make_valid() and st_buffer(., 0).
    # st_make_valid makes GEOMETRYCOLLECTION values which can't be
    # converted to Spatial objects. st_buffer makes MULTIPOLYGONs which
    # are more tractable.
    poly = sf::st_buffer(poly, 0)
    rects = parse_rects(roi)
    parsed = sf::st_sf(tags=tags,
              geometry=sf::st_sfc(poly),
              stringsAsFactors=FALSE)
    parsed$rects = list(rects)
    parsed
  })

  if (nrow(result) == 0)
    result = NULL
  result
}

# Parse the rectangles from an ROI annotation
# @param roi Parsed XML for an ROI annotation
# @return An `sf::sfc` object with geometries for each rectangle in roi
parse_rects = function(roi) {
  fields = xml2::xml_find_all(roi, './/Fields-i[@subtype="RectangleAnnotation"]')
  rects = purrr::map(fields, parse_field)
  sf::st_sfc(rects)
}

parse_field = function(field) {
  x = xml2::xml_find_first(field, './/Origin/X') %>% xml2::xml_double()
  y = xml2::xml_find_first(field, './/Origin/Y') %>% xml2::xml_double()
  width = xml2::xml_find_first(field, './/Size/Width') %>% xml2::xml_double()
  height = xml2::xml_find_first(field, './/Size/Height') %>% xml2::xml_double()
  pts = matrix(c(x, y,
                 x, y+height,
                 x+width, y+height,
                 x+width, y,
                 x, y), ncol=2, byrow=TRUE)
  sf::st_polygon(list(pts))
}

#' Add a geometry column to a cell seg table
#' @param csd A cell seg data table
#' @return The table with a `geometry` column of `sf::st_point` objects
#' @export
add_geometry = function(csd) {
  sf::st_as_sf(csd, coords=c('Cell X Position', 'Cell Y Position'),
               remove=FALSE, dim='XY')
}

#' Compute corrected tissue category areas per annotation by clipping
#' to a provided polygon and counting pixels.
#' @param annotations Names of annotations to process
#' @param roi Clipping polygon, as an `sf::st_sf` or `sf::st_sfc` object.
#' @param export_path Path to a directory containing binary_seg_map files.
#' These files must include TIFF tags for location and resolution.
#' @param plot_path If not `NULL`, the file path at which to save
#'  a check plot showing the given ROI and the trimmed tissue categories.
#' @return A data frame containing computed tissue category areas for
#' all annotations
#' @export
trim_tissue_categories = function(annotations, roi,
                                  export_path, plot_path=NULL) {
  # Sanity check the export path
  expected_maps = paste0(annotations, '_binary_seg_maps.tif')
  existing_maps = list.files(export_path, 'binary_seg_maps.tif')
  if (!all(expected_maps %in% existing_maps))
    stop('The export directory is missing segmentation maps.')

  # raster::mask needs a Spatial* object
  sp_roi = sf::as_Spatial(roi)

  # This will capture the tissue category names and numbers
  tissue_index = NULL

  # Helper to process a single annotation
  pb = dplyr::progress_estimated(length(annotations))
  trim_tissue_categories_single = function(annotation) {
    pb$tick()$print() # Show progress

    # Get the tissue map for annotation and trim by roi
    map_path = get_map_path(annotation, export_path)
    maps = read_maps(map_path)
    if (!'Tissue' %in% names(maps))
      stop('No tissue category map in ', map_path)

    tissue = maps$Tissue
    tissue_ix = get_tissue_category_index(tissue)

    # Merge tissue_ix into tissue_index while preserving names
    # base set operations union and setdiff don't work here
    tissue_index <<- c(tissue_index, tissue_ix[!tissue_ix %in% tissue_index])

    rastr = map_as_raster(tissue)
    rastr = raster::mask(rastr, sp_roi)

    # For each tissue category, get the area within the category in microns^2
    values = raster::getValues(rastr)
    mpp2 = prod(raster::res(rastr)) # microns^2 per pixel
    areas = purrr::map_dbl(tissue_index,
                           ~sum(values==.x, na.rm=TRUE) * mpp2)
    areas['All'] = sum(areas)

    data = tibble::enframe(areas, name='Tissue Category',
                           value='Tissue Category Area (square microns)') %>%
      tibble::add_column(`Annotation ID`=annotation, .before=1)

    # Return the data and optionally the trimmed raster
    ret = list(data=data)
    if (!is.null(plot_path)) ret$raster=rastr
    ret
  }

  # Process data for each field
  cat('Trimming tissue categories\n')
  trimmed = purrr::map(annotations, trim_tissue_categories_single)
  result  = purrr::map_dfr(trimmed, 'data')
  cat('\n')

  if (!is.null(plot_path)) {
    cat('Saving reference plot\n')

    # Make a merged raster from the trimmed annotations
    trimmed_rasters = purrr::map(trimmed, 'raster') %>%
      c(list(tolerance=1)) # Extra parameter for raster::merge

    rm(trimmed) # Free some memory before merging

    # One big raster
    merged_rasters = do.call(raster::merge, trimmed_rasters)

    # Save a plot of the merged rasters
    # One color for each tissue category, skipping black
    # Fill in with grey to get to 256 values
    # Tissue maps may (rarely) include values of 255 when areas can't be classified
    # 255 is not included in tissue_index
    colors = grDevices::palette()[2:(max(tissue_index)+2)]
    colors = c(colors, rep('grey', 256-length(colors)))

   tissue_index = sort(tissue_index) # Make sure colors and legend match

   grDevices::png(plot_path, type='cairo', antialias='gray',
        width=980, height=980)
    raster::plot(merged_rasters, axes=TRUE, legend=FALSE,
         col=colors,
         main=tools::file_path_sans_ext(basename(plot_path)))
    raster::plot(merged_rasters, legend.only=TRUE, col=colors,
         axis.args=list(at=tissue_index,
                        labels=names(tissue_index)))
    plot(roi, add=TRUE, asp=1)
    grDevices::dev.off()
  }

  # Return the trimmed data
  result
}
