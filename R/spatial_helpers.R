#' Make a marked point pattern for cells in a single field
#'
#' Makes a marked point pattern ([`spatstat.geom::ppp`]) for the cells
#' of the specified phenotype in the
#' specified tissue categories and field.
#'
#' @param csd Cell seg data, may contain multiple fields
#' @param export_path Path to a directory containing composite and component
#'   or segmentation map image files from inForm
#' @param pheno Phenotype definition. Either a (possibly named) character
#' vector to be parsed by [parse_phenotypes()] or a named list containing
#' a single phenotype definition as used by [select_rows()].
#' @param field_name Sample Name or Annotation ID for the field of interest.
#' May be omitted if `csd` contains data for only one field.
#' @param tissue_categories Tissue categories of interest. If supplied, the
#' returned `ppp` will contain only cells in these categories and the
#' associated window will be restricted to the extent of these categories.
#' @return Returns a marked point pattern ([`spatstat.geom::ppp`] object)
#' with a single mark value.
#' @examples
#' # ppp for CD8+ cells in the sample data
#' suppressPackageStartupMessages(library(spatstat))
#' pp <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
#'   "CD8+", tissue_categories="Tumor")
#' plot(pp, show.window=FALSE, main='')
#'
#' # To include multiple phenotypes in a single point pattern,
#' # create them separately and join them with [spatstat.geom::superimpose()].
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
            length(pheno) == 1
  )

  # Get data for the correct field
  field_col = field_column(csd)
  if (is.null(field_name)) {
    # If field_name is not provided, csd must be for a single field
    if (dplyr::n_distinct(csd[[field_col]]) != 1)
      stop('Please pass a field_name to make_ppp to select a single field.')

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
  if (is.character(pheno)) {
    if (rlang::is_named(pheno))
      pheno_name = names(pheno)
    else
      pheno_name = pheno

    pheno = parse_phenotypes(as.list(pheno))[[1]]
  } else {
    # Pheno must be a named list containing a single definition
    stopifnot(is.list(pheno), !is.null(names(pheno)))
    pheno_name = names(pheno)
    pheno = pheno[[1]] # Extract the actual phenotype definition
  }

  field_data = field_data[select_rows(field_data, pheno), ]

  if (nrow(field_data) == 0) {
    if (is.null(tissue_categories))
      warning(stringr::str_glue(
        'Field "{field_name}" contains no "{pheno_name}" cells.'))
    else {
      categories = paste0("\"", tissue_categories, "\"", collapse=", ")
      warning(stringr::str_glue(
        'Field "{field_name}" contains no {pheno_name} cells in ',
        'tissue category {categories}'))
    }
  }

  field_info = read_field_info(field_name, export_path)
  stopifnot(!is.null(field_info))

  field_data = correct_for_pixel_data(field_data, field_info)

  # Get the spatstat window
  xrange=c(field_info$location[1], field_info$location[1]+field_info$field_size[1])
  yrange=c(field_info$location[2], field_info$location[2]+field_info$field_size[2])
  if (is.null(tissue_categories)) {
    # Window is the full field
    wind = spatstat.geom::owin(xrange=xrange, yrange=yrange, unitname='micron')
  } else {
    # Window based on tissue categories
    map_path = get_map_path(field_name, export_path)
    maps = read_maps(map_path)
    stopifnot('Tissue' %in% names(maps))

    tissue = maps[['Tissue']]
    layers = parse_tissue_description(tissue)
    layer_nums = layers[tissue_categories]
    if (any(is.na(layer_nums))) {
      missing = tissue_categories[is.na(layer_nums)]
      stop('Tissue category missing from map file: ', basename(map_path), ' - ',
           paste0('"', missing, '"', collapse=', '))
    }

    mask = tissue %in% layer_nums
    dim(mask) = dim(tissue) # Convert back to matrix
    wind = spatstat.geom::owin(mask=mask, xrange=xrange, yrange=yrange,
                          unitname='micron')
  }

  # Finally create the actual point pattern
  pp = spatstat.geom::ppp(field_data$`Cell X Position`, field_data$`Cell Y Position`,
                window=wind,
                marks=factor(rep(pheno_name, nrow(field_data))))
  pp
}


# Operations with ROIs
# These require {sf}

#' Read and combine all tagged ROIS from an annotation file
#'
#' Reads all tagged ROIs from the provided file. ROIs with the same
#' tag are combined with `sf::st_union()`. A single (multi)polygon
#' object is returned for each unique tag.
#'
#' @param annotation_file Path to annotations file
#' @return A named list of (multi)polygons, one for each tag
#' in the annotations file
#' @importFrom rlang .data
#' @export
read_tagged_rois = function(annotation_file) {
  rois = read_phenochart_polygons(annotation_file)
  if (nrow(rois) > 0)
    rois = rois %>% dplyr::filter(.data$tags != '') # Only tagged ROIs

  if (nrow(rois) > 0) {
    # Get a list of all unique tags
    all_roi_names = rois$tags %>%
      stringr::str_split(' ') %>%
      unlist() %>%
      unique() %>%
      sort() %>%
      rlang::set_names()
  } else return(list())

  # Make a single (multi) polygon for each tag
  tagged_rois = purrr::map(all_roi_names, function(tag_name) {
    rois %>%
      dplyr::filter(stringr::str_detect(.data$tags, tag_name)) %>%
      sf::st_union()
  })

  tagged_rois
}


#' Read polygons (ROIs) and tags from a Phenochart annotation.
#'
#' Reads the polygons and (optionally) included rectangles for all
#' ROI annotations in a single file.
#'
#' @param xml_path Path to an annotations file.
#' @param include_rects If `TRUE`, the result will include the rectangles
#' contained in each annotation.
#' @return An `sf::st_sf` object with columns `tags`, `geometry` and
#' (optionally) `rects`.
#' `rects` is a list of `sf::st_sf` objects with columns `tags`, `center_x`,
#' `center_y` and `geometry`.
#' Multiple tags for a single ROI are separated by spaces in a single string.
#' @export
#' @importFrom magrittr %>%
read_phenochart_polygons = function(xml_path, include_rects=TRUE) {
  xml = xml2::read_xml(xml_path)
  ns = xml2::xml_ns(xml)

  rois =
    xml2::xml_find_all(xml, './/Annotations-i[@subtype="ROIAnnotation"]', ns=ns)

  result = purrr::map_dfr(rois, function(roi) {
    tags = xml2::xml_find_all(roi, './Tags/Tags-i', ns=ns) %>%
      xml2::xml_text() %>%
      paste(collapse=' ')
    locs = xml2::xml_find_all(roi, './/Perimeter-i', ns=ns)
    x = xml2::xml_find_all(locs, './X', ns=ns) %>% xml2::xml_double()
    y = xml2::xml_find_all(locs, './Y', ns=ns) %>% xml2::xml_double()

    # A valid polygon has at least four points, making a triangle
    if (length(x)<4) return(NULL)

    poly = sf::st_polygon(list(matrix(c(x, y), ncol=2)))

    # At this point poly may not contain a valid geometry. There are two
    # ways to fix this: st_make_valid() and st_buffer(., 0).
    # st_make_valid makes GEOMETRYCOLLECTION values which can't be
    # converted to Spatial objects. st_buffer makes MULTIPOLYGONs which
    # are more tractable.
    poly = sf::st_buffer(poly, 0)
    parsed = sf::st_sf(tags=tags,
              geometry=sf::st_sfc(poly),
              stringsAsFactors=FALSE)

    if (include_rects) {
      rects = parse_roi_rects(roi, ns)
      parsed$rects = list(rects)
    }
    parsed
  })

  if (nrow(result) == 0)
    result = NULL
  result
}

#' Read the tagged fields from an annotation file
#' @param xml_path Path to an annotations file.
#' @return An `sf::st_sf` object with columns `tags`, `center_x`,
#' `center_y` and `geometry`.
#' Multiple tags for a single field are separated by spaces in a single string.
#' @export
read_phenochart_fields = function(xml_path) {
  xml = xml2::read_xml(xml_path)
  ns = xml2::xml_ns(xml)

  # Find top-level rectangles
  fields =
    xml2::xml_find_all(xml,
                       '//Annotations-i[@subtype="RectangleAnnotation"]', ns=ns)

  # Find nested rectangles
  fields2 =
    xml2::xml_find_all(xml, '//Fields-i[@subtype="RectangleAnnotation"]', ns=ns)
  fields = c(fields, fields2)

  # Read tags first, it is quicker than parsing all the rects
  tags = purrr::map_chr(fields, parse_tags, ns=ns)

  # Just keep the tagged fields
  tagged = tags != ''
  tags = tags[tagged]
  fields = fields[tagged]

  # Now get the rects
  rects = purrr::map_dfr(fields, parse_rect, ns=ns)
  rects$tags = tags
  rects %>%
    dplyr::select(tags, dplyr::everything())
}

#' Parse the rectangles from an ROI annotation
#' @param roi Parsed XML for an ROI annotation
#' @param ns Optional namespace; providing this makes parsing faster.
#' @return An `sf::sf` object with geometries for each rectangle in roi, as
#' returned from `parse_rect()`.
#' @keywords internal
parse_roi_rects = function(roi, ns = xml2::xml_ns(roi)) {
  fields =
    xml2::xml_find_all(roi, './/Fields-i[@subtype="RectangleAnnotation"]', ns=ns)
  purrr::map_dfr(fields, parse_rect, ns=ns)
}

#' Parse tags from a node of an annotations file
#' @param node XML node containing tags
#' @param ns Optional namespace; providing this makes parsing faster.
#' @return A string with space-separated tags
#' @keywords internal
parse_tags = function(node, ns=xml2::xml_ns(node)) {
  xml2::xml_find_all(node, './/Tags-i', ns=ns) %>%
    xml2::xml_text() %>%
    paste(collapse=' ')
}

#' Parse a single rectangle annotation
#' @param field Parsed XML for a single rectangle annotation
#' @param ns Optional namespace; providing this makes parsing many rects
#' much faster.
#' @return An `sf::st_sf` with coordinates of the rectangle
#' and columns `center_x` and `center_y` giving the *rounded*
#' center point that inForm uses to identify the field.
#' @keywords internal
parse_rect = function(field, ns=xml2::xml_ns(field)) {
  x = xml2::xml_find_first(field, './/Origin/X', ns=ns) %>% xml2::xml_double()
  y = xml2::xml_find_first(field, './/Origin/Y', ns=ns) %>% xml2::xml_double()
  width = xml2::xml_find_first(field, './/Size/Width', ns=ns) %>% xml2::xml_double()
  height = xml2::xml_find_first(field, './/Size/Height', ns=ns) %>% xml2::xml_double()
  pts = matrix(c(x, y,
                 x, y+height,
                 x+width, y+height,
                 x+width, y,
                 x, y), ncol=2, byrow=TRUE)
  center_x = round(x+width/2, 0)
  center_y = round(y+height/2, 0)
  sf::st_sf(center_x=center_x, center_y=center_y,
        geometry=sf::st_sfc(sf::st_polygon(list(pts))))
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
  # Sanity check the export path by finding matching maps for each annotation
  expected_maps = paste0(annotations, '_binary_seg_maps.tif')
  existing_maps = list.files(export_path, 'binary_seg_maps.tif',
                             full.names=TRUE, recursive=TRUE)

  # Indices of expected_maps in existing_maps
  matching_maps = match(expected_maps, basename(existing_maps))

  if (any(is.na(matching_maps))) {
    missing_maps = annotations[is.na(matching_maps)]
    stop('The export directory is missing segmentation maps for\n',
         paste(missing_maps, collapse='\n'))
  }

  # Now the actual matching paths
  matching_maps = existing_maps[matching_maps]

  # raster::mask needs a Spatial* object
  sp_roi = sf::as_Spatial(roi)

  # This will capture the tissue category names and numbers
  tissue_index = NULL

  # Helper to process a single annotation
  pb = dplyr::progress_estimated(length(annotations))
  trim_tissue_categories_single = function(annotation, map_path) {
    pb$tick()$print() # Show progress

    # Get the tissue map for annotation and trim by roi
    maps = read_maps(map_path)
    if (!'Tissue' %in% names(maps))
      stop('No tissue category map in ', map_path)

    tissue = maps$Tissue
    tissue_ix = parse_tissue_description(tissue)

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
  trimmed = purrr::map2(annotations, matching_maps, trim_tissue_categories_single)
  result  = purrr::map_dfr(trimmed, 'data')
  cat('\n')

  if (!is.null(plot_path)) {
    cat('Saving reference plot\n')

    # Make a merged raster from the trimmed annotations
    trimmed_rasters = purrr::map(trimmed, 'raster')

    rm(trimmed) # Free some memory before merging

    # To get the plot to come out right we need to invert each raster
    # in Y and negate the bounding box
    invert_raster = function(r) {
      r = raster::flip(r, 'y')
      ext = raster::extent(r)
      raster::setExtent(r, raster::extent(ext@xmin, ext@xmax, -ext@ymax, -ext@ymin))
    }

    trimmed_rasters = purrr::map(trimmed_rasters, invert_raster)

    # One big raster
    if (length(trimmed_rasters) == 1)
      merged_rasters = trimmed_rasters[[1]]
    else
      merged_rasters = rlang::exec(raster::merge, !!!trimmed_rasters, tolerance=1)

    # Save a plot of the merged rasters
    # One color for each tissue category, skipping black
    # If necessary, fill in with grey to get to 256 values
    # Tissue maps may (rarely) include values of 255 when areas can't be classified
    # 255 is not included in tissue_index
    raster_max = merged_rasters@data@max
    colors = grDevices::palette()[2:(max(tissue_index)+2)]
    if (raster_max > max(tissue_index))
      colors = c(colors, rep('grey', raster_max+1-length(colors)))

    tissue_index = sort(tissue_index) # Make sure colors and legend match

    # Prepare to save the plot
    grDevices::png(plot_path, type='cairo', antialias='gray',
        width=980, height=980)

    # Plot the actual raster
    raster::plot(merged_rasters, axes=TRUE, legend=FALSE,
                 col=colors,
                 main=tools::file_path_sans_ext(basename(plot_path)))

    # Legend
    # This gets funky when raster_max != max(tissue_index)
    # raster::plot is expecting a legend for raster_max + 1 colors
    # We are only giving max(tissue_index) + 1 colors
    # The label locations are corrected for the mismatch
    raster::plot(merged_rasters, legend.only=TRUE,
                 col=colors[1:length(tissue_index)],
                 axis.args=list(
                   at=tissue_index * raster_max/max(tissue_index),
                   labels=names(tissue_index)))

    # The ROI that we clipped to; invert it to match the raster
    plot(negate_y(roi), add=TRUE, asp=1)

    grDevices::dev.off()
  }

  # Return the trimmed data
  result
}

#' Parse the description tag of the tissue map image read from an
#' inForm `binary_seg_maps` file
#' @param img The `Tissue` image from \code{\link{read_maps}}.
#' @return A named vector whose names are the tissue classes in `img` and
#' whose values are the mask values for the class in `img`.
#' @md
#' @export
parse_tissue_description = function(img) {
  desc = attr(img, 'description')
  stopifnot(!is.null(desc))
  parsed = xml2::read_xml(desc)
  layers = purrr::map_chr(xml2::xml_find_all(parsed, 'Entry'),
    function(entry) {
      entry %>%
        xml2::xml_find_first('Name') %>%
        xml2::xml_text() %>%
        stringr::str_trim()
    })

  purrr::set_names(seq_along(layers)-1, layers)
}
