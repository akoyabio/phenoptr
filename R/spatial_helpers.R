#' Make a marked point pattern for cells in a single field
#'
#' Makes a marked point pattern ([`spatstat.geom::ppp`]) for the cells
#' of the specified phenotype in the
#' specified tissue categories and field.
#'
#' @param csd Cell seg data, may contain multiple fields
#' @param export_path Path to a directory containing composite and component
#'   image files from inForm
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
    wind = spatstat.geom::owin(xrange=xrange, yrange=yrange, unitname='micron')
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
    wind = spatstat.geom::owin(mask=mask, xrange=xrange, yrange=yrange,
                          unitname='micron')
  }

  # Finally create the actual point pattern
  pp = spatstat.geom::ppp(field_data$`Cell X Position`, field_data$`Cell Y Position`,
                window=wind,
                marks=factor(rep(pheno_name, nrow(field_data))))
  pp
}

