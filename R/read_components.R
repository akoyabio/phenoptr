#' Read a component TIFF file.
#'
#' Read an inForm `component_data.tif` file.
#' If the file contains multiple resolutions or other data,
#' keep only the full-resolution component images,
#' extract the component names from the image descriptions,
#' and return a named list of image matrices.
#'
#' Images are oriented to match the coordinates in a cell seg data file,
#' i.e. (0, 0) at the top left and the row corresponding to Y and
#' column corresponding to X.
#'
#' Limited to 1x1 fields due to
#' limitations in the `tiff` package. Install the
#' [Akoya Biosciences fork](https://github.com/akoyabio/tiff)
#' of the tiff package to remove this limitation.
#' @param path Path to the `component_data.tif` file.
#' @return A named list of matrices, one for each component.
#' @export
#' @family file readers
#' @examples
#' \dontrun{
#' path <- "/path/to/my_component_data.tif"
#' components <- read_components(path)
#' }
#' @md
read_components <- function(path) {
  stopifnot(file.exists(path), endsWith(path, 'component_data.tif'))

  # The readTIFF documentation for `as.is`` is misleading.
  # To read actual float values for components, we want `as.is=FALSE`.
  tif = tiff::readTIFF(path, all=TRUE, info=TRUE)

  # Get the image descriptions and figure out which ones are components
  infos = purrr::map_chr(tif, ~attr(., 'description'))
  images = grepl('FullResolution', infos)

  # Get the component names
  names = stringr::str_match(infos[images], '<Name>(.*)</Name>')[, 2]

  purrr::set_names(tif[images], names)
}

#' Read information about a field from a component image file.
#'
#' Find the location, size and magnification of an inForm field by inspecting a
#' `component_data.tif` file.
#'
#' The field location is determined from the coordinates in the file
#' name. The field size and magnification are read from TIFF tags
#' in the file.
#'
#' Limited to 1x1 fields due to
#' limitations in the `tiff` package. Install the
#' [Akoya Biosciences fork](https://github.com/akoyabio/tiff)
#' of the tiff package to remove this limitation.
#' @param path Path to the `component_data.tif` file.
#' @return A named list with values
#'   \describe{
#'    \item{\code{location}}{The physical (x, y) location of the top-left
#'     corner of the field on the source slide, in microns.}
#'    \item{\code{field_size}}{The physical (width, height)
#'    of the field, in microns.}
#'    \item{\code{image_size}}{The (width, height) of the field image,
#'    in pixels.}
#'    \item{\code{microns_per_pixel}}{The size of each pixel, in microns.}
#'  }
#'
#' @export
#' @family file readers
#'
#' @examples
#'   path <- system.file('extdata', 'sample',
#'     'Set4_1-6plex_[16142,55840]_component_data.tif',
#'     package='phenoptr')
#'   get_field_info(path)
#' @md
get_field_info = function(path) {
  # Use readTIFFDirectory if available, it is faster, more complete
  # and works with tiled images
  if ('readTIFFDirectory' %in% ls(getNamespace('tiff'))) {
    info = tiff::readTIFFDirectory(path, all=FALSE)
    center = NA # Don't need this
   } else {
    tif = try(tiff::readTIFF(path, all=FALSE, info=TRUE))
    if (inherits(tif, 'try-error')) {
      # Give a helpful error message if the problem is missing support
      # for tiled images.
      if (stringr::str_detect(tif, 'tile-based images'))
        stop('Please install akoyabio/tiff ',
             'to support large component data files.\n',
             'Use "devtools::install_github(\'akoyabio/tiff\')"')
      else
        stop()
    }
    info = attributes(tif)
    info$length = info$dim[1]
    info$width = info$dim[2]

    # We have to get the location from the file name; tiff::readTIFF()
    # doesn't read the location tags :-(
    name = basename(path)
    center_pattern = "_\\[([\\d\\.]+),([\\d\\.]+)\\][^\\[]*$"
    center = as.numeric(stringr::str_match(name, center_pattern)[, 2:3])
    if (any(is.na(center)))
      stop("Field location not found in file name.")
  }

  required_attributes = c('width', 'length', 'x.resolution', 'resolution.unit')
  missing_attributes = setdiff(required_attributes, names(info))
  if (length(missing_attributes) > 0) {
    missing = paste(missing_attributes, collapse=', ')
    stop(paste0('Image file is missing required attributes: ', missing))
  }

  if (info$resolution.unit != 'cm')
    stop(paste('Unsupported resolution unit:', info$resolution.unit))

  result = list()
  result$image_size = c(info$width, info$length)
  result$microns_per_pixel = as.numeric(10000/info$x.resolution)
  result$field_size = result$image_size * result$microns_per_pixel
  if (is.na(center)) {
    # Location directly from TIFF info
    result$location = c(info$x.position, info$y.position) * 10000
  } else {
    # Location from file name and resolution
    result$location = center - result$field_size/2
  }
  result
}
