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
#' limitations in the `tiff` package.
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
read_components <- function(path)
{
  stopifnot(file.exists(path), endsWith(path, 'component_data.tif'))

  # The readTIFF documentation for `as.is`` is misleading.
  # To read actual float values for components, we want `as.is=FALSE`.
  tif = tiff::readTIFF(path, all=TRUE, info=TRUE)

  # Get the image descriptions and figure out which ones are components
  infos = purrr::map_chr(tif, ~attr(., 'description'))
  images = grepl('FullResolution', infos)

  # Get the component names
  names = stringr::str_match(infos[images], '<Name>(.*)</Name>')[,2]

  purrr::set_names(tif[images], names)
}
