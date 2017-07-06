#' Read an inForm segmentation map file.
#'
#' Reads an inForm `binary_seg_maps.tif`` file and returns a named
#' list of images. The names reflect the content of the individual images.
#' Possible names are Nucleus, Cytoplasm, Membrane, Object,
#' TissueClassMap, and ProcessRegionImage; not every image file will include all
#' names.
#'
#' Images are oriented to match the coordinates in a cell seg data file,
#' i.e. (0, 0) at the top left and the row corresponding to Y and
#' column corresponding to X.
#'
#' @param map_path Path to the map file or a cell seg data file in the same
#' directory.
#' @return A named list of images, one for each map in the source file.
#' @export
#' @family file readers
#' @examples
#' \dontrun{
#' # Not included at this time
#' path = system.file("extdata",
#'   "TMA/Core[1,5,6,1]_[21302,15107]_binary_seg_maps.tif",
#'   package = "informr")
#' d = read_maps(path)
#' names(d)
#' }
#' @md
read_maps = function(map_path) {
  # Allow a cell seg path to be passed in
  map_path = sub('cell_seg_data.txt', 'binary_seg_maps.tif', map_path)

  stopifnot(file.exists(map_path), endsWith(map_path, 'binary_seg_maps.tif'))

  # Read the mask file and get the image descriptions
  masks = tiff::readTIFF(map_path, all=TRUE, info=TRUE, as.is=TRUE)
  infos = purrr::map_chr(masks, ~attr(., 'description'))

  # All possible maps
  map_keys = c(Nucleus='Nucleus', Cytoplasm='Cytoplasm', Membrane='Membrane',
               Object='Object', Tissue='TissueClassMap',
               ROI='ProcessRegionImage')
  maps = list()
  for (n in names(map_keys)) {
    m = masks[stringr::str_detect(infos, map_keys[n])]
    if (length(m)>0)
      maps[[n]] = m[[1]]
  }

  maps
}
