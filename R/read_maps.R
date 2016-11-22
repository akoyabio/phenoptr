#' Read an inForm segmentation map file
#'
#' Reads an inForm \code{binary_seg_maps} file and returns a list of images.
#' @param map_path Path to the map file.
#' @return A named list of images, one for each map in the source file.
#' @export
#' @examples
#' path = system.file("extdata",
#'   "TMA/Core[1,5,6,1]_[21302,15107]_binary_seg_maps.tif",
#'   package = "informr")
#' d = read_maps(path)
#' names(d)
read_maps = function(map_path) {
  stopifnot(file.exists(map_path))

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
      maps[[n]] = t(m[[1]])
  }

  maps
}
