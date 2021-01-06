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
#' path <- system.file("extdata", "sample",
#'                    "Set4_1-6plex_[16142,55840]_binary_seg_maps.tif",
#'                    package = "phenoptr")
#' maps <- read_maps(path)
#' names(maps)
#' @md
read_maps <- function(map_path) {
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

#' Make a spatially-referenced raster from a map image
#' @param map Map image from `read_maps`
#' @return A spatially-referenced raster
#' @export
map_as_raster = function(map) {
  # Check for required spatial reference in map attributes
  attrs = attributes(map)
  required_attrs = c('x.position', 'y.position',
                     'x.resolution', 'resolution.unit')
  if (!all(required_attrs %in% names(attrs)))
    stop('Map file does not contain required attributes.')

  if (attrs$resolution.unit != 'cm')
    stop('map_as_raster requires centimeter units.')
  microns_per_pixel = 10000 / attrs$x.resolution
  rastr = raster::raster(map,
                    xmn = attrs$x.position,
                    xmx = attrs$x.position + dim(map)[2] * microns_per_pixel,
                    ymn = attrs$y.position,
                    ymx = attrs$y.position + dim(map)[1] * microns_per_pixel
  )

  # Flip to get the correct orientation
  raster::flip(rastr, direction='y')
}

#' Get the mapping from tissue category name to label number for a
#' tissue category map.
#' @param map A tissue category image from read_maps.
#' @return A named vector mapping tissue category names to label numbers
#' for the map.
#' @export
get_tissue_category_index = function(map) {
  desc = attr(map, 'description')
  stopifnot(stringr::str_detect(desc, 'TissueClassMap'))
  xml = xml2::read_xml(desc)
  categories = xml2::xml_find_all(xml, './/Name') %>% xml2::xml_text()

  # Categories are numbered starting from 0
  (seq_along(categories)-1) %>% rlang::set_names(categories)
}

#' Get the path to the segmentation map file for a field
#'
#' @param field_name Name of the field of interest
#' @param export_path Path to the inForm export directory
#' @return The path to the `binary_seg_maps` file for the given field.
#' @export
get_map_path = function(field_name, export_path) {
  field_base = stringr::str_remove(field_name, '\\.im3')
  map_path = file.path(export_path, paste0(field_base,
                                                 '_binary_seg_maps.tif'))
  if(!file.exists(map_path)) {
    warning('File not found: "', map_path, '"')
    return(NULL)
  }

  map_path
}

