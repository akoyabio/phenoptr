#' Read composite images
#'
#' Read images based on inForm composite, including composite with
#' phenotype overlay, and phenotype overlay alone.
#'
#' @details
#' If the selected file is a JPEG file, the returned value will be a list
#' containing a single image matrix and no attributes other than `dim`.
#'
#' If the selected file is a TIFF file, the returned value will be a list
#' possibly containing multiple image matrices with many metadata attributes.
#' In addition to the TIFF metadata provided by `tiff::readTIFF()`,
#' attributes may include
#' - composite_name - the user's name for this composite
#' - components - a named list, where names are component names and values
#'   are hex colors
#' - scheme_name - the user's name for the phenotype scheme
#' - phenotypes -  a named list, where names are phenotype names and values
#'   are hex colors
#' @param path Path to the composite image file
#' @return A list of image matrices, possibly containing image metadata
#' in attributes (see Details).
#' @export
read_composites = function(path) {
  stopifnot(file.exists(path))

  # JPEG is easy
  if (stringr::str_detect(path, '\\.jpe?g$')) {
    img = jpeg::readJPEG(path)
    return(list(img))
  }

  # TIFF is more complex
  imgs = tiff::readTIFF(path, all=TRUE, info=TRUE)

  imgs = purrr::map(imgs, function(img) {
    img_info = parse_composite_info(attr(img, 'description'))
    attributes(img) = c(attributes(img), img_info)
    img
  })

  imgs
}

#' Read just the composite info from images in a TIFF composite image file.
#' @param path Path to a TIFF composite
#' @return A list of lists containing some of
#' - composite_name - the user's name for this composite
#' - components - a named vector, where names are component names and values
#'   are hex colors
#' - scheme_name - the user's name for the phenotype scheme
#' - phenotypes -  a named vector, where names are phenotype names and values
#'   are hex colors
#' @export
read_composite_info = function(path) {
  stopifnot(file.exists(path))
  if (!stringr::str_detect(path, '\\.tiff?$'))
    stop('read_composite_info requires a TIFF image.')

  if (!function_exists('tiff', 'readTIFFDirectory'))
    stop('read_composite_info requires the Akoya tiff package.')

  tiff::readTIFFDirectory(path, all=TRUE) %>%
    purrr::map('description') %>%
    purrr::map(parse_composite_info)
}

#' Parse the ImageDescription of a composite image.
#' @param info The ImageDescription string
#' @return A list, containing none, some or all of these items:
#' - composite_name - the user's name for this composite
#' - components - a named vector, where names are component names and values
#'   are hex colors
#' - scheme_name - the user's name for the phenotype scheme
#' - phenotypes -  a named vector, where names are phenotype names and values
#'   are hex colors
#' @keywords internal
parse_composite_info = function(info) {
  # If info is not XML, return an empty list
  result = list()
  if (!startsWith(info, '<?xml')) return(result)

  xml = xml2::read_xml(info)

  # Look for composite info first
  composite = xml2::xml_find_first(xml, './/Composite')
  if (length(composite) == 2) {
    result$composite_name = xml2::xml_find_first(composite, 'Name') %>%
      xml2::xml_text()
    result$components = parse_names_and_colors(composite, 'Marker')
  }

  # Now phenotypes
  phenotypes = xml2::xml_find_first(xml, './/PhenotypeSchema')
  if (length(phenotypes) == 2) {
    result$scheme_name = xml2::xml_find_first(phenotypes, 'Name') %>%
      xml2::xml_text()
    result$phenotypes = parse_names_and_colors(phenotypes, 'Name')
  }

  result
}

#' Parse the Entry items in a Composite or Phenotypes node
#' @param node The XML node to parse
#' @param name_node_name Name of the node containing item names
#' @return A named vector of colors
#' @keywords internal
parse_names_and_colors = function(node, name_node_name) {
  items = xml2::xml_find_all(node, './/Entry')
  item_names =  xml2::xml_find_all(items, name_node_name) %>%
    xml2::xml_text()
  item_colors = xml2::xml_find_all(items, 'Color') %>%
    xml2::xml_text() %>%             # Get the color text, e.g. "255,0,0,255"
    stringr::str_split(',') %>%      # Split
    purrr::map_chr(. %>% `[`(-1) %>% # Drop alpha
                     as.numeric() %>%
                     {grDevices::rgb(.[1], .[2], .[3], maxColorValue = 255)})
  stats::setNames(item_colors, item_names)
}
