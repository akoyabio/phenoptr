# General utility functions both public and private

# Check if a function is available in a package
# @param package Name of the package
# @param fun Name of the function
# @return TRUE if the package is installed and contains the function.
function_exists =function(package, fun) {
  requireNamespace(package, quietly=TRUE) &&
    (fun %in% getNamespaceExports(package))
}

#' Get the name of the column that distinguishes fields in a merged cell
#' seg data file.
#' @param csd Cell seg data
#' @return A column name (as a string)
#' @export
field_column = function(csd) {
  col = field_column_(csd)
  if (!col %in% names(csd))
    stop('Cell seg table does not include "Sample Name" or "Annotation ID".')
  col
}

#' Stop if `csd` contains multiple fields
#' @keywords internal
stop_if_multiple_fields = function(csd) {
  # Get the field column name
  col = field_column_(csd)

  # Note: If no field column, let it go
  if (col %in% names(csd) && length(unique(csd[[col]]))>1)
    stop('Data contains multiple samples, ',
         'please select one or set whole_slide=TRUE.')
}

# Just get the field column, don't check anything
field_column_ = function(csd) {
  dplyr::if_else('Annotation ID' %in% names(csd),
                 'Annotation ID', 'Sample Name')
}

# Workaround for inForm data that was originally pixels. In that case
# field_data will have the origin at top left; convert to slide origin
# to match the field_info.
# Assume that the data was converted with the default 2 pixels per micron.
# Test for all cells in a field fitting within the field size after
# correcting the location.
correct_for_pixel_data = function(field_data, field_info) {
  if (is.null(field_data) || nrow(field_data)==0) return(field_data)

  correction = 2 / field_info$pixels_per_micron
  if (max(field_data$`Cell X Position`)*correction < field_info$field_size[1]
      && max(field_data$`Cell Y Position`)*correction < field_info$field_size[2]) {
    field_data = field_data %>% dplyr::mutate(
      `Cell Y Position` = `Cell Y Position`*correction+field_info$location[2],
      `Cell X Position` = `Cell X Position`*correction+field_info$location[1]
    )
  }
  field_data
}

#' Read directory info from a TIFF file
#' @param path Path to a TIFF file
#' @param all  If `all=TRUE` then directory info from all images is
#' returned in a list of lists.
#' If `all` is a vector, it gives the (1-based) indices of info to return.
#' Otherwise only info for the first image is returned.
#' @return A list (if `all==FALSE` or length 1) or list of lists
#' (if `all==TRUE`) containing
#' directory info for the images in the file.
#' @export
readTIFFDirectory = function(path, all=FALSE) {
  # This function is a shim that calls either readTIFFDirectory
  # from akoyabio/tiff or readTIFF from the latest s-u/tiff
  if (function_exists('tiff', 'readTIFFDirectory')) {
    # akoyabio/tiff has the same API as this function
    tiff::readTIFFDirectory(path, all=all)
  } else if ('payload' %in% names(formals(tiff::readTIFF))) {
    # s-u/tiff returns a data.frame
    info = tiff::readTIFF(path, all=all, payload=FALSE) %>%
      purrr::transpose()
    if (!all)
        info = info[[1]]
    info
  } else
    stop('Please install a more recent tiff package.')
}
