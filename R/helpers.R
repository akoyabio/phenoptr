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
  col = dplyr::if_else('Annotation ID' %in% names(csd),
                       'Annotation ID', 'Sample Name')
  if (!col %in% names(csd))
    stop('Cell seg table does not include "Sample Name" or "Annotation ID".')
  col
}

# Workaround for inForm data that was originally pixels. In that case
# field_data will have the origin at top left; convert to slide origin
# to match the field_info.
# Assume that the data was converted with the default 2 pixels per micron.
# Test for all cells in a field fitting within the field size after
# correcting the location.
correct_for_pixel_data = function(field_data, field_info) {
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
