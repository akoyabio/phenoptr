# Distance functions

#' Create a distance matrix from cell seg data.
#' @param d A data frame with \code{Cell X Position} and \code{Cell Y Position}
#'          columns.
#' @return A square matrix with both dimensions equal to \code{nrow(d)}
#' @export
distance_matrix = function(d) {
  stopifnot('Cell X Position' %in% names(d), 'Cell Y Position' %in% names(d))
  as.matrix(dist(d[,c('Cell X Position', 'Cell Y Position')]))
}
