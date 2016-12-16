# Distance functions

#' Create a distance matrix from cell seg data.
#' @param csd A data frame with \code{Cell X Position}
#'        and \code{Cell Y Position} columns, such as the result of
#'        \code{\link{read_cell_seg_data}}.
#' @return A square matrix with both dimensions equal to \code{nrow(d)}.
#'         The value at \code{[i, j]} will be the distance from the cell
#'         at row \code{i} of \code{csd} to the cell at row \code{j}.
#'         The returned matrix is symmetric.
#' @export
distance_matrix = function(csd) {
  stopifnot('Cell X Position' %in% names(csd), 'Cell Y Position' %in% names(csd))
  as.matrix(dist(csd[,c('Cell X Position', 'Cell Y Position')]))
}

#' Subset the rows and columns of a distance matrix.
#' @param dist The distance matrix corresponding to \code{csd},
#'        produced by calling \code{\link{distance_matrix}}.
#' @param csd A data frame containing cell segmentation data,
#'        such as the result of
#'        \code{\link{read_cell_seg_data}}.
#' @param row_selection,col_selection Selection criteria for the
#' rows and columns. Accepts all formats accepted by \code{\link{select_rows}}.
#' @return The input matrix \code{dst} subsetted to include only the
#' rows corresponding to \code{row_selection} and columns
#' corresponding to \code{col_selection}.
#' @export
subset_distance_matrix = function(dist, csd, row_selection, col_selection) {
  rows = select_rows(csd, row_selection)
  cols = select_rows(csd, col_selection)
  dist[rows, cols, drop=FALSE]
}
