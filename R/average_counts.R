# Functions to compute average counts

#' Find the number of cells within a radius
#'
#' Find the average number of cells of type \code{to} within \code{radius}
#' of \code{from} cells
#' within tissue category \code{category}.
#' @param csd A data frame with \code{Cell X Position},
#'        \code{Cell Y Position} and \code{Phenotype} columns,
#'        such as the result of calling \code{\link{read_cell_seg_data}}.
#' @param from,to Selection criteria for the
#' rows and columns. Accepts all formats accepted by \code{\link{select_rows}}.
#' @param radius The radius to search within.
#' @param category Optional tissue category to restrict both \code{from} and
#' \code{to}.
#' @param dst Optional distance matrix corresponding to \code{csd},
#'        produced by calling \code{\link{distance_matrix}}.
#' @return A \code{\link{data_frame}} with four columns:
#'   \describe{
#'    \item{\code{from_count}}{The number of \code{from} cells found.}
#'    \item{\code{to_count}}{The number of \code{to} cells found.}
#'    \item{\code{within_count}}{The number of \code{to} cells found within
#'    \code{radius} microns of a \code{from} cell.}
#'    \item{\code{within_mean}}{The average number of \code{to} cells found
#'    within \code{radius} microns of a \code{from} cell.}
#'  }
#' @export
#' @family distance functions
average_count_within = function(csd, from, to, radius, category=NULL, dst=NULL)
{
  if (is.null(dst))
    dst = distance_matrix(csd)

  # Which cells are in the from and to phenotypes?
  if (!is.null(category)) {
    tissue_f = stats::as.formula(paste0('~`Tissue Category`=="', category, '"'))
    if (is.character(from)) from = list(from)
    from = c(from, list(tissue_f))
    if (is.character(to)) to = list(to)
    to = c(to, list(tissue_f))
  }

  dst = subset_distance_matrix(csd, dst, from, to)
  if (prod(dim(dst))>0) {
    within = apply(dst, 1, function(r) sum(r<=radius))
    tibble::data_frame(
      from_count = dim(dst)[1],
      to_count = dim(dst)[2],
      within_count = sum(within),
      within_mean = mean(within)
    )
  } else {
    tibble::data_frame(
      from_count = dim(dst)[1],
      to_count = dim(dst)[2],
      within_count = 0,
      within_mean = 0
    )
  }
}
