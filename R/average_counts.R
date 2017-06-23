# Functions to compute counts and average counts

#' Count cells within a radius
#'
#' Count the number of \code{from} cells having a \code{to} cell within
#' \code{radius} microns in tissue category \code{category}.
#' Compute the average number of \code{to} cells
#' within \code{radius} of \code{from} cells.
#'
#' For each \code{from} cell, count the number of \code{to} cells within
#' \code{radius} microns. Report the number of \code{from} cells containing
#' at least \emph{one} \code{to} cell within \code{radius} as \code{from_with}.
#' Report the \emph{average} number of \code{to} cells per
#' \code{from} cell as \code{within_mean}.
#'
#' There are some subtleties to this calculation. See the example code
#' for a demonstration.
#' \itemize{
#'   \item It is not symmetric in \code{from} and \code{to}.
#'   For example the number of tumor cells with a
#'   macrophage within 25 microns is not the same as the number of macrophages
#'   with a tumor cell within 25 microns.
#'   \item \code{from_count*within_mean} is \emph{not} the number of
#'   \code{to} cells within \code{radius} of a \code{from} cell, it may
#'   count \code{to} cells multiple times.
#'   \item Surprisingly, \code{from_count*within_mean} is symmetric in
#'   \code{from} and \code{to}. The double-counting works out.
#' }
#' To aggregate \code{within_mean} across multiple samples (e.g. by Slide ID)
#' use code such as
#' \preformatted{
#' results \%>\% group_by(`Slide ID`, from, to, radius) \%>\%
#'   summarize(count=sum(from_count),
#'             within=sum(from_count*within_mean),
#'             avg=within/count) \%>\%
#'   ungroup
#' }
#'
#' If \code{category} is specified, all reported values are for cells within
#' the given tissue category. If \code{category} is NULL, values are reported
#' for the entire data set.
#'
#' \code{radius} may be a vector with multiple values.
#'
#' @param csd A data frame with \code{Cell X Position},
#'        \code{Cell Y Position} and \code{Phenotype} columns,
#'        such as the result of calling \code{\link{read_cell_seg_data}}.
#' @param from,to Selection criteria for the
#' rows and columns. Accepts all formats accepted by \code{\link{select_rows}}.
#' @param radius The radius or radii to search within.
#' @param category Optional tissue category to restrict both \code{from} and
#' \code{to}.
#' @param dst Optional distance matrix corresponding to \code{csd},
#'        produced by calling \code{\link{distance_matrix}}.
#'
#' @return A \code{\link{data_frame}} with five columns and one row for each
#'   value in \code{radius}:
#'   \describe{
#'    \item{\code{radius}}{The value of \code{radius} for this row.}
#'    \item{\code{from_count}}{The number of \code{from} cells found in
#'     \code{csd}.}
#'    \item{\code{to_count}}{The number of \code{to} cells found in \code{csd}.}
#'    \item{\code{from_with}}{The number of \code{from} cells with a
#'    \code{to} cell within \code{radius}.}
#'    \item{\code{within_mean}}{The average number of \code{to} cells found
#'    within \code{radius} microns of each \code{from} cell.}
#'  }
#' @export
#' @family distance functions
#' @examples
#' library(tidyverse)
#' s = sample_cell_seg_data
#'
#' # Find the number of macrophages with a tumor cell within 10 or 25 microns
#' count_within(s, from='macrophage CD68', to='tumor', radius=c(10, 25)) %>%
#'   mutate(not_to=from_count*within_mean, to_mean=not_to/to_count)
#'
#' # Find the number of tumor cells with a macrophage within 10 or 25 microns
#' # Show that from_count*within_mean is not the same as from_with in the
#' # previous computation but from_count*within_mean/to_count is the same.
#' count_within(s, from='tumor', to='macrophage CD68', radius=c(10, 25)) %>%
#'   mutate(not_to=from_count*within_mean, to_mean=not_to/to_count)

count_within = function(csd, from, to, radius, category=NULL, dst=NULL) {
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
    purrr::map_df(radius, function(rad) {
      within = apply(dst, 1, function(r) sum(r<=rad))
      tibble::data_frame(
        radius = rad,
        from_count = dim(dst)[1], # Number of from cells
        to_count = dim(dst)[2],   # Number of to cells
        from_with = sum(within>0), # Number of from cells having a
                                   # to cell within radius
        within_mean = mean(within) # Mean number of to cells within
                                   # radius of a from cell
      )}
    )
  } else {
    tibble::data_frame(
      radius = radius,
      from_count = dim(dst)[1],
      to_count = dim(dst)[2],
      from_with = 0,
      within_mean = 0
    )
  }
}
