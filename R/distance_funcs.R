# Distance functions

#' For each phenotype in a single sample,
#' find the distance from
#' each cell to the nearest other cell in the phenotype.
#' @param csd A data frame with \code{Cell X Position},
#'        \code{Cell Y Position} and \code{Phenotype} columns,
#'        such as the result of calling
#'        \code{\link{read_cell_seg_data}}.
#' @param phenotypes Optional list of phenotypes to include. If omitted,
#' \code{unique(csd$Phenotype)} will be used.
#' @return A data_frame containing a 'Distance to <phenotype>' column
#' for each phenotype. Will contain NA values where there is no other cell
#' of the phenotype.
#' @export
#' @examples
#' # Compute distance columns and append them to the source data
#' d = sample_cell_seg_data
#' d = cbind(d, find_nearest_distance(d))
find_nearest_distance = function(csd, phenotypes=NULL)
{
  stopifnot('Phenotype' %in% names(csd))

  # Check for multiple samples, this is probably an error
  if (length(unique((csd$`Sample Name`)))>1)
    stop('Data appears to contain multiple samples.')

  d = distance_matrix(csd)

  if (is.null(phenotypes))
    phenotypes = sort(unique(csd$Phenotype))
  stopifnot(length(phenotypes) > 0)

  result = lapply(phenotypes, FUN=function(phenotype)
  {
    # Which cells are in the target phenotype?
    phenotype_cells = csd$Phenotype==phenotype
    if (sum(phenotype_cells)>0)
    {
      # Find the minimum distance > 0; i.e. from cells to not-self cells
      phenotype_mins = apply(d[,phenotype_cells, drop=FALSE], 1, row_min)
    }
    else
    {
      phenotype_mins = rep(NA, nrow(csd))
    }
    phenotype_mins
  })
  names(result) = paste('Distance to', phenotypes)       # The names for the new columns

  tibble::as_tibble(result)
}

# Find the minimum value > 0 in row
# If none, return NA
row_min = function(row)
{
  row = row[row>0]
  if (length(row) > 0) min(row) else NA
}

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
  as.matrix(stats::dist(csd[,c('Cell X Position', 'Cell Y Position')]))
}

#' Subset the rows and columns of a distance matrix.
#' @param dst The distance matrix corresponding to \code{csd},
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
subset_distance_matrix = function(dst, csd, row_selection, col_selection) {
  rows = select_rows(csd, row_selection)
  cols = select_rows(csd, col_selection)
  dst[rows, cols, drop=FALSE]
}
