# Distance functions

# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c(".", "Sample Name"))

#' Nearest neighbors from a file.
#'
#' Compute nearest distance to each phenotype for each cell in a
#' (possibly merged) inForm cell seg table. Write the result to a new file.
#' @param cell_table_path Path to an inForm cell seg data file, or NULL
#' to prompt for the path.
#' @param out_path Path to the output file, or NULL to create a path from the
#' input file path.
#' @importFrom magrittr "%>%"
#' @export
#' @family distance functions
compute_all_nearest_distance = function(cell_table_path=NULL, out_path=NULL) {
  # Get the path to the cell seg table and check it
  if (is.null(cell_table_path))
    cell_table_path = file.choose()

  # Read the table
  cat('Reading', cell_table_path, '\n')
  csd = read_cell_seg_data(cell_table_path)

  # Compute the distances
  cat('Computing distances\n')
  result = NULL
  phenos = unique(csd$Phenotype)
  result = csd %>%
    dplyr::group_by(`Sample Name`) %>%
    dplyr::do(dplyr::bind_cols(., find_nearest_distance(., phenos)))

  if (is.null(out_path))
    out_path = sub('\\.txt$', '_dist.txt', cell_table_path)
  cat('Writing', out_path, '\n')
  readr::write_tsv(result, out_path, na='#N/A')
}

#' Nearest neighbor distances for each cell and phenotype.
#'
#' For each phenotype in a single sample,
#' find the distance from
#' each cell to the nearest other cell in the phenotype.
#' @param csd A data frame with \code{Cell X Position},
#'        \code{Cell Y Position} and \code{Phenotype} columns,
#'        such as the result of calling
#'        \code{\link{read_cell_seg_data}}.
#' @param phenotypes Optional list of phenotypes to include. If omitted,
#' \code{unique(csd$Phenotype)} will be used.
#'
#' @return A data_frame containing a 'Distance to <phenotype>' column
#' for each phenotype. Will contain NA values where there is no other cell
#' of the phenotype.
#' @export
#' @family distance functions
#' @examples
#' # Compute distance columns and append them to the source data
#' csd = sample_cell_seg_data
#' csd = cbind(csd, find_nearest_distance(csd))
find_nearest_distance = function(csd, phenotypes=NULL) {
  stopifnot('Phenotype' %in% names(csd))

  # Check for multiple samples, this is probably an error
  if (length(unique(csd$`Sample Name`))>1)
    stop('Data appears to contain multiple samples.')

  dst = distance_matrix(csd)

  if (is.null(phenotypes))
    phenotypes = sort(unique(csd$Phenotype))
  stopifnot(length(phenotypes) > 0)

  result = lapply(phenotypes, FUN=function(phenotype) {
    # Which cells are in the target phenotype?
    phenotype_cells = csd$Phenotype==phenotype
    if (sum(phenotype_cells)>0) {
      # Find the minimum distance > 0; i.e. from cells to not-self cells
      phenotype_mins = apply(dst[, phenotype_cells, drop=FALSE], 1, row_min)
    }
    else {
      phenotype_mins = rep(NA, nrow(csd))
    }
    phenotype_mins
  })
  # The names for the new columns
  names(result) = paste('Distance to', phenotypes)

  tibble::as_tibble(result)
}

# Find the minimum value > 0 in row
# If none, return NA
row_min = function(row) {
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
#' @family distance functions
#' @export
distance_matrix = function(csd) {
  stopifnot('Cell X Position' %in% names(csd),
            'Cell Y Position' %in% names(csd))
  as.matrix(stats::dist(csd[, c('Cell X Position', 'Cell Y Position')]))
}

#' Subset the rows and columns of a distance matrix.
#' @param csd A data frame containing cell segmentation data,
#'        such as the result of
#'        \code{\link{read_cell_seg_data}}.
#' @param dst The distance matrix corresponding to \code{csd},
#'        produced by calling \code{\link{distance_matrix}}.
#' @param row_selection,col_selection Selection criteria for the
#' rows and columns. Accepts all formats accepted by
#' \code{\link{select_rows}}.
#' @return The input matrix \code{dst} subsetted to include only the
#' rows corresponding to \code{row_selection} and columns
#' corresponding to \code{col_selection}.
#' @family distance functions
#' @export
subset_distance_matrix = function(csd, dst, row_selection, col_selection) {
  # Check for pre-0.1.0.9002 parameter order
  if (is.matrix(csd) && is.data.frame(dst))
    stop('csd and dst parameters to subset_distance_matrix are in the wrong order')

  rows = select_rows(csd, row_selection)
  cols = select_rows(csd, col_selection)
  dst[rows, cols, drop=FALSE]
}
