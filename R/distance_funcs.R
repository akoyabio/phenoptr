# Distance functions

# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c(".", "Sample Name"))

#' Nearest neighbors from a file.
#'
#' Compute nearest distance to each phenotype for each cell in a
#' (possibly merged) inForm cell seg table. Add `Distance to <phenotype>`
#' columns.
#' Write the result to a new file.
#'
#' NOTE: The input file is read using [read_cell_seg_data] so the conversions
#' and cleanup it does will be applied to the output data.
#'
#' @param cell_table_path Path to an inForm cell seg data file, or NULL
#' to prompt for the path.
#' @param out_path Path to the output file, or NULL to create a path from the
#' input file path.
#' @importFrom magrittr "%>%"
#' @export
#' @seealso [find_nearest_distance] which performs the distance calculation.
#' @family distance functions
#' @md
compute_all_nearest_distance <- function(cell_table_path=NULL, out_path=NULL) {
  # Get the path to the cell seg table and check it
  if (is.null(cell_table_path))
    cell_table_path = file.choose()

  # Read the table
  cat('Reading', cell_table_path, '\n')
  csd = read_cell_seg_data(cell_table_path)

  # Compute the distances
  cat('Computing distances\n')
  result = NULL
  phenos = unique_phenotypes(csd)
  result = csd %>%
    dplyr::group_by(!!rlang::sym(field_column(csd))) %>%
    dplyr::do(dplyr::bind_cols(., find_nearest_distance(., phenos)))

  if (is.null(out_path))
    out_path = sub('\\.txt$', '_dist.txt', cell_table_path)
  cat('Writing', out_path, '\n')
  readr::write_tsv(result, out_path, na='#N/A')
}

#' Nearest neighbor distances for each cell and phenotype.
#'
#' For each cell in a single sample,
#' find the distances from
#' the cell to the nearest neighbor cells in each of the provided phenotypes.
#'
#' If the `rtree` package is available, this will use a fast, memory-efficient
#' algorithm capable of processing fields with many thousand cells. Otherwise,
#' a simple distance matrix algorithm is used. The simple algorithm
#' requires at least 8 * (number of cells)^2 bytes of memory which becomes
#' prohibitive as the number of cells becomes large.
#'
#' Install the `rtree` package from GitHub using the command
#' `devtools::install_github('akoyabio/rtree')`.
#'
#' @param csd A data frame with `Cell X Position`,
#'        `Cell Y Position` and `Phenotype` columns,
#'        such as the result of calling
#'        [read_cell_seg_data].
#' @param phenotypes Optional list of phenotypes to include. If omitted,
#' `unique_phenotypes(csd)` will be used.
#' @param dst Optional distance matrix. If provided, this should be
#' `distance_matrix(csd)`. Not used if `rtree` is available.
#' @return A `tibble` containing a `Distance to <phenotype>` column
#' and `Cell ID <phenotype>` column for each phenotype.
#' Columns will contain `NA` values where there is no other cell
#' of the phenotype.
#' @md
#' @export
#' @seealso [compute_all_nearest_distance] which applies this function to a
#' (possibly merged) data file.
#' @family distance functions
#' @examples
#' # Compute distance columns
#' csd <- sample_cell_seg_data
#' nearest <- find_nearest_distance(csd)
#' dplyr::glimpse(nearest)
#'
#' # Make a combined data frame including original data and distance columns
#' csd <- cbind(csd, find_nearest_distance(csd))
#'
#' \dontrun{
#' # If `merged` is a data frame containing cell seg data from multiple fields,
#' # this code will create a new `tibble` with distance columns computed
#' # for each `Sample Name` in the data.
#' merged_with_distance <- merged %>%
#'   dplyr::group_by(`Sample Name`) %>%
#'   dplyr::do(dplyr::bind_cols(., find_nearest_distance(.)))
#' }
find_nearest_distance <- function(csd, phenotypes=NULL, dst=NULL) {
  stop_if_multiple_fields(csd)

  if (getOption('use.rtree.if.available') &&
      requireNamespace('rtree', quietly=TRUE))
    find_nearest_distance_rtree(csd, phenotypes)
  else
    find_nearest_distance_dist(csd, phenotypes, dst)
}

#' Distance-matrix implementation of `find_nearest_distance`.
#' @param csd A data frame with `Cell X Position`,
#'        `Cell Y Position` and `Phenotype` columns,
#'        such as the result of calling
#'        [read_cell_seg_data].
#' @param phenotypes Optional list of phenotypes to include. If omitted,
#' `unique_phenotypes(csd)` will be used.
#' @param dst Optional distance matrix. If provided, this should be
#' `distance_matrix(csd)`.
#' @return A `tibble` containing a `Distance to <phenotype>` column
#' and `Cell ID <phenotype>` column for each phenotype.
#' Columns will contain `NA` values where there is no other cell
#' of the phenotype.
#' @seealso find_nearest_distance
#' @md
#' @keywords internal
find_nearest_distance_dist = function(csd, phenotypes=NULL, dst=NULL) {

  phenotypes = validate_phenotypes(phenotypes, csd)

  if (is.null(dst))
    dst = distance_matrix(csd)

  # Removing dimnames from dst gives unnamed columns in results
  # The names are clutter, take up space and confuse the tests
  dimnames(dst) = NULL

  result = purrr::map2_dfc(names(phenotypes), phenotypes,
                           function(name, phenotype) {
    # Which cells are in the target phenotype?
    phenotype_cells = select_rows(csd, phenotype)
    if (sum(phenotype_cells)>0) {
      # Subset columns of the distance matrix to just phenotype cells
      phenotype_dist = dst[, phenotype_cells, drop=FALSE]

      # Find the minimum distance > 0; i.e. from cells to not-self cells
      dist_col = apply(phenotype_dist, 1, row_min)

      # Find the index of the minimum distance > 0
      # and use this to index the Cell IDs of the target phenotypes
      which_dist_col = apply(phenotype_dist, 1, which_row_min)
      cell_id_col = csd$`Cell ID`[phenotype_cells][which_dist_col]
      pheno_cols = tibble::tibble(dist_col, cell_id_col)
    }
    else {
      # No cells of the selected phenotype
      na_col = rep(NA_integer_, nrow(csd))
      pheno_cols = tibble::tibble(dist_col=na_col, cell_id_col=na_col)
    }
    pheno_cols  %>%
      rlang::set_names(paste(c('Distance to', 'Cell ID'), name))
  })
}

#' `rtree`-based implementation of `find_nearest_distance`.
#' @param csd A data frame with `Cell X Position`,
#'        `Cell Y Position` and `Phenotype` columns,
#'        such as the result of calling
#'        [read_cell_seg_data].
#' @param phenotypes Optional list of phenotypes to include. If omitted,
#' `unique_phenotypes(csd)` will be used.
#' @return A `tibble` containing a `Distance to <phenotype>` column
#' and `Cell ID <phenotype>` column for each phenotype.
#' Columns will contain `NA` values where there is no other cell
#' of the phenotype.
#' @seealso find_nearest_distance
#' @md
#' @keywords internal
find_nearest_distance_rtree <- function(csd, phenotypes=NULL) {
  phenotypes = validate_phenotypes(phenotypes, csd)
  field_locs = csd %>%
    dplyr::select(X=`Cell X Position`, Y=`Cell Y Position`)
  field_ids = csd$`Cell ID`

  result = purrr::map2_dfc(names(phenotypes), phenotypes,
                           function(name, phenotype) {
    # Which cells are in the target phenotype?
    phenotype_cells = select_rows(csd, phenotype)

    if (sum(phenotype_cells)>0) {
      # Make an rtree of the phenotype cells
      to_cells_locs = field_locs[phenotype_cells, ]
      to_cells_tree = rtree::RTree(as.matrix(to_cells_locs))

      to_cells_ids = field_ids[phenotype_cells]

      # Find nearest neighbor. Get two nearest neighbors so we can
      # discard self when self is a member of the "to" phenotype.
      # Note:
      # It is possible for two distinct cells to have identical positions
      # (when one wraps around the other). This gives two cells with
      # nearest neighbor at distance zero. With n_nn==2 this results in a
      # reported nearest distance of NA. With n_nn==3 the next nearest distance
      # will be reported.
      # This situation is rare, and it's not clear what the "correct"
      # answer is, so we choose n_nn==2 because it is faster.
      n_nn = 2L
      to_cells_nn = rtree::knn(to_cells_tree, as.matrix(field_locs), k=n_nn)

      # rtree::knn() returns a list of vectors. For each "from" cell, it returns
      # a vector of indices of "to" cells. Here we convert the list of vectors
      # to a matrix. This allows fast selection of a single column in the
      # map_dfc loop below.
      # The conversion to matrix will recycle values in any rows that don't have
      # n_nn values; that won't affect the result.
      # Note: If there is only one "to" cell, to_cells_nn will have only one
      # column. Thus the map_dfc below is from 1 to dim(to_cells_nn)[2],
      # not from 1 to n_nn.
      to_cells_nn = do.call(rbind, to_cells_nn)

      # knn gives us the indices of nearest cells, we want distance
      # Look up to_cells_nn in to_cells_locs, combine with field_locs,
      # and compute distance. We have to do this n_nn times and take
      # the minimum > 0.
      # First compute all the distances, giving a list of
      # lists of length dim(to_cells_nn)
      all_distances = purrr::map(1:(dim(to_cells_nn)[2]),
                                      # For each nearest neighbor
        ~{
          # Convert indices of nearest neighbors to locations
          # The neighbors are to_cells.
          nn_indices = to_cells_nn[, .x]
          nn_locs = to_cells_locs[nn_indices, ] %>%
            rlang::set_names(c('X1', 'Y1')) # Avoid duplicate names

          # Combine with original locations and compute distance
          field_locs %>%
            dplyr::bind_cols(nn_locs) %>%
            dplyr::mutate(dist = sqrt( (X-X1)^2 + (Y-Y1)^2 )) %>% # nolint
            dplyr::pull()
        }) %>%
        purrr::transpose() # Get each row as a list

      # The actual minimum distances
      dist_col = all_distances %>%
        purrr::map_dbl(~row_min(unlist(.x)))

      # Get the index to the minimum item > 0 in each list
      which_dist_col = all_distances %>%
        purrr::map_dbl(~which_row_min(unlist(.x))) # Index of Min > 0 or NA

      # Applying which_dist_col to to_cells_nn gives indices into to_cells_ids
      # which then gives us actual cell IDs.
      cell_id_ix = to_cells_nn[cbind(seq_along(which_dist_col), which_dist_col)]
      cell_id_col = to_cells_ids[cell_id_ix]
      pheno_cols = tibble::tibble(dist_col, cell_id_col)
    }
    else {
      # No cells of the selected phenotype
      na_col = rep(NA_integer_, nrow(csd))
      pheno_cols = tibble::tibble(dist_col=na_col, cell_id_col=na_col)
    }
    pheno_cols  %>%
      rlang::set_names(paste(c('Distance to', 'Cell ID'), name))
  })
}

# Find the minimum value > 0 in row
# If none, return NA
row_min = function(row) {
  row = row[!is.na(row) & row>0]
  if (length(row) > 0) min(row) else NA
}

# Find the index of the minimum value > 0 in row
# If none, return NA
which_row_min = function(row) {
  good_values = (!is.na(row) & row>0)
  if (sum(good_values) == 0) return(NA_integer_) # Nothing to see here

  # Replace bad values with Inf so we can find the min of the rest
  row[!good_values] = Inf
  return(which.min(row))
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
distance_matrix <- function(csd) {
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
subset_distance_matrix <- function(csd, dst, row_selection, col_selection) {
  # Check for pre-0.1.0.9002 parameter order
  if (is.matrix(csd) && is.data.frame(dst))
    stop(
      'csd and dst parameters to subset_distance_matrix are in the wrong order')

  rows = select_rows(csd, row_selection)
  cols = select_rows(csd, col_selection)
  dst[rows, cols, drop=FALSE]
}
