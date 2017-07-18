# Functions to compute counts and average counts

# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c(".", "Phenotype", "Tissue Category"))

#' Count cells within a radius for multiple tissue categories, phenotypes
#' and fields.
#'
#' This is a batch version of [count_within]. Given the path to a directory
#' containing cell seg data files, for each given tissue category,
#' 'from' phenotype, 'to' phenotype and radius, it counts the number of
#'  `from` cells
#' having a `to` cell within `radius` microns.
#'
#' `from` and `to` specify the phenotypes of interest.
#' Each one is a list of phenotype names.
#' Each list element is either a single name or a character vector.
#' If the element is a vector,
#' the phenotypes in the vector are treated collectively as a single type.
#' All pairs of `from` and `to` values are included in the output.
#'
#' For example, the parameters
#' `from=list('cd68'), to=list('cd8', 'tumor'), radius=25`
#' would count `cd68` cells having a `cd8` cell within 25 &mu;m and,
#' separately, `cd68` cells with a `tumor` cell within 25 &mu;m.
#'
#' The `category` parameter may be a single category or a list of categories.
#'
#' @param base_path Path to a directory containing at least
#' one `_cell_seg_data.txt` file.
#' @param from Specification of 'from' phenotype(s), see Details.
#' @param to Specification of 'to' phenotype(s), see Details.
#' @param radius The radius or radii to search within.
#' @param category Optional tissue categories to restrict both `from` and
#' `to`.
#' @param verbose If TRUE, display progress.
#' @return A `data_frame` containing these columns:
#'   \describe{
#'    \item{\code{slide_id}}{Slide ID from the data files.}
#'    \item{\code{source}}{Base file name of the source file with
#'    `_cell_seg_data.txt` stripped off for brevity.}
#'    \item{\code{category}}{Tissue category, if provided as a parameter,
#'    or "all".}
#'    \item{\code{from}}{From phenotype.}
#'    \item{\code{to}}{To phenotype.}
#'    \item{\code{radius}, \code{from_count}, \code{to_count},
#'    \code{from_with}, \code{within_mean}}{Results from [count_within]
#'    for this data file and tissue category.}
#'  }
#' @examples
#' base_path = system.file("extdata", "TMA", package = "phenoptr")
#'
#' # Count tumor cells near macrophages and tumor cells near CD8 separately,
#' # in tumor and stroma tissue categories separately.
#' from = list('tumor')
#' to = list('macrophage CD68', 'cytotoxic CD8')
#' radius = c(10, 25)
#' category = list('tumor', 'stroma')
#' count_within_batch(base_path, from, to, radius, category)
#'
#' # Count tumor cells near any T cell in all tissue categories
#' to = list(c('cytotoxic CD8', 'helper CD4', 'T reg Foxp3'))
#' count_within_batch(base_path, from, to, radius)
#' @md
#' @export
#' @family distance functions
#' @importFrom magrittr "%>%"
count_within_batch = function(base_path, from, to, radius,
                              category=NA, verbose=TRUE) {
  files = list_cell_seg_files(base_path)
  if (length(files) == 0)
    stop('No cell seg files found in ', base_path)

  stopifnot(is.list(from), length(from) > 0)
  stopifnot(is.list(to), length(to) > 0)

  all_phenotypes = c(from, to, recursive=TRUE)
  if (any(purrr::map_lgl(all_phenotypes, purrr::is_formula)))
    stop('Formula selection is not supported in count_within.')

  combos = purrr::cross_n(list(from=from, to=to, category=category))

  # Loop through all the cell seg data files
  purrr::map_df(files, function(path) {
    name = basename(path) %>% sub('_cell_seg_data.txt', '', .)
    if (verbose) cat('Processing', name, '\n')

    # Read one file
    csd = read_cell_seg_data(path)

    slide = as.character(csd[1, 'Slide ID'])

    # Subset to what we care about, for faster distance calculation
    csd = csd %>% dplyr::filter(Phenotype %in% all_phenotypes)
    if (!anyNA(category))
      csd = csd %>% dplyr::filter(`Tissue Category` %in% category)

    # Compute the distance matrix for these cells
    dst = distance_matrix(csd)

    # Compute counts for each from, to, and category in combos
      purrr::map_df(combos, function(row) {
      # Call count_within for each item in combos
      # count_within handles multiple radii
      purrr::invoke(.f=count_within, .x=row,
                    csd=csd, radius=radius, dst=dst) %>%
        # Add columns for from, to, category
        tibble::add_column(
          category = ifelse(is.na(row$category), 'all', row$category),
          from=paste(row$from, collapse=' '),
          to=paste(row$to, collapse=' '),
          .before=1)
      }) %>%
      # Add columns for slide and source
      tibble::add_column(slide_id=slide,
                    source=name,
                    .before=1)

  })
}

#' Count cells within a radius.
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
#' the given tissue category. If \code{category} is NA, values are reported
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
#' csd = sample_cell_seg_data
#'
#' # Find the number of macrophages with a tumor cell within 10 or 25 microns
#' count_within(csd, from='macrophage CD68', to='tumor', radius=c(10, 25)) %>%
#'   mutate(not_to=from_count*within_mean, to_mean=not_to/to_count)
#'
#' # Find the number of tumor cells with a macrophage within 10 or 25 microns
#' # Show that from_count*within_mean is not the same as from_with in the
#' # previous computation but from_count*within_mean/to_count is the same.
#' count_within(csd, from='tumor', to='macrophage CD68', radius=c(10, 25)) %>%
#'   mutate(not_to=from_count*within_mean, to_mean=not_to/to_count)

count_within = function(csd, from, to, radius, category=NA, dst=NULL) {
  # Check for multiple samples, this is probably an error
  if (length(unique(csd$`Sample Name`))>1)
    stop('Data appears to contain multiple samples.')

  stopifnot(length(radius) > 0, all(radius>0))

  if (is.null(dst))
    dst = distance_matrix(csd)

  # Which cells are in the from and to phenotypes?
  if (!is.na(category)) {
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
