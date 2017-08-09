#' Find inForm data files.
#'
#' `list_cell_seg_files` finds inForm cell seg data files in a single directory
#' or a directory hierarchy.
#'
#' @param path Path to the base directory to search.
#' @param ... Additional arguments passed to [list.files][base::list.files()].
#' Pass `recursive=TRUE` to search a directory hierarchy.
#' @return A list of  file paths.
#' @export
#' @family file readers
#' @md
list_cell_seg_files <- function(path, ...) {
  list.files(path, pattern='cell_seg_data.txt', full.names=TRUE, ...)
}

#' Read and clean an inForm data file.
#'
#' \code{read_cell_seg_data} makes it easier to use data from PerkinElmer's
#' inForm program. It reads data files written by inForm 2.0 and later and does
#' useful cleanup on the result. Data files written by inForm 2.0 can be read
#' easily using \code{\link[utils]{read.delim}} or
#' \code{\link[readr]{read_tsv}}. However there is still some useful cleanup to
#' be done.
#'
#' \code{read_cell_seg_data} reads both single-image tables and merged tables
#' and does useful cleanup on the data:
#' \itemize{
#' \item Removes columns that are all NA
#'       (these are typically unused summary columns)
#' \item Converts percent columns to numeric fractions
#' \item Converts pixel distances to microns. The conversion factor may be
#' specfied as a parameter or by setting
#' \code{options(phenoptr.pixels.per.micron)}.
#' \item Optionally removes units from expression names
#' \item If the file contains multiple sample names,
#'       a \code{tag} column is created
#'       containing a minimal, unique tag for each sample.
#'       This is useful when a
#'       short name is needed, for example in chart legends.
#' }
#' @param path Path to the file to read, or NA to use a file chooser.
#' @param pixels_per_micron Conversion factor to microns
#'        (default 2 pixels/micron),
#'        set to NA to skip conversion.
#' @param remove_units If TRUE (default),
#'        remove the unit name from expression columns.
#' @return A \code{\link[tibble]{data_frame}}
#'         containing the cleaned-up data set.
#' @export
#' @family file readers
#' @examples
#' path <- sample_cell_seg_path()
#' csd <- read_cell_seg_data(path)
#'
#' # count all the phenotypes in the data
#' table(csd$Phenotype)
#'
#' \dontrun{
#' # Use purrr::map_df to read all cell seg files in a directory
#' # and return a single data_frame.
#' paths <- list_cell_seg_files(path)
#' csd <- purrr::map_df(paths, read_cell_seg_data)
#' }
read_cell_seg_data <- function(
  path=NA,
  pixels_per_micron=getOption('phenoptr.pixels.per.micron'),
  remove_units=TRUE) {
  if (is.na(path)) {
    path <- file.choose() # nocov not going to happen...
    cat('Loading', path)  # nocov
  }
  if (path=='')
    stop("File name is missing.")

  # Read the data. Supplying col_types prevents output of the imputed types
  df <- readr::read_tsv(path, na=c('NA', '#N/A'), col_types=readr::cols())

  sample_name = 'Sample Name'

  # If there are multiple sample names, make 'tag' be an abbreviated
  # Sample.Name column and insert it as the first column
  # Use the 'tag' column when you need a short name for the sample,
  # e.g. in chart legends
  if (length(unique(df[[sample_name]])) > 1 && !('tag' %in% names(df))) {
    tag <- as.factor(remove_extensions(remove_common_prefix(df[[sample_name]])))
    df <- cbind(tag, df)
  }

  # Convert percents if it is not done already
  pcts <- grep('percent|confidence', names(df), ignore.case=TRUE)
  for (i in pcts)
    if (is.character(df[[i]]))
        df[[i]] <- as.numeric(sub('\\s*%$', '', df[[i]]))/100

  # Remove columns that are all NA or all blank
  # The first condition is a preflight that speeds this up a lot
  na_columns <- purrr::map_lgl(df, function(col)
    (is.na(col[1]) | col[1]=='') && all(is.na(col) | col==''))
  df <- df[!na_columns]

  # Convert distance to microns.
  # No way to tell in general if this was done already...
  if (!is.na(pixels_per_micron)) {
    cols = get_pixel_columns(df)
    for (col in cols) {
      # If col has a lot of NAs it may have been read as chr
      if (!is.numeric(df[[col]]))
        df[[col]] = as.numeric(df[[col]])
      df[[col]] = df[[col]] / pixels_per_micron
      names(df)[col] = sub('pixels', 'microns', names(df)[col])
    }

    cols = get_area_columns(df)
    for (col in cols) {
      if (!is.numeric(df[[col]]))
        df[[col]] = as.numeric(df[[col]])
      df[[col]] = df[[col]] / (pixels_per_micron^2)
      names(df)[col] = sub('pixels', 'sq microns', names(df)[col])
    }

    cols = get_density_columns(df)
    for (col in cols) {
      if (!is.numeric(df[[col]]))
        df[[col]] = as.numeric(df[[col]])
      df[[col]] = df[[col]] * (pixels_per_micron^2)
      names(df)[col] = sub('megapixel', 'sq mm', names(df)[col])
    }
  }

  if (remove_units) {
    unit_name = stringr::str_match(names(df), 'Mean( \\(.*\\))$')[, 2]
    unit_name = unit_name[!is.na(unit_name)][1]
    if (!is.na(unit_name))
      names(df) = sub(unit_name, '', names(df), fixed=TRUE)
  }

  dplyr::as_data_frame(df)
}

get_pixel_columns = function(df) {
  rx = 'position|Distance from Process Region Edge|Distance from Tissue Category Edge|axis' # nolint
  grep(rx, names(df), ignore.case=TRUE)
}

get_area_columns = function(df) {
  rx = 'area \\(pixels\\)'
  grep(rx, names(df), ignore.case=TRUE)
}

get_density_columns = function(df) {
  rx = 'megapixel'
  grep(rx, names(df), ignore.case=TRUE)
}


# Remove the common prefix from a vector of strings
remove_common_prefix <- function(x) {
  # Lexicographic min and max
  .min <- min(x)
  .max <- max(x)
  if (.min == .max) return (x)  # All strings are the same

  # Find the first difference by comparing characters
  .split <- strsplit(c(.min, .max), split='')
  suppressWarnings(.match <- .split[[1]] == .split[[2]])
  first_diff <- match(FALSE, .match)

  substring(x, first_diff)
}

# Remove the extensions from a vector of strings
remove_extensions <- function(x) {
  sub('\\.[^.]+$', '', x)
}
