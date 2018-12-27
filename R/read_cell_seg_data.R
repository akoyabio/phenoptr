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
#' \code{read_cell_seg_data} makes it easier to use data from Akoya Biosciences'
#' inForm program. It reads data files written by inForm 2.0 and later and does
#' useful cleanup on the result. Data files written by inForm 2.0 can be read
#' easily using \code{\link[utils]{read.delim}} or
#' \code{\link[readr]{read_tsv}}. However there is still some useful cleanup to
#' be done.
#'
#' \code{read_cell_seg_data} reads both single-image tables and merged tables
#' and does useful cleanup on the data:
#' \itemize{
#' \item Removes columns that are all NA.
#'       These are typically unused summary columns.
#' \item Converts percent columns to numeric fractions.
#' \item Converts pixel distances to microns. The conversion factor may be
#' specified as a parameter, by setting
#' \code{options(phenoptr.pixels.per.micron)}, or by reading an associated
#' \code{component_data.tif} file.
#' \item Optionally removes units from expression names
#' \item If the file contains multiple sample names,
#'       a \code{tag} column is created
#'       containing a minimal, unique tag for each sample.
#'       This is useful when a
#'       short name is needed, for example in chart legends.
#' }
#'
#' If \code{pixels_per_micron='auto'}, \code{read_cell_seg_data} looks for
#' a \code{component_data.tif} file in the same directory as \code{path}.
#' If found, \code{pixels_per_micron} is read from the file \strong{and}
#' the cell coordinates are offset to the correct spatial location.
#'
#' @param path Path to the file to read, or NA to use a file chooser.
#' @param pixels_per_micron Conversion factor to microns
#'        (default 2 pixels/micron, the resolution of 20x MSI fields
#'        taken on Vectra Polaris and Vectra 3.).
#'        Set to NA to skip conversion. Set to \code{'auto'} to read from
#'        an associated \code{component_data.tif} file.
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

  # If any of these fields has missing values, the file may be damaged.
  no_na_cols = c("Path", "Sample Name", "Tissue Category", "Phenotype",
                 "Cell ID", "Cell X Position", "Cell Y Position", "Slide ID")

  bad_na_cols =
    purrr::keep(no_na_cols, ~(.x %in% names(df) && any(is.na(df[[.x]]))))
  if (length(bad_na_cols) > 0)
    warning('Some expected columns have missing data:\n',
            paste(bad_na_cols, collapse=', '), '\n',
            path, ' may be damaged.')
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

  # Convert distance to microns if requested.
  # No way to tell in general if this was done already...
  # Try looking for 'micron' in column names
  if (unit_is_microns(df) && !is.na(pixels_per_micron)) {
    message('Data is already in microns, no conversion performed')
  } else if (!is.na(pixels_per_micron)) {
    if (pixels_per_micron=='auto') {
      # Get pixels_per_micron and field location from component_data.tif
      component_path = sub('_cell_seg_data.txt', '_component_data.tif', path)
      if(!file.exists(component_path))
        stop('To convert cell seg data to microns, a matching component data file is required.')
      info = get_field_info(component_path)
      pixels_per_micron = 1/info$microns_per_pixel
      location = info$location
    } else {
      location = NA
    }

    cols = get_pixel_columns(df)
    for (col in cols) {
      # If col has a lot of NAs it may have been read as chr
      if (!is.numeric(df[[col]]))
        df[[col]] = as.numeric(df[[col]])
      df[[col]] = df[[col]] / pixels_per_micron
      names(df)[col] = sub('pixels', 'microns', names(df)[col])
    }

    # Position columns get an offset, if available
    if (!anyNA(location)) {
      df$`Cell X Position` = df$`Cell X Position` + location[1]
      df$`Cell Y Position` = df$`Cell Y Position` + location[2]
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

# Convert cell locations back to pixels if possible
force_pixel_locations = function(csd, cell_seg_path) {
  if (unit_is_microns(csd)) {
    # The cell seg file's native unit is microns. Get image info
    # from the component data file and convert back to pixels
    component_path = sub('_cell_seg_data.txt', '_component_data.tif', cell_seg_path)
    if(!file.exists(component_path))
      stop('For cell seg data in microns, a matching component data file is required.')
    info = get_field_info(component_path)
    pixels_per_micron = 1/info$microns_per_pixel
    location = info$location
    csd %>% dplyr::mutate(
      `Cell X Position` = (`Cell X Position`-location[1])*pixels_per_micron,
      `Cell Y Position` = (`Cell Y Position`-location[2])*pixels_per_micron)
  } else {
    csd
  }
}

unit_is_microns = function(df) {
  any(stringr::str_detect(names(df), 'micron'))
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
