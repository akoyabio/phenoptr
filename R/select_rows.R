#' Validate a user-specified phenotype definition
#'
#' Check that the specified phenotype
#' can be formed from available phenotypes and data columns.
#' @param pheno Text description of a phenotype,
#' for `phenoptr::parse_phenotypes`.
#' @param available A character vector of available phenotypes
#' @param csd If supplied, any formula arguments will be checked for
#' validity against `csd`.
#' @return An error message or empty string
#' @export
validate_phenotype_definitions = function(pheno, available, csd=NULL) {
  if (is.null(pheno) || pheno==''
      || stringr::str_detect(pheno, 'Total|All'))
    return('')

  phenos = stringr::str_split(pheno, '[,/]')[[1]] %>%
    stringr::str_trim()

  if (!all(stringr::str_detect(phenos, '^~|[+-]$')))
    return('Phenotype definitions must start with ~ or end with + or -.')

  # Check non-formula phenotypes
  pheno_strings = purrr::discard(phenos, ~startsWith(.x, '~'))

  pheno_strings = stringr::str_remove(pheno_strings, '[+-]$')
  missing = !pheno_strings %in% available
  if (any(missing))
    return(paste0('Unknown phenotype(s): ', paste(phenos[missing], sep=', ')))

  # Check formula expressions
  pheno_formulae = purrr::keep(phenos, ~startsWith(.x, '~'))
  if (length(pheno_formulae) > 0 && stringr::str_detect(pheno, ','))
    return("Formula expressions are not allowed in phenotypes combined with ','.")

  for (fmla_str in pheno_formulae) {
    fmla = try(stats::as.formula(fmla_str, globalenv()), silent=TRUE)
    if (class(fmla) == 'try-error')
      return(paste0(fmla_str, ' is not a valid expression.'))

    if (!is.null(csd)) {
      attempt = try(lazyeval::f_eval(fmla, csd), silent=TRUE)
      if (class(attempt) == 'try-error') {
        # This is usually a misspelled variable name
        msg = stringr::str_split(attempt, ' : ')[[1]][[2]]
        return(paste0('Invalid expression: ', stringr::str_trim(msg), '.'))
      } else if (class(attempt) != 'logical' || length(attempt) != nrow(csd)) {
        # This catches other errors such as ~D (returns stats::D)
        # or ~~D (returns the formula ~D)
        return(paste0('Invalid expression: ~', lazyeval::f_text(fmla)))
      }
    }
  }
  return('')
}

#' Parse a vector of phenotype names
#'
#' This helper function takes a user-friendly list of single and
#' multiple phenotype names and converts it to a named list of phenotype
#' selectors for use with [phenoptr::select_rows]. By using `parse_phenotypes`
#' a user does not have to know the (somewhat inscrutable)
#' details of `select_rows`.
#'
#' @param ... Phenotypes to be decoded, or a list of same,
#' optionally with names.
#' @return A named list of phenotype selectors for use with
#'   [phenoptr::select_rows].
#' @section Details:
#' Each phenotype must be either a single phenotype name (e.g. CD3+ or CD8-)
#' or two or more names separated by a slash (/) or comma (,).
#'
#' Phenotypes containing slashes are interpreted as requiring *all* of the
#' individual phenotypes. For example, "CD3+/CD8-" is a CD3+ cell which is
#' also CD8-.
#'
#' Phenotypes containing commas are interpreted as requiring *any* of the
#' individual phenotypes. For example, "CD68+,CD163+" is a cell which is
#' either CD68+ or CD163+ or both.
#'
#' Additionally,
#' - A phenotype without a + or - and containing
#' either "Total" or "All" will be
#' interpreted as meaning "All cells".
#' - A phenotype starting with '~' will be interpreted as a formula expression.
#'   Formulas may be standalone phenotypes or combined with slash (/); they
#'   cannot be combined with comma (,).
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#' # Create selectors for
#' # - All CD3+ cells
#' # - CD3+/CD8+ double-positive cells
#' # - CD3+/CD8- single-positive cells
#' # - CD3+ cells with membrane PDL-1 > 5
#' # - All cells regardless of phenotype
#' # - Macrophages, defined as either CD68+ OR CD163+
#' parse_phenotypes("CD3+", "CD3+/CD8+", "CD3+/CD8-",
#'                  "CD3+/PDL-1+"="CD3+/~`Membrane PDL-1 (Opal 520) Mean`>5",
#'                  "Total Cells", Macrophage="CD68+,CD163+")
#' @md
parse_phenotypes = function(...) {
  phenos = list(...)

  # Allow passing a single list
  if (length(phenos)==1 && is.list(phenos[[1]]))
    phenos = phenos[[1]]

  # Check for non-character parameters
  non_char = !purrr::map_lgl(phenos, is.character)
  if (any(non_char))
    stop('parse_phenotypes only works with text descriptions, not ',
         phenos[non_char])

  # Strip leading/trailing spaces preserving any names
  phenos = purrr::map(phenos, stringr::str_trim)

  # If no names were given, phenos will have names(pheno) == NULL
  # If any names were given, missing names will be ''
  # One way or another, get a named list.
  # Make nicer names for formulas by deleting ~ and `
  clean_names = phenotype_names(phenos)
  if (is.null(names(phenos))) names(phenos)=clean_names else {
    no_names = names(phenos) == ''
    names(phenos)[no_names] = clean_names[no_names]
  }

  # This does the basic decoding
  purrr::map(phenos, function(pheno) {
    if (rlang::is_formula(pheno))
      stop("parse_phenotypes does not support formula definitions.")

    # Multiple AND phenotypes become a list
    if (stringr::str_detect(pheno, '/')) {
      # Can't have comma and slash
      if (stringr::str_detect(pheno, ','))
        stop(paste("Phenotype selectors may not contain both '/' and ',':",
                   pheno))
      ## Split the phenotypes and convert formulae
      purrr::map(split_by_slash(pheno), ~{
        if (startsWith(.x, '~')) stats::as.formula(.x, globalenv()) else .x
      })
    }

    # Multiple OR phenotypes become a character vector
    # Formulae are not supported here, they can't be part of a char vector
    else if (stringr::str_detect(pheno, ','))
      purrr::map_chr(split_and_trim(pheno, ','),
                     ~(if (startsWith(.x, '~'))
                         stop('Invalid phenotype definition: ', .x)
                       else .x))

    # Starts with ~, its a formula
    else if (startsWith(pheno, '~'))
      stats::as.formula(pheno, globalenv())
    # Ends with +- and no '/' or ',' is a single phenotype
    else if (stringr::str_detect(pheno, '[+-]$')) pheno

    # Contains Total or All returns NA which signals "Select All"
    else if (stringr::str_detect(pheno, stringr::regex('Total|All',
                                                       ignore_case=TRUE)))
      NA
    else stop(paste("Unrecognized phenotype selector:", pheno))
  }) %>%
    rlang::set_names(names(phenos))
}

#' Find the columns used by phenotype formulae
#'
#' Given a list of parsed phenotypes, e.g. from `parse_phenotypes`,
#' return a vector containing all the names referenced
#' by formulae in `phenos`.
#'
#' @param phenos A list of parsed phenotypes
#' @return A vector of names or NULL if none found
#' @export
phenotype_columns = function(phenos) {
  result = NULL

  # Recursively get names from a "call"
  call_names = function(cl) {
    rlang::call_args(cl) %>%
      purrr::map(~{
        if (is.name(.x)) rlang::as_string(.x) # Get name as a string
        else if (is.call(.x)) call_names(.x) # Recurse
        else NULL
      }) %>%
      unlist()
  }

  purrr::map(phenos, ~{
    if(is.list(.x)) phenotype_columns(.x) # Recurse
    else if (rlang::is_formula(.x)) call_names(rlang::f_rhs(.x))
    else NULL
  }) %>%
    unlist() %>%
    unname()
}

#' Split a single string and trim white space from the results
#' @param str A single string.
#' @param pattern Pattern to split on.
#' @return A character vector of split components.
#' @keywords internal
split_and_trim = function(str, pattern) {
  stopifnot(is.character(str), length(str)==1)
  stringr::str_trim(stringr::str_split(str, pattern)[[1]])
}

#' Split a single string on / and trim white space from the results.
#' Slashes quoted with backquotes are not split.
#' `str` must not contain commas
#' @param str A single string.
#' @return A character vector of split components.
#' @keywords internal
split_by_slash = function(str) {
  stopifnot(is.character(str), length(str)==1)

  # Split into individual characters
  chars = strsplit(str, split='')[[1]]
  if (',' %in% chars)
    stop('Commas not allowed in split_by_slash().')

  # Find character locations preceded by an even number of `
  even_backticks = cumsum(chars=='`') %% 2 == 0

  # Find slashes that are preceded by an even number of `
  slashes = chars == '/' & even_backticks

  # External constraints prohibit comma in str, use it as a marker
  chars[which(slashes)] = ','
  str2 = paste0(chars, collapse='')

  # Now the / we want to split on have been changed to ,
  result = stringr::str_split(str2, ',')[[1]]
  stringr::str_trim(result)
}

#' Make user-friendly names for phenotypes
#' @param phenos A list or vector of phenotype definitions
#' as used by `parse_phenotypes`.
#' @return A vector of name for `phenos`.
#' @keywords internal
phenotype_names = function(phenos) {
  stringr::str_remove_all(phenos, '[~`]')
}

#' Flexibly select rows of a data frame.
#'
#' Select rows of a data frame based on phenotypes or other
#' expressions.
#'
#' `select_rows` implements a flexible mechanism for selecting cells (rows)
#' from a cell segmentation table. Cells may be selected by single or
#' multiple phenotype, by expression level, or combinations of both.
#'
#' See the tutorial
#' [Selecting cells within a cell segmentation table](https://akoyabio.github.io/phenoptr/articles/selecting_cells.html)
#'for extensive documentation and examples.
#'
#' @param csd A data frame
#' @param sel May be a character vector, a one-sided formula, a list
#'   containing such or `NA`. A character vector is interpreted as
#'   the name(s) of one or
#'   more phenotypes and selects any matching phenotype. A formula is
#'   interpreted as an expression on the columns of `csd`.
#'   Multiple list items are joined with AND. `NA` is interpreted
#'   as "select all". It is convenient for lists of selection criteria.
#' @return A logical vector of length `nrow(csd)` which selects rows
#'   according to `sel`.
#' @export
#' @examples
#' csd <- sample_cell_seg_data
#'
#' # Select tumor cells with PDL1 expression > 3
#' selector <- list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3)
#' pdl1_pos_tumor <- csd[select_rows(csd, selector),]
#' range(pdl1_pos_tumor$`Entire Cell PDL1 (Opal 520) Mean`)
#'
#' # Select all T-cells. Note: Use c() to combine phenotypes, not list()
#' selector <- c('CD8+', 'FoxP3+')
#' tcells <- csd[select_rows(csd, selector),]
#' table(tcells$Phenotype)
#' @md
#' @seealso [parse_phenotypes] for a convenient way to create selectors
#' for most common phenotypes.
select_rows <- function(csd, sel) {
  stopifnot(is.data.frame(csd))

  # Evaluate a single phenotype in a per-marker file
  evaluate_per_marker = function(s) {
    if (!stringr::str_detect(s, '[+-]$'))
      stop(paste0(s, ' is not a valid per-marker phenotype name.'))
    column_name = paste('Phenotype', stringr::str_remove(s, '[+-]$'))
    if (!column_name %in% names(csd))
      stop(paste0("No '", column_name, "' column in data."))
    csd[[column_name]] == s
  }

  # Evaluate a single selector
  select_one = function(s) {
    if (length(s)==1 && is.na(s)) {
      # NA means select all
      rep(TRUE, nrow(csd))
    } else if (is.character(s)) {
      # Selector is one or more phenotype names,
      # look for match with phenotype column
      # Any match qualifies
      if ('Phenotype' %in% names(csd)) {
        csd[['Phenotype']] %in% s
      }
      else {
        # Phenotype per-marker has multiple columns
        col_selections = purrr::map(s, evaluate_per_marker)
        purrr::reduce(col_selections, `|`)
      }
    } else {
      # Selector is a function, evaluate it on csd
      col_selections = lazyeval::f_eval(s, csd)

      # Check for valid result
      if (class(col_selections) != 'logical' ||
          length(col_selections) != nrow(csd))
        stop('Invalid expression in select_rows: ~', lazyeval::f_text(s))
      col_selections
    }
  }

  # Everything is selected by default
  result = rep(TRUE, nrow(csd))
  if (!is.list(sel)) sel = list(sel)
  for (s in sel)
    result = result & select_one(s)

  # Don't return NA values, treat them as false
  result %>% tidyr::replace_na(FALSE)
}

# Helper function to normalize lists of selectors into named lists of selectors,
# so we can give names to the selected items.
normalize_selector = function(sel) {
  if (is.null(sel) || length(sel)==0)
    stop("Empty selector")

  stopifnot(is.list(sel))

  if (!is.null(names(sel)))
    return (sel)

  # Name a single selector
  name_item = function(s) {
    if (is.character(s))
      return (paste(s, collapse='|'))
    else if (lazyeval::is_formula(s))
      return (lazyeval::f_text(s))
    else if (is.list(s))
      return (paste(purrr::map_chr(s, name_item), collapse='&'))
    else
      stop('Unknown selector type')
  }

  names(sel) = purrr::map_chr(sel, name_item)
  sel
}

# Make rules that select phenotypes.
#
# Given a list of phenotype names and a (possibly empty) list of rules
# which create some or all of the phenotypes, return a complete list of
# rules.
# @param phenotypes A list or vector of phenotype names. Values may be
# existing phenotypes or compound phenotypes.
# @param existing_rules A named list of phenotype rules.
# @return A named list of rules containing one entry for each member
# of `phenotypes`.
make_phenotype_rules <- function(phenotypes, existing_rules=NULL) {
  if (is.null(existing_rules))
    existing_rules = list()
  else if (!is.list(existing_rules)
           ||(length(existing_rules)>0 && is.null(names(existing_rules))))
      stop("existing_rules must be a named list.")

  existing_names = names(existing_rules)
  extra_names = setdiff(existing_names, phenotypes)
  if (length(extra_names) > 0)
    stop("A rule was given for an unused phenotype: ",
         paste(extra_names, sep=', '))

  # The default rule is just the phenotype name itself.
  missing_names = setdiff(phenotypes, existing_names)
  new_rules = purrr::set_names(as.list(missing_names))

  c(existing_rules, new_rules)
}

#' Validate a phenotype parameter
#' @param phenotypes Nominally, a list or vector of phenotype definitions,
#' or NULL.
#' @param csd A cell seg table
#' @return A named list of phenotype definitions
#' @export
validate_phenotypes = function(phenotypes, csd) {
  if (is.null(phenotypes))
    phenotypes = unique_phenotypes(csd)
  stopifnot(length(phenotypes) > 0)
  if (!rlang::is_named(phenotypes))
    phenotypes = rlang::set_names(phenotypes)
  as.list(phenotypes)
}

#' Find unique phenotypes in a cell seg table
#'
#' For cell seg tables containing a single `Phenotype` column, this
#' returns a vector containing all the non-blank phenotypes in the table.
#' For cell seg tables containing multiple phenotype columns, it returns
#' a vector with just the positive phenotypes.
#' @param csd A cell seg table such as read by `read_cell_seg_table`.
#' @return A named character vector containing the phenotype names.
#' @export
unique_phenotypes = function(csd) {
  if ('Phenotype' %in% names(csd))
    return(purrr::discard(sort(unique(csd$Phenotype)), ~.x==''))

  phenos = names(csd) %>%
    stringr::str_subset('Phenotype ') %>%
    stringr::str_remove('Phenotype ') %>%
    stringr::str_c('+')

  if (length(phenos)==0)
    stop('Cell seg table does not have a phenotype column.')
  rlang::set_names(phenos)
}

#' Mutate a cell seg table to have a Phenotype column with the
#' desired phenotypes.
#'
#' Note: Cells that satisfy multiple phenotype definitions will appear
#' multiple times in the result.
#'
#' @param csd A cell seg table.
#' @param phenotypes A named vector or list of phenotype rules. If NULL, use
#' `unique_phenotypes(csd)`.
#' @return A new cell seg table
#' @keywords internal
make_phenotype_column = function(csd, phenotypes=NULL) {
  if (is.null(phenotypes)) {
    # If phenotypes==NULL and there is already a Phenotype column,
    # this function is a no-op. Just return the input.
    if ('Phenotype' %in% names(csd))
      return(csd)

    phenotypes = unique_phenotypes(csd) %>% purrr::set_names()
  }

  purrr::map(names(phenotypes),
                   ~(csd[select_rows(csd, phenotypes[[.x]]), ] %>%
                       dplyr::mutate(Phenotype=.x))) %>%
    dplyr::bind_rows()
}
