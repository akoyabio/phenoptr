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
#' [Selecting cells within a cell segmentation table](https://perkinelmer.github.io/phenoptr/articles/selecting_cells.html)
#'for extensive documentation and examples.
#'
#' @param csd A data frame
#' @param sel May be a character vector, a one-sided formula, a list
#'   containing such or `NULL`. A character vector is interpreted as
#'   the name(s) of one or
#'   more phenotypes and selects any matching phenotype. A formula is
#'   interpreted as an expression on the columns of `csd`.
#'   Multiple list items are joined with AND. `NA` is interpreted
#'   as "select all". It is convienent for lists of selection criteria.
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
    if (is.na(s)) {
      # NULL just means selert all
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
      lazyeval::f_eval(s, csd)
    }
  }

  # Everything is selected by default
  result = rep(TRUE, nrow(csd))
  if (!is.list(sel)) sel = list(sel)
  for (s in sel)
    result = result & select_one(s)
  result
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
make_phenotype_rules <- function (phenotypes, existing_rules=NULL) {
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
