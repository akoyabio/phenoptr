#' Flexibly select rows of a data frame.
#'
#' Select rows of a data frame based on phenotypes or other
#' expressions.
#' @param csd A data frame
#' @param sel May be a character vector, a one-sided formula or a list
#'   containing such. A character vector is interpreted as the name(s) of one or
#'   more phenotypes and selects any matching phenotype. A formula is
#'   interpreted as an expression on the columns of \code{csd}. Multiple items are
#'   joined with AND.
#' @param phenotype_column Optional, name of the phenotype column to use for
#'   simple selection
#' @return A logical vector of length \code{nrow(csd)} which selects rows
#'   according to sel.
#' @export
#' @examples
#' csd = sample_cell_seg_data
#'
#' # Select tumor cells with PDL1 expression > 3
#' selector = list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3)
#' pdl1_pos_tumor = csd[select_rows(csd, selector),]
#' range(pdl1_pos_tumor$`Entire Cell PDL1 (Opal 520) Mean`)
#'
#' # Select all T-cells. Note: Use c() to combine phenotypes, not list()
#' selector = c('CD8+', 'FoxP3+')
#' tcells = csd[select_rows(csd, selector),]
#' table(tcells$Phenotype)
select_rows = function(csd, sel, phenotype_column='Phenotype') {
  stopifnot(is.data.frame(csd))

  # Evaluate a single selector
  select_one = function(s) {
    if (is.character(s)) {
      # Selector is one or more phenotype names,
      # look for match with phenotype column
      stopifnot(phenotype_column %in% names(csd))
      csd[[phenotype_column]] %in% s
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
