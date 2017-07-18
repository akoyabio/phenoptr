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
#' # Select tumor cells with PDL1 expression > 4
#' selector = list('tumor', ~`Entire Cell PDL1 (Opal 620) Mean`>4)
#' pdl1_pos_tumor = csd[select_rows(csd, selector),]
#' range(pdl1_pos_tumor$`Entire Cell PDL1 (Opal 620) Mean`)
#'
#' # Select all T-cells. Note: Use c() to combine phenotypes, not list()
#' selector = c('cytotoxic CD8', 'helper CD4', 'T reg Foxp3')
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
