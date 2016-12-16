#' Flexibly select rows of a data frame based on phenotypes or other
#' expressions.
#' @param d A data frame
#' @param sel May be a character vector, a one-sided formula or a list
#'   containing such. A character vector is interpreted as the name(s) of one or
#'   more phenotypes and selects any matching phenotype. A formula is
#'   interpreted as an expression on the columns of \code{d}. Multiple items are
#'   joined with AND.
#' @param phenotype_column Optional, name of the phenotype column to use for
#'   simple selection
#' @return A logical vector of length \code{nrow(d)} which selects rows
#'   according to sel.
#' @export
#' @examples
#' d = sample_cell_seg_data
#'
#' ## Select tumor cells with PDL1 expression > 4
#' selector = list('tumor', ~`Entire Cell PDL1 (Opal 620) Mean`>4)
#' pdl1_pos_tumor = d[select_rows(d, selector),]
#' hist(pdl1_pos_tumor$`Entire Cell PDL1 (Opal 620) Mean`)
#'
#' ## Select all T-cells
#' selector = c('cytotoxic CD8', 'helper CD4', 'T reg Foxp3')
#' tcells = d[select_rows(d, selector),]
#' table(tcells$Phenotype)
select_rows = function(d, sel, phenotype_column='Phenotype') {
  stopifnot(is.data.frame(d))

  # Evaluate a single selector
  select_one = function(s) {
    if (is.character(s)) {
      # Selector is one or more phenotype names, look for match with phenotype column
      stopifnot(phenotype_column %in% names(d))
      d[[phenotype_column]] %in% s
    } else {
      # Selector is a function, evaluate it on d
      lazyeval::f_eval(s, d)
    }
  }

  # Everything is selected by default
  result = rep(TRUE, nrow(d))
  if (!is.list(sel)) sel = list(sel)
  for (s in sel)
    result = result & select_one(s)
  result
}
