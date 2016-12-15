#' Flexibly select rows of a data frame based on phenotypes or other
#' expressions.
#' @param d A data frame
#' @param sel May be a single string, a one-sided formula or a list containing
#'   such. A string is interpreted as the name of a phenotype and selects that
#'   phenotype. A formula is interpreted as an expression on the columns of
#'   \code{d}.
#'   Multiple items are joined with AND.
#' @param phenotype_column Optional, name of the phenotype column to use for
#'   simple selection
#' @return A logical vector of length \code{nrow(d)} which selects rows
#'   according to sel.
#' @export
#' @examples
#' path = system.file("extdata",
#'   "TMA/Core[1,5,6,1]_[21302,15107]_cell_seg_data.txt",
#'   package = "informr")
#' d = read_cell_seg_data(path)
#'
#' ## Select tumor cells with PDL1 expression > 4
#' selector = list('tumor', ~`Entire Cell PDL1 (Opal 620) Mean`>4)
#' pdl1_pos_tumor = d[select_rows(d, selector),]
#' table(pdl1_pos_tumor$Phenotype)
select_rows = function(d, sel, phenotype_column='Phenotype') {
  stopifnot(is.data.frame(d))

  # Evaluate a single selector
  select_one = function(s) {
    if (is.character(s)) {
      # Selector is a phenotype name, look for match with phenotype column
      stopifnot(length(s) == 1, phenotype_column %in% names(d))
      d[[phenotype_column]] == s
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
