#' Helpers for working with inForm data.
#'
#' informr contains functions that make it easier to read and
#'  analyze data tables created by PerkinElmer's inForm software.
#'
#' \subsection{File readers}{
#' [read_cell_seg_data] reads inForm cell segmentation data tables
#' and makes them easier to use in R.
#'
#' [read_components] reads inForm component data and
#' [read_maps] reads inForm segmentation maps.
#' }
#'
#' \subsection{Distance functions}{
#' [distance_matrix] and [subset_distance_matrix]
#' create and subset cell distance matrices from cell segmentation data.
#'
#' [find_nearest_distance] uses a distance matrix to find
#' the distance from each cell to the nearest cell of each phenotype.
#'
#' [count_within] uses a distance matrix to find
#' the number of cells within a radius of other cells.
#'}
#'
#' \subsection{Helpers}{
#' [select_rows] helps select rows in cell segmentation data
#' corresponding to specific phenotypes and expression levels.
#' }
#'
#' @docType package
#' @name informr
#' @aliases informr
#' @md
NULL
