#' Helpers for working with inForm data
#'
#' informr contains functions that make it easier to read and
#'  analyze data tables and images created by PerkinElmer's inForm software.
#'
#' informr is part of the PerkinElmer Phenoptics™ family of
#' Quantitative Pathology Research Solutions. For more information
#' visit the Phenoptics™
#' [home page](http://www.perkinelmer.com/cancer-immunology/index.html).
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
#' [find_nearest_distance] finds
#' the distance from each cell to the nearest cell of each phenotype.
#'
#' [count_within] counts
#' the number of cells within a fixed radius of other cells.
#'
#' [spatial_distribution_report] creates a report visualizing nearest
#' neighbor relationships between cells of two phenotypes.
#'
#' [distance_matrix] and [subset_distance_matrix]
#' create and subset cell distance matrices from cell segmentation data.
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
