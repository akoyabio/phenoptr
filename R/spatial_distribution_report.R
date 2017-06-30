#' Create a spatial distribution report
#'
#' Creates an HTML report showing the location and nearest-neighbor relations
#' of cells in a single field and two chosen phenotypes.
#'
#' A cell seg data file for the field is required. If an image of the field
#' is available, it will be used as a background.
#'
#' The report is saved to the same directory as the input file.
#'
#' @param path Path to a cell seg table.
#' @param phenotypes A list of phenotype specifications in the form used by
#' [select_rows()].
#' @param colors A character vector of color specifications, of the same
#' length as `phenotypes`.
#' @param plot_gx_max Optional, maximum value of radius to show on G
#' function plots.
#' @param micron_line Optional, X location (radius in microns)
#'  of crosshair on G function plots.
#' @export
#' @md
spatial_distribution_report = function(cell_seg_path, phenotypes, colors,
  plot_gx_max = 100, micron_line = 25) {
  stopifnot(grepl('_cell_seg_data.txt', cell_seg_path))
  stopifnot(file.exists(cell_seg_path), length(phenotypes)==length(colors))

  out_path = sub('_cell_seg_data.txt',
                 '_spatial_distribution.html', cell_seg_path)
  rmd_path = system.file("rmd", "SpatialDistributionReport.Rmd")
  rmarkdown::render(rmd_path, 'html_document', out_path)
}
