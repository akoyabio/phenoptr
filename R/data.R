# Documentation for data

#' Sample cell segmentation data.
#'
#' A dataset containing inForm cell segmentation data for a single field.
#'
#' `sample_cell_seg_data` contains the `data_table` returned from
#' ```{r eval=FALSE}
#' read_cell_seg_data(sample_cell_seg_path())
#' ```
#'
#' This table shows the stains used, the epitopes they were bound to,
#' and the colors used to show the stains in the composite and phenotype views.
#'
#' \tabular{llll}{
#'   Stain   \tab Epitope \tab Composite \tab Phenotype \cr
#'   ----------- \tab ----------- \tab ------------- \tab ------------- \cr
#'   DAPI    \tab Nucleus \tab Blue    \tab N/A    \cr
#'   Opal520 \tab PDL1    \tab Red     \tab N/A    \cr
#'   Opal540 \tab CD8     \tab Yellow  \tab Yellow \cr
#'   Opal570 \tab FoxP3   \tab Orange  \tab Orange \cr
#'   Opal620 \tab CD68    \tab Magenta \tab Magenta\cr
#'   Opal650 \tab PD1     \tab Green   \tab N/A    \cr
#'   Opal690 \tab CK      \tab Cyan    \tab Cyan
#' }
#' The sample cell data includes phenotypes `CD8+` (cytotoxic T cell),
#' `CD68+` (macrophage),
#' `FoxP3+` (regulatory T cell), `CK+` (tumor), and `other`.
#'
#' [sample_cell_seg_path] and [sample_cell_seg_folder] are
#' convenience functions which return the path to the sample cell seg file
#' included in this package, or the folder containing
#' sample cell seg data, respectively.
#'
#' To see all the files included for this field, use
#' `list.files(sample_cell_seg_folder())`.
#'
#' For more sample data from this and related samples, see the
#' [phenoptrExamples](https://perkinelmer.github.io/phenoptrExamples) package.
#' @format A data frame with 6072 rows and 199 variables
#' @md
"sample_cell_seg_data"

#' @rdname sample_cell_seg_data
#' @export
sample_cell_seg_path = function()
  system.file("extdata", "sample",
                    "Set4_1-6plex_[16142,55840]_cell_seg_data.txt",
                    package = "phenoptr")

#' @rdname sample_cell_seg_data
#' @export
sample_cell_seg_folder = function()
  system.file("extdata", "sample", package = "phenoptr")
