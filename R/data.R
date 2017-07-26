# Documentation for data

#' Sample cell segmentation data.
#'
#' A dataset containing inForm cell segmentation data for a single field.
#' This is the data from
#' ```{r eval=FALSE}
#' system.file("extdata", "sample",
#'                   "Set4_1-6plex_[16142,55840]_cell_seg_data.txt",
#'                   package = "phenoptr")
#' ```
#'
#' This table shows the Opal stains used, the epitope they were bound to,
#' and the colors used to show the stain in the composite and phenotype views.
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
#' The sample cell data includes phenotypes CD8+ (cytotoxic T cell),
#' CD68+ (macrophage),
#' FoxP3+ (regulatory T cell), CK+ (tumor), and other.
#'
#' To see all the files included for this field, use
#' `list.files(system.file("extdata", "sample", package = "phenoptr"))`.
#'
#' For more sample data from this and related samples, see the
#' [phenoptrExamples](https://perkinelmer.github.io/phenoptrExamples) package.
#' @format A data frame with 6072 rows and 199 variables
#' @md
"sample_cell_seg_data"
