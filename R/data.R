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
#' ```
#'   Stain    Epitope  Composite  Phenotype
#'   -------- -------  ---------  ---------
#'   DAPI     Nucleus  Blue       N/A
#'   Opal520  PDL1     Red        N/A
#'   Opal540  CD8      Yellow     Yellow
#'   Opal570  FoxP3    Orange     Orange
#'   Opal620  CD68     Magenta    Magenta
#'   Opal650  PD1      Green      N/A
#'   Opal690  CK       Cyan       Cyan
#' ```
#'
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
#' [phenoptrExamples](https://akoyabio.github.io/phenoptrExamples) package.
#' @format A data frame with 6072 rows and 199 variables
#' @examples
#' dim(sample_cell_seg_data)
#' table(sample_cell_seg_data$Phenotype)
#' @md
"sample_cell_seg_data"

#' @rdname sample_cell_seg_data
#' @export
sample_cell_seg_path <- function()
  system.file("extdata", "sample",
                    "Set4_1-6plex_[16142,55840]_cell_seg_data.txt",
                    package = "phenoptr")

#' @rdname sample_cell_seg_data
#' @export
sample_cell_seg_folder <- function()
  system.file("extdata", "sample", package = "phenoptr")
