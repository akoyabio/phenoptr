# Documentation for data

#' Sample cell segmentation data.
#'
#' A dataset containing inForm cell segmentation data for a single field.
#' This is the data from
#' ```
#' system.file("extdata", "sample",
#'                   "Set4_1-6plex_[16142,55840]_cell_seg_data.txt",
#'                   package = "phenoptr")
#' ```
#'
#' This table shows the Opal stains used, the epitope they were bound to,
#' and the colors used to show the stain in the composite and phenotype views.
#'
#' Stain   | Epitope | Composite | Phenotype
#' --------|---------|-----------|----------
#' DAPI    | Nucleus | Blue      | N/A
#' Opal520 | PDL1    | Red       | N/A
#' Opal540 | CD8     | Yellow    | Yellow
#' Opal570 | FoxP3   | Orange    | Orange
#' Opal620 | CD68    | Magenta   | Magenta
#' Opal650 | PD1     | Green     | N/A
#' Opal690 | CK      | Cyan      | Cyan
#'
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
