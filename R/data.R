# Documentation for data

#' Sample cell segmentation data.
#'
#' A dataset containing inForm cell segmentation data for a single TMA core.
#' This is the data from
#' \code{system.file("extdata", "TMA",
#'                   "Core[1,5,6,1]_[21302,15107]_cell_seg_data.txt",
#'                   package = "informr")}
#'
#' This table includes columns for components:
#' \itemize{
#'   \item DAPI
#'   \item CD4 (Opal 520)
#'   \item CK (Opal 540)
#'   \item CD8 (Opal 570)
#'   \item PDL1 (Opal 620)
#'   \item CD68 (Opal 650)
#'   \item Foxp3 (Opal 690)
#' }
#'
#' and phenotypes cytotoxic CD8, helper CD4, macrophage CD68, other,
#' T reg Foxp3, and tumor.
#'
#' To see all the included files for this core, use
#' \code{list.files(system.file("extdata", "TMA", package = "informr"))}
#' @format A data frame with 5726 rows and 200 variables
"sample_cell_seg_data"
