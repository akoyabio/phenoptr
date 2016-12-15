# Tests for distance functions
library(testthat)
library(dplyr)

test_that("distance_matrix works", {
  d = sample_cell_seg_data
  m = distance_matrix(d)
  expect_equal(dim(m), c(nrow(d), nrow(d)))
})
