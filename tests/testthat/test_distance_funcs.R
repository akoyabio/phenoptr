# Tests for distance functions
library(testthat)
library(dplyr)

test_that("distance_matrix works", {
  csd = sample_cell_seg_data %>% filter(Phenotype != 'other')
  dst = distance_matrix(csd)
  expect_equal(dim(dst), c(nrow(csd), nrow(csd)))
  s = subset_distance_matrix(dst, csd, 'tumor', 'cytotoxic CD8')
  expect_equal(dim(s)[1], sum(csd$Phenotype=='tumor'))
  expect_equal(dim(s)[2], sum(csd$Phenotype=='cytotoxic CD8'))
})

test_that("find_nearest_distance works", {
  csd = sample_cell_seg_data %>% filter(Phenotype != 'other')
  phenos = sort(unique(csd$Phenotype))
  nearest = find_nearest_distance(csd)
  expect_equal(ncol(nearest), length(phenos))
  expect_equal(nrow(nearest), nrow(csd))
  expect_equal(names(nearest), paste('Distance to', phenos))
})
