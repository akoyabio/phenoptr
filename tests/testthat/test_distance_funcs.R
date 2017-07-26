# Tests for distance functions
# Mostly smoke tests to make sure the code runs...
library(testthat)

test_that("distance_matrix works", {
  csd = sample_cell_seg_data %>% dplyr::filter(Phenotype != 'other')
  dst = distance_matrix(csd)
  expect_equal(dim(dst), c(nrow(csd), nrow(csd)))

  s = subset_distance_matrix(csd, dst, 'CK+', 'CD8+')
  expect_equal(dim(s)[1], sum(csd$Phenotype=='CK+'))
  expect_equal(dim(s)[2], sum(csd$Phenotype=='CD8+'))
  expect_error(subset_distance_matrix(dst, csd, 'CK+', 'CD8+'), 'wrong order')
})

test_that("find_nearest_distance works", {
  csd = sample_cell_seg_data %>% dplyr::filter(Phenotype != 'other')
  phenos = sort(unique(csd$Phenotype))
  nearest = find_nearest_distance(csd)
  expect_equal(ncol(nearest), length(phenos))
  expect_equal(nrow(nearest), nrow(csd))
  expect_equal(names(nearest), paste('Distance to', phenos))

  phenos = c(phenos, 'other') # A mising phenotype
  nearest = find_nearest_distance(csd, phenos)
  expect_equal(ncol(nearest), length(phenos))
  expect_equal(nearest$`Distance to other`, rep(NA_real_, nrow(csd)))

  # Force an error
  csd$`Sample Name`[1] = 'foo.im3'
  expect_error(find_nearest_distance(csd), 'multiple')
})

test_that('compute_all_nearest_distance works', {
  skip_on_cran()
  path = sample_cell_seg_path()
  out_path = tempfile()
  compute_all_nearest_distance(path, out_path)
  expect_true(file.exists(out_path))
  file.remove(out_path)
})
