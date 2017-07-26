# Tests for read_components and read_maps
# Uses data from phenoptrExamples so don't run on CRAN
library(testthat)

test_that("read_components works", {
  skip_on_cran()
  skip_if_not_installed('phenoptrExamples')
  path = system.file('extdata', 'samples',
                     'Set4_1-6plex_[16142,55840]_component_data.tif',
                     package='phenoptrExamples')

  images = read_components(path)

  expected_names = c("PDL1 (Opal 520)", "CD8 (Opal 540)", "FoxP3 (Opal 570)",
                     "CD68 (Opal 620)", "PD1 (Opal 650)", "CK (Opal 690)",
                     "DAPI", "Autofluorescence")
  expect_equal(names(images), expected_names)
  expect_equivalent(purrr::map(images, dim), rep(list(c(1400, 1868)), 8))
})

test_that('read_maps works', {
  skip_on_cran()
  skip_if_not_installed('phenoptrExamples')
  path = system.file('extdata', 'samples',
                     'Set4_1-6plex_[16142,55840]_binary_seg_maps.tif',
                     package='phenoptrExamples')

  maps = read_maps(path)

  expected_names = c("Nucleus", "Cytoplasm", "Membrane", "Tissue")
  expect_equal(names(maps), expected_names)
  expect_equivalent(purrr::map(maps, dim), rep(list(c(1400, 1868)), 4))
})
