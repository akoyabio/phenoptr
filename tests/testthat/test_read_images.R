# Tests for read_components and read_maps
library(testthat)

test_that("read_components works", {
  # component_data files are big, even for our toy test data.
  # Use phenoptrExamples data here to keep the package size down a bit,
  # though it makes the test take longer...
  skip_if_not_installed('phenoptrExamples')
  path = system.file('extdata', 'samples',
                     'Set4_1-6plex_[16142,55840]_component_data.tif',
                     package='phenoptrExamples')

  images = read_components(path)

  expected_names = c("PDL1 (Opal 520)", "CD8 (Opal 540)", "FoxP3 (Opal 570)",
                     "CD68 (Opal 620)", "PD1 (Opal 650)", "CK (Opal 690)",
                     "DAPI", "Autofluorescence")
  expect_equal(names(images), expected_names)

  # Check values only (ignore names)
  expect_equal(purrr::map(images, dim), rep(list(c(1400, 1868)), 8),
               ignore_attr=TRUE)
  expect_equal(purrr::map_dbl(images, min), rep(0, 8),
               ignore_attr=TRUE)
  expect_equal(purrr::map_int(images, ~as.integer(max(.))),
                    c(97, 30, 15, 6, 55, 41, 35, 31),
               ignore_attr=TRUE)
})

test_that('read_maps works', {
  path =
    file.path('test_data',
              'FIHC4__0929309_HP_IM3_2_binary_seg_maps.tif')
  maps = read_maps(path)

  expected_names = c("Nucleus", "Cytoplasm", "Membrane", "Tissue")
  expect_equal(names(maps), expected_names)
  expect_equal(purrr::map(maps, dim), rep(list(c(260, 348)), 4),
               ignore_attr=TRUE)
})

test_that('get_field_info works', {
  skip_if_not_installed('phenoptrExamples')
  path = system.file('extdata', 'samples',
                     'Set4_1-6plex_[16142,55840]_component_data.tif',
                     package='phenoptrExamples')
  info = get_field_info(path)
  expect_equal(info$image_size, c(1868, 1400))
  expect_equal(info$microns_per_pixel, 0.4981226, tolerance=0.0000001)
  expect_equal(info$field_size, c(930.4929, 697.3716), tolerance=0.0001)
  expect_equal(info$location, c(15676.75, 55491.31), tolerance=0.01)
})
