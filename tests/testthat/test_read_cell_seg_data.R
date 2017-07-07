# Tests for read_cell_seg_data
library(testthat)

path = system.file("extdata", "TMA",
   "Core[1,5,6,1]_[21302,15107]_cell_seg_data.txt",
   package = "informr")

expect_contains = function(container, items) {
  for (item in items)
    expect_true(item %in% container, info=item)
}

expect_does_not_contain = function(container, items) {
  for (item in items)
    expect_false(item %in% container, info=item)
}


test_that("read_cell_seg_data works", {
  d = read_cell_seg_data(path)
  # Check size
  expect_equal(nrow(d), 5726)
  expect_equal(ncol(d), 200)

  # Check column names
  expect_contains(names(d), c('Sample Name', 'Cell Y Position',
                              'Nucleus Area (sq microns)',
                              'Nucleus DAPI Mean', 'Phenotype'))
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Nucleus DAPI Mean (Normalized Counts, Total Weighting)', # shortened
          'Nucleus Area (pixels)' # Converted
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 1211/2)
  expect_equal(as.numeric(d[1, 'Nucleus Area (percent)']), 0.0001)
})

test_that('Skipping pixels_per_micron works', {
  d = read_cell_seg_data(path, pixels_per_micron=NA)
  # Check size
  expect_equal(nrow(d), 5726)
  expect_equal(ncol(d), 200)

    # Check column names
  expect_contains(names(d), 'Nucleus Area (pixels)')
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Nucleus Area (sq microns)' # Converted
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 1211)
})

test_that('Setting pixels_per_micron works', {
  d = read_cell_seg_data(path, pixels_per_micron=4)
  # Check size
  expect_equal(nrow(d), 5726)
  expect_equal(ncol(d), 200)

    # Check column names
  expect_contains(names(d), 'Nucleus Area (sq microns)')
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Nucleus Area (pixels)' # Converted
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 1211/4)
})

test_that('remove_units=FALSE works', {
  d = read_cell_seg_data(path, remove_units=FALSE)
  # Check size
  expect_equal(nrow(d), 5726)
  expect_equal(ncol(d), 200)

    # Check column names
  expect_contains(names(d), c('Nucleus Area (sq microns)',
                  'Nucleus DAPI Mean (Normalized Counts, Total Weighting)'))
  expect_does_not_contain(names(d),
        c('Nucleus DAPI Mean' # Not shortened
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 1211/2)
})
