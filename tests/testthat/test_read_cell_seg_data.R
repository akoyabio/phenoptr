# Tests for read_cell_seg_data
library(testthat)

path = system.file("extdata", "sample",
   "Set4_1-6plex_[16142,55840]_cell_seg_data.txt",
   package = "phenoptr")

summary_path = system.file("extdata", "sample",
   "Set4_1-6plex_[16142,55840]_cell_seg_data_summary.txt",
   package = "phenoptr")

expect_contains = function(container, items) {
  for (item in items)
    expect_true(item %in% container, info=item)
}

expect_does_not_contain = function(container, items) {
  for (item in items)
    expect_false(item %in% container, info=item)
}

test_that('list_cell_seg_files works', {
  files = list_cell_seg_files(sample_cell_seg_folder())
  expect_equal(length(files), 1)
})

test_that("Blank file name is an error", {
  expect_error(read_cell_seg_data(''), 'missing')
})

test_that("read_cell_seg_data works", {
  d = read_cell_seg_data(path)
  # Check size
  expect_equal(nrow(d), 6072)
  expect_equal(ncol(d), 199)

  # Check column names
  expect_contains(names(d),
                  c('Sample Name', 'Cell Y Position',
                    'Nucleus Area (sq microns)',
                    'Nucleus DAPI Mean', 'Phenotype'))
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Tissue Category Area (sq microns)',
          'Nucleus DAPI Mean (Normalized Counts, Total Weighting)', # shortened
          'Nucleus Area (pixels)', # Converted
          "Cell Density (per megapixel)", # Not in regular cell table
          "Cell Density (per sq mm)"
          ))

  # Check values
  # Calling as.numeric converts from a tibble to an atomic value
  expect_equal(as.numeric(d[1, 'Cell X Position']), 515/2)
  expect_equal(as.numeric(d[1, 'Nucleus Area (percent)']), 0.0001)

  # Summary file
  d = read_cell_seg_data(summary_path)
  expect_contains(names(d),
                  c('Sample Name',
                  'Phenotype',
                  'Tissue Category',
                  'Tissue Category Area (sq microns)',
                  'Nucleus Area (sq microns)',
                  'Nucleus DAPI Mean',
                  "Cell Density (per sq mm)"
                  ))
  expect_does_not_contain(names(d),
        c('Cell Y Position',
          'Tissue Category Area (pixels)', # Converted
          'Nucleus DAPI Mean (Normalized Counts, Total Weighting)', # shortened
          'Nucleus Area (pixels)', # Converted
          "Cell Density (per megapixel)" # Converted
          ))

  expect_equal(as.character(d[1, 'Cell ID']), 'all')
  expect_equal(as.numeric(d[1, 'Tissue Category Area (sq microns)']), 1208520/4)
  expect_equal(as.numeric(d[1, 'Cell Density (per sq mm)']), 83.57*4)
  expect_equal(as.numeric(d[1, 'Nucleus Area (percent)']), 0.0173)
})

test_that('Skipping pixels_per_micron works', {
  d = read_cell_seg_data(path, pixels_per_micron=NA)
  # Check size
  expect_equal(nrow(d), 6072)
  expect_equal(ncol(d), 199)

  # Check column names
  expect_contains(names(d), 'Nucleus Area (pixels)')
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Nucleus Area (sq microns)' # Not converted
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 515)
})

test_that('Setting pixels_per_micron works', {
  d = read_cell_seg_data(path, pixels_per_micron=4)
  # Check size
  expect_equal(nrow(d), 6072)
  expect_equal(ncol(d), 199)

  # Check column names
  expect_contains(names(d), 'Nucleus Area (sq microns)')
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Nucleus Area (pixels)' # Converted
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 515/4)
})

test_that('remove_units=FALSE works', {
  d = read_cell_seg_data(path, remove_units=FALSE)
  # Check size
  expect_equal(nrow(d), 6072)
  expect_equal(ncol(d), 199)

    # Check column names
  expect_contains(names(d), c('Nucleus Area (sq microns)',
                  'Nucleus DAPI Mean (Normalized Counts, Total Weighting)'))
  expect_does_not_contain(names(d),
        c('Nucleus DAPI Mean' # Not shortened
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 515/2)
})

test_that('Making tags works', {
  names = c("Set12_20-6plex_[14146,53503].im3",
    "Set4_1-6plex_[15206,60541].im3",
    "Set8_11-6plex_[17130,56449].im3")
  expected = c("12_20-6plex_[14146,53503].im3", "4_1-6plex_[15206,60541].im3",
    "8_11-6plex_[17130,56449].im3")
  expect_equal(remove_common_prefix(names), expected)

  names = expected
  expected = c("12_20-6plex_[14146,53503]", "4_1-6plex_[15206,60541]",
    "8_11-6plex_[17130,56449]")
  expect_equal(remove_extensions(names), expected)

  names = rep("Set12_20-6plex_[14146,53503].im3", 3)
  expect_equal(remove_common_prefix(names), names)
})
