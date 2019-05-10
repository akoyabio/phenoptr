# Tests for read_cell_seg_data
context('read_cell_seg_data')

library(testthat)

path = test_path('test_data',
              'FIHC4__0929309_HP_IM3_2_cell_seg_data.txt')

summary_path =
  test_path('test_data',
          'FIHC4__0929309_HP_IM3_2_cell_seg_data_summary.txt')

expect_contains = function(container, items) {
  for (item in items)
    expect_true(item %in% container, info=item)
}

expect_does_not_contain = function(container, items) {
  for (item in items)
    expect_false(item %in% container, info=item)
}

test_that('list_cell_seg_files works', {
  files = list_cell_seg_files(test_path('test_data'))
  expect_equal(length(files), 5)
})

test_that('sample_cell_seg_path is correct', {
  # We don't use this anywhere else in the tests, spot check it here
  p = sample_cell_seg_path()
  expect_true(!is.na(p))
  expect_true(file.exists(p))
  expect_true(endsWith(p, '_cell_seg_data.txt'))
})

test_that("Blank file name is an error", {
  expect_error(read_cell_seg_data(''), 'missing')
})

test_that("read_cell_seg_data works", {
  d = read_cell_seg_data(path)
  # Check size
  expect_equal(dim(d), c(270, 40))

  # Check column names
  expect_contains(names(d),
                  c('Sample Name', 'Cell Y Position',
                    'Nucleus Area (square microns)',
                    'Nucleus DAPI Mean', 'Phenotype'))
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Tissue Category Area (square microns)',
          'Nucleus DAPI Mean (Normalized Counts, Total Weighting)', # shortened
          'Nucleus Area (pixels)', # Converted
          "Cell Density (per megapixel)", # Not in regular cell table
          "Cell Density (per square mm)"
          ))

  # Check values
  # Calling as.numeric converts from a tibble to an atomic value
  expect_equal(as.numeric(d[1, 'Cell X Position']), 268/2)
  expect_equal(as.numeric(d[1, 'Nucleus Area (percent)']), 0.0017)

  # Summary file
  d = read_cell_seg_data(summary_path)
  expect_contains(names(d),
                  c('Sample Name',
                  'Phenotype',
                  'Tissue Category',
                  'Tissue Category Area (square microns)',
                  'Nucleus Area (square microns)',
                  'Nucleus DAPI Mean',
                  "Cell Density (per square mm)"
                  ))
  expect_does_not_contain(names(d),
        c('Cell Y Position',
          'Tissue Category Area (pixels)', # Converted
          'Nucleus DAPI Mean (Normalized Counts, Total Weighting)', # shortened
          'Nucleus Area (pixels)', # Converted
          "Cell Density (per megapixel)" # Converted
          ))

  expect_equal(as.character(d[1, 'Cell ID']), 'all')
  expect_equal(as.numeric(d[1, 'Tissue Category Area (square microns)']),
               89904/4)
  expect_equal(as.numeric(d[1, 'Cell Density (per square mm)']), 2770*4)
  expect_equal(as.numeric(d[1, 'Nucleus Area (percent)']), 0.5961)
})

test_that('Comma as decimal separator works', {
  comma_path = test_path('test_data',
                   'FIHC4__0929309_HP_IM3_2_comma_cell_seg_data.txt')
  comma_summary_path =
    test_path('test_data',
              'FIHC4__0929309_HP_IM3_2_comma_cell_seg_data_summary.txt')

  # These should match except for the Sample Name which has commas
  d = read_cell_seg_data(path)
  d_comma = expect_message(read_cell_seg_data(comma_path), 'comma')
  expect_equal(d %>% dplyr::select(-`Sample Name`),
               d_comma %>% dplyr::select(-`Sample Name`))

  d = read_cell_seg_data(summary_path)
  d_comma = expect_message(read_cell_seg_data(comma_summary_path), 'comma')
  expect_equal(d %>% dplyr::select(-`Sample Name`),
               d_comma %>% dplyr::select(-`Sample Name`))

})

test_that('Skipping pixels_per_micron works', {
  d = read_cell_seg_data(path, pixels_per_micron=NA)
  # Check size
  expect_equal(dim(d), c(270, 40))

  # Check column names
  expect_contains(names(d), 'Nucleus Area (pixels)')
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Nucleus Area (square microns)' # Not converted
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 268)
})

test_that('Setting pixels_per_micron works', {
  d = read_cell_seg_data(path, pixels_per_micron=4)
  # Check size
  expect_equal(dim(d), c(270, 40))

  # Check column names
  expect_contains(names(d), 'Nucleus Area (square microns)')
  expect_does_not_contain(names(d),
        c('Tissue Category Area (pixels)', # Empty column
          'Nucleus Area (pixels)' # Converted
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 268/4)
})

test_that('auto pixels_per_micron works', {
  skip_if_not_installed('phenoptrExamples')
  path = system.file('extdata', 'samples',
                     'Set4_1-6plex_[16142,55840]_cell_seg_data.txt',
                     package='phenoptrExamples')
  d = read_cell_seg_data(path, pixels_per_micron='auto')

  # Check for X, Y in known range
  component_path = sub('_cell_seg_data.txt', '_component_data.tif', path)
  info = get_field_info(component_path)
  min_x = info$location[1]
  max_x = info$location[1] + info$image_size[1]*info$microns_per_pixel
  expect_equal(sum(d$`Cell X Position`>=min_x), nrow(d))
  expect_equal(sum(d$`Cell X Position`<=max_x), nrow(d))

  min_y = info$location[2]
  max_y = info$location[2] + info$image_size[2]*info$microns_per_pixel
  expect_equal(sum(d$`Cell Y Position`>=min_y), nrow(d))
  expect_equal(sum(d$`Cell Y Position`<=max_y), nrow(d))
})

test_that('remove_units=FALSE works', {
  d = read_cell_seg_data(path, remove_units=FALSE)
  # Check size
  expect_equal(dim(d), c(270, 40))

    # Check column names
  expect_contains(names(d), c('Nucleus Area (square microns)',
                  'Nucleus DAPI Mean (Normalized Counts, Total Weighting)'))
  expect_does_not_contain(names(d),
        c('Nucleus DAPI Mean' # Not shortened
          ))

  # Check values
  expect_equal(as.numeric(d[1, 'Cell X Position']), 268/2)
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
