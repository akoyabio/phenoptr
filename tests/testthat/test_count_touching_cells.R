library(testthat)

test_that('count_touching_cells works', {
  skip_on_cran()

  cell_seg_path = sample_cell_seg_path()

  pairs = list(c("CD8+", "CK+"),
               c("CD8+", "CD68+"),
               c("CD8+", "FOO+"))

  # Colors for all the phenotypes mentioned in pairs
  colors = list(
    'CD8+' = 'yellow',
    'CK+' = 'cyan',
    'CD68+' = 'magenta',
    'FOO+' = 'pink'
  )

  output_base = tempdir()

  # Error checking
  expect_error(count_touching_cells(cell_seg_path, pairs, write_images=TRUE),
               'required')
  expect_error(count_touching_cells(cell_seg_path, pairs,
                                    colors=list('CD8+' = 'yellow')),
               'names\\(colors\\)')

  counts = count_touching_cells(cell_seg_path, pairs, colors,
    output_base=output_base)

  expect_equal(dim(counts), c(3, 9))
  expect_equal(counts$phenotype1, rep('CD8+', 3))
  expect_equal(counts$phenotype2, c('CK+', 'CD68+', 'FOO+'))
  expect_equal(counts$total1, rep(228, 3))
  expect_equal(counts$total2, c(2257, 417, 0))

  image_names = list.files(output_base, 'touching')
  expect_equal(length(image_names), 2)

  test_base = 'test_data'
  for (image_name in image_names) {
    # Read the result images and compare to a reference
    actual = tiff::readTIFF(file.path(output_base, image_name), as.is=TRUE)
    expected = tiff::readTIFF(file.path(test_base, image_name), as.is=TRUE)
    expect_equal(actual, expected, info=image_name)
  }
})
