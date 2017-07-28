library(testthat)

test_that('count_touching_cells works', {
  cell_seg_path =
    file.path('test_data',
              'FIHC4__0929309_HP_IM3_2_cell_seg_data.txt')

  pairs = list(c("Helper T", "B"),
               c("Helper T", "Cytotoxic T"),
               c("Helper T", "FOO+"))

  # Colors for all the phenotypes mentioned in pairs
  colors = list(
    'Cytotoxic T' = 'yellow',
    'Helper T' = 'green',
    'B' = 'red',
    'FOO+' = 'pink'
  )

  output_base = tempdir()

  # Error checking
  expect_error(count_touching_cells(cell_seg_path, pairs, write_images=TRUE),
               'required')
  expect_error(count_touching_cells(cell_seg_path, pairs,
                                    colors=list('CD8+' = 'yellow')),
               'names\\(colors\\)')

  # This one should work
  expect_warning(counts <- count_touching_cells(cell_seg_path, pairs, colors,
                                                categories='Tumor',
                                                output_base=output_base),
                 "No image for .* Helper T touching FOO+")

  expect_equal(dim(counts), c(3, 9))
  expect_equal(counts$phenotype1, rep('Helper T', 3))
  expect_equal(counts$phenotype2, c('B', 'Cytotoxic T', 'FOO+'))
  expect_equal(counts$total1, rep(15, 3))
  expect_equal(counts$total2, c(249, 6, 0))
  expect_equal(counts$p1_touch_p2, c(14, 3, 0))
  expect_equal(counts$p2_touch_p1, c(23, 4, 0))
  expect_equal(counts$touch_pairs, c(31, 5, 0))

  image_names = list.files(output_base, '_touch_')
  expect_equal(length(image_names), 2)

  test_base = 'test_results'
  for (image_name in image_names) {
    # Read the result images and compare to a reference
    actual = tiff::readTIFF(file.path(output_base, image_name), as.is=TRUE)
    expected = tiff::readTIFF(file.path(test_base, image_name), as.is=TRUE)
    expect_equal(actual, expected, info=image_name)
  }
})

test_that('replace_invalid_path_characters works', {
  bad = "<>:\"/\\|?* ";

  expect_equal(replace_invalid_path_characters(bad, "_"), "__________")

  bad = paste0("a", bad, "b")
  expect_equal(replace_invalid_path_characters(bad, "_"), "a__________b")

  bad = "a<b>c:d\"e/f\\g|h?i*j k"
  expect_equal(replace_invalid_path_characters(bad, "_"),
    "a_b_c_d_e_f_g_h_i_j_k")

  # Replacement can be multiple chars
  expect_equal(replace_invalid_path_characters(bad, "_="),
    "a_=b_=c_=d_=e_=f_=g_=h_=i_=j_=k")

  # or nothing
  expect_equal(replace_invalid_path_characters(bad, ""), "abcdefghijk")

  # Replacement can't contain invalid characters
  expect_error(replace_invalid_path_characters(bad, "*"), 'invalid')
})
