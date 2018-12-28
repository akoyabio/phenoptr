# Create a spatial distribution report
context('spatial_distribution_report')
library(testthat)

test_that('spatial_distribution_report runs', {
  cell_seg_path =
    file.path('test_data',
              'FIHC4__0929309_HP_IM3_2_cell_seg_data.txt') %>%
    normalizePath

  pairs = list(
    c("B", 'Cytotoxic T'),
    c("B", 'Helper T'))
  colors = c(B='red', 'Helper T'='green', 'Cytotoxic T'="yellow")
  out_path = tempfile()

  expect_warning(spatial_distribution_report(cell_seg_path, pairs, colors,
                                             output_path=out_path),
                 'missing values', all=TRUE)

  expect_true(file.exists(out_path))
  file.remove(out_path)
})
