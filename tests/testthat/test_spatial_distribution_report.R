# Create a spatial distribution report
library(testthat)

test_that('spatial_distribution_report runs', {
  skip_on_cran()

  cell_seg_path = sample_cell_seg_path()

  pairs = list(
    c("CK+", "CD8+"),
    c("CK+", "CD68+"))
  colors = c('CK+'="cyan", "CD68+"="magenta", "CD8+"="yellow")
  out_path = tempfile()

  spatial_distribution_report(cell_seg_path, pairs, colors,
    output_path=out_path)

  expect_true(file.exists(out_path))
  file.remove(out_path)
})
