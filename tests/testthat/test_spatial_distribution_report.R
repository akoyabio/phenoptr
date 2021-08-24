# Create a spatial distribution report
library(testthat)

test_that('spatial_distribution_report runs', {
  cell_seg_path =
    test_path('test_data',
              'FIHC4__0929309_HP_IM3_2_cell_seg_data.txt') %>%
    normalizePath

  pairs = list(
    c("B", 'Cytotoxic T'),
    c("B", 'Helper T'))
  colors = c(B='red', 'Helper T'='green', 'Cytotoxic T'="yellow")

  # Make sure tempdir exists
  tempdir(check=TRUE)
  out_path = tempfile(fileext='.html')

  spatial_distribution_report(cell_seg_path, pairs, colors,
                                             output_path=out_path)

  expect_true(file.exists(out_path))
  file.remove(out_path)
})

test_that('spatial_distribution_report works with consolidated data', {
  # Hack some test data to avoid keeping more copies of stuff
  raw_path =
    test_path('test_data', 'consolidated',
              'FIHC4_consolidated_merge_cell_seg_data.txt') %>%
    normalizePath

  csd = read_cell_seg_data(raw_path) %>%
    dplyr::filter(`Sample Name`=="FIHC4__0929309_HP_IM3_2.im3")

  temp_dir = tempdir()
  cell_seg_path =  file.path(temp_dir,
                             "FIHC4__0929309_HP_IM3_2_cell_seg_data.txt")
  readr::write_tsv(csd, cell_seg_path)
  file.copy(file.path('test_data',
                      'FIHC4__0929309_HP_IM3_2_component_data.tif'),
            file.path(temp_dir,
                      'FIHC4__0929309_HP_IM3_2_component_data.tif'))

  pairs = list(
    c("B+", 'Cytotoxic_T+'),
    c("B+", 'Helper_T+'))
  colors = c('B+'='red', 'Helper_T+'='green', 'Cytotoxic_T+'="yellow")
  out_path = tempfile(fileext='.html')

  spatial_distribution_report(cell_seg_path, pairs, colors,
                                             output_path=out_path)

  expect_true(file.exists(out_path))
  file.remove(out_path)
})
