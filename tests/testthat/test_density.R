# Test for density_bands
library(testthat)
library(dplyr)

check_density_bands <- function(values, pixels_per_micron) {
  expect_equal(names(values), c("densities", "cells", "distance"))

  densities = values$densities
  expect_equal(range(densities$midpoint), c(-112.5, 212.5))

  dens_summary = densities %>%
    group_by(phenotype) %>%
    summarize(count=sum(count), pixels=sum(area)*pixels_per_micron^2)

  # Sum of pixels is image size
  expect_equal(dens_summary$pixels, rep(1868*1400, 3), tolerance=1)

  # Cell counts should match cell_seg_data
  csd_summary = sample_cell_seg_data %>%
    group_by(Phenotype) %>%
    summarize(count=n()) %>%
    right_join(dens_summary, by=c(Phenotype='phenotype'))

  expect_equal(nrow(csd_summary), 3)
  expect_equal(csd_summary$count.x, csd_summary$count.y)

  # Now by tissue
  dens_summary = densities %>%
    mutate(tissue = ifelse(midpoint<0, 'Tumor', 'Stroma')) %>%
    group_by(tissue, phenotype) %>%
    summarize(count=sum(count))

  # Cell counts should be close to cell_seg_data
  # There are differences at the tissue boundary...
  csd_summary = sample_cell_seg_data %>%
    group_by(`Tissue Category`, Phenotype) %>%
    summarize(count=n()) %>%
    right_join(dens_summary,
               by=c('Tissue Category'='tissue', Phenotype='phenotype'))

  expect_equal(nrow(csd_summary), 6)
  expect_equal(csd_summary$count.x, csd_summary$count.y, tolerance=4)
}

test_that('density_bands works', {
  # Compute density for the sample data
  values <- density_bands(sample_cell_seg_path(),
    list("CD8+", "CD68+", "FoxP3+"),
    positive="Stroma", negative="Tumor")

  check_density_bands(values, pixels_per_micron=2.007538)

  # Try with consolidated format
  cell_seg_path = testthat::test_path('test_data', 'consolidated',
                  'Set4_1-6plex_[16142,55840]_cons_cell_seg_data.txt')

  # Get the path to the segmentation map & component data for the sample data,
  # it is the same image as the consolidated data.
  map_path = sub('_cell_seg_data.txt', '_binary_seg_maps.tif',
                 sample_cell_seg_path())
  component_path = sub('_cell_seg_data.txt', '_component_data.tif',
                 sample_cell_seg_path())
  values <- density_bands(cell_seg_path,
                          list("CD8+", "CD68+", "FoxP3+"),
                          positive="Stroma", negative="Tumor",
                          map_path=map_path, component_path=component_path)

  check_density_bands(values, pixels_per_micron=2.007538)


})
