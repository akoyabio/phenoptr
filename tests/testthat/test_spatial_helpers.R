# Tests for spatial_helpers

test_that('make_ppp works', {
  # Try a couple of different ways to make a ppp
  pp <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
                 "CD8+", tissue_categories="Tumor")
  expect_equal(pp$n, 45)
  expect_equal(levels(spatstat::marks(pp)), "CD8+")

  pp <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
                 pheno=list(CD8="CD8+"),
                 field_name="Set4_1-6plex_[16142,55840].im3",
                 tissue_categories="Tumor")
  expect_equal(pp$n, 45)
  expect_equal(levels(spatstat::marks(pp)), "CD8")

  pp <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
                 pheno=c(CD8="CD8+"))
  expect_equal(pp$n, sum(sample_cell_seg_data$Phenotype=='CD8+'))
  expect_equal(levels(spatstat::marks(pp)), "CD8")

  # Test compound phenotypes with some nonsense values
  # There are no actual multiple positive cells in the dataset so
  # this just checks that the syntax is handled without error.
  pp <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
                 pheno="CD8+/CD8+")
  expect_equal(pp$n, sum(sample_cell_seg_data$Phenotype=='CD8+'))

  pp <- make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
                 pheno=list(CD8=list('CD8+', 'CD8+')))
  expect_equal(pp$n, sum(sample_cell_seg_data$Phenotype=='CD8+'))
  expect_equal(levels(spatstat::marks(pp)), "CD8")

  # More than one phenotype is an error
  expect_error(make_ppp(sample_cell_seg_data, sample_cell_seg_folder(),
                        pheno=c("CD8+", "CK+")),
               'length\\(pheno)')
})
