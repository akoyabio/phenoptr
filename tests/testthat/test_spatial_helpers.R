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


test_that('read_phenochart_polygons works', {
  xml_path = test_path('test_data/test_annotations.xml')
  rois = read_phenochart_polygons(xml_path)

  # There are two ROIs. The first is not tagged and contains rectangles,
  # the second is tagged and has no rectangles
  expect_equal(rois$tags, c('', '#IncludeInResults #Test'))
  expect_equal(purrr::map_int(rois$rects, length), c(10, 0))

  expect_equal(as.numeric(sf::st_bbox(rois$rects[[1]])),
               c(11547.55, 42073.73, 17093.25, 47617.44),
               tolerance=0.01)
})
