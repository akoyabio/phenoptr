# Tests for distance functions
# Mostly smoke tests to make sure the code runs...
context('distance_functions')
library(testthat)

path = test_path('test_data',
              'FIHC4__0929309_HP_IM3_2_cell_seg_data.txt')

# Compute the Euclidean distance between two cells given their row numbers
cell_dist = function(csd, ix1, ix2) {
  d = csd[c(ix1, ix2), c('Cell X Position', 'Cell Y Position')]
  purrr::map_dbl(d, ~diff(.)^2) %>% sum %>% sqrt
}

test_that("distance_matrix works", {
  csd = read_cell_seg_data(path)
  dst = distance_matrix(csd)
  expect_equal(dim(dst), c(nrow(csd), nrow(csd)))

  # Spot check a few values by hand-calculation
  pairs = list(c(1, 2), c(10, 90), c(200, 250))
  for (pair in pairs) {
    ix1 = pair[1]
    ix2 = pair[2]
    info = paste(pair, collapse=', ')
    expect_equal(dst[ix1, ix2], dst[ix2, ix1], info=info) # Should be symmetric
    expect_equal(dst[ix1, ix2], cell_dist(csd, ix1, ix2), info=info)
  }

  s = subset_distance_matrix(csd, dst, 'B', 'Helper T')
  expect_equal(sum(csd$Phenotype=='B'), 249) # Make sure this is a real test
  expect_equal(sum(csd$Phenotype=='Helper T'), 15)
  expect_equal(dim(s)[1], sum(csd$Phenotype=='B'))
  expect_equal(dim(s)[2], sum(csd$Phenotype=='Helper T'))
  expect_error(subset_distance_matrix(dst, csd, 'B', 'Helper T'), 'wrong order')
})

test_that("find_nearest_distance works", {
  csd = read_cell_seg_data(path)
  phenos = sort(unique(csd$Phenotype))
  nearest = find_nearest_distance_dist(csd)
  expect_equal(ncol(nearest), 2*length(phenos))
  expect_equal(nrow(nearest), nrow(csd))
  expect_equal(names(nearest)[seq(1, 2*length(phenos), 2)],
               paste('Distance to', phenos))
  expect_equal(names(nearest)[seq(2, 2*length(phenos), 2)],
               paste('Cell ID', phenos))

  # Spot check some cells where nearest neighbors are clear from inspection
  dst = distance_matrix(csd)
  pairs = list(
    c(64, 118),
    c(16, 30),
    c(187, 214),
    c(209, 215)
  )
  for (pair in pairs) {
    ix1 = pair[1]
    ix2 = pair[2]
    info = paste(ix1, ix2, sep=', ')
    pheno1 = paste('Distance to', csd$Phenotype[ix1])
    pheno2 = paste('Distance to', csd$Phenotype[ix2])
    dist = as.numeric(dst[ix1, ix2])
    expect_equal(as.numeric(nearest[ix1, pheno2]), dist, info=info)
    expect_equal(as.numeric(nearest[ix2, pheno1]), dist, info=info)
  }

  # Compare with rtree version
  nearest_rtree = find_nearest_distance_rtree(csd)

  # There are three places where Cell ID B may be different
  # because there are two cells at the same distance and the two methods
  # find different ones.
  # Just check the distances for these
  diff_b = c(76, 154, 269)
  distance_b = c(9.055385138, 6.800735254, 5.385164807)
  expect_equal(nearest$`Distance to B`[diff_b], distance_b)
  expect_equal(nearest_rtree$`Distance to B`[diff_b], distance_b)

  # The rest should match
  # all.equal.tbl_df is broken, coerce to data.frame.
  expect_equal(as.data.frame(nearest)[-diff_b, ],
               as.data.frame(nearest_rtree)[-diff_b, ])

  phenos = c(phenos, 'other') # A mising phenotype
  nearest = find_nearest_distance_dist(csd, phenos)
  expect_equal(ncol(nearest), 2*length(phenos))
  expect_equal(nearest$`Distance to other`, rep(NA_real_, nrow(csd)))

  nearest_rtree = find_nearest_distance_rtree(csd, phenos)
  expect_equal(nearest$`Distance to other`, nearest_rtree$`Distance to other`)

  # Force an error
  csd$`Sample Name`[1] = 'foo.im3'
  expect_error(find_nearest_distance(csd), 'multiple')
})

test_that('compute_all_nearest_distance works', {
  merge_path = test_path('test_data',
              'merge/FIHC4_merge_cell_seg_data.txt')
  out_path = tempfile()
  compute_all_nearest_distance(merge_path, out_path)
  expect_true(file.exists(out_path))
  all_d = readr::read_tsv(out_path, na=c('NA', '#N/A'), col_types=readr::cols())
  expect_equal(nrow(all_d), 66+68+270+215)

  # The result for the sample data should be the same
  csd = read_cell_seg_data(path)
  nearby = find_nearest_distance(csd)
  all_d_2 = all_d %>%
    dplyr::filter(`Sample Name`=="FIHC4__0929309_HP_IM3_2.im3") %>%
    dplyr::arrange(`Cell ID`) %>%
    dplyr::select(`Distance to B`,
                  `Cell ID B`,
                  `Distance to Cytotoxic T`,
                  `Cell ID Cytotoxic T`,
                  `Distance to Helper T`,
                  `Cell ID Helper T`)
  expect_equivalent(as.data.frame(all_d_2), as.data.frame(nearby))

  file.remove(out_path)

  # Repeat with a consolidated file
  merge_path =
    test_path('test_data',
              'consolidated/FIHC4_consolidated_merge_cell_seg_data.txt')
  out_path = tempfile()
  compute_all_nearest_distance(merge_path, out_path)
  expect_true(file.exists(out_path))
  all_d = readr::read_tsv(out_path, na=c('NA', '#N/A'), col_types=readr::cols())
  expect_equal(nrow(all_d), 66+68+270+215)

  # The result for the sample data should be the same
  # Column names are different
  all_d_2 = all_d %>%
    dplyr::filter(`Sample Name`=="FIHC4__0929309_HP_IM3_2.im3") %>%
    dplyr::arrange(`Cell ID`) %>%
    dplyr::select(`Distance to B`=`Distance to B+`,
                  `Cell ID B` = `Cell ID B+`,
                  `Distance to Cytotoxic T`=`Distance to Cytotoxic_T+`,
                  `Cell ID Cytotoxic T` = `Cell ID Cytotoxic_T+`,
                  `Distance to Helper T`=`Distance to Helper_T+`,
                  `Cell ID Helper_T` = `Cell ID Helper_T+`)
  expect_equivalent(as.data.frame(all_d_2), as.data.frame(nearby))

  file.remove(out_path)
})

test_that("which_row_min works", {
  tests = list(
    # Check all possibilities with rows of length 2
    c(1, 2),
    c(2, 1),
    c(0, 2),
    c(2, 0),
    c(NA, 2),
    c(2, NA),
    c(0, NA),
    c(NA, 0),
    c(NA, NA),
    c(0, 0),

    # Check a few longer ones too
    c(5, 4, 0),
    c(5, NA, 0),
    c(0, 4, 0),
    c(0, NA, 0),
    c(0, 0, 0),
    c(NA, NA, NA)
  )

  # Indices in the vectors of `test` which have the smallest
  # non-zero, non-NA value.
  expected = c(1, 2, 2, 1, 2, 1, NA, NA, NA, NA,
               2, 1, 2, NA, NA, NA)

  purrr::walk2(tests, expected,
               ~expect_equal(which_row_min(.x), .y, info=.x))
})
