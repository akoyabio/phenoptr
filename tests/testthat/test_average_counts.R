# Tests for average count functions
library(testthat)
library(dplyr)

check_within = function(within, radius, from, to) {
  expect_equal(nrow(within), length(radius))
  expect_equal(names(within),
               c("radius", "from_count", "to_count",
                 "from_with", "within_mean"))
  expect_equal(within$from_count, from)
  expect_equal(within$to_count, to)
}

# Smoke test of average_counts. Mostly checks that it doesn't barf and
# that the correct cells are selected.
test_that("count_within works", {
  csd = sample_cell_seg_data %>% filter(Phenotype != 'other')
  dst = distance_matrix(csd)

  within15 = count_within(csd, 'tumor', 'cytotoxic CD8', 15, dst=dst)
  check_within(within15, 15, 3303, 293)

  within30 = count_within(csd, 'tumor', 'cytotoxic CD8', 30, dst=dst)
  check_within(within30, 30, 3303, 293)
  expect_gt(within30$within_mean, within15$within_mean)

  within15tumor = count_within(csd, 'tumor', 'cytotoxic CD8', 15,
                                       'tumor', dst=dst)
  check_within(within15tumor, 15, 3221, 129)

  within = count_within(csd, 'tumor',
                        c('cytotoxic CD8', 'helper CD4'), 15, 'tumor', dst=dst)
  check_within(within, 15, 3221, 129+6)

  within = count_within(csd,
                        c('cytotoxic CD8', 'helper CD4'),
                        'tumor',  15, 'tumor', dst=dst)
  check_within(within, 15, 129+6, 3221)
})

test_that("count_within works with no data", {
  csd = sample_cell_seg_data %>% filter(Phenotype != 'other')
  dst = distance_matrix(csd)

  within = count_within(csd, 'other', 'cytotoxic CD8', 15, dst=dst)
  check_within(within, 15, 0, 293)
  expect_equal(within$from_with, 0)
  expect_equal(within$within_mean, 0)

  within = count_within(csd, 'tumor', 'other', 15, dst=dst)
  check_within(within, 15, 3303, 0)
  expect_equal(within$from_with, 0)
  expect_equal(within$within_mean, 0)

  within = count_within(csd, 'other', 'other', 15, dst=dst)
  check_within(within, 15, 0, 0)
  expect_equal(within$from_with, 0)
  expect_equal(within$within_mean, 0)
})

test_that("count_within works with multiple radii", {
  csd = sample_cell_seg_data %>% filter(Phenotype != 'other')
  dst = distance_matrix(csd)

  within = count_within(csd, 'tumor', 'cytotoxic CD8', c(15, 30), dst=dst)
  check_within(within, c(15, 30), c(3303, 3303), c(293, 293))

  within = count_within(csd, 'other', 'other', c(15, 30), dst=dst)
  check_within(within, c(15, 30), c(0, 0), c(0, 0))
  expect_equal(within$from_with, c(0, 0))
  expect_equal(within$within_mean, c(0, 0))
})
