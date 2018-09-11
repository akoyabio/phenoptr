library(testthat)

check_results = function(sels) {
  expect_equal(sels[[1]], 'CD3+')
  expect_equal(sels[[2]], list('CD3+', 'CD8-'))
  expect_equal(sels[[3]], NA)
  expect_equal(sels[[4]], c('CD68+', 'CD163+'))
}

test_that('parse_phenotypes works with unnamed args', {
  # Unnamed args get self-named
  vals = c("CD3+", "CD3+/CD8-", "Total Cells", "CD68+,CD163+")
  sels = do.call(parse_phenotypes, as.list(vals))
  expect_equal(names(sels), vals)
  check_results(sels)
})

test_that('parse_phenotypes works with named args', {
  # Unnamed args get self-named
  sels = parse_phenotypes("CD3+", "CD3+/CD8-", All="Total Cells", Macrophage="CD68+,CD163+")
  expect_equal(names(sels), c("CD3+", "CD3+/CD8-", 'All', 'Macrophage'))
  check_results(sels)
})

test_that('parse_phenotypes works with spaces in args', {
  vals = c(" CD3+ ", " CD3+ / CD8- ", " Total Cells ", " CD68+ , CD163+ ")
  sels = do.call(parse_phenotypes, as.list(vals))
  expect_equal(names(sels), stringr::str_trim(vals))
  check_results(sels)
})

test_that('parse_phenotypes works with a single list arg', {
  args = list("CD3+", "CD3+/CD8-", All="Total Cells", Macrophage="CD68+,CD163+")
  sels = parse_phenotypes(args)
  expect_equal(names(sels), c("CD3+", "CD3+/CD8-", 'All', 'Macrophage'))
  check_results(sels)
})

test_that('parse_phenotypes error checking works', {
  expect_error(parse_phenotypes('CD3+/CD8+,CD68+'))
  expect_error(parse_phenotypes('CD3'))
  expect_error(parse_phenotypes(PDL1=~`Membrane PDL1 (Opal 520) Mean`>1))
})

test_that('split_and_trim works', {
  expect_equal(split_and_trim(' xx / yy ', '/'), c('xx', 'yy'))
  expect_equal(split_and_trim(' xx , yy ', ','), c('xx', 'yy'))
  expect_error(split_and_trim(c('xx', 'yy'), '/'))
})
