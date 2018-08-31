library(testthat)

test_that('parse_phenotypes works with unnamed args', {
  # Unnamed args get self-named
  vals = c("CD3+", "CD3+/CD8-", "Total Cells", "CD68+,CD163+")
  sels = do.call(parse_phenotypes, as.list(vals))
  expect_equal(names(sels), vals)
  expect_equal(sels[[1]], 'CD3+')
  expect_equal(sels[[2]], list('CD3+', 'CD8-'))
  expect_equal(sels[[3]], NA)
  expect_equal(sels[[4]], c('CD68+', 'CD163+'))
})

test_that('parse_phenotypes works with named args', {
  # Unnamed args get self-named
  sels = parse_phenotypes("CD3+", "CD3+/CD8-", All="Total Cells", Macrophage="CD68+,CD163+")
  expect_equal(names(sels), c("CD3+", "CD3+/CD8-", 'All', 'Macrophage'))
  expect_equal(sels[[1]], 'CD3+')
  expect_equal(sels[[2]], list('CD3+', 'CD8-'))
  expect_equal(sels[[3]], NA)
  expect_equal(sels[[4]], c('CD68+', 'CD163+'))
})

test_that('parse_phenotypes works with a single list arg', {
  args = list("CD3+", "CD3+/CD8-", All="Total Cells", Macrophage="CD68+,CD163+")
  sels = parse_phenotypes(args)
  expect_equal(names(sels), c("CD3+", "CD3+/CD8-", 'All', 'Macrophage'))
  expect_equal(sels[[1]], 'CD3+')
  expect_equal(sels[[2]], list('CD3+', 'CD8-'))
  expect_equal(sels[[3]], NA)
  expect_equal(sels[[4]], c('CD68+', 'CD163+'))
})

test_that('parse_phenotypes error checking works', {
  expect_error(parse_phenotypes('CD3+/CD8+,CD68+'))
  expect_error(parse_phenotypes('CD3'))
  expect_error(parse_phenotypes(PDL1=~`Membrane PDL1 (Opal 520) Mean`>1),
               regexp='formula')
})
