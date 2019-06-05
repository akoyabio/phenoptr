context('parse_phenotypes')
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
  sels = parse_phenotypes("CD3+", "CD3+/CD8-",
                          All="Total Cells", Macrophage="CD68+,CD163+")
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

test_that('parse_phenotypes works with formulae', {
  expect_equal(parse_phenotypes(PDL1='~`Membrane PDL1`>1'),
               list(PDL1=~`Membrane PDL1`>1))
  expect_equal(parse_phenotypes('CD8+', PDL1='~`Membrane PDL1`>1'),
               list(`CD8+`='CD8+', PDL1=~`Membrane PDL1`>1))
  expect_equal(parse_phenotypes(Mixed='CD8+/~`Membrane PDL1`>1'),
               list(Mixed=list('CD8+', ~`Membrane PDL1`>1)))
  expect_equal(parse_phenotypes('CD3+', Mixed='CD8+/~`Membrane PDL1`>1'),
               list(`CD3+`='CD3+', Mixed=list('CD8+', ~`Membrane PDL1`>1)))
})

test_that('parse_phenotypes error checking works', {
  expect_error(parse_phenotypes('CD3+/CD8+,CD68+'))
  expect_error(parse_phenotypes('CD3'))

  # Parse phenotypes expects string arguments
  expect_error(parse_phenotypes(PDL1=~`Membrane PDL1 (Opal 520) Mean`>1))

  # ORing of ordinary phenotype and formula is not supported
  # by the select_rows syntax and thus not here.
  expect_error(parse_phenotypes('CD3+,~`Membrane PDL1`>1'))
})

test_that('split_and_trim works', {
  expect_equal(split_and_trim(' xx / yy ', '/'), c('xx', 'yy'))
  expect_equal(split_and_trim(' xx , yy ', ','), c('xx', 'yy'))
  expect_error(split_and_trim(c('xx', 'yy'), '/'))
})

test_that('validate_phenotype_definitions works', {
  # These are all valid regardless of the second argument
  expect_equal(validate_phenotype_definitions(NULL, ''), '')
  expect_equal(validate_phenotype_definitions('', ''), '')
  expect_equal(validate_phenotype_definitions('Total', ''), '')
  expect_equal(validate_phenotype_definitions('All', ''), '')
  expect_equal(validate_phenotype_definitions('~`Mean x`>3', ''), '')

  expect_equal(validate_phenotype_definitions('CD3+', 'CD3'), '')
  expect_equal(validate_phenotype_definitions('CD3+/~`Mean x`>3', 'CD3'), '')

  expect_match(validate_phenotype_definitions('CD3', 'CD3'), 'must start')
  expect_match(validate_phenotype_definitions('CD3+', 'CD8'), 'Unknown')
  expect_match(validate_phenotype_definitions('CD3+,~`Mean x`>3', 'CD3'),
               'not allowed')

  df = tibble::tibble(x=1:2)
  expect_equal(validate_phenotype_definitions('~x==1', '', df), '')
  expect_match(validate_phenotype_definitions('~x==', ''),
               'not a valid expression')
  expect_match(validate_phenotype_definitions('~y==1', '', df),
               'not found')
})
