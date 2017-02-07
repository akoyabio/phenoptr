# Tests for select_rows
library(dplyr)
library(testthat)

test_data = data_frame(
  Phenotype = c(rep('tumor', 3), rep('cd8', 2)),
  Expr = c(1:2, 1:3),
  E2 = c(1, 1, 2, 2, 1),
  E3 = c(1, 2, 1, 2, 1)
)

test_that("Select phenotype works", {
  expect_equal(select_rows(test_data, 'tumor'), c(T, T, T, F, F))
  expect_equal(select_rows(test_data, 'cd8'), c(F, F, F, T, T))
})

test_that("Alternate name works", {
  td = test_data %>% rename(pheno=Phenotype)
  expect_equal(select_rows(td, 'tumor', 'pheno'), c(T, T, T, F, F))
  expect_equal(select_rows(td, 'cd8', 'pheno'), c(F, F, F, T, T))
})

test_that('Single expression works', {
  expect_equal(select_rows(test_data, ~(Expr==1)), c(T, F, T, F, F))
  expect_equal(select_rows(test_data, ~E2==1), c(T, T, F, F, T))
})

test_that('Phenotype + expression works', {
  expect_equal(select_rows(test_data, list('tumor', ~Expr==1)), c(T, F, T, F, F))
  expect_equal(select_rows(test_data, list(~Expr==1, 'tumor')), c(T, F, T, F, F))
  expect_equal(select_rows(test_data, list('cd8', ~E2==1)), c(F, F, F, F, T))
  expect_equal(select_rows(test_data, list(~E2==1, 'cd8')), c(F, F, F, F, T))
})

test_that('Multiple phenotypes work', {
  expect_equal(select_rows(test_data, c('tumor', 'cd8')), c(T, T, T, T, T))
  expect_equal(select_rows(test_data, list(c('tumor', 'cd8'))), c(T, T, T, T, T))
  expect_equal(select_rows(test_data, list(~E3==1, c('tumor', 'cd8'))),
               c(T, F, T, F, T))
})

test_that('Phenotype in as list are ANDed', {
  expect_equal(select_rows(test_data, list('tumor', 'cd8')), c(F, F, F, F, F))
})

test_that('Multiple expressions work', {
  expect_equal(select_rows(test_data, list(~Expr==1, ~E2==2)), c(F, F, T, F, F))
  expect_equal(select_rows(test_data, list(~E2==1, ~Expr==2)), c(F, T, F, F, F))
})

test_that('All together now', {
  expect_equal(select_rows(test_data, list('tumor', ~E2==1, ~E3==1)), c(T, F, F, F, F))
  expect_equal(select_rows(test_data, list(~E2==1, 'tumor', ~E3==1)), c(T, F, F, F, F))
})

test_that('Normalize selector works', {
  expect_error(normalize_selector(NULL))
  expect_error(normalize_selector(c()))
  expect_error(normalize_selector('cd8'))
  expect_error(normalize_selector(list()))

  expect_equal(normalize_selector(list(a='b')), list(a='b'))

  expect_equal(normalize_selector(list('cd8')), list(cd8='cd8'))
  expect_equal(normalize_selector(list(~Expr==1)), list(`Expr == 1`= ~Expr==1))

  # Note these are all different!
  # Three selectors, each an individual phenotype
  expect_equal(normalize_selector(list('a', 'b', 'c')), list(a='a', b='b', c='c'))

  # One selector, OR of three phenotypes
  expect_equal(normalize_selector(list(c('a', 'b', 'c'))),
               list(`a|b|c`=c('a', 'b', 'c')))

  # One selector, AND of three phenotypes
  expect_equal(normalize_selector(list(list('a', 'b', 'c'))),
               list(`a&b&c`=list('a', 'b', 'c')))

  # Mixed phenotype and threshold
  expect_equal(normalize_selector(list(list('cd8', ~Expr==1))),
               list(`cd8&Expr == 1`=list('cd8', ~Expr==1)))
  expect_equal(normalize_selector(list(list('a', ~Expr==1, 'c'))),
               list(`a&Expr == 1&c`=list('a', ~Expr==1, 'c')))

  # Multiple selectors
  expect_equal(normalize_selector(list(
      list('cd8', ~Expr==1),
      c('a', 'b', 'c'),
      'tumor'
    )),
    list(
      `cd8&Expr == 1`=list('cd8', ~Expr==1),
      `a|b|c`=c('a', 'b', 'c'),
      tumor='tumor'
    ))
})
