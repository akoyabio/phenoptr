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

test_that('Multiple expressions work', {
  expect_equal(select_rows(test_data, list(~Expr==1, ~E2==2)), c(F, F, T, F, F))
  expect_equal(select_rows(test_data, list(~E2==1, ~Expr==2)), c(F, T, F, F, F))
})

test_that('All together now', {
  expect_equal(select_rows(test_data, list('tumor', ~E2==1, ~E3==1)), c(T, F, F, F, F))
  expect_equal(select_rows(test_data, list(~E2==1, 'tumor', ~E3==1)), c(T, F, F, F, F))
})
