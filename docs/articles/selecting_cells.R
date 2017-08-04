## ----setup,echo=FALSE----------------------------------------------------
#knitr::opts_chunk$set(eval=FALSE)
library(phenoptr)

## ------------------------------------------------------------------------
csd = sample_cell_seg_data
rows = select_rows(csd, 'CK+')
sum(rows) # The number of selected rows

# Select just the desired rows by subsetting
ck = csd[rows, ]
dim(ck)

## ------------------------------------------------------------------------
dst = distance_matrix(csd) # Compute this just once and re-use it
count_within(csd, from='CK+', to='CD8+', radius=15, dst=dst)

## ------------------------------------------------------------------------
tcells = csd[select_rows(csd, c('CD8+', 'FoxP3+')), ]
dim(tcells)

count_within(csd, from='CK+', to=c('CD8+', 'FoxP3+'), radius=15, dst=dst)

## ------------------------------------------------------------------------
rows = select_rows(csd, list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3))
ck_pdl1 = csd[rows, ]
dim(ck_pdl1)

count_within(csd, from=list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3), 
             to='CD8+', radius=15, dst=dst)

## ------------------------------------------------------------------------
pairs = list(c('CK+', 'CD8+'))

## ------------------------------------------------------------------------
pairs = c('CK+', 'CD8+')

## ------------------------------------------------------------------------
pairs = list(c('CK+', 'CD8+'),
             c('CK+', 'CD68+'))

## ------------------------------------------------------------------------
pairs = c('PDL1+ CK+', 'T Cell')
phenotype_rules = list(
  'PDL1+ tumor'=list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3),
  'T Cell'=c('CD8+', 'FoxP3+'))

## ------------------------------------------------------------------------
pairs = list(
  c('PDL1+ CK+', 'T Cell'),
  c('PDL1+ CK+', 'CD68+'))
phenotype_rules = list(
  'PDL1+ CK+'=list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3),
  'T Cell'=c('CD8+', 'FoxP3+')
)

