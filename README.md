# phenoptr - Helpers for working with inForm data

`phenoptr` contains functions that make it easier to read and analyze data tables
and images created by PerkinElmer's inForm<sup>&reg;</sup> software.

`phenoptr` is part of the PerkinElmer Phenoptics&trade; family of
Quantitative Pathology Research Solutions. For more information
visit the Phenoptics&trade;
[home page](http://www.perkinelmer.com/cancer-immunology/index.html).

## Installation

Installation is from GitHub:

```
# install.packages("devtools")
devtools::install_github("PerkinElmer/phenoptr", build_vignettes=TRUE)
```

## Usage

`phenoptr` is designed to work in harmony with packages in the 
[tidyverse](http://tidyverse.org/). 

- [readr](http://readr.tidyverse.org/) is used to read data files.
- A [tibble](http://tibble.tidyverse.org/) (also known as `data_frame`) 
  is the preferred representation of tabular data.
- [dplyr](http://dplyr.tidyverse.org/), [purrr](http://purrr.tidyverse.org/) 
  and the pipe operator ([`%>%`](http://magrittr.tidyverse.org/)) 
  are used extensively in 
  package code and examples.


If you'd like to learn how to use the tidyverse, 
a good place to start is Garrett Grolemund and Hadley Wickham's book,
available free online at
[R for data science](http://r4ds.had.co.nz/).

### Package summary

#### Read and preprocess inForm data and image files

[list_cell_seg_files](https://perkinelmer.github.io/phenoptr/reference/list_cell_seg_files.html)
finds all inForm cell segmentation files in a directory.

[read_cell_seg_data](https://perkinelmer.github.io/phenoptr/reference/read_cell_seg_data.html)
reads inForm cell segmentation
data tables and makes them easier to use in R.

[read_components](https://perkinelmer.github.io/phenoptr/reference/read_components.html) 
reads inForm component data and [read_maps](https://perkinelmer.github.io/phenoptr/reference/read_maps.html) 
reads inForm
segmentation maps.

#### Compute and visualize spatial relationships between cells

[find_nearest_distance](https://perkinelmer.github.io/phenoptr/reference/find_nearest_distance.html)
finds the distance
from each cell to the nearest cell of each phenotype.

[count_within](https://perkinelmer.github.io/phenoptr/reference/count_within.html)
and 
[count_within_batch](https://perkinelmer.github.io/phenoptr/reference/count_within_batch.html) 
count the number of cells within a
fixed radius of other cells.

[spatial_distribution_report](https://perkinelmer.github.io/phenoptr/reference/spatial_distribution_report.html) 
creates a report visualizing nearest neighbor
relationships between cells of two phenotypes.

[count_touching_cells](https://perkinelmer.github.io/phenoptr/reference/count_touching_cells.html)
uses morphological
analysis to find, count, and visualize touching cells in paired
phenotypes.

#### Helpers

[distance_matrix](https://perkinelmer.github.io/phenoptr/reference/distance_matrix.html)
and 
[subset_distance_matrix](https://perkinelmer.github.io/phenoptr/reference/subset_distance_matrix.html)
create and subset cell distance
matrices from cell segmentation data.

[select_rows](https://perkinelmer.github.io/phenoptr/reference/select_rows.html)
helps select rows in cell segmentation
data corresponding to specific phenotypes and expression levels.

### Full documentation

See the 
[Articles](https://perkinelmer.github.io/phenoptr/articles/index.html) and 
[Reference](https://perkinelmer.github.io/phenoptr/reference/index.html)
sections of the full documentation for details.

For extended examples and more sample data see the 
[phenoptrExamples](https://perkinelmer.github.io/phenoptrExamples) package.

[![Travis-CI Build Status](https://travis-ci.org/PerkinElmer/phenoptr.svg?branch=master)](https://travis-ci.org/PerkinElmer/phenoptr)
[![Coverage Status](https://img.shields.io/codecov/c/github/PerkinElmer/phenoptr/master.svg)](https://codecov.io/github/PerkinElmer/phenoptr?branch=master)
