# phenoptr - Helpers for working with inForm data

[![Travis-CI Build Status](https://travis-ci.org/PerkinElmer/phenoptr.svg?branch=master)](https://travis-ci.org/PerkinElmer/phenoptr)
[![Coverage Status](https://img.shields.io/codecov/c/github/PerkinElmer/phenoptr/master.svg)](https://codecov.io/github/PerkinElmer/phenoptr?branch=master)

`phenoptr` contains functions that make it easier to read and analyze data tables
and images created by PerkinElmer's inForm<sup>&reg;</sup> software.

`phenoptr` is part of the PerkinElmer Phenoptics&trade; family of
Quantitative Pathology Research Solutions. For more information
visit the Phenoptics&trade;
[home page](http://www.perkinelmer.com/cancer-immunology/index.html).

----

## Installation

`phenoptr` requires the R environment for statistical computing. To install R,
see the [R home page](https://www.r-project.org/).
The [RStudio IDE](https://www.rstudio.com/products/rstudio/)
is highly recommended as well.

Installation of `phenoptr` is from GitHub:

```
install.packages("devtools")
devtools::install_github("PerkinElmer/phenoptr", build_vignettes=TRUE)
```

----

## Getting Started

The Tutorials walk through most of `phenoptr`.

- [Reading and exploring inForm tables](https://perkinelmer.github.io/phenoptr/articles/reading_tables.html)
demonstrates reading and processing inForm cell segmentation tables. This is a 
good place to start.
- [Computing inter-cellular distances](https://perkinelmer.github.io/phenoptr/articles/computing_distances.html)
introduces most of `phenoptr`'s spatial processing capabilities.
- [Find and count touching cells](https://perkinelmer.github.io/phenoptr/articles/find_and_count_touching_cells.html)
covers the remaining spatial processing functions.

<div class="panel panel-default"><div class="panel-body">
For extended examples and more sample data see the Tutorials in the
<a href="https://perkinelmer.github.io/phenoptrExamples">phenoptrExamples</a>
package.</div></div>

### Learning R

R is a powerful and popular environment for data manipulation. 
Search for [learn R](https://www.google.com/search?q=learn+r) to find many 
resources for beginners.

`phenoptr` is designed to work in harmony with packages in the 
[tidyverse](http://tidyverse.org/). 

- [readr](http://readr.tidyverse.org/) is used to read data files.
- A [tibble](http://tibble.tidyverse.org/) (also known as `data_frame`) 
  is the preferred representation of tabular data.
- [dplyr](http://dplyr.tidyverse.org/), [purrr](http://purrr.tidyverse.org/) 
  and the pipe operator ([%>%](http://magrittr.tidyverse.org/)) 
  are used extensively in 
  package code and examples.

If you'd like to learn more about the tidyverse packages, 
a good place to start is Garrett Grolemund and Hadley Wickham's book,
available free online at
[R for data science](http://r4ds.had.co.nz/).
If you are new to R, the book's
[Introduction](http://r4ds.had.co.nz/introduction.html)
will help you get started.

----

## Full documentation

See the [Reference](https://perkinelmer.github.io/phenoptr/reference/index.html)
section of the documentation for details on individual functions.

[<i class='fa fa-smile-o'></i>](articles/README1000.html)
