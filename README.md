# phenoptr <img src="man/figures/Akoya.png" align="right" width="175" />

## Helpers for working with inForm data

`phenoptr` contains functions that make it easier to read and analyze data tables
and images created by Akoya Biosciences' inForm<sup>&reg;</sup> software.

`phenoptr` is part of the Akoya Biosciences Phenoptics&trade; family of
Quantitative Pathology Research Solutions. For more information visit the 
[Phenoptics&trade; home page](https://www.akoyabio.com/phenopticstm/technology/quantitative-pathology-research).

----

## Installation

`phenoptr` requires the [R environment](https://www.r-project.org/) 
for statistical computing, version 4.0.0 or higher. To install R,
visit the [R download](https://cloud.r-project.org/) page.
The [RStudio IDE](https://www.rstudio.com/products/rstudio/)
is highly recommended as well.

1. Install R. Download the most recent version from  https://cloud.r-project.org/.
1. Install RStudio. Download the desktop version from https://www.rstudio.com/products/rstudio/.
1. Start RStudio.
1. Install `phenoptr` from GitHub. In the RStudio console, copy and paste or type these commands (press Enter after each line):
```
install.packages("devtools")
devtools::install_github("akoyabio/phenoptr")
```
1. When requested, enter `1` (Yes) to install BiocInstaller.

### Optional packages

- **Spatial metrics**
  The Akoya Biosciences `rtree` package 
  dramatically speeds calculation and reduces memory requirements of
  spatial metrics such as nearest neighbors and count within.
  See the [installation instructions](https://github.com/akoyabio/rtree/tree/master#installation)
  in the package README file.

----

## Getting Started

These [Tutorials](https://akoyabio.github.io/phenoptr/articles/index.html)
introduce the most important features of `phenoptr`.

- [Reading and exploring inForm tables](https://akoyabio.github.io/phenoptr/articles/reading_tables.html)
demonstrates reading and processing inForm cell segmentation tables. This is a 
good place to start.
- [Computing inter-cellular distances](https://akoyabio.github.io/phenoptr/articles/computing_distances.html)
introduces most of `phenoptr`'s spatial processing capabilities---finding
nearest
neighbor distances, counting cells within a radius, and visualizing nearest
neighbors.
- [Find and count touching cells](https://akoyabio.github.io/phenoptr/articles/find_and_count_touching_cells.html)
shows how to count touching cells of paired phenotypes.

<div class="panel panel-default"><div class="panel-body">
For examples showing aggregation across multiple inForm fields and multiple
slides, see the Tutorials in the
<a href="https://akoyabio.github.io/phenoptrExamples">phenoptrExamples</a>
package.
</div></div>

### Learning R

R is a powerful and popular environment for data manipulation and 
statistical analysis. Many learning resources are available online.

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

See the [Reference](https://akoyabio.github.io/phenoptr/reference/index.html)
section of the documentation for details on individual functions.

To cite package `phenoptr` in publications use:

&nbsp;&nbsp;Kent S Johnson (2021). phenoptr: inForm Helper Functions. R package version 0.2.10.  
&nbsp;&nbsp;https://akoyabio.github.io/phenoptr/

<!-- badges: start -->
[![R build status](https://github.com/akoyabio/phenoptr/workflows/R-CMD-check/badge.svg)](https://github.com/akoyabio/phenoptr/actions)
<!-- badges: end -->

