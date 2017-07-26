# Package development help
# Hat tip to http://www.masalmon.eu/2017/06/17/automatictools/

# Notes:
# Adding a new function requires several touches:
# - The code for the function
# - Add to _pkgdowr.yml
# - Add to README.md

# The tip about including biocViews: in DESCRIPTION is from
# https://stackoverflow.com/questions/14343817/cran-package-depends-on-bioconductor-package-installing-error

# Build the package reference. Output will be in a temp dir, look at
# the console output to see where it is.
devtools::check(manual = TRUE)

# Linting
devtools::install_github('jimhester/lintr')
lintr::lint_package()

# Good practices
# This takes a while to run and requires working R CMD CHECK
devtools::install_github('mangothecat/goodpractice')
checks = goodpractice::all_checks()
checks_to_omit = c("lintr_assignment_linter")
checks = setdiff(checks, checks_to_omit)
(gp=goodpractice::gp(checks=checks))

# Run code coverage on a single file and show a nice report
Sys.setenv("NOT_CRAN"="true") # Make not-CRAN tests run
covr::report(covr::file_coverage('R/read_cell_seg_data.R',
                                 'tests/testthat/test_read_cell_seg_data.R'))

# Run coverage for the entire package, tests only
covr::report(covr::package_coverage())

# Run coverage for the entire package, including examples and vignettes
covr::report(covr::package_coverage(type='all'))

devtools::spell_check()

# Build documentation site. `build_site` seems to work better than
# the RStudio menu which calls `build_site_rstudio`.
devtools::install_github("hadley/pkgdown")
pkgdown::build_site()

