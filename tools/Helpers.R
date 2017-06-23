# Package development help
# Hat tip to http://www.masalmon.eu/2017/06/17/automatictools/

# Build the package reference. Output will be in a temp dir, look at
# the console output to see where it is.
devtools::check(manual = TRUE)

# Linting
devtools::install_github('jimhester/lintr')
lintr::lint_package()

# Good practices
# This takes a while to run and requires working R CMD CHECK
devtools::install_github('mangothecat/goodpractice')
goodpractice::gp()

devtools::spell_check()

# Build documentation site - not working for NEWS.md?
devtools::install_github("hadley/pkgdown")
pkgdown::build_site()

# Hack for NEWS
rmarkdown::render("NEWS.md", output_file='docs/news/index.html')
