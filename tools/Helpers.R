# Package development help

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
