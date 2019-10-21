# Package load function
# Don't report code coverage of this file
# nocov start

.onLoad <- function(libname, pkgname) { # nolint
  # Set default options() related to functionality in 'phenoptr' pkg
  # if they are not already set
  op <- options()
  my_opts = list(
    phenoptr.pixels.per.micron=2,
    use.rtree.if.available=TRUE)
  toset = !(names(my_opts) %in% names(op))
  if(any(toset)) options(my_opts[toset])
}
# nocov end
