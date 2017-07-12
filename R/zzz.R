.onLoad <- function(libname, pkgname)
{
    # Set default options() related to functionality in 'informr' pkg
    op <- options()
    my_opts = list(
      informr.pixels.per.micron=2)
    toset = !(names(my_opts) %in% names(op))
    if(any(toset)) options(my_opts[toset])
}
