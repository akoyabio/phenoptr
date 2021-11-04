#' Negate the y-axis of a simple features object
#'
#' Use `geom_sf_invert()` with `scale_sf_invert()` to match the orientation
#' of Polaris images.
#' @export
#' @rdname geom_sf_invert
#' @param stat The statistical transformation to use on the data for this layer,
#' as a string or Stat object.
#' @param ... Passed on to `ggplot2::geom_sf()` or `ggplot2::scale_y_continuous()`
#' @export
geom_sf_invert = function(stat=StatSfInvert, ...) {
  ggplot2::geom_sf(stat=stat, ...)
}

#' @export
#' @rdname geom_sf_invert
#' @usage NULL
#' @format NULL
StatSfInvert <- ggplot2::ggproto("StatSfInvert", ggplot2::StatSf,
  compute_group = function(data, scales, coord) {
    # Flip the data
    data[[ ggplot2:::geom_column(data) ]] =
      negate_y(data[[ ggplot2:::geom_column(data) ]])

    parent = ggplot2::ggproto_parent(ggplot2::StatSf, NULL)
    parent$compute_group(data, scales, coord)
  },

  required_aes = c("geometry")
)

#' @export
#' @rdname geom_sf_invert
#' @usage NULL
stat_sf_invert <- function(mapping = NULL, data = NULL, geom = "rect",
                    position = "identity", na.rm = FALSE, show.legend = NA,
                    inherit.aes = TRUE, ...) {
  ggplot2::layer_sf(
    stat = StatSfInvert,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_sf_invert
#' @export
scale_sf_invert = function(...) {
  ggplot2::scale_y_continuous(labels = function(br) as.character(-br), ...)
}

#' Negate the y-axis
#' @param sf_obj A simple features object
#' @return `negate_y()` returns `sf_obj` with the y-axis negated.
#' @export
#' @rdname geom_sf_invert
negate_y = function(sf_obj) {
  flip = matrix(c(1, 0, 0, -1), ncol=2) # Matrix to negate Y
  sf_obj * flip

}
