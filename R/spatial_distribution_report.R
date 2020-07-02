# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c(
    "Cell ID",
    "Cell X Position", "Cell Y Position",
    "To X Position", "To Y Position",
    "x", "y"
    ))

#' Create a spatial distribution report.
#'
#' Creates an HTML report visualizing the location
#' and nearest-neighbor relations for cells of pairs of phenotypes
#' in a single field.
#'
#' A cell seg data file and component data file for the field are required.
#' If a tissue segmentation
#' or composite image of the field
#' is available, it will be used as a background.
#'
#' The report is saved to the same directory as the input file.
#'
#' See the tutorial
#' [Selecting cells within a cell segmentation table](https://akoyabio.github.io/phenoptr/articles/selecting_cells.html)
#'for more on
#' the use of `pairs` and `phenotype_rules`.
#'
#' @param cell_seg_path Path to a cell seg data file.
#' @param pairs A list of pairs of phenotypes. Each entry is a two-element
#'   vector. The report will contain one section for each pair showing the
#'   nearest neighbor relations between the two phenotypes.
#' @param colors A named list of phenotype colors to use when drawing
#'   the output.
#' @param phenotype_rules (Optional) A named list.
#'   Item names are phenotype names and must match entries in `pairs`.
#'   Item values are selectors for [select_rows].
#' @param output_path Optional, path to the output HTML file. If omitted,
#' output will be written to the directory containing `cell_seg_path`.
#' @export
#' @family distance functions
#' @examples
#' \dontrun{
#' # This example creates a report in the current user's
#' # home or Document directory. It requires the phenoptrExamples
#' # package, which includes extended sample data.
#'
#' library(phenoptrExamples)
#' cell_seg_path = system.file("extdata", "samples",
#'                             "Set4_1-6plex_[16142,55840]_cell_seg_data.txt",
#'                             package = "phenoptrExamples")
#'
#' pairs <- list(
#'   c("CK+", "CD8+"),
#'   c("CK+", "CD68+"))
#' colors <- c('CK+'="cyan", "CD68+"="magenta", "CD8+"="yellow")
#' out_path <- path.expand('~/spatial_distribution_report.html')
#'
#' spatial_distribution_report(cell_seg_path, pairs, colors,
#'   output_path=out_path)
#'
#' # This example extends the previous example to restrict tumor cells to
#' # only PDL1+ cells.
#' # Note that spatial_distribution_report uses read_cell_seg_table
#' # so units are removed from the expression column names.
#'
#' pairs <- list(
#'   c("CK+ PDL1+", "CD8+"),
#'   c("CK+ PDL1+", "CD68+"))
#' colors <- c("CK+ PDL1+"="cyan", "CD68+"="magenta", "CD8+"="yellow")
#' phenotype_rules <- list(
#'   "CK+ PDL1+"=list("CK+", ~`Entire Cell PDL1 (Opal 520) Mean`>3))
#'
#' spatial_distribution_report(cell_seg_path, pairs, colors, phenotype_rules,
#'   output_path=out_path)
#' }
#' @md
spatial_distribution_report <- function(cell_seg_path, pairs, colors,
                      phenotype_rules=NULL, output_path=NULL
                      ) {
  stopifnot(endsWith(cell_seg_path, '_cell_seg_data.txt'))
  if (!file.exists(cell_seg_path))
    stop(paste('File not found:', cell_seg_path))

  # Make phenotype_rules for any not already specified
  phenotypes = unique(do.call(c, pairs))
  phenotype_rules = make_phenotype_rules(phenotypes, phenotype_rules)

  stopifnot(length(colors)==length(phenotypes),
            all(phenotypes %in% names(colors)))

  if (is.null(output_path))
    output_path = sub('_cell_seg_data.txt',
                      '_spatial_distribution.html', cell_seg_path)
  rmd_path = system.file("rmd", "SpatialDistributionReport.Rmd",
                         package="phenoptr")

  rmarkdown::render(rmd_path, 'html_document', output_path)
}

## Helpers
## Note these are not exported and the doc comments are not roxygen
### Spatial functions ###

# Get a point pattern for a single phenotype
# @param data Cell seg data
# @param pheno List of (name, color, select)
# @param window spatstat::owin containing the points
# @return List of (data subset, ppp, pheno)
phenotype_as_ppp = function(data, pheno, window) {
  # Keep just the rows corresponding to the desired expression and phenotypes
  # for our data and columns containing fields of interest
  fixed_columns = c('Cell ID', 'Cell X Position', 'Cell Y Position')
  data = data[select_rows(data, pheno$select), fixed_columns]

  data = stats::na.omit(data)

  # Convert to point pattern objects. Guess the dimensions...
  pp = with(data, spatstat::ppp(`Cell X Position`, `Cell Y Position`,
                              window=window,
                              marks=factor(rep(pheno$name, nrow(data)))))
  list(data=data, pp=pp, pheno=pheno)
}

# Compute nearest neighbor distance between positively expressing phenotypes.
# @param pheno_data1 Data for the 'from' phenotype.
# @param pheno_data2 Data for the 'from' phenotype.
# @return A data.frame showing ID and location of from and to points
# with the distance between them.
# @importFrom magrittr "%>%"
find_nearest_neighbor = function(pheno_data1, pheno_data2) {
  d = spatstat::nncross(pheno_data1$pp, pheno_data2$pp)
  names(d)[1] = 'Distance'
  from_data = pheno_data1$data %>%
    dplyr::select(`Cell ID`, `Cell X Position`, `Cell Y Position`)
  to_data = pheno_data2$data[d$which, ] %>%
    dplyr::select('To Cell ID'=`Cell ID`,
                  'To X Position'=`Cell X Position`,
                  'To Y Position'=`Cell Y Position`)
  result = cbind(from_data, d, to_data)
  result$which = NULL
  result
}

### Point plots ###

# Make a plot of all points in pheno_data with lines to nearest neighbors
# as defined by nn_dist
# @param pheno_data result of a call to getPhenotypes
# @param nn_dist result of a call to find_nearest_neighbor
# @param lineColor color for the connecting lines
# @return a ggplot object
nn_plot <- function(pheno_data1, pheno_data2, nn_dist,
                    background, xlim, ylim,
                    line_color='gray90', scale_color='white') {
  title = bquote('Nearest neighbor from'
                 ~ italic(.(pheno_data1$pheno$name))
                 ~ 'to'
                 ~ italic(.(pheno_data2$pheno$name)))
  nn_plot_impl(nn_dist, pheno_data1, pheno_data2, title,
               background, xlim, ylim,
               line_color, scale_color)
}

# Make a plot of mutual nearest neighbors in phenoData -
# points for which each is the NN of the other
nn_plot_mutual = function(pheno_data1, pheno_data2, nn_mutual,
                          background, xlim, ylim,
                          line_color='gray90', scale_color='white') {
  title = bquote('Mutual nearest neighbors for'
                 ~ italic(.(pheno_data1$pheno$name))
                 ~ 'and'
                 ~ italic(.(pheno_data2$pheno$name)))
  nn_plot_impl(nn_mutual, pheno_data1, pheno_data2, title,
               background, xlim, ylim, line_color, scale_color)
}

# Shared implementation for nearest-neighbor plots
nn_plot_impl <- function(nn_dist, pheno_data1, pheno_data2, title,
                         background, xlim, ylim,
                         lineColor, scale_color) {
  p = nn_plot_base(title, background, xlim, ylim, scale_color)
  p = add_dist_data(p, nn_dist, lineColor)

  # Draw non-paired cells fainter than paired cells
  alphas = c(0.3, 1)
  from_alpha = alphas[(pheno_data1$data$`Cell ID` %in% nn_dist$`Cell ID`) + 1]
  p = p + ggplot2::geom_point(data=pheno_data1$data,
                     ggplot2::aes(x=`Cell X Position`, y=`Cell Y Position`),
                       color=pheno_data1$pheno$color, alpha=from_alpha)

  to_alpha = alphas[(pheno_data2$data$`Cell ID` %in% nn_dist$`To Cell ID`) + 1]
  p = p + ggplot2::geom_point(data=pheno_data2$data,
                     ggplot2::aes(`Cell X Position`, `Cell Y Position`),
                       color=pheno_data2$pheno$color, alpha=to_alpha)
  p
}

nn_plot_base = function(title, background, xlim, ylim, scale_color) {
  # Fake d.f needed to get background to draw...
  if (length(xlim)==1) xlim = c(0, xlim)
  if (length(ylim)==1) ylim = c(0, ylim)
  p = ggplot2::ggplot(data=data.frame(x=xlim, y=ylim), ggplot2::aes(x=x, y=y))
  p = p + ggplot2::labs(x='Cell X Position', y='Cell Y Position', title=title)
  add_scales_and_background(p, background, xlim, ylim, scale_color)
}

# Add line segments according to nn_dist
add_dist_data = function(p, nn_dist, lineColor) {
  if (nrow(nn_dist) > 0) {
    p = p + ggplot2::geom_segment(data = nn_dist,
                     ggplot2::aes(x=`Cell X Position`,
                                  y=`Cell Y Position`,
                                  xend=`To X Position`,
                                  yend=`To Y Position`),
                     color=lineColor)
  }
  p
}

# Add scales, scale line and background image to a ggplot object
add_scales_and_background = function(p, background, xlim, ylim,
                                     scale_color='black') {
  if (length(xlim)==1) xlim = c(0, xlim)
  if (length(ylim)==1) ylim = c(0, ylim)

  # Add scales at the image limits. Reverse the y scale to match the image
  p = p + ggplot2::scale_x_continuous(limits=xlim,
                                      expand=expansion()) +
    ggplot2::scale_y_reverse(limits=rev(ylim),
                             expand=expansion())

  # Force square aspect ratio
  p = p + ggplot2::coord_fixed()

  # Add background image if we have one
  if (length(background) > 1) {
    if (inherits(background, 'nativeRaster')) {
    p = p + annotation_raster_native(background,
                              xmin=xlim[1], xmax=xlim[2],
                              ymin=-ylim[1], ymax=-ylim[2])
    } else {
      p = p + ggplot2::annotation_raster(background,
                                         xmin=xlim[1], xmax=xlim[2],
                                         ymin=-ylim[1], ymax=-ylim[2])
    }
  }

  p = add_scale_line(p, xlim, ylim, scale_color, scale_alpha=1)
  p
}

# Add a 200-micron line segment for scale reference
add_scale_line <- function(p, xlim, ylim, scale_color, scale_alpha) {
  p = p + ggplot2::geom_segment(ggplot2::aes(x=xlim[2]-50-200, xend=xlim[2]-50,
                                             y=ylim[2]-100, yend=ylim[2]-100),
                                color=scale_color, alpha=scale_alpha, size=1)
  p = p + ggplot2::geom_text(ggplot2::aes(x=xlim[2]-50-200/2, y=ylim[2]-90,
                                          label=paste(200, '~mu*m')),
                             size=3, hjust=0.5, vjust=1,
                             color=scale_color, alpha=scale_alpha,
                             parse=TRUE)
  p
}

# A version of ggplot2::annotation_raster which takes a nativeRaster
# image. This is much faster than standard raster.
annotation_raster_native <- function(raster, xmin, xmax, ymin, ymax,
                                     interpolate = FALSE) {
  stopifnot(inherits(raster, 'nativeRaster'))
  ggplot2::layer(
    data = ggplot2:::dummy_data(),
    mapping = NULL,
    stat = ggplot2::StatIdentity,
    position = ggplot2::PositionIdentity,
    geom = ggplot2::GeomRasterAnn,
    inherit.aes = FALSE,
    params = list(
      raster = raster,
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      interpolate = interpolate
    )
  )
}

# ggplot2 3.3.0 deprecated expand_scale in favor of expansion.
# Using this function supports the old version and
# avoids deprecation warnings with the new version
expansion = function(mult=0, add=0) {
  if (function_exists('ggplot2', 'expansion'))
    ggplot2::expansion(mult, add)
  else
    ggplot2::expand_scale(mult, add)
}
