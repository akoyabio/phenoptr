# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c(
    "Cell X Position", "Cell Y Position",
    "To X Position", "To Y Position",
    "x", "y"
    ))

#' Create a spatial distribution report.
#'
#' Creates an HTML report showing the location and nearest-neighbor relations
#' of cells in a single field and two chosen phenotypes.
#'
#' A cell seg data file for the field is required. If a tissue segmentation
#' or composite image of the field
#' is available, it will be used as a background.
#'
#' The report is saved to the same directory as the input file.
#'
#' @param cell_seg_path Path to a cell seg data file.
#' @param phenotypes A character vector containing two phenotype names.
#' @param colors A character vector of two color specifications.
#' @param output_path Optional, path to the output HTML file. If omitted,
#' output will be written to the directory containing `cell_seg_path`.
#' @param pixels_per_micron Conversion factor to microns
#'        (default 2 pixels/micron),
#' @param plot_gx_max Optional, maximum value of radius to show on G
#' function plots.
#' @param micron_line Optional, X location (radius in microns)
#'  of crosshair on G function plots.
#' @export
#' @examples
#' \dontrun{
#' cell_seg_path = system.file("extdata", "TMA",
#'                        "Core[1,5,6,1]_[21302,15107]_cell_seg_data.txt",
#'                        package = "informr")
#'
#' phenotypes = c("macrophage CD68", "cytotoxic CD8")
#' colors = c('red', 'blue')
#' out_path = path.expand('~/spatial_distribution_report.html')
#'
#' spatial_distribution_report(cell_seg_path, phenotypes, colors, out_path)
#' }
#' @md
spatial_distribution_report = function(cell_seg_path, phenotypes, colors,
                        output_path=NULL,
                        pixels_per_micron=getOption('informr.pixels.per.micron'),
                        plot_gx_max = 100,
                        micron_line = 25) {
  stopifnot(grepl('_cell_seg_data.txt', cell_seg_path))
  stopifnot(file.exists(cell_seg_path),
            length(phenotypes)==2,
            length(colors)==2)

  if (is.null(output_path))
    output_path = sub('_cell_seg_data.txt',
                      '_spatial_distribution.html', cell_seg_path)
  rmd_path = system.file("rmd", "SpatialDistributionReport.Rmd",
                         package="informr")

  rmarkdown::render(rmd_path, 'html_document', output_path)
}

## Helpers
## Note these are not exported and the doc comments are not roxygen
### Spatial functions ###

# Get a point pattern for a single phenotype
# @param data Cell seg data
# @param pheno List of (name, color)
# @param window spatstat::owin containing the points
# @return List of (data subset, ppp, pheno)
phenotype_as_ppp = function(data, pheno, window)
{
  stopifnot(pheno$name %in% unique(data$Phenotype))

  # Keep just the rows corresponding to the desired expression and phenotypes
  # for our data and columns containing fields of interest
  fixedColumns = c('Cell ID', 'Cell X Position', 'Cell Y Position', 'Phenotype')
  data = droplevels(data[data$Phenotype==pheno$name,fixedColumns])

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
find_nearest_neighbor = function(pheno_data1, pheno_data2)
{
  d = spatstat::nncross(pheno_data1$pp, pheno_data2$pp)
  names(d)[1] = 'Distance'
  from_data = pheno_data1$data %>%
    dplyr::select(`Cell X Position`, `Cell Y Position`)
  to_data = pheno_data2$data[d$which,] %>%
    dplyr::select('To X Position'=`Cell X Position`,
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
nn_plot <- function (pheno_data1, pheno_data2, nn_dist,
                     background, xlim, ylim,
                     lineColor='gray40') {
  title = bquote('Nearest neighbor from'
                 ~ italic(.(pheno_data1$pheno$name))
                 ~ 'to'
                 ~ italic(.(pheno_data2$pheno$name)))
  nn_plot_impl(nn_dist, pheno_data1, pheno_data2, title,
               background, xlim, ylim,
               lineColor)
}

# Shared implementation for nearest-neighbor plots
nn_plot_impl <- function (nn_dist, pheno_data1, pheno_data2, title,
                          background, xlim, ylim,
                          lineColor) {
  p = nn_plot_base(title, background, xlim, ylim)
  p = add_dist_data(p, nn_dist, lineColor)

  p = p + geom_point(data=pheno_data1$data,
                     aes(x=`Cell X Position`, y=`Cell Y Position`),
                     color=pheno_data1$pheno$color)
  p = p + geom_point(data=pheno_data2$data,
                     aes(`Cell X Position`, `Cell Y Position`),
                     color=pheno_data2$pheno$color)
  p
}

#' @import ggplot2
nn_plot_base = function(title, background, xlim, ylim) {
  # Fake d.f needed to get background to draw...
  p = ggplot(data=data.frame(x=0, y=0), aes(x=x, y=y))
  p = p + labs(x='Cell X Position', y='Cell Y Position', title=title)
  add_scales_and_background(p, background, xlim, ylim)
}

# Add line segments according to nn_dist
#' @import ggplot2
add_dist_data = function(p, nn_dist, lineColor) {
  if (nrow(nn_dist) > 0) {
    p = p + geom_segment(data = nn_dist,
                         aes(x=`Cell X Position`,
                             y=`Cell Y Position`,
                             xend=`To X Position`,
                             yend=`To Y Position`),
                         color=lineColor)
  }
  p
}

# Add scales, scale line and background image to a ggplot object
#' @import ggplot2
add_scales_and_background = function(p, background, xlim, ylim) {
  # Add scales at the image limits. Reverse the y scale to match the image
  p = p + scale_x_continuous(limits=c(0, xlim)) +
    scale_y_reverse(limits=c(ylim, 0))

  # Force square aspect ratio
  p = p + coord_fixed()

  # Add background image if we have one
  if (length(background) > 1)
  {
    p = p + annotation_raster(background, xmin=0, xmax=xlim, ymin=0, ymax=-ylim)
  }

  # Add a 200-micron line segment for scale reference
  p = p + geom_segment(aes(x=xlim-50-200, xend=xlim-50,
                           y=ylim-100, yend=ylim-100),
                       color='black', size=1)
  p = p + geom_text(aes(x=xlim-50-200/2, y=ylim-90, label=paste(200, '~mu*m')),
                    size=3, hjust=0.5, vjust=1, color='black', parse=TRUE)
  p
}
