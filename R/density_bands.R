# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c("counts", "mids"))

#' Estimate cell density in bands from a tissue boundary.
#'
#' Given a cell seg table and an image containing masks for two tissue
#' classes, estimate the density of cells of each specified
#' phenotype in bands from the boundary between the two
#' tissue classes.
#'
#' `density_bands` uses a counting approach similar to a histogram.
#' First the image is divided into bands based on distance from the
#' specified boundary. Next, the number of cells of each phenotype
#' within each distance band is counted and the area of each band
#' is estimated. The density estimates are the ratios of the cell
#' counts to the area estimates.
#'
#' Density estimates are
#' in cells per square micron; multiply by 1,000,000 for cells per square
#' millimeter.
#'
#' The returned value includes the cell counts and area of each band,
#' making it straightforward to aggregate across multiple fields from a
#' single sample. The aggregate density is computed by summing the
#' cell counts and areas across all fields from a sample, then dividing
#' to compute density.
#'
#' @param cell_seg_path Path to a cell segmentation data file.
#' @param phenotypes Optional named list of phenotypes to process.
#'   \code{names(phenotypes)} are the names of the resulting phenotypes.
#'   The values are in any format accepted by \code{\link{select_rows}}.
#'   If omitted, will use all phenotypes in the cell seg data.
#' @param positive Name of the tissue category used as positive distance,
#' e.g. "stroma".
#' @param negative Name of the tissue category used as negative distance,
#' e.g. "tumor".
#' @param width Width of the bands, in microns
#' @param pixels_per_micron Conversion factor to microns.
#' @param map_path Path to the segmentation map file. If NULL, look for the
#'  map in the same directory as `cell_seg_path`.
#' @return Returns a `list` with three items:
#' \tabular{ll}{
#'   `densities` \tab A `data_frame` with five columns (see below).\cr
#'  `cells` \tab Cell seg data with phenotypes updated per the `phenotypes`
#'  parameter and an additional `distance` column.\cr
#'  `distance` \tab The distance map, a pixel image
#'      (\code{\link[spatstat]{im.object}}).\cr
#'  }
#'
#' The `densities` item contains five columns:
#' \tabular{ll}{
#'   `phenotype` \tab The supplied phenotypes.\cr
#'   `midpoint` \tab The midpoint of the distance band.\cr
#'   `count` \tab The number of cells of the phenotype found
#'      within the band.\cr
#'   `area` \tab The area of the band, in square microns.\cr
#'   `density` \tab The density of cells of the phenotype in the band,
#'      in cells per square micron.\cr
#'    }
#' @examples
#' # Compute density for the sample data
#' values <- density_bands(sample_cell_seg_path(),
#'   list("CD8+", "CD68+", "FoxP3+"),
#'   positive="Stroma", negative="Tumor")
#'
#' # Plot the densities in a single plot
#' library(ggplot2)
#' ggplot(values$densities, aes(midpoint, density*1000000, color=phenotype)) +
#'   geom_line(size=2) +
#'   labs(x='Distance from tumor boundary (microns)',
#'        y='Estimated cell density (cells per sq mm)')
#' @family density estimation
#' @export
#' @md
#' @importFrom magrittr "%>%"
#' @importFrom foreach "%dopar%"
#' @import dplyr
density_bands = function(cell_seg_path, phenotypes, positive, negative,
   width=25, pixels_per_micron=getOption('phenoptr.pixels.per.micron'),
   map_path=NULL)
{
  if (!file.exists(cell_seg_path))
    stop(paste('File not found:', cell_seg_path))

  if (is.null(map_path))
    map_path = sub('_cell_seg_data.txt', '_binary_seg_maps.tif', cell_seg_path)
  if (!file.exists(map_path))
    stop(paste('File not found:', map_path))

  csd = read_cell_seg_data(cell_seg_path, pixels_per_micron)

  # Check for multiple samples, this is probably an error
  if (length(unique(csd$`Sample Name`))>1)
    stop('Data appears to contain multiple samples.')

  if (is.null(phenotypes)) {
    # Use phenotypes from file
    phenotypes = unique_phenotypes(csd) %>% purrr::set_names()
  } else {
    # Phenotypes were provided
    stopifnot(length(phenotypes) > 0)
    # Allow an unnamed list if they are all single character vectors
    # Otherwise it is an error, user must supply names for compound phenotypes
    if (all(purrr::map_lgl(phenotypes, ~is.character(.) && length(.)==1)))
      phenotypes = purrr::set_names(phenotypes)
    if (is.null(names(phenotypes)))
      stop('phenotypes parameter must be a named list.')
  }

  # Mutate csd to have the desired phenotypes
  csd = purrr::map(names(phenotypes),
      ~(csd[select_rows(csd, phenotypes[[.x]]),] %>% mutate(Phenotype=.x))) %>%
    bind_rows()

  # Read the mask and create separate masks for positive and negative regions
  maps = read_maps(map_path)
  if (!'Tissue' %in% names(maps))
    stop('No tissue segmentation in segmentation map file.')
  tissue = maps[['Tissue']]
  layers = parse_tissue_description(tissue)
  stopifnot(positive %in% names(layers), negative %in% names(layers))

  pos_mask = tissue==layers[positive]
  neg_mask = tissue==layers[negative]

  xrange=c(0, dim(tissue)[2]/pixels_per_micron)
  yrange=c(0, dim(tissue)[1]/pixels_per_micron)

  # Make distance maps for distance from positive and negative
  pos_win = spatstat::owin(mask=pos_mask, xrange=xrange, yrange=yrange)
  dist_from_pos = spatstat::distmap(pos_win)

  neg_win = spatstat::owin(mask=neg_mask, xrange=xrange, yrange=yrange)
  dist_from_neg = spatstat::distmap(neg_win)

  # Positive distance is into the positive mask == away from negative
  distance = dist_from_neg - dist_from_pos

  # For each cell, look up its distance
  # Use the nearest point in the distance matrix
  # Index is [row, col] so [y, x]!
  ix = as.matrix(csd[,c('Cell Y Position', 'Cell X Position')])
  ix = ix * pixels_per_micron %>% round() %>% as.integer()
  csd$distance = distance$v[ix]

  # Compute cut-points for distance
  min_cut = floor(min(distance$v)/width)*width
  max_cut = ceiling(max(distance$v)/width) * width
  cut_points = seq(min_cut, max_cut, width)

  # Compute areas by counting pixels in the distance matrix
  area = graphics::hist(distance$v, breaks=cut_points, plot=FALSE)

  # Normalize by pixel size ^ 2
  areas = as_data_frame(area[c('mids', 'counts')]) %>%
    mutate(area = counts / pixels_per_micron^2) %>%
    select(-counts)

  # Count cells for each phenotype separately
  cell_counts = csd %>% group_by(Phenotype) %>%
    summarize(count =
            list(graphics::hist(distance, breaks=cut_points, plot=FALSE))) %>%
    mutate(count = purrr::map(count,
                              ~as_data_frame(.x[c('mids', 'counts')]))) %>%
    tidyr::unnest() %>%
    rename(count=counts, phenotype=Phenotype)

  densities = cell_counts %>% inner_join(areas, by='mids') %>%
    mutate(density=count/area) %>%
    rename(midpoint=mids)

  list(densities=densities, cells=csd, distance=distance)
}
