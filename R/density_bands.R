# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c("count", "counts", "mids"))

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
#' @param map_path Path to the segmentation map file. If NULL, look for the
#'  map in the same directory as `cell_seg_path`.
#' @param component_path Path to the component_data file corresponding to
#' `cell_seg_path`; if omitted, look for the component data file in the
#' same directory as `cell_seg_path`.
#' @return Returns a `list` with three items:
#' \describe{
#'  \item{`densities`}{A `tibble` with five columns (see below).}
#'  \item{`cells`}{Cell seg data with phenotypes updated per the `phenotypes`
#'  parameter and an additional `distance` column.}
#'  \item{`distance`}{The distance map, a pixel image
#'      (\code{\link[spatstat.geom]{im.object}}).}
#'  }
#'
#' The `densities` item contains five columns:
#' \describe{
#'   \item{`phenotype`}{The supplied phenotypes.}
#'   \item{`midpoint`}{The midpoint of the distance band.}
#'   \item{`count`}{The number of cells of the phenotype found
#'      within the band.}
#'   \item{`area`}{The area of the band, in square microns.}
#'   \item{`density`}{The density of cells of the phenotype in the band,
#'      in cells per square micron.}
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
density_bands = function(cell_seg_path, phenotypes, positive, negative,
   width=25, map_path=NULL, component_path=NULL)
{
  if (!file.exists(cell_seg_path))
    stop(paste('File not found:', cell_seg_path))

  if (is.null(map_path))
    map_path = sub('_cell_seg_data.txt', '_binary_seg_maps.tif', cell_seg_path)
  if (!file.exists(map_path))
    stop(paste('File not found:', map_path))

  csd = read_cell_seg_data(cell_seg_path, pixels_per_micron='auto')

  # Get field metadata
  if (is.null(component_path)) component_path = cell_seg_path
  field_info = get_field_info(component_path)
  if(is.null(field_info))
    stop('density_bands() requires a matching component data or ',
    'segmentation map file to access the spatial reference for the field.')

  pixels_per_micron = 1/field_info$microns_per_pixel

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
  csd = make_phenotype_column(csd, phenotypes)

  # Read the mask and create separate masks for positive and negative regions
  maps = read_maps(map_path)
  if (!'Tissue' %in% names(maps))
    stop('No tissue segmentation in segmentation map file.')
  tissue = maps[['Tissue']]
  layers = parse_tissue_description(tissue)
  stopifnot(positive %in% names(layers), negative %in% names(layers))

  pos_mask = tissue==layers[positive]
  neg_mask = tissue==layers[negative]

  xrange=c(0, field_info$field_size[1])
  yrange=c(0, field_info$field_size[2])

  # Make distance maps for distance from positive and negative
  pos_win = spatstat.geom::owin(mask=pos_mask, xrange=xrange, yrange=yrange)
  dist_from_pos = spatstat.geom::distmap(pos_win)

  neg_win = spatstat.geom::owin(mask=neg_mask, xrange=xrange, yrange=yrange)
  dist_from_neg = spatstat.geom::distmap(neg_win)

  # Positive distance is into the positive mask == away from negative
  distance = dist_from_neg - dist_from_pos

  # For each cell, find its coordinates in the distance matrix
  # Mask index is [row, col] so [y, x]!
  ix = csd %>% dplyr::select(`Cell Y Position`, `Cell X Position`)

  # Adjust ix to have the origin at top left if needed
  if (max(ix$`Cell X Position`) > field_info$field_size[1]
          || max(ix$`Cell Y Position`) > field_info$field_size[2]) {
    ix = ix %>% dplyr::mutate(
      `Cell Y Position` = `Cell Y Position`-field_info$location[2],
      `Cell X Position` = `Cell X Position`-field_info$location[1]
    )
  }

  # Use the nearest point in the distance matrix
  ix = ix %>% as.matrix()
  ix = ix * pixels_per_micron %>% round() %>% as.integer()

  # For each cell, look up its distance from the boundary
  csd$distance = distance$v[ix]

  # Compute cut-points for distance
  min_cut = floor(min(distance$v)/width)*width
  max_cut = ceiling(max(distance$v)/width) * width
  cut_points = seq(min_cut, max_cut, width)

  # Compute areas by counting pixels in the distance matrix
  area = graphics::hist(distance$v, breaks=cut_points, plot=FALSE)

  # Normalize by pixel size ^ 2
  areas = tibble::as_tibble(area[c('mids', 'counts')]) %>%
    dplyr::mutate(area = counts / pixels_per_micron^2) %>%
    dplyr::select(-counts)

  # Count cells for each phenotype separately
  cell_counts = csd %>%
    dplyr::group_by(Phenotype) %>%
    dplyr::summarize(count =
            list(graphics::hist(distance, breaks=cut_points, plot=FALSE))) %>%
    dplyr::mutate(count = purrr::map(count,
                              ~tibble::as_tibble(.x[c('mids', 'counts')]))) %>%
    tidyr::unnest(cols=c(count)) %>%
    dplyr::rename(count=counts, phenotype=Phenotype)

  densities = cell_counts %>%
    dplyr::inner_join(areas, by='mids') %>%
    dplyr::mutate(density=count/area) %>%
    dplyr::rename(midpoint=mids)

  list(densities=densities, cells=csd, distance=distance)
}

#' Compute density of each phenotype in each given region
#' @param csd Cell seg table
#' @param regions An `sf::st_sf` or `sf::st_sfc` object containing the regions
#' of interest, for example the output of `create_buffers()`.
#' @param phenotypes Phenotypes of interest; if omitted, will use
#' `phenoptr::unique_phenotypes(csd)`.
#' @return A data frame, `regions` annotated with an area column
#' (in square microns) and two columns per phenotype
#' giving the count and density of that phenotype in that region.
#' @export
density_by_region = function(csd, regions, phenotypes=NULL) {
  stopifnot(!is.null(csd))

  if (is.null(phenotypes))
    phenotypes = phenoptr::unique_phenotypes(csd)
  if (!rlang::is_named(phenotypes))
    phenotypes = rlang::set_names(phenotypes)

  # Make an st_sf (data frame) if we don't already have one
  if (inherits(regions, 'sfc'))
    regions = sf::st_sf(regions)
  if (!inherits(regions, 'sf'))
    stop('regions must be an sf::st_sf or sf::st_sfc object.')

  # Compute areas once
  regions$area = purrr::map_dbl(sf::st_geometry(regions), sf::st_area)

  densities = purrr::imap_dfc(phenotypes, function(phenotype, name) {
    # Which cells are in the target phenotype?
    pheno_cells = csd[phenoptr::select_rows(csd, phenotype),
                      c('Cell X Position', 'Cell Y Position')]
    pheno_pts = sf::st_as_sf(pheno_cells, coords=1:2)
    cells_in_region = sf::st_contains(regions, pheno_pts) %>%
      purrr::map_int(length)
    tibble::tibble(count=cells_in_region,
                   density=cells_in_region/regions$area) %>%
      rlang::set_names(paste0(name, ' count'), paste0(name, ' density'))
  })

  dplyr::bind_cols(regions, densities)
}

#' Create buffer regions around a boundary and within a region of interest
#'
#' For use with `density_by_region()`.
#' @param boundary Polygon for boundary line. Must be at least partially
#' within `roi`.
#' @param roi Polygon for ROI
#' @param n Number of buffers (on each side of boundary)
#' @param width Buffer width
#' @return A list with two items
#' - buffers - a data frame with index from -n to n and geometry column
#' - divider - the dividing line for `boundary` within `roi`
#' @importFrom magrittr %>%
#' @export
create_buffer_bands = function(boundary, roi, n, width) {
  # Get just the dividing line
  divider = boundary %>%
    sf::st_cast('LINESTRING') %>% # So the intersection is a line, not a polygon
    sf::st_intersection(roi) %>%
    sf::st_cast('LINESTRING')

  if (sum(sf::st_length(divider)) == 0)
    stop('boundary and roi do not intersect.')

  # divider may have two segments, depending on where it broke into
  # a linestring. If so, bind them back together.
  if (length(sf::st_length(divider)) > 1) {
    pts = unclass(sf::st_geometry(divider))
    divider = sf::st_linestring(do.call(rbind, rev(pts))) %>% sf::st_sfc()
  }

  # Make buffers
  # First make inclusive buffers that span both sides of the divider
  # and subsets for inside and outside boundary
  fat_buffers = list()
  inside_buffers = list()
  outside_buffers = list()
  seed = divider

  1:n %>% purrr::walk(function(i) {
    seed <<- sf::st_buffer(seed, width) %>% sf::st_intersection(roi)
    fat_buffers <<- c(fat_buffers, list(seed))
    inside_buffers <<- c(inside_buffers, list(sf::st_intersection(seed, boundary)))
    outside_buffers <<- c(outside_buffers, list(sf::st_difference(seed, boundary)))
  })

  # If the next-to-last fat buffer is the same as the roi, the buffers are
  # are too many / too large
  if (length(sf::st_difference(roi, fat_buffers[[n-1]])) == 0)
    warning('Largest buffers exceed the roi, use smaller n or width.')

  # Now sort it out to n distinct buffers on each side of the divide
  buffers = purrr::map_dfr(1:n, function(i) {
    # Don't use ifelse here, it loses the class attribute
    if (i==1) {
      inside = inside_buffers[[1]]
      outside = outside_buffers[[1]]
    } else {
      inside = sf::st_difference(inside_buffers[[i]], inside_buffers[[i-1]])
      outside = sf::st_difference(outside_buffers[[i]], outside_buffers[[i-1]])
    }

    inside = if (length(inside) > 0) inside else NA
    outside = if (length(outside) > 0) outside else NA
    tibble::tibble(
      inner=c(-i+1, i-1) * width,
      outer=c(-i, i) * width,
      poly=list(inside, outside))
  }) %>%
    dplyr::filter(!is.na(.data$poly)) %>%
    dplyr::mutate(poly=sf::st_sfc(.data$poly))

  buffers=sf::st_sf(buffers, sf_column_name='poly') %>%
    dplyr::arrange(.data$inner)
  list(divider=divider, buffers=buffers)
}

#' Plot density bands with divider and ROI
#' @param buffers The band buffers
#' @param divider The tumor boundary line
#' @param roi The region of interest boundary
#' @return A `ggplot` object
#' @export
plot_buffers = function(buffers, divider, roi)
{
  ggplot2::ggplot() +
    geom_sf_invert(data=buffers,
                           color=scales::alpha('lightblue'), fill=NA) +
    geom_sf_invert(data=divider, color='red', fill=NA) +
    geom_sf_invert(data=roi, color='blue', fill=NA) +
    scale_sf_invert() +
    ggplot2::labs(title='Area of analysis showing distance bands') +
    ggplot2::theme_minimal()
}

#' Plot density by band
#' @param densities Result of calling `density_by_region()`
#' @param colors Named vector of phenotype colors
#' @return A `ggplot` object
#' @export
plot_density_by_band = function(densities, colors) {
  # To plot densities, make a tall data frame
  tall = densities %>% sf::st_drop_geometry() %>%
    dplyr::select(outer, dplyr::ends_with('density')) %>%
    tidyr::gather('Phenotype', 'Density', -outer) %>%
    dplyr::mutate(Phenotype = stringr::str_remove(Phenotype, ' density'))

  # Our density values are cells/micron^2 so multiply by 1e6
  ggplot2::ggplot(tall, ggplot2::aes(outer, .data$Density*1e6, color=Phenotype)) +
    ggplot2::geom_vline(xintercept=0, linetype=2, color='gray') +
    ggplot2::geom_line(size=2) +
    ggplot2::annotate('text', x=20, y=800, angle=-90, vjust=0,
                      label='Tumor margin') +
    ggplot2::scale_color_manual(values=colors) +
    ggplot2::scale_y_continuous(trans='sqrt') +
    ggplot2::labs(x='Distance from tumor margin (negative is inside tumor)',
                  y=expression(paste('Cell Density (', cells/mm^2, ') (non-linear scale)')),
                  title='Cell density in bands from tumor margin') +
    ggplot2::theme_minimal()
}
