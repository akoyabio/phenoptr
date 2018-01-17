#' Estimate cell density at distance from a tissue boundary.
#'
#' Given a cell seg table and an image containing masks for two tissue
#' classes, estimate the dependence of the density of cells of each
#' phenotype on the distance from the boundary between the two tissue classes.
#'
#' The density estimate is computed by \code{\link[spatstat]{rhohat}}.
#' The signed distance from the boundary between the two tissue classes
#' is used as the covariate and density is estimated separately for each
#' phenotype.
#'
#' The `rho_hat` element of the returned list is a
#' \code{\link[tibble]{data_frame}} containing the results of the density
#' estimation for each phenotype. It has six columns:
#'   \tabular{ll}{
#'    `phenotype` \tab The phenotype name.\cr
#'    `distance` \tab Distance from tissue boundary (positive or negative),
#'    i.e. the covariate for density estimation.\cr
#'    `rho` \tab The density estimate from \code{\link[spatstat]{rhohat}},
#'      in cells per square micron.\cr
#'    `var`, `hi`, `lo` \tab The variance and high and low estimates of `rho`.
#' }
#'
#' @section Edge correction:
#'
#' The \code{\link[spatstat]{rhohat}} function does not have any built-in
#' edge correction. This may lead to incorrect density estimates
#' because it does not account for cells at the edge of the image which may
#' be near a tissue boundary which is not part of the image.
#'
#' The `mask` parameter, if set to `TRUE`, restricts the density estimation
#' to cells which are closer to the tissue boundary in the image than they are
#' to the edge of the image or any other tissue class.
#'
#' @param cell_seg_path Path to cell segmentation data.
#' @param phenotypes Optional named list of phenotypes to process.
#'   \code{names(phenotypes)} are the names of the resulting phenotypes.
#'   The values are in any format accepted by \code{\link{select_rows}}.
#'   If omitted, will use all phenotypes in the cell seg data.
#' @param positive Name of the tissue category used as positive distance,
#' e.g. "stroma".
#' @param negative Name of the tissue category used as negative distance,
#' e.g. "tumor".
#' @param mask If true, compensate for edge effects by masking the image
#' to exclude cells which are closer to the edge of the image than to the
#' tissue boundary.
#' @param pixels_per_micron Conversion factor to microns.
#' @param smoother,method,... Additional arguments passed to
#' \code{\link[spatstat]{rhohat}}.
#' @return Returns a `list` containing two or three items:
#'   \tabular{ll}{
#'    `rho_hat` \tab The density estimate (see Details).\cr
#'    `distance` \tab The distance map, a pixel image of class
#'      \code{\link[spatstat]{im}}.\cr
#'    `mask` \tab If `mask` is `TRUE`,
#'    a matrix showing the locations used.\cr
#'  }
#' @references A. Baddeley, E. Rubak and R.Turner.
#' Spatial Point Patterns: Methodology and Applications with R.
#' Chapman and Hall/CRC Press, 2015. Section 6.6.3.
#' @examples
#' \dontrun{
#' # Compute density for the sample data
#' values <- density_at_distance(sample_cell_seg_path(),
#'   list("CD8+", "CD68+", "FoxP3+"),
#'   positive="Stroma", negative="Tumor")
#'
#' # Plot the densities
#' library(ggplot2)
#' ggplot(values[['rho_hat']], aes(distance, rho, color=phenotype)) +
#'   geom_line()
#'
#' # Show the distance map with the mask superimposed and CD8+ cells drawn
#' plot_diverging(values$distance, title='Distance')
#' plot(values$mask, add=TRUE, col=c('black', 'transparent'))
#' cells = subset(sample_cell_seg_data, Phenotype=='CD8+')
#' points(cells$`Cell X Position`, cells$`Cell Y Position`)
#' }
#' @export
#' @md
#' @importFrom magrittr "%>%"
density_at_distance = function(cell_seg_path, phenotypes, positive, negative,
   mask=TRUE, pixels_per_micron=getOption('phenoptr.pixels.per.micron'),
   smoother='local', method='reweight', ...)
{
  map_path = sub('_cell_seg_data.txt', '_binary_seg_maps.tif', cell_seg_path)
  stopifnot(file.exists(cell_seg_path), file.exists(map_path))

  csd = read_cell_seg_data(cell_seg_path, pixels_per_micron)
  stopifnot('Phenotype' %in% names(csd))

  # Check for multiple samples, this is probably an error
  if (length(unique(csd$`Sample Name`))>1)
    stop('Data appears to contain multiple samples.')

  if (is.null(phenotypes)) {
    # Use phenotypes from file
    phenotypes = sort(unique(csd$Phenotype)) %>% purrr::set_names()
  } else {
    # Phenotypes were provided
    stopifnot(length(phenotypes) > 0)
    # Allow an unnamed list if they are all single character vectors
    # Otherwise it is an error
    if (all(purrr::map_lgl(phenotypes, ~is.character(.) && length(.)==1)))
      phenotypes = purrr::set_names(phenotypes)
    stopifnot(!is.null(names(phenotypes)))

    # Mutate csd to have the desired phenotypes
    csd = purrr::map(names(phenotypes),
        ~(csd[select_rows(csd, phenotypes[[.x]]),] %>%
          dplyr::mutate(Phenotype=.x))) %>%
      dplyr::bind_rows()
  }

  # Read the mask and create separate masks for positive and negative regions
  maps = read_maps(map_path)
  stopifnot('Tissue' %in% names(maps))
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

  if (mask) {
    # We want a distance map from anything that is not positive or negative
    boundary_mask = (!pos_mask) & (!neg_mask)

    # Put a 1-pixel border around boundary_mask so we get distance from the edge
    boundary_mask[c(1, nrow(boundary_mask)),]= TRUE
    boundary_mask[,c(1, ncol(boundary_mask))]= TRUE

    boundary_win =
      spatstat::owin(mask=boundary_mask, xrange=xrange, yrange=yrange)
    dist_from_boundary = spatstat::distmap(boundary_win)

    # Valid region is anything closer to positive or negative than to boundary
    valid_mask =
      dist_from_pos<dist_from_boundary & dist_from_neg<dist_from_boundary
    valid_win = spatstat::owin(mask=as.matrix(valid_mask),
                               xrange=xrange, yrange=yrange)
  } else {
    valid_win = spatstat::owin(xrange=xrange, yrange=yrange)
  }

  pp = spatstat::ppp(csd$`Cell X Position`, csd$`Cell Y Position`,
                    window=valid_win, marks=factor(csd$Phenotype))

  # Distance estimation for each phenotype separately
  all_rho = purrr::map_dfr(split(pp), function(points)
    spatstat::rhohat(points, covariate=distance,
                    smoother=smoother, method=method, ...),
    .id='phenotype')

  all_rho = tibble::as_tibble(all_rho)

  value = list(rho_hat=all_rho, distance=distance)
  if (mask)
    value[['mask']] = valid_mask
  value
}

#' Plot a distance map using a diverging color scale with white (grey) at 0
#'
#' @param im A pixel image of class \code{\link[spatstat]{im}}.
#' @param title Title for the plot.
#' @md
#' @export
plot_diverging = function(im, title) {
  if (missing(title))
    title <- deparse(substitute(im))

  # create color ramp from blue to red
  col5 <- grDevices::colorRampPalette(c('blue', 'gray96', 'red'))
  max_absolute_value=max(abs(im)) # what is the maximum absolute value of im?
  color_levels=100 # the number of colors to use
  color_sequence=seq(-max_absolute_value,max_absolute_value,length.out=color_levels+1)
  plot(im, main=title, col=col5(n=color_levels), breaks=color_sequence)
}

#' Parse the description tag of the tissue map image read from an
#' inForm `binary_seg_maps` file
#' @param img The `Tissue` image from \code{\link{read_maps}}.
#' @return A named vector whose names are the tissue classes in `img` and
#' whose values are the mask values for the class in `img`.
#' @md
#' @export
parse_tissue_description = function(img) {
  desc = attr(img, 'description')
  stopifnot(!is.null(desc))
  parsed = xml2::read_xml(desc)
  layers = purrr::map_chr(xml2::xml_find_all(parsed, 'Entry'),
                          function(entry) {
    entry %>% xml2::xml_find_first('Name') %>% xml2::xml_text()
  })

  setNames(1:length(layers)-1, layers)
}
