# Suppress CMD CHECK notes for things that look like global vars
if (getRversion() >= "2.15.1")
  utils::globalVariables(c("points"))

#' Estimate cell density at distance from a tissue boundary.
#'
#' Given a cell seg table and an image containing masks for two tissue
#' classes, estimate the dependence of the density of cells of each
#' phenotype on the distance from the boundary between the two tissue classes.
#'
#' The cell density estimate is computed by \code{\link[spatstat]{rhohat}}.
#' The signed distance from the boundary between the two tissue classes
#' is used as the covariate and density is estimated separately for each
#' phenotype. `rhohat` uses kernel density estimation to
#' estimate the dependence of cell count on distance and the dependence
#' of area on distance. The ratio of these two estimates is the estimate
#' of density vs distance.
#'
#' The `rhohat` element of the returned list is a
#' `list` containing the results of the cell density
#' estimation for each phenotype. Each list value is a `rhohat` object,
#' see \code{\link[spatstat]{methods.rhohat}}.
#'
#' Density estimates are
#' in cells per square micron; multiply by 1,000,000 for cells per square
#' millimeter.
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
#' For both `mask==TRUE` and `mask==FALSE`, the results tend to be unreliable
#' near the distance extremes because the tissue area at that distance
#' will be relatively small so random variation in cell counts is magnified.
#'
#' @section Parallel computation:
#'
#' `density_at_distance` supports parallel computation using the
#' \code{\link[foreach]{foreach-package}}. When this is enabled, the densities
#' for each phenotype will be computed in parallel. To use this feature,
#' you must enable a parallel backend, for example using
#' \code{\link[doParallel]{registerDoParallel}}.
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
#' @param mask If true, compensate for edge effects by masking the image
#' to exclude cells which are closer to the edge of the image than to the
#' tissue boundary.
#' @param pixels_per_micron Conversion factor to microns.
#' @param ... Additional arguments passed to \code{\link[spatstat]{rhohat}}.
#' Default parameters are `method="ratio", smoother="kernel", bw="nrd"`.
#' @return Returns a `list` containing four items:
#'   \tabular{ll}{
#'    `points` \tab The points used, marked with their phenotype,
#'    a \code{\link[spatstat]{ppp.object}}.\cr
#'    `rhohat` \tab The density estimates (see Details).\cr
#'    `distance` \tab The distance map, a pixel image
#'      (\code{\link[spatstat]{im.object}}).\cr
#'    `mask` \tab A mask matrix showing the locations which are closer to
#'    the tissue boundary than to the border or other regions.\cr
#'  }
#' @references A. Baddeley, E. Rubak and R.Turner.
#' Spatial Point Patterns: Methodology and Applications with R.
#' Chapman and Hall/CRC Press, 2015. Sections 6.6.3-6.6.4.
#' @examples
#' # Compute density for the sample data
#' values <- density_at_distance(sample_cell_seg_path(),
#'   list("CD8+", "CD68+", "FoxP3+"),
#'   positive="Stroma", negative="Tumor")
#'
#' # Combine all the densities into a single data_frame
#' # and filter out the extremes
#' library(dplyr)
#' all_rho <- purrr::map_dfr(values$rhohat, ~., .id='phenotype') %>%
#'   as_data_frame %>% filter(X>=-50, X<=100)
#'
#' # Plot the densities in a single plot
#' library(ggplot2)
#' ggplot(all_rho, aes(X, rho*1000000, color=phenotype)) +
#'   geom_line(size=2) +
#'   labs(x='Distance from tumor boundary (microns)',
#'        y='Estimated cell density (cells per sq mm)')
#'
#' # Show the distance map with CD68+ cells superimposed
#' plot_diverging(values$distance, show_boundary=TRUE,
#'   title=paste('Distance from tumor for CD68+cells'),
#'   sub='Positive (blue) distances are away from tumor')
#' plot(values$points[values$points$marks=='CD68+', drop=TRUE],
#'   add=TRUE, use.marks=FALSE, cex=0.5, col=rgb(0,0,0,0.5))
#' @family density estimation
#' @export
#' @md
#' @importFrom magrittr "%>%"
#' @importFrom foreach "%dopar%"
density_at_distance = function(cell_seg_path, phenotypes, positive, negative,
   mask=FALSE, pixels_per_micron=getOption('phenoptr.pixels.per.micron'),
   ...)
{
  if (!file.exists(cell_seg_path))
    stop(paste('File not found:', cell_seg_path))

  map_path = sub('_cell_seg_data.txt', '_binary_seg_maps.tif', cell_seg_path)
  if (!file.exists(map_path))
    stop(paste('File not found:', map_path))

  csd = read_cell_seg_data(cell_seg_path, pixels_per_micron)
  if (!'Phenotype' %in% names(csd))
    stop('Cell seg data does not contain a Phenotype column.')

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
    if (is.null(names(phenotypes)))
      stop('phenotypes parameter must be a named list.')
  }

  # Mutate csd to have the desired phenotypes
  csd = purrr::map(names(phenotypes),
      ~(csd[select_rows(csd, phenotypes[[.x]]),] %>%
        dplyr::mutate(Phenotype=.x))) %>%
    dplyr::bind_rows()

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

  # Make the mask
  # We want a distance map from anything that is not positive or negative
  boundary_mask = (!pos_mask) & (!neg_mask)

  # Put a 1-pixel border around boundary_mask so we get distance from the edge
  boundary_mask[c(1, nrow(boundary_mask)),]= TRUE
  boundary_mask[,c(1, ncol(boundary_mask))]= TRUE

  boundary_win =
    spatstat::owin(mask=boundary_mask, xrange=xrange, yrange=yrange)
  dist_from_boundary = spatstat::distmap(boundary_win)

  valid_mask =
    dist_from_pos<dist_from_boundary & dist_from_neg<dist_from_boundary

  if (mask) {
    # Valid region is anything closer to positive or negative than to boundary
    valid_win = spatstat::owin(mask=as.matrix(valid_mask),
                               xrange=xrange, yrange=yrange)
  } else {
    valid_win = spatstat::owin(xrange=xrange, yrange=yrange)
  }

  pp = spatstat::ppp(csd$`Cell X Position`, csd$`Cell Y Position`,
                    window=valid_win, marks=factor(csd$Phenotype))

  # Marshal args for rhohat
  params = list(...)
  params$covariate = distance

  # Add defaults
  if (!'method' %in% names(params)) params$method = 'ratio'
  if (!'smoother' %in% names(params)) params$smoother = 'kernel'
  if (!'bw' %in% names(params)) params$bw = 'nrd'

  # Distance estimation for each phenotype separately
  split_pp = split(pp)
  all_rho = foreach::foreach(points=split_pp, .packages='spatstat') %dopar% {
    local_params = params
    local_params$object = points
    do.call(spatstat::rhohat, local_params)
   }
  names(all_rho) = names(split_pp)

  list(points=pp, rhohat=all_rho, distance=distance, mask=valid_mask)
}

#' Plot a distance map using a diverging color scale with white at 0
#'
#' @param im A pixel image of class \code{\link[spatstat]{im}}.
#' @param title Title for the plot.
#' @param show_boundary Should the boundary be highlighted?
#' @param ... Additional arguments passed to \code{\link[spatstat]{plot.im}}
#' @md
#' @export
plot_diverging = function(im, title, show_boundary=FALSE, ...) {
  if (missing(title))
    title <- deparse(substitute(im))

  # create color ramp from blue to red
  color_levels=101 # the number of colors to use
  colors =
    grDevices::colorRampPalette(c('blue', 'gray96', 'red'))(n=color_levels)
  if (show_boundary) colors[(color_levels+1)/2] = 'gray90'
  max_absolute_value=max(abs(im)) # what is the maximum absolute value of im?
  color_sequence =
    seq(-max_absolute_value,max_absolute_value,length.out=color_levels+1)
  graphics::plot(im, main=title, col=colors, breaks=color_sequence, ...)
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

  purrr::set_names(seq_along(layers)-1, layers)
}
