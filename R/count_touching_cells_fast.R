#' Count touching cells for a single pair of phenotypes
#'
#' Fast implementation of count touching cells.
#' This function requires a membrane label image such as
#' created by adaptive cell segmentation in inForm 2.4.3 or newer.
#' It will not work with
#' membrane masks created by older segmentation algorithms.
#'
#' @param csd Cell seg data, may include multiple fields
#' @param field_name Sample Name or Annotation ID to process
#' @param export_path Path to a directory containing composite
#'   and segmentation image files from inForm
#' @param phenos Named list of phenotype definitions. Must have length 2.
#' @param color1,color2 Colors for cells matching `phenos`.
#' @param discard_dups If `TRUE`, cells matching both phenotypes will be
#' discarded from the output. If `FALSE`, overlapping phenotypes is an error.
#' @return Returns a `list` containing two items:
#' \describe{
#'   \item{image}{An `EBImage::Image` of touching cells.}
#'   \item{data}{A \code{\link[tibble]{tibble}} containing
#'   data about the touching cells.}
#' }
#' @export
#' @md
count_touching_cells_fast = function(csd, field_name, export_path,
                                       phenos, color1, color2,
                                     discard_dups=FALSE) {
  stopifnot(is.list(phenos) && length(phenos) == 2)

  # Get the phenotype definitions
  phenos = validate_phenotypes(phenos, csd)
  pheno1 = phenos[[1]]
  pheno_name1 = names(phenos)[[1]]

  pheno2 = phenos[[2]]
  pheno_name2 = names(phenos)[[2]]

  # Make a fake cell_seg_path that will be used to find related images
  name = stringr::str_remove(field_name, '.im3')
  fake_cell_seg_path = file.path(export_path,
                                 paste0(name, '_cell_seg_data.txt'))

  # Subset csd to just the cells of interest and check that the phenotype
  # definitions don't overlap
  csd = csd[csd[[field_column(csd)]]==field_name, ]
  csd = force_pixel_locations(csd, fake_cell_seg_path)

  # First just the boolean selections
  selected_cells1 = select_rows(csd, pheno1)
  selected_cells2 = select_rows(csd, pheno2)
  selected_overlap = selected_cells1 & selected_cells2
  if (any(selected_overlap)) {
    if (discard_dups) {
      message('Discarding ', sum(selected_overlap),
              ' cells which match both phenotypes.')
      selected_cells1 = selected_cells1 & !selected_overlap
      selected_cells2 = selected_cells2 & !selected_overlap
    } else {
      stop('Phenotypes for count_touching_cells_fast ',
           'must define distinct cells. Found ', sum(selected_overlap),
           ' cells which match both phenotypes.')
    }
  }

  # Now the actual cells
  selected_cells1 = csd[selected_cells1, ]
  selected_cells2 = csd[selected_cells2, ]

  # Read membrane and nuclear masks
  masks = read_masks_label(fake_cell_seg_path)

  # Make label images for the membranes in each phenotype
  cell_image1 = masks$membrane
  cell_image1[!cell_image1 %in% selected_cells1$`Cell ID`] = 0

  cell_image2 = masks$membrane
  cell_image2[!cell_image2 %in% selected_cells2$`Cell ID`] = 0

  # Find the touching pairs
  touch_pairs = find_touching_cell_pairs(cell_image1, cell_image2, extra_size=0)

  # Need unique IDs for imaging
  touching_ids = list(unique(touch_pairs[, 1]), unique(touch_pairs[, 2]))

  # Make an image of the touching (and non-touching) cells
  colors = list(color1, color2) %>%
    rlang::set_names(c(pheno_name1, pheno_name2))
  composite_path = find_composite_image(fake_cell_seg_path)
  touching_image = make_touching_image_fast(pheno_name1, pheno_name2,
                                  cell_image1, cell_image2,
                                  touching_ids, masks,
                                  composite_path, colors)

  # Make a nice data frame of the touching pairs.
  # First get columns of interest from each selection and rename for clarity.
  field_col = field_column(selected_cells1)
  touching1 = selected_cells1 %>%
    # Select columns of interest
    dplyr::select(
      dplyr::contains('Slide ID'),
      !!field_col,
      `Cell ID`,
      dplyr::matches('Cell . Position'),
      dplyr::contains('Tissue Category')) %>%
    # Suffix column names with the phenotype
    dplyr::rename_at(.vars=dplyr::vars(-dplyr::contains('Slide ID'),
                                       -!!field_col),
                     .funs = ~paste0(.x, '.', pheno_name1))

  touching2 = selected_cells2 %>%
    dplyr::select(
      `Cell ID`,
      dplyr::matches('Cell . Position'),
      dplyr::contains('Tissue Category')) %>%
    dplyr::rename_all(.funs = ~paste0(.x, '.', pheno_name2))

  # Join tables using touch_pairs to correctly handle duplicates - one cell
  # may be touching multiple neighbors.
  id_names = paste0('Cell ID.', c(pheno_name1, pheno_name2))
  colnames(touch_pairs) = id_names
  join_table = touch_pairs %>% dplyr::as_tibble()
  data = touching1 %>%
    dplyr::right_join(join_table, by=id_names[[1]]) %>%
    dplyr::left_join(touching2, by=id_names[[2]])

  list(image=touching_image, data=data)
}


# Read the membrane and nuclear masks as label images.
# Returns a list with `nuclei` and `membrane`  members.
read_masks_label = function(cell_seg_path) {
  mask_path = sub('cell_seg_data.txt', 'binary_seg_maps.tif', cell_seg_path)
  if (!file.exists(mask_path))
    stop('count_touching_cells_fast requires a membrane label image ',
         'from inForm 2.4.3 or newer.')

  masks = read_maps(mask_path)
  if (!all(c('Membrane', 'Nucleus') %in% names(masks)))
    stop('binary_seg_maps file must contain ',
         'nuclear and membrane segmentation.')

  membrane = masks[['Membrane']]

  # We require a membrane label image; this will have a max > 0
  if(max(membrane) == 1)
    stop('count_touching_cells_fast requires a membrane label image.')

  nuclei = masks[['Nucleus']]
  membrane = t(membrane)
  nuclei = t(nuclei)
  list(nuclei=nuclei, membrane=membrane)
}

# Make a pretty picture showing the touch points.
make_touching_image_fast <- function(pheno_name1, pheno_name2,
                                cell_image1, cell_image2,
                                touching_ids, masks,
                                composite_path, colors) {
  composite = EBImage::readImage(composite_path, all=FALSE)

  # Fill touching cells
  i1_membrane = cell_image1
  i1_membrane[!i1_membrane %in% touching_ids[[1]]] = 0
  i1_membrane_filled = EBImage::fillHull(i1_membrane)
  composite = EBImage::paintObjects(i1_membrane_filled, composite,
                       col=c(colors[[pheno_name1]], colors[[pheno_name1]]))

  i2_membrane = cell_image2
  i2_membrane[!i2_membrane %in% touching_ids[[2]]] = 0
  i2_membrane_filled = EBImage::fillHull(i2_membrane)
  composite = EBImage::paintObjects(i2_membrane_filled, composite,
                       col=c(colors[[pheno_name2]], colors[[pheno_name2]]))

  # fillHull will not fill open objects, i.e. membrane of cells at the edge
  # of the image. These may be touching cells so fill the nuclei as well,
  # ensuring that these cells are marked in some way.
  i1_nuclei = masks$nuclei
  i1_nuclei[!i1_nuclei %in% touching_ids[[1]]] = 0
  composite = EBImage::paintObjects(i1_nuclei, composite, thick=TRUE,
                       col=c(colors[[pheno_name1]], colors[[pheno_name1]]))

  i2_nuclei = masks$nuclei
  i2_nuclei[!i2_nuclei %in% touching_ids[[2]]] = 0
  composite = EBImage::paintObjects(i2_nuclei, composite, thick=TRUE,
                       col=c(colors[[pheno_name2]], colors[[pheno_name2]]))

  # Outline all cells of each type by painting the membranes
  composite = EBImage::paintObjects(cell_image1, composite,
                                    col=c(colors[[pheno_name1]], NA))
  composite = EBImage::paintObjects(cell_image2, composite,
                                    col=c(colors[[pheno_name2]], NA))

  # Outline touching cells in white so cell boundaries are visible
  composite = EBImage::paintObjects(i1_membrane, composite,
                                    col=c('white', NA))
  composite = EBImage::paintObjects(i2_membrane, composite,
                                    col=c('white', NA))

  names(composite) = NULL # No need for names and they confuse tests
  composite
}
