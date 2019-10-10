#' Count touching cells for a single pair of phenotypes
#'
#' @param csd Cell seg data
#' @param field_name Sample Name or Annotation ID to process
#' @param export_path Path to a directory containing composite and component
#'   image files from inForm
#' @param phenos Named list of phenotype definitions. Must have length 2.
#' @param color1,color2 Colors to draw the phenotype dots
#' @return Returns a `list` containing two items:
#' \describe{
#'   \item{`image`}{An `EBImage::Image` of touching cells.}
#'   \item{`data`}{A \code{\link[tibble]{tibble}} containing
#'   data about the touching cells.}
#' }
#' @export
count_touching_cells_single = function(csd, field_name, export_path,
                                       phenos, color1, color2) {
  stopifnot(is.list(phenos) && length(phenos) == 2)

  # Get the phenotype definitions
  phenos = validate_phenotypes(phenos, csd)
  pheno1 = phenos[[1]]
  pheno_name1 = names(phenos)[[1]]

  pheno2 = phenos[[2]]
  pheno_name2 = names(phenos)[[2]]

  # Make a fake cell_seg_path that will be used to find related images
  name = stringr::str_remove(field_name, '.im3')
  fake_cell_seg_path = file.path(export_path, paste0(name, '_cell_seg_data.txt'))

  # Subset csd to just the cells of interest and check that the phenotype
  # definitions don't overlap
  csd = csd[csd[[field_column(csd)]]==field_name,]
  csd = force_pixel_locations(csd, fake_cell_seg_path)

  # First just the boolean selections
  selected_cells1 = select_rows(csd, pheno1)
  selected_cells2 = select_rows(csd, pheno2)
  selected_overlap = any(selected_cells1 & selected_cells2)
  if (selected_overlap)
    stop('Phenotypes for count_touching_cells_single must define distinct cells.')

  # Now the actual cells
  selected_cells1 = csd[selected_cells1,]
  selected_cells2 = csd[selected_cells2,]

  # Read membrane and nuclear masks
  masks = read_masks(fake_cell_seg_path)

  # Make label images for the cells in each phenotype
  cell_image1 = make_cell_image(selected_cells1,
                                masks$nuclei, masks$membrane)
  cell_image2 = make_cell_image(selected_cells2,
                                masks$nuclei, masks$membrane)

  # Find the touching pairs
  extra_size = ifelse(masks$membrane_width==1, 0, 2)
  touch_pairs = find_touching_cell_pairs(cell_image1, cell_image2, extra_size)

  # Need unique IDs for imaging
  touching_ids = list(unique(touch_pairs[,1]), unique(touch_pairs[,2]))

  # Make an image of the touching (and non-touching) cells
  colors = list(color1, color2) %>%
    rlang::set_names(c(pheno_name1, pheno_name2))
  composite_path = find_composite_image(fake_cell_seg_path)
  touching_image = make_touching_image(pheno_name1, pheno_name2,
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
      dplyr::starts_with('Phenotype '),
      dplyr::contains('Tissue Category')) %>%
    # Suffix column names with the phenotype
    dplyr::rename_at(.vars=dplyr::vars(-dplyr::contains('Slide ID'), -!!field_col),
                     .funs = ~paste0(.x, '.', pheno_name1))

  touching2 = selected_cells2 %>%
    dplyr::select(
      `Cell ID`,
      dplyr::matches('Cell . Position'),
      dplyr::starts_with('Phenotype '),
      dplyr::contains('Tissue Category')) %>%
    dplyr::rename_all(.funs = ~paste0(.x, '.', pheno_name2))

  # Use touch_pairs to join to correctly handle duplicates - one cell
  # may be touching multiple neighbors.
  colnames(touch_pairs) = paste0('Cell ID.', c(pheno_name1, pheno_name2))
  join_table = touch_pairs %>% dplyr::as_tibble()
  data = touching1 %>% dplyr::right_join(join_table) %>% dplyr::left_join(touching2)

  list(image=touching_image, data=data)
}
