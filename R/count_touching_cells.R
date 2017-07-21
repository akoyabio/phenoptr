# Find, count and image touching cells by morphological analysis
# of the membrane mask.


#' Find and count touching cells for pairs of phenotypes.
#'
#' `count_touching_cells` uses morphological analysis of nuclear and
#' membrane segmentation maps to find touching cells of paired phenotypes.
#' It reports the number of touching cells found and, optionally,
#' writes image files showing the touching cells.
#'
#' The number of touching cells is reported in three ways. For a pair of
#' phenotypes A and B, this function reports the number of A touching a B,
#' the number of B touching an A, and the number of mutually touching pairs.
#' Note that the number of mutual touches is often larger than the count of
#' either "A touching B" or "B touching A" because a single touching cell
#' may be part of multiple pairs.
#'
#' The image files written show cells of the selected phenotypes on a background
#' of the composite. Touching cells are filled in the provided color;
#' cells which are not touching the other phenotype are outlined.
#' Image files are written as
#' TIFF files to preserve the fine detail of cell boundaries.
#'
#' Images are only written when both phenotypes of the pair are represented.
#'
#' See `vignette('selecting_cells', package='phenoptr')` for more on
#' the use of `pairs` and `phenotype_rules`.
#'
#' @param cell_seg_path The path to the cell seg data file. The same directory
#'   must also contain `_memb_seg_map.tif` or `_binary_seg_maps.tif` and, if
#'   `write_images` is true, a  TIFF or JPEG composite image.
#' @param pairs A list of pairs of phenotypes. Each entry is a two-element
#'   vector. The result will contain one line for each pair showing the
#'   number of cells and number of touches.
#' @param colors A named list of phenotype colors to use when drawing
#'   the output. Only used when `write_images` is `TRUE`.
#' @param phenotype_rules (Optional) A named list.
#'   Item names are phenotype names and must match entries in `pairs`.
#'   Item values are selectors for [select_rows].
#' @param categories If given, a vector or list of tissue category names.
#' Categories not in the list will be excluded from the analysis.
#' @param write_images If `TRUE`, for each pair, write an image showing the
#' touching pairs. Requires `colors` and a composite image in the same
#' directory as the cell seg table.
#' @param output_base Base path for image output.
#' If `NULL`, output will be to the same
#' directory as the cell table.
#' @return Returns a `data_frame` with one row for each pair in
#' `pairs`, containing these columns:
#'   \describe{
#'    \item{\code{slide_id}}{Slide ID from the data file.}
#'    \item{\code{source}}{Base file name of the source file with
#'    `_cell_seg_data.txt` stripped off for brevity.}
#'    \item{\code{phenotype1}}{The name of the first phenotype in
#'    the touching pair.}
#'    \item{\code{phenotype2}}{The name of the second phenotype in
#'    the touching pair.}
#'    \item{\code{total1}}{The total number of `phenotype1` cells
#'    in the image.}
#'    \item{\code{total2}}{The total number of `phenotype2` cells
#'    in the image.}
#'    \item{\code{p1_touch_p2}}{The number of `phenotype1` cells
#'    touching a `phenotype2` cell.}
#'    \item{\code{p2_touch_p1}}{The number of `phenotype2` cells
#'    touching a `phenotype1` cell.}
#'    \item{\code{touch_mutual}}{The number of mutually touching pairs.}
#'  }
#' @examples
#' \dontrun{
#' # This example creates an image in a subdirectory of the
#' # current user's directory.
#' cell_seg_path = system.file("extdata", "TMA",
#'                        "Core[1,5,6,1]_[21302,15107]_cell_seg_data.txt",
#'                        package = "phenoptr")
#'
#' pairs = list(c("macrophage CD68", "cytotoxic CD8"))
#' colors = c("macrophage CD68"='red', "cytotoxic CD8"='blue')
#' output_base = path.expand('~/touches')
#'
#' count_touching_cells(cell_seg_path, pairs, color,
#'   output_base=output_base)
#'
#' # This example will count and image all files in the `base_path` directory.
#' base_path = '/path/to/data'
#' output_base = file.path(base_path, 'touches')
#' files = list_cell_seg_files(base_path)
#'
#' # The phenotype pairs to locate. This will find CD8 cells touching
#' # tumor cells, and, separately, CD8 cells touching CD68 cells.
#' pairs = list(c("CD8", "Tumor"),
#'              c("CD8", "CD68"))
#'
#' # Colors for all the phenotypes mentioned in pairs
#' colors = list(
#'   CD8 = 'yellow',
#'   Tumor = 'cyan',
#'   CD68 = 'pink'
#' )
#'
#' # Count and visualize touching cells
#' touch_counts = purrr::map_df(files, function(path) {
#'   cat('Processing', path, '\n')
#'   count_touching_cells(path, pairs, colors, output_base=output_base)
#' })
#'
#' # Save the result
#' touches_path = file.path(output_base, 'TouchCounts.csv')
#' readr::write_csv(touch_counts, touches_path)
#'
#' # The phenotype definitions can be more complex. The default is to use
#' # the names in `pairs`. Using `phenotype_rules`, the definition can be
#' # anything allowed by select_rows().
#'
#' # You can also limit the tissue category.
#'
#' # For example, find all touches between lymphocytes and tumor cells
#' # within the tumor:
#' pairs = list(c('Tumor', 'Lymphocyte'))
#' colors = list(Tumor='cyan', Lymphocyte='yellow')
#' phenotype_rules = list(
#'   Lymphocyte=c('CD4', 'CD8', 'TReg')
#' )
#'
#' touch_counts = map_df(files, function(path) {
#'   cat('Processing', path, '\n')
#'   count_touching_cells(path, pairs, colors, phenotype_rules,
#'                        categories='tumor',
#'                        output_base=output_base)
#' })
#'
#' # Then write the results as above.
#' }
#' @md
#' @export
#' @family distance functions
#' @importFrom magrittr %>%
count_touching_cells = function(cell_seg_path, pairs, colors=NULL,
                                phenotype_rules=NULL, categories=NULL,
                                write_images=!is.null(colors), output_base=NULL)
{
  # Allow a single pair to be specified as a plain vector
  if (is.character(pairs) && length(pairs)==2)
    pairs = list(pairs)

    # Make phenotype_rules for any not already specified
  phenotypes = unique(do.call(c, pairs))
  phenotype_rules = make_phenotype_rules(phenotypes, phenotype_rules)

  # Check the requirements for writing images, if requested
  if (write_images) {
    if (is.null(colors))
      stop('Colors are required when write_images is TRUE.')
    stopifnot(is.null(colors) || all(phenotypes %in% names(colors)))

    # Look for composite as TIFF or JPEG
    composite_path =
      sub('cell_seg_data.txt', 'composite_image.tif', cell_seg_path)
    if (!file.exists(composite_path))
      composite_path = sub('tif', 'jpg', composite_path)

    if (!file.exists(composite_path))
      stop('write_images requires a TIFF or JPEG composite image.')
  }

  # Make the output directory
  if (!is.null(output_base) && !file.exists(output_base))
    dir.create(output_base, showWarnings=FALSE, recursive=TRUE)

  # Read the data. We don't want to convert to microns here, we need
  # image coordinates, so don't call `read_cell_seg_table`
  name = basename(cell_seg_path)
  name = sub('_cell_seg_data.txt', '', name)
  csd = readr::read_tsv(cell_seg_path, na=c('NA', '#N/A'),
                        col_types=readr::cols())

  # Filter out unwanted tissue categories
  if (!is.null(categories))
  {
    stopifnot('Tissue Category' %in% names(csd))
    csd = csd %>% dplyr::filter(`Tissue Category` %in% categories)
  }

  slide = as.character(csd[1, 'Slide ID'])

  # Read the membrane and nuclear masks.
  # Convert the membrane to single values. Use 0.5 to avoid conflict with cell labeling
  mask_path = sub('cell_seg_data.txt', 'memb_seg_map.tif', cell_seg_path)
  if (file.exists(mask_path))
  {
    # Old-style membrane mask
    membrane = EBImage::readImage(mask_path)
    membrane[membrane>0] = 0.5

    # Don't use readImage to read nuclear mask, it converts to 0-1 scale!
    nuc_path = sub('cell_seg_data.txt', 'nuc_seg_map.tif', cell_seg_path)
    nuclei = tiff::readTIFF(nuc_path, as.is=TRUE)
    nuclei = t(nuclei)
  }
  else
  {
    mask_path = sub('cell_seg_data.txt', 'binary_seg_maps.tif', cell_seg_path)
    masks = read_maps(mask_path)
    membrane = masks[['Membrane']]
    nuclei = masks[['Nucleus']]
    rm(masks)
    membrane[membrane>0] = 0.5
    membrane = t(membrane)
    nuclei = t(nuclei)
  }
  stopifnot(exists('membrane'), exists('nuclei'))

  # Make images for the cells in each phenotype by filling in the membrane mask
  # at each cell, then removing the membrane.
  # Fill with a unique ID value for each cell (here, the cell ID)
  cell_images = lapply(phenotypes, function(phenotype)
  {
    rule = phenotype_rules[[phenotype]]
    d = csd[select_rows(csd, rule),]
    make_cell_image(d, nuclei, membrane)
  })
  names(cell_images) = phenotypes

  # Will be a data frame with counts etc
  result = NULL

  # Process each pair of phenotypes
  for (pair in pairs)
  {
    # Get the names and images for each phenotype
    p1 = pair[1]
    i1 = cell_images[[p1]]

    p2 = pair[2]
    i2 = cell_images[[p2]]

    p1_count = sum(select_rows(csd, phenotype_rules[[p1]]))
    p2_count = sum(select_rows(csd, phenotype_rules[[p2]]))

    touches_found = 0

    if (p1_count == 0 || p2_count == 0)
    {
      # No data for one of the phenotypes in this pair
      # Report empty result and go on
      result = rbind(result,
                   tibble::data_frame(
                      slide_id=slide,
                      source=name,
                      phenotype1=p1,
                      phenotype2=p2,
                      total1=p1_count,
                      total2=p2_count,
                      p1_touch_p2=0,
                      p2_touch_p1=0,
                      touch_pairs=0))
      if (write_images)
        warning('No image for ', name, ', ', p1, ' touching ', p2)
      next
    }

    touch_pairs = find_touching_cell_pairs(i1, i2)

    # Need individual IDs for imaging and counting
    touching_ids = list(unique(touch_pairs[,1]), unique(touch_pairs[,2]))
    p1_touching_count = length(touching_ids[[1]])
    p2_touching_count = length(touching_ids[[2]])
    touches_found = nrow(touch_pairs)

    result = rbind(result,
                 tibble::data_frame(
                    slide_id=slide,
                    source=name,
                    phenotype1=p1,
                    phenotype2=p2,
                    total1=p1_count,
                    total2=p2_count,
                    p1_touch_p2=p1_touching_count,
                    p2_touch_p1=p2_touching_count,
                    touch_pairs=touches_found))

    if (write_images)
    {
      tag = paste0(p1, '_', p2, '_touching.tif')
      composite_out = sub('composite_image.(tif|jpg)', tag, composite_path)
      if (!is.null(output_base))
        composite_out = file.path(output_base, basename(composite_out))

      # Make a pretty picture showing the touch points.
      # First make a mask containing just the touching cells by searching for
      # the touching IDs in the original cell images. Then dilate to overlap the
      # membrane mask and mask out anything not in the membrane mask.
      both = EBImage::Image(dim=dim(i1))
      both[i1 %in% touching_ids[[1]]] = 1
      both[i2 %in% touching_ids[[2]]] = 1

      kern3 = EBImage::makeBrush(3, shape='diamond')
      both = EBImage::dilate(both, kern3)
      both[membrane==0] = 0

      composite = EBImage::readImage(composite_path)
      if (!is.null(colors))
      {
        # If we have colors, outline all cells of a type; fill the touching cells
        i1_touching = i1
        i1_touching[!i1 %in% touching_ids[[1]]] = 0

        composite = EBImage::paintObjects(i1, composite, col=c(colors[[p1]], NA))
        composite = EBImage::paintObjects(i1_touching, composite, col=c(colors[[p1]], colors[[p1]]))

        i2_touching = i2
        i2_touching[!i2 %in% touching_ids[[2]]] = 0

        composite = EBImage::paintObjects(i2, composite, col=c(colors[[p2]], NA))
        composite = EBImage::paintObjects(i2_touching, composite, col=c(colors[[p2]], colors[[p2]]))
      }

      # Draw the cell outlines onto the composite image and save
      composite = EBImage::paintObjects(both, composite, col='white')
      EBImage::writeImage(composite, composite_out, compression='LZW')
    }
  }

  # Return the data
  result
}

# Given a data frame of cells and membrane and nuclear masks, make an image with
# a region for each cell. Returns NULL if d is empty
make_cell_image <- function (d, nuclei, membrane) {
  stopifnot('Cell X Position' %in% names(d),
            'Cell Y Position' %in% names(d),
            'Cell ID' %in% names(d))

  if (nrow(d)==0) return(NULL) # No cells in this data

  image = membrane
  for (i in 1:nrow(d))
  {
    cell_id = d$`Cell ID`[i]

    nuc_locations = find_interior_point(nuclei, cell_id)

    if (is.null(nuc_locations)) next # Didn't find this cell

    image = EBImage::floodFill(image, nuc_locations, cell_id)
  }
  # Remove the membrane outlines
  image[image==0.5] = 0
  image
}

# This uses a distance map to find the interior-most point in a nucleus
find_interior_point = function(nuclei, cell_id)
{
  # Locate the nucleus in the overall map and extract it as a patch.
  # We could mask out the nucleus in the full image and compute the distance
  # map on that; this is faster. Extracting the patch is fast and the distance
  # map calculation is O(M*N*log(max(M,N))) so smaller M, N helps a lot.
  # We do have to mask out the nucleus because it may be touching other nuclei.
  # Alternatively we could mask out the membrane and compute the full distance
  # map. The approach taken here seems safer.

  # Where is the nucleus?
  nuc_locations = which(nuclei==cell_id, arr.ind=TRUE)
  if (nrow(nuc_locations)==0) return(NULL) # Didn't find this cell

  # Figure out the bounds of the patch. Include a 1-pixel border if possible
  row_min = max(1, range(nuc_locations[,1])[1]-1)
  row_max = min(dim(nuclei)[1], range(nuc_locations[,1])[2]+1)
  row_range =  seq.int(row_min, row_max)

  col_min = max(1, range(nuc_locations[,2])[1]-1)
  col_max = min(dim(nuclei)[2], range(nuc_locations[,2])[2]+1)
  col_range =  seq.int(col_min, col_max)

  # Extract the actual patch
  n = nuclei[row_range, col_range]

  # Clear out extra stuff and compute a distance map on the patch
  n[n!=cell_id] = 0
  dm = EBImage::distmap(n)

  # Where is the interior-most point? We are happy to take the first one.
  ix_raw = which.max(dm)

  # Arrays are stored column-wise in R
  col_length = dim(dm)[1]
  row = ix_raw %% col_length # modulus
  col = ix_raw %/% col_length + 1 # integer division

  # Offset by the patch origin. Subtract one because we are adding
  # two one-based indices
  c(row+row_min-1, col+col_min-1)
}

# Given two cell images (from make_cell_image), find the IDs of the cells in
# each image that touch cells in the other image.
find_touching_cell_ids <- function (i1, i2) {
  # Make binary masks for the cells
  i1_mask = i1
  i1_mask[i1_mask>0] = 1

  i2_mask = i2
  i2_mask[i2_mask>0] = 1

  # Dilate and look for intersections.
  # Order matters here. If we want to know how many p1s are touching a p2, we
  # should dilate the p2s. To count p2's touching a p1, dilate the p2s.
  # We are doing both...
  kern5 = EBImage::makeBrush(5, shape='diamond')
  i1_big = EBImage::dilate(i1_mask, kern5)
  i2_big = EBImage::dilate(i2_mask, kern5)

  # Find p1s touching p2s. overlap will have a non-zero strip within each p1
  # that touches a p2. The value in the strip will be the p1 ID. The number of
  # unique, non-zero IDs is the count of touching cells.
  overlap = i1
  overlap[i2_big==0] = 0
  p1_touching_ids = unique(as.numeric(overlap))
  p1_touching_ids = p1_touching_ids[p1_touching_ids>0]

  # Same thing for p2s touching p1s.
  overlap2 = i2
  overlap2[i1_big==0] = 0
  p2_touching_ids = unique(as.numeric(overlap2))
  p2_touching_ids = p2_touching_ids[p2_touching_ids>0]

  list(p1_touching_ids, p2_touching_ids)
}

# Given two cell images (from make_cell_image), find the IDs of the cells in
# each image that touch each other. Returns a matrix with two columns,
# the cell numbers in i1 and i2
find_touching_cell_pairs <- function (i1, i2) {
  # Dilate i1 and look for intersections with i2
  # Order doesn't matter here, we are looking for kissing pairs.
  kern3 = EBImage::makeBrush(3, shape='diamond')
  i1_big = EBImage::dilate(i1, kern3)
  i2_big = EBImage::dilate(i2, kern3)

  # Find p1s touching p2s as pairs.
  overlap = cbind(as.numeric(i1_big), as.numeric(i2_big))
  overlap = overlap[overlap[,1]>0.1 & overlap[,2]>0.1,]
  overlap = unique(overlap)
  overlap
}
