# phenoptr 0.1.1.9000

## Backwards incompatible changes

- Renamed the package from `informr` to `phenoptr`.
- This version changed the order of arguments to `subset_distance_matrix`
  in a backwards-incompatible way. This was done to put the `csd` parameter
  first, matching other functions with a `csd` parameter.
  
## New features

- `spatial_distribution_report` creates an HTML report showing the 
  location and nearest-neighbor relations of cells in a single field.
- `count_touching_cells` uses morphological analysis of nuclear and
  membrane segmentation maps to find touching cells of paired phenotypes.
- `read_components` reads component image files.
- `count_within` and `count_within_batch` count the number of `from` 
   cells having a `to` cell within a given radius.
- `list_cell_seg_files` lists all cell seg data files in a folder.
- Added vignettes
    - "Reading and Exploring inForm Tables"
    - "Reading and Displaying inForm Image Files"
    - "Computing Inter-cellular Distances"
    - "Find and count touching cells"
    - "Selecting cells within a cell segmentation table"

## Bug fixes

- Better handling of `NA` values in distance columns of cell seg tables.
  Previously `NA` values could cause the column to be read as character data.
  
## Other changes

- `read_maps` will find the correct path when given a cell seg table path.
- Many documentation improvements.
- Internal cleanup using `lintr` and `goodpractices`.
- zlib license

# phenoptr 0.1.0.9001

## New features

- `compute_all_nearest_distance` is a convenience functien which reads a
  cell seg table, adds `Distance to <phenotype>` columns, and writes it out
  again.
