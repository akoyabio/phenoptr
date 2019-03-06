# phenoptr 0.1.5.9000

- Fix `get_field_info` to work with the standard (CRAN) `tiff` package.

# phenoptr 0.1.5
**2019-03-03**

- Add `count_within_many` in support of phenoptrReports.
- `count_within` does not count self when 'from' and 'to' phenotypes are the same.
- Add `field_column` (moved from phenoptrReports) and `validate_phenotypes`.
- Fix `compute_all_nearest_distance` to recognize "Annotation ID" as a field 
  name.
  
# phenoptr 0.1.4
**2019-02-08**

- Update `find_nearest_distance`, `compute_all_nearest_distance`,
  `density_at_distance`, `density_bands`, `count_within` and
  `spatial_distribution_report` to work
  with consolidated data from `phenoptrReports` (#1).
  
Bug fixes

- Change `sq microns` and `sq mm` to `square microns` and `square mm` to match
  inForm output.
  
# phenoptr 0.1.3
**2018-12-27**

- `spatial_distribution_report` works correctly with cell seg data in microns.
  To do so, it requires a component data file for the 
  target field. Image dimensions are taken from that file.
- `count_touching_cells` works correctly when cell seg data is in microns.
- Require readr version 1.2.0 or higher to avoid incorrectly reading numeric
  data columns as integer.
- Add a simple data quality check to `read_cell_seg_data`
  [#9](https://github.com/PerkinElmer/phenoptr/issues/9)
- `count_within` won't complain if the data doesn't have a `Sample Name` column.
- More accurate image location from component file if PerkinElmer tiff library
  is available.
  
# phenoptr 0.1.2
**2018-10-29**

- Update link to Phenoptics home page.

# phenoptr 0.1.1.9007
**2018-09-28**

- Remove import dependency on tidyverse. Examples and vignettes use dplyr, 
  ggplot2 and magrittr instead. #6

# phenoptr 0.1.1.9006
**2018-09-28**

- New `parse_phenotypes` function simplifies creation of selectors for 
  `select_rows`.
- `select_rows` recognizes the column format of phenotype-per-marker data.
- `select_rows` treats a `NA` selector as "select all". This is helpful
  in lists of selectors.
- `read_cell_seg_data` attempts to skip microns conversion if it was already 
  done.

# phenoptr 0.1.1.9005
**2018-07-11**

- Fix `count_touching_cells` to work with two-pixel-wide membrane map. #8
  
# phenoptr 0.1.1.9004
**2018-04-24**

- Fix `count_within_batch` to work with cell seg files which don't 
  have `Slide ID` fields.
  
# phenoptr 0.1.1.9003

- Add `density_at_distance` and `density_bands` to estimate cell density
  at a distance from a boundary.
  
# phenoptr 0.1.1.9002

- Add `pixels_per_micron="auto"` option to `read_cell_seg_data()`.
- Supports tiled component files using the PerkinElmer
  [fork of the tiff package](https://github.com/PerkinElmer/tiff).

# phenoptr 0.1.1.9001

- Faster count_touching_cells when EBImage 4.19.9 is available.

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

- `compute_all_nearest_distance` is a convenience function which reads a
  cell seg table, adds `Distance to <phenotype>` columns, and writes it out
  again.
