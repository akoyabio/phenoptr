# phenoptr 0.3.2.9000

# phenoptr 0.3.2
**2022-01-03**

- Fix `parse_composite_info` to work with older composite images
  that are missing ImageDescription (akoyabio/phenoptrReports#52)
- Correctly parse phenotype expressions containing slashes (#18)
- Tighten up the column exclusion criteria for 
  `read_cell_seg_data` with `col_select="phenoptrReports"` (#16)
- Update tests to work with testthat 3.1.

# phenoptr 0.3.1
**2021-08-24**

- Add `read_composites` and `read_composite_info` to read and parse
  metadata in inForm composite images.

# phenoptr 0.3.0
**2021-08-03**

Misc:
- Update dependencies to modern versions
- Fix for [bug in readr 2.0.0](https://github.com/tidyverse/readr/issues/1258)

# phenoptr 0.2.10
**2021-06-01**

Changes to `read_cell_seg_data`:
- Add a `col_select` parameter to support reading only
  selected fields, with a `'phenoptrReports'` option to read only fields
  needed by `phenoptrReports`.
- Use `vroom::vroom()` instead of `readr::read_tsv` for speed and support
  of `col_select`.
- For an export file containing all inForm fields,
  `read_cell_seg_data(csd_path, col_select='phenoptrReports')` can be
  4-5x faster and use 75% less memory for the results,
  compared to the previous version.

Spatial distribution report:
- Fix report generation to work when paths contain spaces.

Misc:
- Fix `select_rows` to always return `FALSE` instead of `NA`.
- Update to work with CRAN release of tiff 0.1-8
- Update to work with (and require) spatstat 2.0.0
- Require R 4.0.0 or greater for compatibility with latest rtree

# phenoptr 0.2.9
**2020-11-11**

- Improved the way `read_cell_seg_data` determines column types to prevent
  errors in some edge cases.
- Fixed `make_ppp` to work correctly when given compound phenotypes 
  such as "CD8+/PD-L1+"
  
# phenoptr 0.2.8
**2020-08-22**

- Add `make_ppp` function to create marked point patterns.
- `parse_phenotypes` makes nicer names for expression-based phenotypes.
- Add `phenotype_columns` in support of akoyabio/phenoptrReports#36.
- Add `read_field_info` and `get_map_path` functions.

# phenoptr 0.2.7
**2020-06-02**

Fix to work with dplyr 1.0.0.

# phenoptr 0.2.6
**2020-04-21**

- Avoid deprecation warning for `expand_scale` with ggplot2 version 3.3.0

# phenoptr 0.2.5
**2020-02-28**

- Make `read_cell_data` more robust against data with comma as the decimal separator.
- Better documentation and error messages in `count_touching_cells_fast`.

# phenoptr 0.2.4
**2019-10-22**

- Add `count_touching_cells_fast` for phenoptrReports spatial map viewer.
- Don't return NA values from `select_rows()`; return `FALSE` instead.
- Now requires tidyr >= 1.0.0.

# phenoptr 0.2.3
**2019-08-07**

- Add `annotation_raster_native` for faster plotting with background images.

# phenoptr 0.2.2
**2019-06-11**

- Remove the explicit dependency on `phenoptrReports` to avoid installation pain.
  `phenoptrReports` is still needed to build the vignettes and to run a few
  of the examples.
- `parse_phenotypes` and `validate_phenotype_definitions` support
  expressions in phenotype definitions.

# phenoptr 0.2.1
**2019-05-31**

- The spatial analysis report now uses the full-intensity composite image as 
  its background image. Connecting lines and the scale bar are white to be 
  visible against the dark background of the composite image.

Bug fixes:
- Fix `density_bands` to work with cell seg data in microns with slide
  origin (#10).
  
# phenoptr 0.2.0
**2019-05-10**

New feature (and breaking change):
- The nearest neighbor functions `compute_all_nearest_distance` and
  `find_nearest_distance` now create columns containing the Cell ID
  of the nearest cell, as well as the columns with the actual distance.
  This column can be used to find the locations of nearest cells
  and to find mutual nearest neighbors. The "Computing inter-cell
  distances" tutorial shows some uses of this data.
  
Bug fixes:
- `read_cell_seg_data` recognizes and correctly reads inForm data 
  which uses comma as the decimal separator (#8).
- Fix SpatialDistributionReport to work with more than 8 phenotypes (#6).
- `count_touching_cells` ignores phenotype pairs of self to self with a warning
  because the touching algorithm does not handle this case (#7).

# phenoptr 0.1.6.1
**2019-04-22**

- Remove the scary `Failed with error: ‘there is no package called ‘rtree’’`
  message that `requireNamespace('rtree')` generates (#5).

# phenoptr 0.1.6
**2019-04-12**

New features:
- Add fast implementation of `count_within_many` and `count_within_batch` 
  using `akoyabio/rtree` package.
- Add fast `count_within_detail` to give per-cell counts.

Bug fixes:
- Fix `get_field_info` to work with the standard (CRAN) `tiff` package.
- `count_within` and related functions return `NA` for `within_mean`
  when there are no `from` cells (rather than a mean of 0).
- Update `spatial_distribution_report` vignette and example to use
  `phenoptrExamples`, which contains the required component data file (#3).
- Fix problem in `count_touching_cells` when there is exactly 
  one touching cell pair (#4)
  
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
