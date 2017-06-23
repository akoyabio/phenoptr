# informr 0.1.0.9002

## Backwards incompatible change

- This version changed the order of arguments to `subset_distance_matrix`
  in a backwards-incompatible way. This was done to put the `csd` parameter
  first, matching other functions with a `csd` parameter.
  
## New features

- `count_within` counts the number of `from` cells having a `to` cell
   within a given radius.

## Bug fixes

- Better handling of `NA` values in distance columns of cell seg tables.
  Previously `NA` values could cause the column to be read as character data.
  
## Other changes

- Many documentation improvements.
- Internal cleanup using `lintr` and `goodpractices`.

# informr 0.1.0.9001

## New features

- `compute_all_nearest_distance` is a convenience functien which reads a
  cell seg table, adds `Distance to <phenotype>` columns, and writes it out
  again.
