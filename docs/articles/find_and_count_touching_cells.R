## ----setup,echo=FALSE----------------------------------------------------
knitr::opts_chunk$set(eval=FALSE)

## ----pairs_only----------------------------------------------------------
#  cell_seg_path = sample_cell_seg_path()
#  pairs = list(
#    c('CK+', 'CD8+'),
#    c('CK+', 'CD68+')
#  )
#  colors=list('CK+'='cyan', 'CD8+'='yellow', 'CD68+'='magenta')
#  count_touching_cells(cell_seg_path, pairs, colors)

## ----echo=FALSE, eval=TRUE-----------------------------------------------
# Show cached touch counts
readr::read_csv('touch_counts.csv', col_types=readr::cols())

## ----phenotype_rules-----------------------------------------------------
#  pairs = list(
#    c('CK+ PLD1+', 'CD8+'),
#    c('CK+ PLD1+', 'CD68+')
#  )
#  
#  phenotype_rules = list(
#    'CK+ PLD1+'=list('CK+', ~`Entire Cell PDL1 (Opal 520) Mean`>3)
#  )
#  colors=list('CK+ PLD1+'='cyan', 'CD8+'='yellow', 'CD68+'='magenta')
#  count_touching_cells(cell_seg_path, pairs, colors, phenotype_rules)

## ----process_directory---------------------------------------------------
#  # Directory containing data files
#  base_path = '/path/to/data'
#  
#  # A subdirectory for the results
#  output_base = file.path(base_path, 'touches')
#  
#  # All cell seg data files in base_path
#  files = list_cell_seg_files(base_path)
#  
#  # Count and visualize touching cells
#  touch_counts = purrr::map_df(files, function(path) {
#    cat('Processing', path, '\n')
#    count_touching_cells(path, pairs, colors, phenotype_rules,
#                         output_base=output_base)
#  })

## ----write_csv-----------------------------------------------------------
#  touches_path = file.path(output_base, 'TouchCounts.csv')
#  readr::write_csv(touch_counts, touches_path)

