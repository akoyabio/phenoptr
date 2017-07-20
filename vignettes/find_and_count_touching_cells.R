## ----setup,echo=FALSE----------------------------------------------------
knitr::opts_chunk$set(eval=FALSE)

## ----pairs_only----------------------------------------------------------
#  cell_seg_path = '/path/to/cell_seg_data/'
#  pairs = list(
#    c('tumor', 'CD8'),
#    c('tumor', 'CD68')
#  )
#  colors=list(tumor='cyan', CD8='yellow', CD68='magenta')
#  count_touching_cells(cell_seg_path, pairs, colors=colors)

## ----phenotype_rules-----------------------------------------------------
#  pairs = list(
#    c('tumor PLD1+', 'CD8'),
#    c('tumor PLD1+', 'CD68')
#  )
#  
#  phenotype_rules = list(
#    'tumor PLD1+'=list('tumor', ~`Entire Cell PDL1 (Opal 620) Mean`>0.5)
#  )
#  colors=list('tumor PLD1+'='cyan', CD8='yellow', CD68='magenta')
#  count_touching_cells(cell_seg_path, pairs, phenotype_rules, colors=colors)

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
#    count_touching_cells(path, pairs, colors=colors, phenotype_rules,
#                         output_base=output_base)
#  })

## ----write_csv-----------------------------------------------------------
#  touches_path = file.path(output_base, 'TouchCounts.csv')
#  readr::write_csv(touch_counts, touches_path)

