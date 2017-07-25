### Configuration
# Input directory must contain cell seg data, membrane seg image and composite for each field to be processed
base_path = 'F:/SpatialAnalysis/Fox melanoma 2x2s/Export2/'

# Output directory gets the results
output_base = 'F:/SpatialAnalysis/Fox melanoma 2x2s/Spatial2Test/'

filenames = c(
  "PCC_FIHC6__13651-10_INT-NA-2_HP_IM3_8_[27172,16862]_cell_seg_data.txt",
  "PCC_FIHC6__11243-11_INT-NA-7_HP_IM3_6_[24003,6072]_cell_seg_data.txt",
  "PCC_FIHC6__5523-09_INT-NA-L10_HP_IM3_3_[28784,13460]_cell_seg_data.txt",
  "PCC_FIHC6__3805-07_INT-NA-2_HP_IM3_3_[27733,20059]_cell_seg_data.txt",
  "PCC_FIHC6__3332-10_INT-NA_HP_IM3_0_[30922,5468]_cell_seg_data.txt",
  "PCC_FIHC6__2475-12_INT-NA-1_HP_IM3_4_[29538,6893]_cell_seg_data.txt"
)

# The phenotype pairs to locate
pairs = list(c("Cytotoxic T Cells", "PDL1+ Tumor Cells"),
             c("Cytotoxic T Cells", "macrophages"),
             c("T Regs", "macrophages"))

# Colors for all the phenotypes mentioned in pairs
colors = list(
  `Cytotoxic T Cells` = 'yellow',
  `PDL1+ Tumor Cells` = 'red',
  `PDL1- Tumor Cells` = 'cyan',
  macrophages = 'pink',
  `T Regs` = 'green'
)

### End of configuration

library(tidyverse)
library(phenoptr)

# Start with an empty output directory
if (dir.exists(output_base))
  file.remove(list.files(output_base, full.names=TRUE))

full_paths = map_chr(filenames, ~file.path(base_path, .))

# Three pairs, not mutual touches
r = map_df(full_paths, function(path) {
  cat('Processing', path, '\n')
  count_touching_cells(path, pairs, colors, output_base=output_base)
})

r$fraction = round(r$fraction, 4)

outPath = file.path(output_base, 'TouchCounts.csv')
write.csv(r, outPath, row.names=FALSE)

# Test combined phenotype, tumor only
pairs = list(c('tumor', 'lymphocyte'))
colors = list(tumor='cyan', lymphocyte='yellow')
phenotype_rules = list(
  tumor=c('PDL1- Tumor Cells', 'PDL1+ Tumor Cells'),
  lymphocyte=c('B Cells', 'Cytotoxic T Cells', 'T Regs')
)

r = map_df(full_paths, function(path) {
  cat('Processing', path, '\n')
  count_touching_cells(path, pairs, colors, phenotype_rules,
                       categories='tumor',
                       output_base=output_base)
})

r$fraction = round(r$fraction, 4)

outPath = file.path(output_base, 'TouchCounts2.csv')
write.csv(r, outPath, row.names=FALSE)
