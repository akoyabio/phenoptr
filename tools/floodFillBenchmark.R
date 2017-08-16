devtools::load_all(".")
cell_seg_path <- sample_cell_seg_path()
pairs <- list(c("CD68+", "CD8+"))
phenotype_rules=NULL
phenotypes = unique(do.call(c, pairs))
phenotype_rules = make_phenotype_rules(phenotypes, phenotype_rules)
csd = read_cell_seg_data(cell_seg_path, pixels_per_micron=NA)
mask_path = sub('cell_seg_data.txt', 'binary_seg_maps.tif', cell_seg_path)
masks = read_maps(mask_path)
membrane = masks[['Membrane']]
nuclei = masks[['Nucleus']]
rm(masks)
membrane[membrane>0] = 0.5
membrane = t(membrane)
nuclei = t(nuclei)

phenotype='CD68+'
rule = phenotype_rules[[phenotype]]
d = csd[select_rows(csd, rule),]
image = membrane
cell_ids = d$`Cell ID`
nuc_locations = purrr::map(cell_ids, ~find_interior_point(nuclei, .x))
cell_ids = cell_ids[!is.null(nuc_locations)]
nuc_locations = nuc_locations[!is.null(nuc_locations)]
x = EBImage::floodFill(image, list(nuc_locations),
                           list(as.list(cell_ids)))
display(x)

y = image
for (i in seq_len(length(nuc_locations)))
{
  y = EBImage::floodFill(y, nuc_locations[i], cell_ids[i])
}
display(y)
all(y==x)

library(microbenchmark)
dim(image)
length(nuc_locations)
microbenchmark(EBImage::floodFill(image, list(nuc_locations),
                                  list(as.list(cell_ids))), times=10)

microbenchmark({y=image
  for (i in seq_len(length(nuc_locations)))
    y = EBImage::floodFill(y, nuc_locations[i], cell_ids[i])},
  times=10)

