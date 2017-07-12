## ----setup, echo=FALSE,include=FALSE,message=FALSE-----------------------
suppressPackageStartupMessages(library(tidyverse))
knitr::opts_chunk$set(fig.width=7, fig.height=5, 
                      comment=NA, warning=FALSE, message=FALSE)
theme_set(theme_bw())
# No figure margins
par(mar=rep(0, 4))

## ------------------------------------------------------------------------
library(tidyverse)
library(informr)

csd = sample_cell_seg_data
csd %>% count(Phenotype)

## ------------------------------------------------------------------------
csd = csd %>% filter(!Phenotype %in% c('T reg Foxp3', 'other'))

## ------------------------------------------------------------------------
distances = find_nearest_distance(csd)
glimpse(distances)
nrow(csd)

## ------------------------------------------------------------------------
csd_with_distance = bind_cols(csd, distances)

## ----eval=FALSE----------------------------------------------------------
#  merged_with_distance = merged %>%
#    dplyr::group_by(`Sample Name`) %>%
#    dplyr::do(dplyr::bind_cols(., find_nearest_distance(.)))

## ------------------------------------------------------------------------
csd_with_distance %>% group_by(Phenotype) %>% 
  select(starts_with('Distance to')) %>% 
  summarize_all(~round(mean(.), 1))

## ------------------------------------------------------------------------
ggplot(csd_with_distance, aes(`Distance to helper CD4`, color=Phenotype)) +
  geom_density()

## ------------------------------------------------------------------------
count_within(csd, from='macrophage CD68', to='tumor', radius=25)

## ------------------------------------------------------------------------
count_within(csd, from='tumor', to='macrophage CD68', radius=25)

## ----eval=FALSE----------------------------------------------------------
#  base_path = "/path/to/my_directory"
#  
#  from = list('T Reg')
#  to = list('Cytotoxic T Cell', 'PDL1+ Tumor Cell')
#  radii = c(10, 25)
#  categories = c('stroma', 'tumor')
#  
#  count_within_batch(base_path, from, to, radii, categories)

## ----eval=FALSE----------------------------------------------------------
#  cell_seg_path = system.file("extdata", "TMA",
#                         "Core[1,5,6,1]_[21302,15107]_cell_seg_data.txt",
#                         package = "informr")
#  
#  phenotypes = c("macrophage CD68", "cytotoxic CD8")
#  colors = c('red', 'blue')
#  out_path = path.expand('~/spatial_distribution_report.html')
#  
#  spatial_distribution_report(cell_seg_path, phenotypes, colors, out_path)
#  

