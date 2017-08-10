## ----setup, echo=FALSE,include=FALSE,message=FALSE-----------------------
suppressPackageStartupMessages(library(tidyverse))
knitr::opts_chunk$set(fig.width=7, fig.height=5, 
                      comment=NA, warning=FALSE, message=FALSE)
theme_set(theme_bw())
# No figure margins
par(mar=rep(0, 4))

## ------------------------------------------------------------------------
library(tidyverse)
library(phenoptr)

csd <- sample_cell_seg_data
csd %>% count(Phenotype)

## ------------------------------------------------------------------------
csd <- csd %>% filter(Phenotype!='other')

## ------------------------------------------------------------------------
distances <- find_nearest_distance(csd)
glimpse(distances)
nrow(csd)

## ------------------------------------------------------------------------
csd_with_distance <- bind_cols(csd, distances)

## ----eval=FALSE----------------------------------------------------------
#  merged_with_distance <- merged %>%
#    dplyr::group_by(`Sample Name`) %>%
#    dplyr::do(dplyr::bind_cols(., find_nearest_distance(.)))

## ------------------------------------------------------------------------
csd_with_distance %>% group_by(Phenotype) %>% 
  select(starts_with('Distance to')) %>% 
  summarize_all(~round(mean(.), 1))

## ------------------------------------------------------------------------
ggplot(csd_with_distance, aes(`Distance to CD8+`, color=Phenotype)) +
  geom_density(size=1) + theme_minimal()

## ------------------------------------------------------------------------
count_within(csd, from='CD68+', to='CK+', radius=25)

## ------------------------------------------------------------------------
count_within(csd, from='CK+', to='CD68+', radius=25)

## ----eval=FALSE----------------------------------------------------------
#  base_path <- "/path/to/my_directory"
#  
#  pairs <- list(c('FoxP3+', 'CD8+'),
#                c('FoxP3+', 'CK+'))
#  radii <- c(10, 25)
#  categories <- c('Stroma', 'Tumor')
#  
#  count_within_batch(base_path, pairs, radii, categories)

## ----eval=FALSE----------------------------------------------------------
#  cell_seg_path <- sample_cell_seg_path()
#  
#  pairs <- list(c('CD68+', 'CD8+'))
#  colors <- c('CD68+'='magenta', 'CD8+'='yellow')
#  out_path <- path.expand('~/spatial_distribution_report.html')
#  
#  spatial_distribution_report(cell_seg_path, pairs, colors, output_path=out_path)
#  

## ----eval=FALSE----------------------------------------------------------
#  base_path <- '/path/to/data/'
#  paths <- list_cell_seg_files(base_path)
#  for (path in paths)
#    spatial_distribution_report(path, pairs, colors)

