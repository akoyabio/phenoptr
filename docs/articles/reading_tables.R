## ----setup, echo=FALSE,include=FALSE,message=FALSE-----------------------
suppressPackageStartupMessages(library(tidyverse))
knitr::opts_chunk$set(fig.width=6, fig.height=4, 
                      comment=NA, warning=FALSE, message=FALSE)
theme_set(theme_bw())
# No figure margins
par(mar=rep(0, 4))

## ----read----------------------------------------------------------------
# Load libraries
library(tidyverse)
library(phenoptr)

# sample_cell_seg_path gives the path to a sample file included with phenoptr.
# Change this to be the path to your data. For example you might use
# path <- 'C:/data/my_experiment/my_image_cell_seg_data.txt'
path <- sample_cell_seg_path()

# Read the data file
csd <- read_cell_seg_data(path)

# Show some nicely shortened names
# The suffix "(Normalized Counts, Total Weighting)" has been removed.
grep('Nucleus.*Mean', names(csd), value=TRUE)

## ----inspect-------------------------------------------------------------
# How many cells did we read?
nrow(csd)

# How many cells of each phenotype are in each tissue category?
table(csd$`Tissue Category`, csd$Phenotype)

## ----mutate--------------------------------------------------------------
csd <- csd %>% mutate(pdl1_plus=`Entire Cell PDL1 (Opal 520) Mean`>3)
table(csd$pdl1_plus, csd$Phenotype)

## ----aggregate-----------------------------------------------------------
csd %>% 
  filter(Phenotype!='other') %>% 
  group_by(Phenotype) %>% 
  summarize(mean_pdl1=mean(`Entire Cell PDL1 (Opal 520) Mean`))

## ----ggplot2-------------------------------------------------------------
ggplot(csd, aes(Phenotype, `Entire Cell PDL1 (Opal 520) Mean`, color=Phenotype)) +
  geom_boxplot() + 
  scale_color_brewer(palette='Set1') + 
  labs(y='PDL1 Expression', title='PDL1 Expression per Phenotype')

ggplot(csd %>% filter(Phenotype!='other'), 
       aes(`Entire Cell PDL1 (Opal 520) Mean`, 
                `Entire Cell PD1 (Opal 650) Mean`,
                color=Phenotype)) +
  geom_point(size=1, alpha=0.2) + 
  facet_wrap(~Phenotype) +
  scale_x_log10() + scale_y_log10() + scale_color_brewer(palette='Set1') +
  labs(x='PDL1 Expression', y='PD1 Expression', 
       title='Comparison of PD1 and PDL1 Expression per Phenotype')

