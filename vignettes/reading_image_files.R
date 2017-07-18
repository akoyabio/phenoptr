## ----setup, echo=FALSE,include=FALSE,message=FALSE-----------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5, 
                      comment=NA, warning=FALSE, message=FALSE)
# No figure margins
par(mar=rep(0, 4))

## ----read_composite------------------------------------------------------
path = system.file("extdata", "TMA", 
                   "Core[1,5,6,1]_[21302,15107]_composite_image.jpg", 
                   package = "phenoptr")
img = jpeg::readJPEG(path)
dim(img)

## ----show_composite------------------------------------------------------
snippet = as.raster(img[800:1000, 1200:1500,])
plot(snippet)

## ----show_composite_ebimage,eval=FALSE-----------------------------------
#  img_transposed = aperm(img, c(2, 1, 3))
#  EBImage::display(img_transposed)
#  EBImage::display(EBImage::Image(img_transposed, colormode='Color'))

## ----read_map------------------------------------------------------------
map_path = system.file("extdata",
   "TMA/Core[1,5,6,1]_[21302,15107]_binary_seg_maps.tif", package = "phenoptr")
maps = phenoptr::read_maps(map_path)
names(maps)

## ----nucleus-------------------------------------------------------------
nucleus = maps[['Nucleus']]
dim(nucleus)
range(nucleus)
range(phenoptr::sample_cell_seg_data$`Cell ID`)

## ----display_map---------------------------------------------------------
nucleus_snippet = nucleus[800:1000, 1200:1500]
plot(as.raster(nucleus_snippet, max=max(nucleus_snippet)))

## ----display_map_ebimage,eval=FALSE--------------------------------------
#  EBImage::display(t(nucleus/max(nucleus)))
#  EBImage::display(t(maps[['Membrane']]))

## ----plotting------------------------------------------------------------
library(ggplot2)
d = phenoptr::sample_cell_seg_data
d = subset(d, `Cell X Position`>=600 & `Cell X Position`<=750 &
              `Cell Y Position`>=400 & `Cell Y Position`<=500)
ggplot(data=d, aes(`Cell X Position`, `Cell Y Position`, color=Phenotype)) +
  scale_x_continuous(limits=c(600, 750)) + 
  scale_y_reverse(limits=c(500, 400)) +
  annotation_raster(snippet, 600, 750, -400, -500) +
  geom_point(size=2) + coord_equal() +
  scale_color_manual(values=c("tumor"="cyan", "cytotoxic CD8"="yellow",
                              "other"="blue", "macrophage CD68"="red", 
                              "helper CD4"="green"))

