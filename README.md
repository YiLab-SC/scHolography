# scHolography: a workflow for single-cell 3D spatial reconstruction <a href='https://github.com/YiLab-SC/scHolography'><img src="img/HexSticker.png" align="right" height="200" /></a>


scHolography is an neural network-based computational toolkit for integrative scRNA-seq and ST data analysis. Our pipeline enables 3D spatial inference at a high resolution. Instead of mapping cells to a spot on the fixed ST slice, our reconstruction result is orientation-free, and the inferred structure will be varied for different scRNA-seq data inputs. Together with our downstream analytical functions, we aim to bring new perspectives on scRNA-seq and ST data for researchers.

## 0. Installation
### Dependencies

The deep learning functionalities of our package powers by the Keras API. 
To install Keras:
```r
install.packages("remotes")
remotes::install_github(sprintf("rstudio/%s", c("reticulate", "tensorflow", "keras")))
reticulate::miniconda_uninstall() 
reticulate::install_miniconda()
keras::install_keras()
```

To confirm if the installation is successful:
```r
library(tensorflow)
tf$constant("Hello Tensorflow!")
```

Successful installation will give the output:
```
Loaded TensorFlow version 2.10.0
tf.Tensor(b'Hello Tensorflow!', shape=(), dtype=string)
```

DO NOT PROCEED TO THE NEXT STEP IF THE KERAS INSTALLATION FAILS. Please refer to the Keras page for FAQ: https://github.com/rstudio/keras 


### Install scHolography

```r
install.packages("devtools")
devtools::install_github("YiLab-SC/scHolography")
```
Load scHolography for use:
```r
library(scHolography)
```

## 1. scHolography 3D Reconstruction

For demonstration, we are using mouse brain data. This dataset has been widely used for ST computational method devlopment and pocessed Seurat objects can be downloaded from CellTrek site:

The mouse brain scRNA-seq data: https://www.dropbox.com/s/ruseq3necn176c7/brain_sc.rds?dl=0

The mouse brain ST data: https://www.dropbox.com/s/azjysbt7lbpmbew/brain_st_cortex.rds?dl=0

Load data into our workspace:
```r
brain_st_cortex <-readRDS("~/Downloads/brain_st_cortex.rds")
brain_sc <- readRDS("~/Downloads/brain_sc.rds")
```

Load packages:
```r
library(scHolography)
library(dplyr)
```

scHolography can take in objects directly after quality control and automatically perform normalization before integration. `dataAlign` integrations expression data from scRNA-seq and ST modalities. `trainHolography` trains neural network models and infers a SMN graph.
```r
options(future.globals.maxSize = 3000 * 1024^2)
sp.integrated <- dataAlign(low.res.sp =  brain_st_cortex,high.res.sp =  brain_sc,nPCtoUse = 32)
brain.obj<-trainHolography(sp.integrated,n.repeat = 30)
```

The reconstruction results can be visualized in 3D with `scHolographyPlot`:
```r
scHolographyPlot(brain.obj,color.by = "cell_type")
```

```
scene = list(camera = list(eye = list(x = 0., y = 0.01, z = -2)))
fig <- scHolographyPlot(brain.obj,color.by = "cell_type")%>% plotly::layout(scene = scene) 
fig
```
![](img/mouse.brain.3D.svg)








## X. Running Time
```r
system.time(brain.obj<-trainHolography(sp.integrated,n.repeat = 30))
```
```
user  system elapsed 
830.224 172.764 901.495 
```








