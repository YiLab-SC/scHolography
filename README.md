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

**DO NOT PROCEED TO THE NEXT STEP IF KERAS INSTALLATION FAILS**. Please refer to the Keras page for FAQ: https://github.com/rstudio/keras 


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

## 2. SMN Distance and First-Degree Neighbors

### clusterDistanceBoxplot
Based on the reconstructed graph, we can investigate spatial relationship among cell clusters of interest by calculating their SMN distances. Here we demonstrate a boxplot for distances between the L2/3 IT layer to all glutamatergic cells.

```r
clusterDistanceBoxplot(brain.obj,annotationToUse = "cell_type",query.cluster.list = c( "L2/3 IT", "L4","L5 IT","L5 PT","NP", "L6 IT", "L6 CT","L6b" ),reference.cluster="L2/3 IT")
```
![](img/Brain.clusterDist.png)

### scHolographyNeighborCompPlot
For each cell cluster, we can further examine its inferred microenvironment by dissecting their first-degree neighbor cell type composition on the SMN graph. scHolography enables this query with `scHolographyNeighborCompPlot`. This function output **1.** first-degree neighbor composition plot of all neighbors; **2.** first-degree neighbor composition plot for only enriched neighboring cell types; **3.** significance levels for each enriched neighboring cell type.
```{r}
neighbor.comp <- scHolographyNeighborCompPlot(brain.obj,annotationToUse = "cell_type")
neighbor.comp$neighbor.comp
neighbor.comp$neighbor.comp.sig
neighbor.comp$significance
```
![](img/Brain.first.degree.png)
![](img/Brain.first.degree.sig.png)

```
$Astro
        Astro          Endo       L2/3 IT           Vip 
2.179335e-245  3.020038e-68  5.175522e-20  3.232111e-47 

$Endo
        Astro          Endo       L2/3 IT           Vip 
 9.082251e-48 1.499351e-204  1.916450e-02  4.714201e-10 

$`L2/3 IT`
       Astro      L2/3 IT           L4        Lamp5         Sncg          Vip 
3.894444e-11 0.000000e+00 1.186417e-10 5.794239e-04 7.579714e-14 5.953635e-28 

$L4
     L2/3 IT           L4        L5 IT 
7.940044e-29 0.000000e+00 1.504830e-81 

$`L5 IT`
           L4         L5 IT         L5 PT           Sst 
 1.141931e-78 2.030590e-248  2.883319e-83  7.886836e-13 

$`L5 PT`
           L4         L5 IT         L5 PT         Pvalb           Sst 
 2.303106e-02  4.908108e-89 1.772265e-228  4.370237e-20  1.243123e-12 

$`L6 CT`
       L6 CT        L6 IT          L6b        Pvalb 
0.000000e+00 1.378626e-47 1.159740e-15 7.248710e-43 

$`L6 IT`
       L6 CT        L6 IT          L6b 
4.922746e-37 0.000000e+00 6.988967e-22 

$L6b
        Astro         L6 CT         L6 IT           L6b           Vip 
 7.002880e-03  3.995450e-24  2.264252e-05 3.157963e-265  5.500931e-06 

$Lamp5
   L2/3 IT      Lamp5       Sncg 
3.2255e-07 0.0000e+00 2.4025e-15 

$NP
NP 
 0 

$Pvalb
       L5 IT        L5 PT        L6 CT        L6 IT        Pvalb 
6.860777e-04 1.434897e-20 8.129352e-40 3.895256e-04 0.000000e+00 

$Sncg
        Astro       L2/3 IT           L6b         Lamp5          Sncg           Vip 
 7.999573e-11  7.587078e-16  4.789452e-02  9.099911e-26 3.767464e-113  2.104632e-17 

$Sst
       L5 IT        L5 PT           NP          Sst 
7.141141e-09 1.587150e-14 9.446149e-09 0.000000e+00 

$Vip
        Astro          Endo           L6b          Sncg           Vip 
3.857367e-113  3.341504e-17  2.485397e-03  5.728714e-20  0.000000e+00 

```

## 2. SMN Distance

Based on the reconstructed graph, we can investigate spatial relationship among cell clusters of interest by calculating their SMN distances. Here we demonstrate a boxplot for distances between the L2/3 IT layer to all glutamatergic cells.

```r
clusterDistanceBoxplot(brain.obj,annotationToUse = "cell_type",query.cluster.list = c( "L2/3 IT", "L4","L5 IT","L5 PT","NP", "L6 IT", "L6 CT","L6b" ),reference.cluster="L2/3 IT")
```
![](img/Brain.clusterDist.png)




## X. Running Time
```r
system.time(brain.obj<-trainHolography(sp.integrated,n.repeat = 30))
```
```
user  system elapsed 
830.224 172.764 901.495 
```








