library(scHolography)
library(tidyverse)
setwd("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human Skin/")

scHolography.obj <- readRDS("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human Skin/scHolography.obj_new.rds")
high.res.sp<-readRDS("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Spatial Transcriptomics/VAE/f2.rds")
low.res.sp<-Seurat::Load10X_Spatial("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Data/human/st_Foreskin/20220712/A1/")
sp.integrated <- dataAlign(low.res.sp,high.res.sp,nPCtoUse = 32,scProcessed = T)


out <- scHolography:::getData(sp.integrated, nPCtoOut=32, is.fov = F, fov=NULL) 
low.res.sp[["SIC1"]] <- out[[3]][,1]
low.res.sp[["SIC2"]] <- out[[3]][,2]
low.res.sp[["SIC3"]] <- out[[3]][,3]
library(RColorBrewer)
tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human Skin/SIC1_sp.tiff", units="in", width=10, height=8, res=300)
Seurat::SpatialFeaturePlot(low.res.sp,"SIC1")+ggplot2::scale_fill_gradientn(colors=rev(brewer.pal(11, "RdBu")))
dev.off()

tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human Skin/SIC2_sp.tiff", units="in", width=10, height=8, res=300)
Seurat::SpatialFeaturePlot(low.res.sp,"SIC2")+ggplot2::scale_fill_gradientn(colors=rev(brewer.pal(11, "RdBu")))
dev.off()

tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human Skin/SIC3_sp.tiff", units="in", width=10, height=8, res=300)
Seurat::SpatialFeaturePlot(low.res.sp,"SIC3")+ggplot2::scale_fill_gradientn(colors=rev(brewer.pal(11, "RdBu")))
dev.off()
  
scene = list(camera = list(eye = list(x = -1, z = 0, y = 2)))


SIC1 <- scHolography.obj$raw.prediction[,1]
names(SIC1) <- colnames(high.res.sp)
scHolography.obj$scHolography.sc[["pre_SIC1"]] <- SIC1[colnames(scHolography.obj$scHolography.sc)]
fig1 <- scHolographyPlot(scHolography.obj,feature = "pre_SIC1",feature.pal = "rdbu")%>% plotly::layout(scene = scene)
plotly::orca(fig1, "SIC1_sc.svg",width = 7,height = 5)
fig1

SIC2 <- scHolography.obj$raw.prediction[,2]
names(SIC2) <- colnames(high.res.sp)
scHolography.obj$scHolography.sc[["pre_SIC2"]] <- SIC2[colnames(scHolography.obj$scHolography.sc)]
fig2 <- scHolographyPlot(scHolography.obj,feature = "pre_SIC2",feature.pal = "rdbu")%>% plotly::layout(scene = scene)
plotly::orca(fig2, "SIC2_sc.svg",width = 7,height = 5)
fig2

SIC3 <- scHolography.obj$raw.prediction[,3]
names(SIC3) <- colnames(high.res.sp)
scHolography.obj$scHolography.sc[["pre_SIC3"]] <- SIC3[colnames(scHolography.obj$scHolography.sc)]
fig3 <- scHolographyPlot(scHolography.obj,feature = "pre_SIC3",feature.pal = "rdbu")%>% plotly::layout(scene = scene)
plotly::orca(fig3, "SIC3_sc.svg",width = 7,height = 5)
fig3


