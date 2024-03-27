load("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/Enriched pathways for SPneighborhoods_part2.RData")


genes_Squamous  = unlist(strsplit(a@result[1,"geneID" ],"/"))
genes_Epithelial = unlist(strsplit(a@result[2,"geneID" ],"/"))
genes_Stress = unlist(strsplit(b@result[1,"geneID" ],"/"))
genes_EMT = unlist(strsplit(b@result[2,"geneID" ],"/"))
genes_Cycle = unlist(strsplit(c@result[1,"geneID" ],"/"))
genes_Interferon = unlist(strsplit(c@result[4,"geneID" ],"/"))




# use gsva in Seurat object scHolography for each cluster
library(scHolography)
library(Seurat)
library(ggplot2)
SCC.obj <- readRDS("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/scholography.obj.rds")
scc.sp.neighbor.tumor_new <- readRDS("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/scc.sp.neighbor.tumor_new.rds")
P6<-readRDS("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Spatial Transcriptomics/PK Data//P6.rds")
scSCC.meta <- read.table("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Data/human/sc_squamous cell carcinoma/GSE144236_patient_metadata_new.txt",sep="\t")

celltype <- names(table(P6$level2_celltype))
P6.sub <- subset(P6, cells = which(P6$level2_celltype%in%celltype[-which(celltype%in%c("Keratinocyte","Multiplet"))]))
P6<-P6.sub 
P6$level2_celltype<-factor(P6$level2_celltype)

P6@active.ident<-P6$level2_celltype


P6$my.level <- as.character(P6$level2_celltype)
P6$my.level[which(P6$my.level%in%c("Mac","ASDC","CD1C","CLEC9A","LC","MDSC","PDC"))] <- "Myeloid"





og.color <- colorRampPalette( RColorBrewer::brewer.pal(12,"Set1"))(length(unique(P6$my.level)))
my.col <- (og.color)[-c(2,5)]

kc1.neighbor <- which(colSums(scc.sp.neighbor.tumor_new$scHolography.obj$adj.mtx[which(scc.sp.neighbor.tumor_new$scHolography.obj$scHolography.sc$spatial.neighborhood%in%"1"),])>0)
kc2.neighbor <- which(colSums(scc.sp.neighbor.tumor_new$scHolography.obj$adj.mtx[which(scc.sp.neighbor.tumor_new$scHolography.obj$scHolography.sc$spatial.neighborhood%in%"2"),])>0)
kc3.neighbor <- which(colSums(scc.sp.neighbor.tumor_new$scHolography.obj$adj.mtx[which(scc.sp.neighbor.tumor_new$scHolography.obj$scHolography.sc$spatial.neighborhood%in%"3"),])>0)

DefaultAssay(SCC.obj$scHolography.sc) <- "RNA"
SCC.obj$scHolography.sc <- NormalizeData(SCC.obj$scHolography.sc)

SCC.obj$scHolography.sc$my.level <- SCC.obj$scHolography.sc$level2_celltype
SCC.obj$scHolography.sc$my.level[which(SCC.obj$scHolography.sc$my.level%in%c("Mac","ASDC","CD1C","CLEC9A","LC","MDSC","PDC"))] <- "Myeloid"

seurat.sub.kc1 <- subset(SCC.obj$scHolography.sc, cells = kc1.neighbor)
seurat.sub.kc2 <- subset(SCC.obj$scHolography.sc, cells = kc2.neighbor)
seurat.sub.kc3 <- subset(SCC.obj$scHolography.sc, cells = kc3.neighbor)

seurat.sub.kc1 <- AddModuleScore(seurat.sub.kc1, features  = list(genes_Squamous), name = "Squamous", assay = "RNA")
tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/seurat.sub.kc1_Squamous1.tiff", units="in", width=8, height=6, res=300)
VlnPlot(seurat.sub.kc1,features = "Squamous1",group.by = "my.level",cols = my.col[sort(unique(SCC.obj$scHolography.sc$my.level))%in%c(
  "Fibroblast","Myeloid",
  "Normal_KC_Basal","Normal_KC_Cyc",
  "Normal_KC_Diff","Pilosebaceous","Tcell",
  "TSK","Tumor_KC_Basal","Tumor_KC_Cyc",
  "Tumor_KC_Diff")])+ 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 15,face = "bold"),axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  ggtitle("Tumor KC1 -- Squamous (Barkley 2022)")
dev.off()

seurat.sub.kc1 <- AddModuleScore(seurat.sub.kc1, features  = list(genes_Epithelial), name = "Epithelial", assay = "RNA")
tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/seurat.sub.kc1_Epithelial1.tiff", units="in", width=8, height=6, res=300)
VlnPlot(seurat.sub.kc1,features = "Epithelial1",group.by = "my.level",cols = my.col[sort(unique(SCC.obj$scHolography.sc$my.level))%in%c(
  "Fibroblast","Myeloid",
  "Normal_KC_Basal","Normal_KC_Cyc",
  "Normal_KC_Diff","Pilosebaceous","Tcell",
  "TSK","Tumor_KC_Basal","Tumor_KC_Cyc",
  "Tumor_KC_Diff")])+ 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 15,face = "bold"),axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  ggtitle("Tumor KC1 -- Epithelial (Gavish 2023)")
dev.off()

seurat.sub.kc2 <- AddModuleScore(seurat.sub.kc2, features  = list(genes_Stress), name = "Stress", assay = "RNA")
tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/seurat.sub.kc2_Stress1.tiff", units="in", width=8, height=6, res=300)
VlnPlot(seurat.sub.kc2,features = "Stress1",group.by = "my.level",cols = my.col[sort(unique(SCC.obj$scHolography.sc$my.level))%in%c("Endothelial Cell",
                                                                                                                                    "Fibroblast","Myeloid",
                                                                                                                                    "Normal_KC_Basal","Normal_KC_Cyc",
                                                                                                                                    "Normal_KC_Diff","Pilosebaceous","Tcell",
                                                                                                                                    "TSK","Tumor_KC_Basal","Tumor_KC_Cyc",
                                                                                                                                    "Tumor_KC_Diff")])+ 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 15,face = "bold"),axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  ggtitle("Tumor KC2 -- Stress (Barkley 2022)")
dev.off()

seurat.sub.kc2 <- AddModuleScore(seurat.sub.kc2, features  = list(genes_EMT), name = "EMT", assay = "RNA")
tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/seurat.sub.kc2_EMT1.tiff", units="in", width=8, height=6, res=300)
VlnPlot(seurat.sub.kc2,features = "EMT1",group.by = "my.level",cols = my.col[sort(unique(SCC.obj$scHolography.sc$my.level))%in%c("Endothelial Cell",
                                                                                                                                 "Fibroblast","Myeloid",
                                                                                                                                 "Normal_KC_Basal","Normal_KC_Cyc",
                                                                                                                                 "Normal_KC_Diff","Pilosebaceous","Tcell",
                                                                                                                                 "TSK","Tumor_KC_Basal","Tumor_KC_Cyc",
                                                                                                                                 "Tumor_KC_Diff")])+ 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 15,face = "bold"),axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  ggtitle("Tumor KC2 -- EMT IV (Gavish 2023)")
dev.off()

seurat.sub.kc3 <- AddModuleScore(seurat.sub.kc3, features  = list(genes_Cycle), name = "Cycle", assay = "RNA")
tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/seurat.sub.kc3_Cycle1.tiff", units="in", width=8, height=6, res=300)
VlnPlot(seurat.sub.kc3,features = "Cycle1",group.by = "my.level",cols = my.col[sort(unique(SCC.obj$scHolography.sc$my.level))%in%c(
  "Fibroblast","Myeloid","NK", 
  "Normal_KC_Basal","Normal_KC_Cyc",
  "Normal_KC_Diff","Pilosebaceous","Tcell",
  "TSK","Tumor_KC_Basal","Tumor_KC_Cyc",
  "Tumor_KC_Diff")])+ 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 15,face = "bold"),axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  ggtitle("Tumor KC3 -- Cycle (Barkley 2022)")
dev.off()

seurat.sub.kc3 <- AddModuleScore(seurat.sub.kc3, features  = list(genes_Interferon), name = "Interferon", assay = "RNA")
tiff("~/OneDrive - Northwestern University/Northwestern/Yi Lab/Manuscript/Code_2ndDraft/Human SCC/seurat.sub.kc3_Interferon1.tiff", units="in", width=8, height=6, res=300)
VlnPlot(seurat.sub.kc3,features = "Interferon1",group.by = "my.level",cols = my.col[sort(unique(SCC.obj$scHolography.sc$my.level))%in%c(
  "Fibroblast","Myeloid","NK", 
  "Normal_KC_Basal","Normal_KC_Cyc",
  "Normal_KC_Diff","Pilosebaceous","Tcell",
  "TSK","Tumor_KC_Basal","Tumor_KC_Cyc",
  "Tumor_KC_Diff")])+ 
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 15,face = "bold"),axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  ggtitle("Tumor KC3 -- Interferon (Barkley 2022)")
dev.off()




