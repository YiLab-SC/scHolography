library(scHolography)
library(ggplot2)
library(dplyr)

scHolo.obj <- readRDS("/projects/b1042/YiLab/Cady/Code_3rdDraft/Mouse_Hippo_sim/hippo.vizgen.visium.rds")
real.coord <- scHolo.obj$scHolography.sc@images$s2r1@boundaries$centroids@coords
data <- data.frame(x = real.coord[,1], y = real.coord[,2])



n=30
# Create the plot
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point() +theme_minimal()
# Add vertical lines
x_breaks <- seq(min(data$x), max(data$x), length.out = (n+1))
p <- p + geom_vline(xintercept = x_breaks, linetype = "dashed", color = "red")

# Add horizontal lines
y_breaks <- seq(min(data$y), max(data$y), length.out = (n+1))
p <- p + geom_hline(yintercept = y_breaks, linetype = "dashed", color = "blue")

tiff("/projects/b1042/YiLab/Cady/Code_3rdDraft/Mouse_Hippo_sim/density.plot.tiff", units="in", width=8, height=6, res=300)
p +theme(axis.text.x = element_blank(),axis.text.y = element_blank())
dev.off()


xbin <- cbind(x_breaks[1:n],x_breaks[2:(n+1)])
ybin <- cbind(y_breaks[1:n],y_breaks[2:(n+1)])
count.mat <- matrix(0, ncol = n, nrow = n)
cor.mat <- matrix(0, ncol = n, nrow = n)
kl.mat <- matrix(0, ncol = n, nrow = n)

real<- as.matrix(dist(real.coord))
graph <- igraph::graph_from_adjacency_matrix(scHolo.obj$adj.mtx, mode = "undirected")
dist <- igraph::distances(graph, mode = "out")
colnames(dist) = as.character(1:ncol(dist))
rownames(dist) = as.character(1:ncol(dist))

findKL<-function(realmat,mat){
  realmat.sub <- realmat[rownames(mat),colnames(mat)]
  
  out <- lapply(rownames(mat),function(nam){
    eps <- 1e-20
    X_ <- cbind(realmat.sub[nam,],mat[nam,])
    X_[X_ < eps] <- eps
    X_ <- scale(X_, center = F, scale = colSums(X_)) %>% data.frame
    KL_D <- suppressMessages(philentropy::KL(t(X_), unit = "log"))
    KL_D
  })
  out <- unlist(out)
  names(out) <- rownames(mat)
  out
}


scHolo.kl<-findKL(realmat =real ,mat = dist)

for (i in 1:n) {
  if (i==1){
    xind= intersect(which(data$x>=xbin[i,1]),which(data$x<=xbin[i,2]))
  }else{
    xind= intersect(which(data$x>xbin[i,1]),which(data$x<=xbin[i,2]))
  }
  for (j in 1:n) {
    if(j==1){
      yind= intersect(which(data$y>=ybin[j,1]),which(data$y<=ybin[j,2]))
    }else{
      yind= intersect(which(data$y>ybin[j,1]),which(data$y<=ybin[j,2]))
    }
    count.mat[i,j]=length(intersect(xind, yind))
    cor.mat[i,j]=mean(unlist(lapply(intersect(xind, yind),
                                    function(cell) cor(real[cell,],dist[cell,],method="spearman") )))
    kl.mat[i,j] = mean(scHolo.kl[intersect(xind, yind)])
  }
  
}

count.mat.melt <-reshape2::melt(count.mat)
corr.mat.melt <-reshape2::melt(cor.mat)
kl.mat.melt <-reshape2::melt(kl.mat)

tiff("/projects/b1042/YiLab/Cady/Code_3rdDraft/Mouse_Hippo_sim/binCount.plot.tiff", units="in", width=8, height=6, res=300)
ggplot(count.mat.melt, aes(Var1, Var2)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = value))+viridis::scale_fill_viridis()+               
  labs(x = "x", y = "y", fill = "Bin Count") + theme_classic() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())
dev.off()


tiff("/projects/b1042/YiLab/Cady/Code_3rdDraft/Mouse_Hippo_sim/binCorr.plot.tiff", units="in", width=8, height=6, res=300)
ggplot(corr.mat.melt, aes(Var1, Var2)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = value))+viridis::scale_fill_viridis()+               
  labs(x = "x", y = "y", fill = "Bin Correlation") + theme_classic() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())
dev.off()

tiff("/projects/b1042/YiLab/Cady/Code_3rdDraft/Mouse_Hippo_sim/binKL.plot.tiff", units="in", width=8, height=6, res=300)
ggplot(kl.mat.melt, aes(Var1, Var2)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = -value))+viridis::scale_fill_viridis()+               
  labs(x = "x", y = "y", fill = "Neg. Bin K-L") + theme_classic() +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())
dev.off()


