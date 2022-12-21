#' Layer Composition Plot from Given Celltypes
#' @export
#' @import Seurat
#' @import RColorBrewer
#' @param  scHolography.obj scHolography object list
#' @param  annotationToUse Which annotation to call identities from. Default is orig.cluster
#' @param  query.cluster A vector of query identity types
#' @param  palette Color palette to use for coloring. Default is Paired
neighborCompositionPlot<-function(scHolography.obj, annotationToUse="orig.cluster", query.cluster,palette = "Paired"){
  scHolography.sc<-scHolography.obj$scHolography.sc
  adj.mtx <- scHolography.obj$adj.mtx
  if( is.null(levels(scHolography.sc[[annotationToUse]][[1]]))){
    scHolography.sc[[annotationToUse]][[1]]<-factor(scHolography.sc[[annotationToUse]][[1]],levels = stringr::str_sort(unique(scHolography.sc[[annotationToUse]][[1]]),numeric = T))
  }
  graph <- igraph::graph_from_adjacency_matrix(adj.mtx,mode = "undirected")

  dist <- igraph::distances(graph, mode="out")
  anno<-scHolography.sc[[annotationToUse]][[1]]
  cell.ind <- which(anno%in%query.cluster)
  sum.mat<-matrix(0, ncol =length(levels(anno)),nrow =  (1+max(dist)))
  rownames(sum.mat)<-as.character(0:max(dist))
  colnames(sum.mat)<-levels(anno)
  for (i in cell.ind) {
    sig.mtx <- as.matrix(table(dist[i,],anno))
    sig.mtx <- sig.mtx[,stringr::str_sort(colnames(sig.mtx),numeric = T)]
    sig.mtx.perc <- sig.mtx/rowSums(sig.mtx)
    m3<-sum.mat
    mcol<-match(rownames(sig.mtx.perc),rownames(sum.mat))
    m3[mcol,]<-m3[mcol,]+sig.mtx.perc
    sum.mat <- m3

  }
  sum.mat.perc <- sum.mat/rowSums(sum.mat)
  mycolors <- colorRampPalette(brewer.pal(12, palette))(length(levels(anno)))
  df <- reshape2::melt(sum.mat.perc)
  colnames(df) <- c("Layer","Identity","Perc")
  df$Layer <- as.factor(df$Layer)
  df$Identity <- as.factor(df$Identity)
  plot<-ggplot(df, aes(fill=Identity, y=Perc, x=Layer)) +
    geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = mycolors)+
    ggplot2::theme_classic()

  plot

}



#' Cluster Distance Boxplot
#' @export
#' @import ggplot2
#' @import igraph
#' @import RColorBrewer
#' @import plotly
#' @param  scHolography.obj scHolography object list
#' @param  annotationToUse Which annotation to call identities from. Default is orig.cluster
#' @param  query.cluster.list A vector of query identity types to be measured SMN distance to reference
#' @param  reference.cluster A vector of reference identity types
#' @param  palette Color palette to use for coloring. Default is Paired
#' @param  n.neighbor Number of nearest cells to use to define distance. Default is 30
clusterDistanceBoxplot <- function(scHolography.obj,annotationToUse = "orig.cluster",query.cluster.list,
                                   reference.cluster,palette = "Paired", n.neighbor = 30){
  scHolography.sc<-scHolography.obj$scHolography.sc
  adj.mtx <- scHolography.obj$adj.mtx
  if( is.null(levels(scHolography.sc[[annotationToUse]][[1]]))){
    scHolography.sc[[annotationToUse]][[1]]<-factor(scHolography.sc[[annotationToUse]][[1]],levels = stringr::str_sort(unique(scHolography.sc[[annotationToUse]][[1]]),numeric = T))
  }
  graph <- igraph::graph_from_adjacency_matrix(adj.mtx,mode = "undirected")

  dist <- igraph::distances(graph, mode="out")
  anno<-scHolography.sc[[annotationToUse]][[1]]
  ref.ind <- which(anno%in%reference.cluster)
  dist.list <- lapply(query.cluster.list, function(query.cluster){
    query.ind <- which(anno%in%query.cluster)
    dist.mtx <- dist[ref.ind,query.ind]
    colMeans(apply(dist.mtx,2,sort)[1:n.neighbor,])
  })
  cell_count <- unlist(lapply(dist.list, function(x) length(x)))
  dist.dat <- data.frame(Distance = unlist(dist.list), Celltype = unlist(lapply(1:length(cell_count), function(x) rep(query.cluster.list[x],cell_count[x]))))
  dist.dat$Celltype <- factor(dist.dat$Celltype, levels = query.cluster.list)

  my.color <- colorRampPalette(brewer.pal(12, palette))(length(levels(scHolography.sc[[annotationToUse]][[1]])))
  color.used <- unlist(lapply(query.cluster.list, function(x){which(levels(scHolography.sc[[annotationToUse]][[1]]) %in%x)}))
  fig <- ggplot(dist.dat, aes(x=Celltype, y=Distance,fill=Celltype))+geom_boxplot()+theme_classic()+scale_fill_manual(values=my.color[color.used])
  fig
}



