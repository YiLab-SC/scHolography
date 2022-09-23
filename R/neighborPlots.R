#' Layer Composition Plot from Given Celltypes
#' @export
#' @import Seurat
#' @import RColorBrewer

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


#' Find Neighbor Determining Markers
#' @export
#' @import Seurat
#' @import RColorBrewer
#' @import dplyr


neighborMarkerPlot <- function(scHolography.obj, annotationToUse="orig.cluster", query.cluster, target.cluster, palette = "Paired",assayToUse="SCT",cutoff=0){
  library(dplyr)
  scHolography.sc<-scHolography.obj$scHolography.sc
  adj.mtx <- scHolography.obj$adj.mtx
  if(is.null(target.cluster)){
    target.cluster=levels(scHolography.sc[[annotationToUse]][[1]])
  }
  anno<-as.character(scHolography.sc[[annotationToUse]][[1]])
  target.cluster <- stringr::str_sort(target.cluster,numeric = T)
  query.cluster <-stringr::str_sort(query.cluster,numeric = T)
  cell.ind <- which(anno%in%query.cluster)
  dist.sub <-adj.mtx[cell.ind,]
  ind.ls<-lapply(target.cluster, function(x){
    tg.cell.ind <- which(anno%in%x)
    which(rowSums(dist.sub[,tg.cell.ind])>cutoff)
  })
  sub.obj.ls <- lapply(ind.ls, function(x){
    subset(scHolography.sc,cells = cell.ind[x])
  })
  sub.obj.ls <- lapply(1:length(target.cluster), function(x){
    sub.obj.ls[[x]][["target.clus"]]<- target.cluster[[x]]
    sub.obj.ls[[x]]
  })

  names(sub.obj.ls)<-paste("c",1:length(target.cluster),sep = "_")
  combined <- merge (
    x = sub.obj.ls[[1]],
    y = within(sub.obj.ls, rm(c_1)),
    add.cell.ids = target.cluster
  )
  combined<-SetIdent(combined,value = factor(combined$target.clus,levels = target.cluster))
  DefaultAssay(combined) <- "RNA"
  combined <- NormalizeData(combined)
  all.genes <- rownames(combined)
  combined <- ScaleData(combined, features = all.genes)
  markers<-FindAllMarkers(combined,assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  num.clus <- (levels(scHolography.sc[[annotationToUse]][[1]]))
  mycolors <- colorRampPalette(brewer.pal(12, palette))(length(num.clus))[which(num.clus%in%levels(combined@active.ident))]

  heatmap<-DoHeatmap(combined,group.colors = mycolors,features = top10$gene,assay = assayToUse)+ggplot2::scale_fill_gradientn(colors = c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"))
  show(heatmap)
  list(heatmap=heatmap, DEG=top10)
}





#' Neighbor Sankey Plot
#' @export
#' @import ggplot2
#' @import plyr
#' @import RColorBrewer
#' @import plotly
neighborSankeyPlot<- function(scHolography.obj, annotationToUse="orig.cluster",query.cluster,palette="Paired",font.size=12){
  scHolography.sc<-scHolography.obj$scHolography.sc
  adj.mtx <- scHolography.obj$adj.mtx
  if( is.null(levels(scHolography.sc[[annotationToUse]][[1]]))){
    scHolography.sc[[annotationToUse]][[1]]<-factor(scHolography.sc[[annotationToUse]][[1]],levels = stringr::str_sort(unique(scHolography.sc[[annotationToUse]][[1]]),numeric = T))
  }
  graph <- igraph::graph_from_adjacency_matrix(adj.mtx,mode = "undirected")

  dist <- igraph::distances(graph, mode="out")


  anno<-as.character(scHolography.sc[[annotationToUse]][[1]])
  cell.ind <- which(anno%in%query.cluster)
  nb.ls <- lapply(cell.ind,function(x){
    which(dist[x,]==1)
  })

  Sankey.ls<-lapply(1:length(cell.ind), function(x){
    cbind(rep(cell.ind[x],length(nb.ls[[x]])),nb.ls[[x]])

  })

  Sankey.mat<-do.call(rbind, Sankey.ls)

  Sankey.dat <- data.frame(x=as.factor(rep("source",nrow(Sankey.mat))),next_x=as.factor(rep("target",nrow(Sankey.mat)) ))
  Sankey.dat$node <- factor(paste("s",anno[Sankey.mat[,1]],sep = "_"),levels = paste("s",levels(scHolography.sc[[annotationToUse]][[1]]),sep = "_"))
  Sankey.dat$next_node <- factor(paste("t",anno[Sankey.mat[,2]],sep = "_"),levels = paste("t",levels(scHolography.sc[[annotationToUse]][[1]]),sep = "_"))

  Sankey.net <-plyr::ddply(Sankey.dat,.(node, next_node), nrow)
  nodes <- data.frame(
    name=c(as.character(Sankey.net$node),
           as.character(Sankey.net$next_node)) %>% unique()
  )

  nodes.clus<-unlist(lapply(nodes$name, function(x) strsplit(x,split = "_")[[1]][2]))
  mycolors <- colorRampPalette(brewer.pal(12, palette))(length(levels(scHolography.sc[[annotationToUse]][[1]])))


  fig <- plot_ly(
    type = "sankey",
    orientation = "h",node = list(
      label = nodes$name,
      color = mycolors[unlist(lapply(nodes.clus,function(x) which(levels(scHolography.sc[[annotationToUse]][[1]])%in%x)))],
      pad = 15,
      thickness = 20,
      line = list(
        color = "white"
      )
    ),

    link = list(
      source = unlist(lapply(as.vector(Sankey.net$node), function(x) which(nodes$name%in%x)))-1,
      target = unlist(lapply(as.vector(Sankey.net$next_node), function(x) which(nodes$name%in%x)))-1,
      value =  as.vector(Sankey.net$V1),
      color = mycolors[unlist(lapply(unlist(lapply(as.character(Sankey.net$node), function(x) strsplit(x,split = "_")[[1]][2])),function(x) which(levels(scHolography.sc[[annotationToUse]][[1]])%in%x)))]
    )
  )
  fig <- fig %>% layout(
    title = "Stable Neighbor Sankey Diagram",
    font = list(
      size = font.size
    )
  )

  fig
}





