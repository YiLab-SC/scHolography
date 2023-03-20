
#' scHolography 3D Tissue Reconstruction Visualization
#' @export
#' @import RColorBrewer
#' @import plotly
#' @import Seurat
#' @param  scHolography.obj scHolography object list
#' @param  cells Which cells are used for plotting. Default is NULL (all single cells from scHolography.obj will be plotted)
#' @param  feature Which feature is used for plotting. Default is NULL, where dots will be colored by annotation. When a feature from expression or metadata (e.g., a gene name) is given, dots will be colored by the intensity of the feature.
#' @param  cutoff Cutoff for learning variance of prediction. Only cells with learning variance no greater than the cutoff will be plotted. Default is NULL (all cells will be plotted)
#' @param  color.by Which annotation to use for coloring the dot plot. Default is orig.cluster
#' @param  dot.size Dot size. Default is 5
#' @param  color A vector of user defined color to use for coloring. Default is NA (pre-fix color palette will be used)
#' @param  palette Color palette to use for coloring. Default is Paired
#' @param  highlight Index of cells to highlight on the plot. Default is NULL (No cells will be highlighted)
#' @param  assayToUse Which assay to get feature from. Default is SCT
#' @param  feature.pal Color palette to use for feature plot coloring. Default is viridis. Other options are "rdbu", "magma", and "rdbu_1"
#' @param  show.grid If to show the grid. Default is False
#' @param  g.width Width of the grid. Default is 5
#' @param  dim Number of space dimensions to visualize the reconstructed SMN graph

scHolographyPlot<-function(scHolography.obj, cells=NULL,feature=NULL,cutoff=NULL, color.by="orig.cluster", dot.size = 5 , color=NA,palette="Paired", highlight=NULL,assayToUse="SCT",feature.pal="viridis",show.grid=F,g.width=5, dim=3){
  scHolography.sc<-scHolography.obj$scHolography.sc
  if(is.null(levels(scHolography.sc[[color.by]][[1]]))){
    scHolography.sc[[color.by]][[1]]<-as.factor(scHolography.sc[[color.by]][[1]])
  }
  if(sum(is.na(color))>0){
    getPalette =colorRampPalette(brewer.pal(12, palette))
    my.scheme = getPalette(length(levels(scHolography.sc[[color.by]][[1]])))
  }else{
    my.scheme =color
  }

  if (feature.pal == "rdbu") {
    fea.col <- rev(RColorBrewer::brewer.pal(n = 25, name = "RdBu"))
  }else if (feature.pal == "magma") {
    fea.col <- viridis::viridis(25, option = "A")
  }else if (feature.pal == "rdbu_1") {
    fea.col <- colorspace::diverge_hsv(25)
  }else {
    fea.col <- viridis::viridis(25)
  }

  if(is.null(cells)){
    cellToPlot=1:ncol(scHolography.sc)
  }else if(is.numeric(cells)){
    cellToPlot=cells
  }else{
    cellToPlot<-which(scHolography.sc[[color.by]][[1]]%in%cells)
  }
  if(is.null(cutoff)==F){
    cellToPlot<-intersect(cellToPlot,which(scHolography.sc$motility<=cutoff))
  }

  if(dim==2){
    if(is.null(highlight)){
      if(is.null(feature)){
        out.plot<-plotly::plot_ly(x=scHolography.sc$x_sp[cellToPlot], y=scHolography.sc$y_sp[cellToPlot], type="scatter", mode="markers", color=scHolography.sc[[color.by]][[1]][cellToPlot],colors = my.scheme,marker=list(size=dot.size))
      }else{
        if(feature%in%rownames(Seurat::GetAssayData(scHolography.sc,assay =as.character(assayToUse)))){
          out.plot<-plotly::plot_ly( x=scHolography.sc$x_sp[cellToPlot], y=scHolography.sc$y_sp[cellToPlot], colors=c(fea.col),type="scatter", mode="markers", color=as.vector(Seurat::GetAssayData(scHolography.sc,assay = assayToUse)[feature,])[cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

        }else{
          out.plot<-plotly::plot_ly( x=scHolography.sc$x_sp[cellToPlot], y=scHolography.sc$y_sp[cellToPlot], colors=c(fea.col),type="scatter", mode="markers", color=scHolography.sc[[feature]][[1]][cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

        }
      }
    }else{
      new.scheme<-c(my.scheme,"gray80")
      if(is.numeric(highlight)){
        cellToHighlight=highlight
      }else{
        cellToHighlight<-which(scHolography.sc[[color.by]][[1]]%in%highlight)
      }
      color.highlight <- scHolography.sc[[color.by]][[1]]
      levels(color.highlight)<-c(levels(color.highlight),"")
      color.highlight[which((1:length(color.highlight)%in%cellToHighlight)==F)] <- as.factor("")
      out.plot<-plotly::plot_ly(x=scHolography.sc$x_sp[cellToPlot], y=scHolography.sc$y_sp[cellToPlot], type="scatter", mode="markers", color=color.highlight[cellToPlot],colors = new.scheme,marker=list(size=dot.size))
    }
  }else if(dim==3){
    if(is.null(highlight)){
      if(is.null(feature)){
        out.plot<-plotly::plot_ly(x=scHolography.sc$x3d_sp[cellToPlot], y=scHolography.sc$y3d_sp[cellToPlot], z=scHolography.sc$z3d_sp[cellToPlot],  type="scatter3d", mode="markers", color=scHolography.sc[[color.by]][[1]][cellToPlot],colors = my.scheme,marker=list(size=dot.size))%>%plotly::layout(scene=list(xaxis = list(showticklabels = show.grid,gridwidth=g.width), yaxis = list( showticklabels = show.grid,gridwidth=g.width),  zaxis = list(showticklabels = show.grid,gridwidth=g.width)))
      }else{
        if(feature%in%rownames(Seurat::GetAssayData(scHolography.sc,assay =as.character(assayToUse)))){
          out.plot<-plotly::plot_ly( x=scHolography.sc$x3d_sp[cellToPlot], y=scHolography.sc$y3d_sp[cellToPlot], z=scHolography.sc$z3d_sp[cellToPlot], colors=c(fea.col),type="scatter3d", mode="markers", color=as.vector(Seurat::GetAssayData(scHolography.sc,assay = assayToUse)[feature,])[cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature),scene=list(xaxis = list(showticklabels = show.grid,gridwidth=g.width), yaxis = list( showticklabels = show.grid,gridwidth=g.width),  zaxis = list(showticklabels = show.grid,gridwidth=g.width)))

        }else{
          out.plot<-plotly::plot_ly( x=scHolography.sc$x3d_sp[cellToPlot], y=scHolography.sc$y3d_sp[cellToPlot],z=scHolography.sc$z3d_sp[cellToPlot],  colors=c(fea.col),type="scatter3d", mode="markers", color=scHolography.sc[[feature]][[1]][cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature),scene = list(xaxis = list(showticklabels = show.grid,gridwidth=g.width), yaxis = list( showticklabels = show.grid,gridwidth=g.width),  zaxis = list(showticklabels = show.grid,gridwidth=g.width)))

        }
      }

    }else{
      new.scheme<-c(my.scheme,"gray80")
      if(is.numeric(highlight)){
        cellToHighlight=highlight
      }else{
        cellToHighlight<-which(scHolography.sc[[color.by]][[1]]%in%highlight)
      }
      color.highlight <- scHolography.sc[[color.by]][[1]]
      levels(color.highlight)<-c(levels(color.highlight),"")
      color.highlight[which((1:length(color.highlight)%in%cellToHighlight)==F)] <- as.factor("")
      out.plot<-plotly::plot_ly(x=scHolography.sc$x3d_sp[cellToPlot], y=scHolography.sc$y3d_sp[cellToPlot], z=scHolography.sc$z3d_sp[cellToPlot], type="scatter3d", mode="markers",color=color.highlight[cellToPlot],colors = new.scheme,marker=list(size=dot.size))
    }

  }else{
    print("ERROR: Wrong dimension input")
  }
  out.plot
}



#' Spatial Neighborhood Composition Plot
#' @export
#' @import RColorBrewer
#' @import ggplot2
#' @import stringr
#' @param  scHolography.obj scHolography object list
#' @param  annotationToUse Which annotation to call identities from. Default is orig.cluster
#' @param  query.cluster A vector of query identity types
#' @param  pal Color palette to use for coloring. Default is Paired
#' @param  color A vector of user defined color to use for coloring. Default is NA (pre-fix color palette will be used)
#' @param  row.group.level Level order of annotation to plot. Default is NULL (alphabetical order will be adopted)

scHolographyNeighborCompPlot <- function (scHolography.obj, annotationToUse, query.cluster = NULL,
                                          pal = "Paired", color = NULL, row.group.level = NULL)
{
  clus <- scHolography.obj$scHolography.sc[[annotationToUse]][[1]]
  if (is.null(query.cluster)) {
    query.cluster <- stringr::str_sort(unique(clus), numeric = T)
  }
  matrix <- matrix(ncol = length(unique(clus)), nrow = length(query.cluster),
                   0)
  colnames(matrix) <- stringr::str_sort(unique(clus), numeric = T)
  rownames(matrix) <- query.cluster
  uni.clus <- stringr::str_sort(unique(clus), numeric = T)
  for (i in query.cluster) {
    ind <- which(clus == i)
    if(length(ind)>1){
      tab <- table(clus[which(colSums(scHolography.obj$adj.mtx[ind, ]) > 0)])}else{
        tab <- table(clus[which((scHolography.obj$adj.mtx[ind, ]) > 0)])
      }
    matrix[i, names(tab)] <- tab
  }
  data.fra <- reshape2::melt(matrix/rowSums(matrix))
  if (is.null(row.group.level) == F) {
    data.fra$Var1 <- factor(data.fra$Var1, row.group.level)
  }
  getPalette = colorRampPalette(brewer.pal(12, pal))
  sig.class.ls <- lapply(query.cluster, function(ind.clus) {
    query.cluster.ls <- list(c(ind.clus), uni.clus[-which(uni.clus %in%
                                                         c(ind.clus))])
    matrix.ls <- vector("list", 2)
    for (i in 1:2) {
      ind <- which(clus %in% query.cluster.ls[[i]])
      matrix.this <- matrix(unlist(lapply(uni.clus, function(x) {
        this.clus <- which(clus == x)
        if(length(this.clus)>1){
          rowSums(scHolography.obj$adj.mtx[, this.clus])
        }else{
          scHolography.obj$adj.mtx[, this.clus]
        }
      })), ncol = length(uni.clus), byrow = F)
      colnames(matrix.this) <- uni.clus
      matrix.ls[[i]] <- matrix.this[ind, ]
    }
    pval.ls <- lapply(uni.clus, function(x) {
      if(is.null(dim(matrix.ls[[1]]))){
        m1 <- matrix.ls[[1]][x]

      }else{
        m1 <- matrix.ls[[1]][, x]
      }
      if(is.null(dim(matrix.ls[[2]]))){
        m2 <- matrix.ls[[2]][x]

      }else{
        m2 <- matrix.ls[[2]][, x]
      }

      wilcox.test(m1, m2, alternative = "greater")$p.value
    })
    names(pval.ls) <- uni.clus
    unlist(pval.ls)[names(which(unlist(pval.ls) >= 0.05))]
  })
  names(sig.class.ls) <- query.cluster
  norm.mat <- matrix/rowSums(matrix)
  norm.mat.og <-norm.mat
  for (i in rownames(norm.mat)) {
    norm.mat[i, names(sig.class.ls[[i]])] <- 0
  }
  data.fra.sig <- reshape2::melt(norm.mat)
  col.use <- getPalette(length(uni.clus))
  if (is.null(color) == F) {
    col.use <- color
  }
  if (is.null(row.group.level) == F) {
    data.fra.sig$Var1 <- factor(data.fra.sig$Var1, row.group.level)
  }
  neighbor.comp <- ggplot(data.fra, aes(x = Var1, y = value,
                                        fill = Var2)) + geom_bar(position = "stack", stat = "identity") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45,
                                                       hjust = 1)) + scale_fill_manual(values = col.use) + labs(fill = "",
                                                                                                                x = "", y = "Composition")
  neighbor.comp.sig <- ggplot(data.fra.sig, aes(x = Var1, y = value,
                                                fill = Var2)) + geom_bar(position = "stack", stat = "identity") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45,
                                                       hjust = 1)) + scale_fill_manual(values = col.use) + labs(fill = "",
                                                                                                                x = "", y = "Composition")
  sig.class.ls <- lapply(query.cluster, function(ind.clus) {
    query.cluster.ls <- list(c(ind.clus), uni.clus[-which(uni.clus %in%
                                                            c(ind.clus))])
    matrix.ls <- vector("list", 2)
    for (i in 1:2) {
      ind <- which(clus %in% query.cluster.ls[[i]])
      matrix.this <- matrix(unlist(lapply(uni.clus, function(x) {
        this.clus <- which(clus == x)
        if(length(this.clus)>1){
          rowSums(scHolography.obj$adj.mtx[, this.clus])
        }else{
          scHolography.obj$adj.mtx[, this.clus]
        }
      })), ncol = length(uni.clus), byrow = F)
      colnames(matrix.this) <- uni.clus
      matrix.ls[[i]] <- matrix.this[ind, ]
    }
    pval.ls <- lapply(uni.clus, function(x) {
      if(is.null(dim(matrix.ls[[1]]))){
        m1 <- matrix.ls[[1]][x]

      }else{
        m1 <- matrix.ls[[1]][, x]
      }
      if(is.null(dim(matrix.ls[[2]]))){
        m2 <- matrix.ls[[2]][x]

      }else{
        m2 <- matrix.ls[[2]][, x]
      }

      wilcox.test(m1, m2, alternative = "greater")$p.value
    })
    names(pval.ls) <- uni.clus
    unlist(pval.ls)[names(which(unlist(pval.ls) < 0.05))]
  })
  names(sig.class.ls) <- query.cluster


  neighbor.comp.sig.ls <- lapply(1:nrow(norm.mat),function(r){
    sort(norm.mat[r,which(norm.mat[r,]>0)],decreasing = T)
  })
  names(neighbor.comp.sig.ls) <- rownames(norm.mat)

  list(neighbor.comp.plot = neighbor.comp, neighbor.comp.sig.plot = neighbor.comp.sig,
       significance = sig.class.ls, neighbor.comp.sig=neighbor.comp.sig.ls,neighbor.comp=norm.mat.og)
}


