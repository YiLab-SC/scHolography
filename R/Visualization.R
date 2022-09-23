
#' scHolography 3D Tissue Reconstruction Visualization
#' @export
#' @import RColorBrewer
#' @import plotly
#' @import Seurat

scHolographyPlot<-function(scHolography.obj, dim=3, cells=NULL,feature=NULL,cutoff=NULL, color.by="orig.cluster", dot.size = 5 , palette="Paired", highlight=NULL,assayToUse="SCT"){
  scHolography.sc<-scHolography.obj$scHolography.sc
  if(is.null(levels(scHolography.sc[[color.by]][[1]]))){
    scHolography.sc[[color.by]][[1]]<-as.factor(scHolography.sc[[color.by]][[1]])
  }
  getPalette =colorRampPalette(brewer.pal(12, palette))
  my.scheme = getPalette(length(levels(scHolography.sc[[color.by]][[1]])))

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
          out.plot<-plotly::plot_ly( x=scHolography.sc$x_sp[cellToPlot], y=scHolography.sc$y_sp[cellToPlot], colors=c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"),type="scatter", mode="markers", color=as.vector(Seurat::GetAssayData(scHolography.sc,assay = assayToUse)[feature,])[cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

        }else{
          out.plot<-plotly::plot_ly( x=scHolography.sc$x_sp[cellToPlot], y=scHolography.sc$y_sp[cellToPlot], colors=c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"),type="scatter", mode="markers", color=scHolography.sc[[feature]][[1]][cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

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
        out.plot<-plotly::plot_ly(x=scHolography.sc$x3d_sp[cellToPlot], y=scHolography.sc$y3d_sp[cellToPlot], z=scHolography.sc$z3d_sp[cellToPlot],  type="scatter3d", mode="markers", color=scHolography.sc[[color.by]][[1]][cellToPlot],colors = my.scheme,marker=list(size=dot.size))
      }else{
        if(feature%in%rownames(Seurat::GetAssayData(scHolography.sc,assay =as.character(assayToUse)))){
          out.plot<-plotly::plot_ly( x=scHolography.sc$x3d_sp[cellToPlot], y=scHolography.sc$y3d_sp[cellToPlot], z=scHolography.sc$z3d_sp[cellToPlot], colors=c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"),type="scatter3d", mode="markers", color=as.vector(Seurat::GetAssayData(scHolography.sc,assay = assayToUse)[feature,])[cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

        }else{
          out.plot<-plotly::plot_ly( x=scHolography.sc$x3d_sp[cellToPlot], y=scHolography.sc$y3d_sp[cellToPlot],z=scHolography.sc$z3d_sp[cellToPlot],  colors=c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"),type="scatter3d", mode="markers", color=scHolography.sc[[feature]][[1]][cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

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



#' Cell State Connection Plot in 3D
#' @export
#' @import RColorBrewer
#' @import pracma
#' @import igraph
scHolographyConnectionPlot<-function(scHolography.obj, color.by="orig.cluster" , palette="Paired",seed=60611, n.closeNeighbor=30,n.mutualNeighbor=NA){
  scHolography.sc<-scHolography.obj$scHolography.sc
  #scHolography.sc<-scHolography.obj$starPoints.sc
  coord <- cbind(scHolography.sc$x3d_sp,scHolography.sc$y3d_sp,scHolography.sc$z3d_sp)
  affinity <- pracma::distmat(coord,coord)
  perc.vec<-c()
  if(is.null(levels(scHolography.sc[[color.by]][[1]]))){
    scHolography.sc[[color.by]][[1]]<-as.factor(scHolography.sc[[color.by]][[1]])
  }

  if(is.na(n.mutualNeighbor)){
    n.mutualNeighbor=5
  }
  for (i in levels(scHolography.sc[[color.by]][[1]])) {
    for (j in levels(scHolography.sc[[color.by]][[1]])) {
      rInd<-which(scHolography.sc[[color.by]][[1]]%in%c(i))
      cInd<-which(scHolography.sc[[color.by]][[1]]%in%c(j))
      sub.adj<-affinity[rInd,cInd]
      perc<-sum(apply(sub.adj,1,function(x) sort(x)[1:min(n.closeNeighbor,length(cInd))]))/(length(rInd)*min(n.closeNeighbor,length(cInd)))
      perc.vec<-c(perc.vec,perc)
    }
  }
  clus.con<-matrix(perc.vec,ncol = length(levels(scHolography.sc[[color.by]][[1]])),byrow = T)

  rslt.ls <- interpretDecomp(clus.con,nslot = n.mutualNeighbor,is3D = F)
  adj.mtx<-rslt.ls[[1]]
  utls<-rslt.ls[[2]]
  diag(adj.mtx)<-0
  set.seed(seed)
  adj.mtx.toUse <- adj.mtx*utls
  colnames(adj.mtx.toUse)<-levels(scHolography.sc[[color.by]][[1]])
  rownames(adj.mtx.toUse)<-levels(scHolography.sc[[color.by]][[1]])
  set.seed(seed)
  g.clus.con<- igraph::graph.adjacency(adj.mtx.toUse, mode="undirected",weighted = T)
  getPalette =colorRampPalette(brewer.pal(12, palette))
  my.scheme = getPalette(length(levels(scHolography.sc[[color.by]][[1]])))
  igraph::V(g.clus.con)$color <-my.scheme
  igraph::E(g.clus.con)$edge.color <- "gray80"
  igraph::V(g.clus.con)$size <- (200/(length(levels(scHolography.sc[[color.by]][[1]]))*2))+15*(table(scHolography.sc[[color.by]][[1]])-min(table(scHolography.sc[[color.by]][[1]])))/(max(table(scHolography.sc[[color.by]][[1]]))-min(table(scHolography.sc[[color.by]][[1]])))
  igraph::E(g.clus.con)$width <- (3/length(levels(scHolography.sc[[color.by]][[1]])))+5*(igraph::E(g.clus.con)$weight-min(igraph::E(g.clus.con)$weight))/(max(igraph::E(g.clus.con)$weight)-min(igraph::E(g.clus.con)$weight))

  #plot(g.clus.con, vertex.label.color="black", vertex.frame.color="#ffffff",vertex.label.cex=0.8+igraph::V(g.clus.con)$size/100 )
  #igraph::V(g.clus.con)$name <- as.character(levels(scHolography.sc[[color.by]][[1]]))
  plot(g.clus.con, vertex.frame.color="#ffffff",vertex.label.font=2, vertex.label.color="black",vertex.label.cex=0.8+igraph::V(g.clus.con)$size/100 )
}



#' 3D Cell State Distribution Plot
#' @export
#' @import RColorBrewer
#' @import ggplot2
#' @import ggridges
scHolography3dDistributionPlot<-function(scHolography.obj,dimToPlot,color.by="celltype",palette="Paired",transparency=0.7){
  scHolography.sc<-scHolography.obj$scHolography.sc
  if(is.null(levels(scHolography.sc[[color.by]][[1]]))){
    scHolography.sc[[color.by]][[1]]<-as.factor(scHolography.sc[[color.by]][[1]])
  }
  getPalette =colorRampPalette(brewer.pal(12, palette))
  my.scheme = getPalette(length(levels(scHolography.sc[[color.by]][[1]])))


  ggplot2::ggplot(scHolography.sc@meta.data, ggplot2::aes(x=scHolography.sc[[dimToPlot]][[1]],y=scHolography.sc[[color.by]][[1]], fill=scHolography.sc[[color.by]][[1]])) +
    ggridges::geom_density_ridges (alpha=transparency, color=NA)+ggplot2::scale_fill_manual(values=my.scheme)+ ggplot2::theme_classic()+ggplot2::xlab(as.character(dimToPlot))+ggplot2::ylab("Group")+ ggplot2::theme(legend.position = "none")
}
