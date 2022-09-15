#' scHolography 2D Single Cell Spatial Mapping to ST Data
#' @export
#' @import dplyr
#' @import matchingR
scHolographyMapToPixel<-function(scHolography.obj,nCperP=NULL,seed=60611){

  scHolography.sc<-scHolography.obj$scHolography.sc

  if(is.null(nCperP)){
    nCperP = ceiling(ncol(scHolography.obj$scHolography.sc)/ncol(scHolography.obj$scHolography.sp))
  }

  iter.mean<-scHolography.obj$est.array
  diag.mean<-(iter.mean+t(iter.mean))/2
  cell.to.pix.dis<-diag.mean[(ncol(scHolography.obj$scHolography.sp)+1):(ncol(scHolography.obj$scHolography.sp)+ncol(scHolography.obj$scHolography.sc)),1:ncol(scHolography.obj$scHolography.sp)]
  cell.to.pix.dis.M<-mean(cell.to.pix.dis)
  cell.to.pix.dis.SD<-sqrt(var(as.vector(cell.to.pix.dis)))
  cell.to.pix.dis.norm<-cell.to.pix.dis-cell.to.pix.dis.M
  cell.to.pix.dis.norm<-cell.to.pix.dis.norm/cell.to.pix.dis.SD
  cell.to.pix.ut<- -cell.to.pix.dis.norm
  cell.to.pix_results <- matchingR::galeShapley.collegeAdmissions(collegeUtils =cell.to.pix.ut, studentUtils =t(cell.to.pix.ut),slot=nCperP, studentOptimal = F)

  amp.vec<-as.vector(apply(scHolography.obj$scHolography.sp@images[[1]]@coordinates[,c("row","col")],1, function(x){
    t(rbind(x)[rep(1,nCperP), ])
  }))
  set.seed(seed)
  amp.vec<-amp.vec+runif(length(amp.vec))-0.5
  amp.mat<-matrix(amp.vec,ncol = 2,byrow = T)
  amp.mat<-cbind(amp.mat,as.vector(t(cell.to.pix_results$matched.colleges)))
  final.mat<-amp.mat[complete.cases(amp.mat), ]
  colnames(final.mat)<-c("row","col","ind")
  rownames(final.mat)<-final.mat[,"ind"]
  final.mat<-final.mat[as.character(1:ncol(scHolography.obj$scHolography.sc)),]
  scHolography.obj$scHolography.sc[["toPixel_x"]]<-final.mat[,"row"]
  scHolography.obj$scHolography.sc[["toPixel_y"]]<-final.mat[,"col"]
  scHolography.obj

}


#' Visualization of scHolography 2D Single Cell Spatial Mapping
#' @export
#' @import RColorBrewer
#' @import plotly
#' @import Seurat
scHolographyToPixelPlot<-function(scHolography.obj,color.by="orig.cluster", cutoff=NULL,assayToUse="SCT",feature=NULL,palette="Paired",dot.size=5, cells=NULL, highlight=NULL){
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
  if(is.null(highlight)){
    out.plot<-plotly::plot_ly(x=scHolography.sc$toPixel_x[cellToPlot], y=scHolography.sc$toPixel_y[cellToPlot], type="scatter", mode="markers", color=scHolography.sc[[color.by]][[1]][cellToPlot],colors = my.scheme,marker=list(size=dot.size)) }else{
      new.scheme<-c(my.scheme,"gray80")
      if(is.numeric(highlight)){
        cellToHighlight=highlight
      }else{
        cellToHighlight<-which(scHolography.sc[[color.by]][[1]]%in%highlight)
      }
      color.highlight <- scHolography.sc[[color.by]][[1]]
      levels(color.highlight)<-c(levels(color.highlight),"")
      color.highlight[which((1:length(color.highlight)%in%cellToHighlight)==F)] <- as.factor("")
      out.plot<-plotly::plot_ly(x=scHolography.sc$toPixel_x[cellToPlot], y=scHolography.sc$toPixel_y[cellToPlot], type="scatter", mode="markers", color=color.highlight[cellToPlot],colors = new.scheme,marker=list(size=dot.size))
    }
  if(is.null(feature)==F){
    if(feature%in%rownames(Seurat::GetAssayData(scHolography.obj$scHolography.sc,assay =as.character(assayToUse)))){
      out.plot<-plotly::plot_ly( x=scHolography.sc$toPixel_x[cellToPlot], y=scHolography.sc$toPixel_y[cellToPlot], colors=c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"),type="scatter", mode="markers", color=as.vector(Seurat::GetAssayData(scHolography.sc,assay = assayToUse)[feature,])[cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

    }else{
      out.plot<-plotly::plot_ly( x=scHolography.sc$toPixel_x[cellToPlot], y=scHolography.sc$toPixel_y[cellToPlot], colors=c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"),type="scatter", mode="markers", color=scHolography.sc[[feature]][[1]][cellToPlot],marker=list(size=dot.size))%>% plotly::layout(title=as.character(feature))

    }
  }
  out.plot
}
