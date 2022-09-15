#' Do MultiBar Heatmap
#' @import Seurat
#' @import dplyr
#' @import RColorBrewer
#' @import rlang
#' @import ggplot2
#' @import grid
DoMultiBarHeatmap <- function (object,
                               features = NULL,
                               cells = NULL,
                               group.by = "ident",
                               additional.group.by = NULL,
                               additional.group.sort.by = NULL,
                               cols.use = NULL,
                               group.bar = TRUE,
                               disp.min = -2.5,
                               disp.max = NULL,
                               slot = "scale.data",
                               assay = NULL,
                               label = TRUE,
                               size = 5.5,
                               hjust = 0,
                               angle = 45,
                               raster = TRUE,
                               draw.lines = TRUE,
                               lines.width = NULL,
                               group.bar.height = 0.02,
                               combine = TRUE,
                               palette="Paired",
                               disPlot=F
)
{

  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% Seurat::DefaultAssay(object = object)
  Seurat::DefaultAssay(object = object) <- assay
  features <- features %||% Seurat::VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data",
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = Seurat::GetAssayData(object = object,
                                                         slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot,
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ",
            slot, " slot for the ", assay, " assay: ", paste(bad.features,
                                                             collapse = ", "))
  }

  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ",
                paste(bad.sorts, collapse = ", "))
      }
    }
  }

  data <- as.data.frame(x = as.matrix(x = t(x = Seurat::GetAssayData(object = object,
                                                                     slot = slot)[features, cells, drop = FALSE])))

  object <- suppressMessages(expr = Seurat::StashIdent(object = object,
                                                       save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }

    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]

    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }
    }

    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) *
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) *
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }

      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells

      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells

      group.use <- rbind(group.use, placeholder.groups)

      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }

      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells),
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells,
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }

    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])

    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster,
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features,
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])

    if (group.bar) {
      pbuild <- ggplot2::ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {disPlot
        if (disPlot) {
          colid = ""
        }else if(colname == i){
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }

        # Default
        cols[[colname]] <- colorRampPalette(brewer.pal(12, palette))(length(x = levels(x = group.use[[colname]])))
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }

        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])

        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)

        plot <- suppressMessages(plot +
                                   ggplot2::annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) +
                                   ggplot2::annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = grid::gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   ggplot2::coord_cartesian(ylim = c(0, y.max), clip = "off"))

        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos),
                                    label.x.pos)
          plot <- plot + ggplot2::geom_text(stat = "identity",
                                            data = label.x.pos, ggplot2::aes_string(label = "group",
                                                                                    x = "label.x.pos"), y = y.max + y.max *
                                              0.03 * 0.5, angle = angle, hjust = hjust,
                                            size = size)
          plot <- suppressMessages(plot + ggplot2::coord_cartesian(ylim = c(0,
                                                                            y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) *
                                                                              size), clip = "off"))
        }
      }
    }
    plot <- plot + ggplot2::theme(line = ggplot2::element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <-Seurat::CombinePlots(plots = plots)
  }
  return(plots)
}



#' Relative Spatial Analysis
#' @export
#' @import dplyr
#' @import RColorBrewer
#' @import pracma
#' @import ggplot2
#' @import stringr
#' @import Seurat
#' @import grid
#' @import cowplot
relativeSpatialAnalysis <- function(scHolography.obj, query.cluster, ref.cluster,plotByDist=F,geneOI=NULL,
                                    annotationToUse="orig.cluster",assayToUse="SCT", quant=0.1,
                                    nCperL=NULL,nL=NULL,heatmapGp=NULL,pal="Paired",n.neighbor=30){

  # query.cluster=c("c3","c10","c17")
  # ref.cluster="c0"
  # plotByDist=T
  # geneOI=c("KRT5","KRT14")
  # annotationToUse="orig.cluster"
  # assayToUse="SCT"
  # nCperL=NULL
  # nL=NULL
  # quant=0.1
  # heatmapGp=NULL
  # n.neighbor=30
  # pal="Paired"


  scHolography.sc<-scHolography.obj$scHolography.sc
  scHolography.sc[["orig.cluster"]][[1]]<-factor(paste("c",scHolography.sc[["orig.cluster"]][[1]],sep = ""),levels = paste("c",levels(scHolography.sc[[annotationToUse]][[1]]),sep = ""))
  query.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%query.cluster)
  ref.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%ref.cluster)
  coord <- cbind(scHolography.sc$x3d_sp,scHolography.sc$y3d_sp,scHolography.sc$z3d_sp)
  clus.dist <- pracma::distmat(coord[query.cluster.ind,],coord[ref.cluster.ind,])
  query.to.ref.dis <- (colMeans(apply(clus.dist, 1, sort)[1:n.neighbor,]))


  if(plotByDist){
    nCperL=1
    quant = NULL
  }
  if(is.null(nCperL)){
    if(is.null(nL)){
      "Need at leat one value for quant, nCperL, or nL"
    }
    else{nCperL <- floor(length(query.cluster.ind)/nL)}
  }else{
    nL <- ceiling(length(query.cluster.ind)/nCperL)

  }


  if(is.null(quant)==F){
    close.ind <-which(query.to.ref.dis<quantile(query.to.ref.dis,quant))
    far.ind <- which(query.to.ref.dis>quantile(query.to.ref.dis,(1-quant)))
    layer.seq <- rep(NA,length(query.cluster.ind))
    layer.seq[close.ind]<-paste(paste(query.cluster,collapse = "_"),"1Close",sep = "_")
    layer.seq[far.ind]<-paste(paste(query.cluster,collapse = "_"),"3Far",sep = "_")
    layer.seq[which(is.na(layer.seq))]<-paste(paste(query.cluster,collapse = "_"),"2Mid",sep = "_")
    df <- data.frame(dist=query.to.ref.dis)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=dist)) +
      ggplot2::geom_density()+ ggplot2::geom_vline(ggplot2::aes(xintercept=quantile(dist,quant)),color="blue", linetype="dashed", size=1)+ ggplot2::geom_vline(ggplot2::aes(xintercept=quantile(dist,(1-quant))),color="blue", linetype="dashed", size=1)
    show(p)
  }else{
    layer.seq <- rep(NA,length(query.cluster.ind))
    for (i in 1:(nL-1)) {
      l.ind<-which(query.to.ref.dis%in%(sort(query.to.ref.dis)[((i-1)*nCperL+1):(i*nCperL)])==T)
      layer.seq[l.ind]<-paste(paste(query.cluster,collapse = "_"),i,sep = "_")
    }
    layer.seq[which(is.na(layer.seq))]<-paste(paste(query.cluster,collapse = "_"),nL,sep = "_")
  }

  query.cluster.sub<-subset(scHolography.sc,cells = c(query.cluster.ind,ref.cluster.ind))
  query.cluster.sub$cal.dist <- (c(query.to.ref.dis,rep(NA,length(ref.cluster.ind))))
  query.cluster.sub[[paste(paste(query.cluster,collapse = "_"),paste(ref.cluster,collapse = "_"),sep = "To")]]<-factor(c(layer.seq,scHolography.sc[[annotationToUse]][[1]][ref.cluster.ind]),levels=c(stringr::str_sort(unique(layer.seq),numeric = T),ref.cluster))
  query.cluster.sub$sublayer<-query.cluster.sub[[paste(paste(query.cluster,collapse = "_"),paste(ref.cluster,collapse = "_"),sep = "To")]]
  query.cluster.sub <- Seurat::SetIdent(query.cluster.sub,
                                        value = paste(paste(query.cluster,collapse = "_"),paste(ref.cluster,collapse = "_"),sep = "To"),sep = "To")

  if(plotByDist==F){
    if(is_null(geneOI)){
      markers<-Seurat::FindAllMarkers(subset(query.cluster.sub,cells = which((query.cluster.sub@active.ident %in% ref.cluster)==F)),
                                      only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      markers %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = -p_val_adj) -> top10
    }}

  if(is.null(geneOI)){
    if(plotByDist){
      show("Need the Gene of Interest list if ploting by ordered distance")
    }else{
      p1 <- DoMultiBarHeatmap(subset(query.cluster.sub,cells = which((query.cluster.sub@active.ident %in% ref.cluster)==F)),
                              features = top10$gene,additional.group.by = if(is.null(heatmapGp)==F){heatmapGp},label = T,palette = pal) +
        Seurat::NoLegend()+ggplot2::scale_fill_gradientn(colors = c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"))
      show(p1)
    }
  }else{
    obj.sub <- subset(query.cluster.sub,cells = which((query.cluster.sub@active.ident %in% ref.cluster)==F))

    p1 <- DoMultiBarHeatmap(obj.sub,features = geneOI, draw.lines = if(plotByDist){F}else{T},additional.group.by = if(is.null(heatmapGp)==F){heatmapGp},label = if(plotByDist){F}else{T},palette = pal) + Seurat::NoLegend()+ggplot2::scale_fill_gradientn(colors = c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F"))
    if(plotByDist==F){
      show(p1)
    }else{

      p1 <- DoMultiBarHeatmap(obj.sub,features = geneOI, draw.lines = if(plotByDist){F}else{T},size = 3,disPlot = T,additional.group.by = annotationToUse,label = F,palette = pal) + Seurat::NoLegend()+ggplot2::scale_fill_gradientn(colors = c("#053061","#3784BB","#A7CFE4","#F7F7F7","#F7B698","#CA4741","#67001F")) +ggplot2::theme(legend.title=ggplot2::element_blank())
      fea <- Seurat::VariableFeatures(obj.sub)

      obj.sub[[annotationToUse]][[1]] <-droplevels(obj.sub[[annotationToUse]][[1]])

      labal.col<-colorRampPalette(brewer.pal(12, pal))(length(levels(scHolography.sc[[annotationToUse]][[1]])))
      p<- Seurat::DoHeatmap(obj.sub,features = fea[1],group.by = annotationToUse,group.colors =labal.col[unlist(lapply(levels(obj.sub[[annotationToUse]][[1]] ),function(x){which(levels(scHolography.sc[[annotationToUse]][[1]])%in%x)}))] )+ggplot2::scale_fill_gradientn(colors = colorRampPalette(brewer.pal(12, pal))(100),name="Distance",labels=c("close","far"),n.breaks=2)
      # show(p)
      #p<- Seurat::DoHeatmap(obj.sub,features = "KRT5",group.by = annotationToUse)
      grid::grid.newpage()
      legend <- cowplot::get_legend(p)
      show( cowplot::plot_grid( p1, NULL, legend, rel_widths = c(1, -0.1, 1), align = "hv",nrow = 1 ))

    }
  }

  scHolography.obj$scHolography.sc <- query.cluster.sub
  if(exists("top10")){list(scHolography.obj=scHolography.obj,DEG=top10)}else{list(scHolography.obj=scHolography.obj,DEG=NULL)}

}


# sub.scHolography.obj<-relativeSpatialAnalysis(scHolography.obj = scHolography.obj,query.cluster = "Fibroblast",ref.cluster = "Basal",annotationToUse = "celltype",nL = 3)
# scHolographyPlot(sub.scHolography.obj,color.by = "FibroblastToBasal")
# Seurat::VlnPlot(sub.scHolography.obj$scHolography.sc,"SPARC")
# table(sub.scHolography.obj$scHolography.sc$FibroblastToBasal,sub.scHolography.obj$scHolography.sc$orig.cluster)
#
# for (i in levels(sub.scHolography.obj$scHolography.sc@active.ident)) {
#   neighbor.cluster.i<-colSums(scHolography.obj$adj.mtx[which(query.cluster.sub@active.ident==i),])
#   neighbor.type.cluster.i<-unlist(lapply(levels(scHolography.sc[[annotationToUse]][[1]]),function(x) sum(neighbor.cluster.i[which(scHolography.sc[[annotationToUse]][[1]]==x)])))
#   #show(neighbor.type.cluster.i/as.vector(table(scHolography.sc[[annotationToUse]][[1]])))
#   show(neighbor.type.cluster.i/sum(neighbor.type.cluster.i))
# }



#' Find Space Driver Gene
#' @export
#' @import dplyr
#' @import pracma
#' @import stats
#' @import enrichR

findDriverGene <- function(scHolography.obj,query.cluster,ref.cluster,k1=1.96,k2=1.96,annotationToUse="orig.cluster",n.closeNeighbor=30, assayToUse="SCT",bandwidth=1){

  scHolography.sc<-scHolography.obj$scHolography.sc

  query.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%query.cluster)
  ref.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%ref.cluster)
  coord <- cbind(scHolography.sc$x3d_sp,scHolography.sc$y3d_sp,scHolography.sc$z3d_sp)
  clus.dist <- pracma::distmat(coord[query.cluster.ind,],coord[ref.cluster.ind,])
  query.to.ref.dis <- (colMeans(apply(clus.dist, 1, sort)[1:n.closeNeighbor,]))

  sub.query.sc.obj<-subset(scHolography.sc, cells = query.cluster.ind)
  p<-hist(query.to.ref.dis,breaks =seq(from=0,to=max(ceiling(query.to.ref.dis/bandwidth)*bandwidth),by=bandwidth))
  count<-p$count
  mid <- rowMeans(embed(p$breaks,2))

  avg.expr <- apply(embed(p$breaks,2),1,function(x){
    if(length(intersect(which(query.to.ref.dis>x[2]),which(query.to.ref.dis<=x[1])))>0){
      sub.query.sc.int <- subset(sub.query.sc.obj,cells=intersect(which(query.to.ref.dis>x[2]),which(query.to.ref.dis<=x[1])))
      apply(as.matrix(sub.query.sc.int@assays$SCT@data),1,median)}else{
        rep(NA, nrow(sub.query.sc.obj@assays$SCT@data))
      }
  })
  colnames(avg.expr) <- mid
  rownames(avg.expr) <-  rownames(sub.query.sc.obj@assays$SCT@data)
  avg.expr.sub <- avg.expr[which(apply(avg.expr,1,function(x) sum(x,na.rm = T))>0),]
  regression <- apply(avg.expr.sub, 1,function(x){
    pois <- (stats::glm( x~ as.numeric(colnames(avg.expr.sub)),"poisson"))
    summary(pois)$coefficients[2,]
  })
  regression<-t(regression)
  regression<-as.data.frame(regression)
  genes <- regression$`z value`
  setEnrichrSite("Enrichr") # Human genes
  dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
  enriched <- enrichr(rownames(regression)[which(regression$`z value`<(-k1))], dbs)
  enriched[["GO_Biological_Process_2021"]]
  show(plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "Proximal Gene Ontology Analysis")
  )
  enriched <- enrichr(rownames(regression)[which(regression$`z value`>(k2))], dbs)
  enriched[["GO_Biological_Process_2021"]]
  show(plotEnrich(enriched[[3]], showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", title = "Distal Gene Ontology Analysis")
  )
  list(regression=regression,avg.expr.sub=avg.expr.sub,count=count,bandwidth=bandwidth)

}


#' Visualization of Space Driver Genes
#' @export
#' @import ggplot2

driverGenePlot <- function(findDG.boj,gene){
  avg.expr.sub <- findDG.boj[[2]]
  count <- findDG.boj[[3]]
  bandwidth <- findDG.boj[[4]]
  line.dat<-data.frame(dist=as.numeric(names(avg.expr.sub[gene,])),expr=avg.expr.sub[gene,],ct=count)
  coeff <- max(line.dat$ct,na.rm = T)/max(line.dat$expr,na.rm = T)
  ggplot2::ggplot(line.dat, ggplot2::aes(x=dist))+
    ggplot2::geom_bar(stat="identity",fill="white",color="grey70", alpha=0.5, width=bandwidth,ggplot2::aes(y=ct))+
    ggplot2::geom_point(ggplot2::aes(y=expr*coeff))+
    ggplot2::geom_smooth(ggplot2::aes(y=expr*coeff),method="loess",se=F)+
    ggplot2::scale_y_continuous(
      name = "Frequency",
      sec.axis = ggplot2::sec_axis(~./coeff, name="Expression"))+ggplot2::ggtitle(as.character(gene))+
    ggplot2::theme_classic()

}

