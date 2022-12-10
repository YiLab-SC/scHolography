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
    if (!is.null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]
      if (!is.null(additional.group.sort.by)){
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
        if (!is.null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is.null(names(cols.use[[colname]]))) {
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

#' Compare Spatial Gene Dynamics
#' @export

findDiffGeneDynamics <- function(gene.dyn1, gene.dyn2){
  inter.gene <- intersect(rownames(gene.dyn1), rownames(gene.dyn2))
  gene.dyn1.sub <- gene.dyn1[rownames(gene.dyn1) %in% inter.gene,
  ]
  gene.dyn2.sub <- gene.dyn2[rownames(gene.dyn2) %in% inter.gene,
  ]
  gene.dyn1.sub$rank <- 1:nrow(gene.dyn1.sub)
  rank1 <- gene.dyn1.sub$rank
  names(rank1) <- gene.dyn1.sub$gene
  gene.dyn2.sub$rank <- 1:nrow(gene.dyn2.sub)
  rank2 <- gene.dyn2.sub$rank
  names(rank2) <- gene.dyn2.sub$gene
  rank.dif <- sort(rank2[names(rank1)]-rank1)
  rank1.ord <- rank1[names(rank.dif)]
  rank2.ord <- rank2[names(rank.dif)]
  sig <-unlist(lapply(1:length(rank.dif), function(x){
    intermediate <- pnorm(rank.dif[x],mean = mean(rank.dif),sd = sd(rank.dif))
    if(intermediate>=.5){
      intermediate <- 1-intermediate
    }
    intermediate
  }))
  sig <-2*sig
  z.sc <- ((rank.dif-mean(rank.dif))/sd(rank.dif))

  out <- cbind(rank.dif,sig,z.sc)
  rownames(out) <-names( rank.dif )
  colnames(out) <- c("Rank.Diff","P.Val","Z.Score")
  out<-as.data.frame(out)
  out
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
#' @import igraph
#' @import viridis
#' @import colorspace

relativeSpatialAnalysis <- function(scHolography.obj, query.cluster, ref.cluster,plotByDist=F,geneOI=NULL,
                                    annotationToUse="orig.cluster",assayToUse="SCT", quant.left=0.1, quant.right=0.9,
                                    nCperL=NULL,nL=NULL,heatmapGp=NULL,pal="Paired",n.neighbor=30,extreme.comp=F,heatmap.pal="viridis"){

  scHolography.sc<-scHolography.obj$scHolography.sc
  scHolography.sc[[annotationToUse]][[1]]  <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]), function(x) paste0(strsplit(x,split = " ")[[1]],collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]]  <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]), function(x) paste0(strsplit(x,split = "/")[[1]],collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]] <- factor(scHolography.sc[[annotationToUse]][[1]],levels = stringr::str_sort(unique(scHolography.sc[[annotationToUse]][[1]]),numeric = T))
  query.cluster  <- unlist(lapply(as.character(query.cluster), function(x) paste0(strsplit(x,split = " ")[[1]],collapse = "_")))
  query.cluster  <- unlist(lapply(as.character(query.cluster), function(x) paste0(strsplit(x,split = "/")[[1]],collapse = "_")))
  ref.cluster  <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,split = " ")[[1]],collapse = "_")))
  ref.cluster  <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,split = "/")[[1]],collapse = "_")))

  if(heatmap.pal=="rdbu"){
    heat.col <- rev(RColorBrewer::brewer.pal(n = 25, name = "RdBu"))
  }else if(heatmap.pal=="magma"){
    heat.col <-viridis::viridis(25,option = "A")
  }else if (heatmap.pal=="rdbu_1"){
    heat.col<-colorspace::diverge_hsv(25)
  }else{
    heat.col <-viridis::viridis(25)
  }
  if(sum(is.na(as.numeric(as.character(levels(scHolography.sc[["orig.cluster"]][[1]])))))==0){
    scHolography.sc[[annotationToUse]][[1]]<-factor(paste("c",scHolography.sc[[annotationToUse]][[1]],sep = ""),levels = paste("c",levels(scHolography.sc[[annotationToUse]][[1]]),sep = ""))
    query.cluster<-paste("c",query.cluster,sep = "")

    ref.cluster<-paste("c",ref.cluster,sep = "")

  }

  query.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%query.cluster)
  ref.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%ref.cluster)
  graph <- igraph::graph_from_adjacency_matrix(scHolography.obj$adj.mtx,mode = "undirected")
  dist <- igraph::distances(graph, mode="out")
  clus.dist <- dist[query.cluster.ind,ref.cluster.ind]
  query.to.ref.dis <- (colMeans(apply(clus.dist, 1, sort)[1:n.neighbor,]))

  if(plotByDist){
    nCperL=1
    quant.left = NULL
    quant.right = NULL
  }
  if(is.null(nCperL)){
    if(is.null(nL)){
      "Need at leat one value for quant, nCperL, or nL"
    }
    else{nCperL <- floor(length(query.cluster.ind)/nL)
    quant.left = NULL
    quant.right = NULL}
  }else{
    nL <- ceiling(length(query.cluster.ind)/nCperL)
    quant.left = NULL
    quant.right = NULL
  }
  if(is.null(quant.left)==F){
    close.ind <-which(query.to.ref.dis<quantile(query.to.ref.dis,quant.left))
    far.ind <- which(query.to.ref.dis>quantile(query.to.ref.dis,quant.right))
    layer.seq <- rep(NA,length(query.cluster.ind))
    layer.seq[close.ind]<-paste(paste(query.cluster,collapse = "_"),"1Proximal",sep = "_")
    layer.seq[far.ind]<-paste(paste(query.cluster,collapse = "_"),"3Distal",sep = "_")
    layer.seq[which(is.na(layer.seq))]<-paste(paste(query.cluster,collapse = "_"),"2Intermediate",sep = "_")
    df <- data.frame(dist=query.to.ref.dis)
    p <- ggplot2::ggplot(df, ggplot2::aes(x=dist)) +
      ggplot2::geom_density()+ ggplot2::geom_vline(ggplot2::aes(xintercept=quantile(dist,quant.left)),color="blue", linetype="dashed", size=1)+ ggplot2::geom_vline(ggplot2::aes(xintercept=quantile(dist,(quant.right))),color="blue", linetype="dashed", size=1)
    dist.plot <- p
    show(p)
  }else{
    layer.seq <- rep(NA,length(query.cluster.ind))
    if(plotByDist){
      for (i in 1:(nL)) {
        l.ind<-which(query.to.ref.dis%in%(sort(query.to.ref.dis)[((i-1)*nCperL+1):(i*nCperL)])==T)
        layer.seq[l.ind]<-paste(paste(query.cluster,collapse = "_"),i,sep = "_")
      }
    }else{
      for (i in 1:(nL-1)) {
        l.ind<-which(query.to.ref.dis%in%(sort(query.to.ref.dis)[((i-1)*nCperL+1):(i*nCperL)])==T)
        layer.seq[l.ind]<-paste(paste(query.cluster,collapse = "_"),i,sep = "_")
      }
      layer.seq[which(is.na(layer.seq))]<-paste(paste(query.cluster,collapse = "_"),nL,sep = "_")
    }
  }

  query.cluster.sub<-subset(scHolography.sc,cells = c(query.cluster.ind,ref.cluster.ind))
  query.cluster.sub$cal.dist <- (c(query.to.ref.dis,rep(NA,length(ref.cluster.ind))))

  query.cluster.sub[[paste(paste(query.cluster,collapse = "_"),paste(ref.cluster,collapse = "_"),sep = "To")]]<-factor(c(layer.seq,as.character(scHolography.sc[[annotationToUse]][[1]][ref.cluster.ind])),levels=c(ref.cluster,stringr::str_sort(unique(layer.seq),numeric = T)))
  if(length(intersect(query.cluster,ref.cluster))>0){
    names.ls <-colnames(query.cluster.sub)
    ind.rg <-((length(layer.seq)+1):ncol(query.cluster.sub))
    names.ls[ind.rg]<-paste(names.ls[ind.rg],"_ref",sep = "")
    query.cluster.sub <- Seurat::RenameCells(query.cluster.sub,new.names =names.ls)
  }
  if(extreme.comp){
    query.cluster.sub <- subset(query.cluster.sub, cells=which(query.cluster.sub[[paste(paste(query.cluster,collapse = "_"),paste(ref.cluster,collapse = "_"),sep = "To")]]!=paste(paste(query.cluster,collapse = "_"),"2Intermediate",sep = "_")))
  }

  query.cluster.sub <- Seurat::SetIdent(query.cluster.sub,
                                        value = paste(paste(query.cluster,collapse = "_"),paste(ref.cluster,collapse = "_"),sep = "To"),sep = "To")

  if(plotByDist==F){
    if(is.null(geneOI)){
      marker.obj <- subset(query.cluster.sub,cells = which((query.cluster.sub@active.ident %in% ref.cluster)==F))
      Seurat::DefaultAssay(marker.obj) <- "RNA"
      marker.obj <- Seurat::NormalizeData(marker.obj)
      all.genes <- rownames(marker.obj)
      marker.obj <- Seurat::ScaleData(marker.obj, features = all.genes)
      markers<-Seurat::FindAllMarkers(marker.obj,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      top10 <- markers %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC)
    }}

  if(is.null(geneOI)){
    if(plotByDist){
      show("Need the Gene of Interest list if ploting by ordered distance")
    }else{
      p1 <- DoMultiBarHeatmap(subset(query.cluster.sub,cells = which((query.cluster.sub@active.ident %in% ref.cluster)==F)),assay = assayToUse,
                              features = top10$gene,additional.group.by = if(is.null(heatmapGp)==F){heatmapGp},label = T,palette = pal) +
        Seurat::NoLegend()+ggplot2::scale_fill_gradientn(colors = heat.col)
      expr.plot <- p1
      show(p1)
    }
  }else{
    obj.sub <- subset(query.cluster.sub,cells = which((query.cluster.sub@active.ident %in% ref.cluster)==F))
    p1 <- DoMultiBarHeatmap(obj.sub,features = geneOI, assay = assayToUse,draw.lines = if(plotByDist){F}else{T},additional.group.by = if(is.null(heatmapGp)==F){heatmapGp},label = if(plotByDist){F}else{T},palette = pal) + Seurat::NoLegend()+ggplot2::scale_fill_gradientn(colors = heat.col)
    if(plotByDist==F){
      expr.plot <- p1
      show(p1)
    }else{

      p1 <- DoMultiBarHeatmap(obj.sub,features = geneOI, assay = assayToUse,draw.lines = if(plotByDist){F}else{T},size = 3,disPlot = T,additional.group.by = annotationToUse,label = F,palette = pal) + Seurat::NoLegend()+ggplot2::scale_fill_gradientn(colors = heat.col) +ggplot2::theme(legend.title=ggplot2::element_blank())
      fea <- Seurat::VariableFeatures(obj.sub)

      obj.sub[[annotationToUse]][[1]] <-droplevels(obj.sub[[annotationToUse]][[1]])
      labal.col<-colorRampPalette(brewer.pal(12, pal))(length(levels(scHolography.sc[[annotationToUse]][[1]])))
      p<- Seurat::DoHeatmap(obj.sub,label = F,features = fea[1],group.by = annotationToUse,group.colors =labal.col[unlist(lapply(levels(obj.sub[[annotationToUse]][[1]] ),function(x){which(levels(scHolography.sc[[annotationToUse]][[1]])%in%x)}))] )+ggplot2::scale_fill_gradientn(colors = colorRampPalette(brewer.pal(12, pal))(100),name="Distance",labels=c("Proximal","Distal"),n.breaks=2)
      grid::grid.newpage()
      legend <- cowplot::get_legend(p)
      expr.plot <- cowplot::plot_grid( p1, NULL, legend, rel_widths = c(1, -0.1, 1), align = "hv",nrow = 1 )
      show( expr.plot)

    }
  }

  scHolography.obj$scHolography.sc <- query.cluster.sub
  if(exists("top10")){
    if(exists("dist.plot")){
      list(scHolography.obj=scHolography.obj,DEG=top10,dist.plot=dist.plot,expr.plot=expr.plot)
    }else{
      list(scHolography.obj=scHolography.obj,DEG=top10,expr.plot=expr.plot)
    }
  }else{
    if(exists("dist.plot")){
      list(scHolography.obj=scHolography.obj,DEG=NULL,dist.plot=dist.plot,expr.plot=expr.plot)
    }else{
      list(scHolography.obj=scHolography.obj,DEG=NULL,expr.plot=expr.plot)
    }
  }
}

#' Find Gene Spatial Dynamics
#' @import Seurat
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import MASS
#' @import igraph
#' @import viridis
#' @import colorspace
#' @import stringr
#' @export
findGeneSpatialDynamics  <- function (scHolography.obj, query.cluster, ref.cluster, annotationToUse = "orig.cluster",
                                      assayToUse = "SCT", n.neighbor = 30)
{

  #plotByDist = F; geneOI = NULL; annotationToUse = "celltype"; assayToUse = "SCT";quant.left = 0.1; quant.right = 0.9; nCperL = NULL ;nL = NULL;heatmapGp = NULL; pal = "Paired"; n.neighbor = 30;extreme.comp = F;heatmap.pal = "viridis"
  scHolography.sc <- scHolography.obj$scHolography.sc
  scHolography.sc[[annotationToUse]][[1]] <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]),
                                                           function(x) paste0(strsplit(x, split = " ")[[1]], collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]] <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]),
                                                           function(x) paste0(strsplit(x, split = "/")[[1]], collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]] <- factor(scHolography.sc[[annotationToUse]][[1]],
                                                    levels = stringr::str_sort(unique(scHolography.sc[[annotationToUse]][[1]]),
                                                                               numeric = T))
  query.cluster <- unlist(lapply(as.character(query.cluster),
                                 function(x) paste0(strsplit(x, split = " ")[[1]], collapse = "_")))
  query.cluster <- unlist(lapply(as.character(query.cluster),
                                 function(x) paste0(strsplit(x, split = "/")[[1]], collapse = "_")))
  ref.cluster <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,
                                                                                      split = " ")[[1]], collapse = "_")))
  ref.cluster <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,
                                                                                      split = "/")[[1]], collapse = "_")))
  if (sum(is.na(as.numeric(as.character(levels(scHolography.sc[["orig.cluster"]][[1]]))))) ==
      0) {
    scHolography.sc[[annotationToUse]][[1]] <- factor(paste("c",
                                                            scHolography.sc[[annotationToUse]][[1]], sep = ""),
                                                      levels = paste("c", levels(scHolography.sc[[annotationToUse]][[1]]),
                                                                     sep = ""))
    query.cluster <- paste("c", query.cluster, sep = "")
    ref.cluster <- paste("c", ref.cluster, sep = "")
  }
  query.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]] %in%
                               query.cluster)
  ref.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]] %in%
                             ref.cluster)
  graph <- igraph::graph_from_adjacency_matrix(scHolography.obj$adj.mtx,
                                               mode = "undirected")
  dist <- igraph::distances(graph, mode = "out")
  clus.dist <- dist[query.cluster.ind, ref.cluster.ind]
  query.to.ref.dis <- (colMeans(apply(clus.dist, 1, sort)[1:n.neighbor,
  ]))

  names(query.to.ref.dis) <- colnames(scHolography.sc)[query.cluster.ind]
  feature.obj <-subset(scHolography.sc,cells = query.cluster.ind)
  feature.obj <- SCTransform(feature.obj)
  glm.reg <- lapply(VariableFeatures(feature.obj), function(x){
    dat <-data.frame(dist=query.to.ref.dis,expr=as.numeric(scHolography.sc@assays[[assayToUse]]@data[x,query.cluster.ind]))
    suppressWarnings( mod <- glm(expr ~ dist, dat, family = "poisson"))
    names(x) <- "gene"
    c(x,summary(mod)$coefficients[2,])

  })

  reg.mat <-do.call(rbind,glm.reg)
  reg.mat <-as.data.frame(reg.mat)

  reg.mat$Estimate <- as.numeric(reg.mat$Estimate)
  reg.mat$`Std. Error` <- as.numeric(reg.mat$`Std. Error` )
  reg.mat$`Pr(>|z|)` <- as.numeric(reg.mat$`Pr(>|z|)`)
  reg.mat$`z value` <- as.numeric(reg.mat$`z value`)
  out <- arrange(reg.mat,`z value`)
  rownames(out) <- out$gene
  out
}


#' Expression by Distance Plot
#' @import Seurat
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import igraph
#' @import viridis
#' @import colorspace
#' @import stringr
#' @export
expressionByDistPlot <- function(scHolography.obj, query.cluster, ref.cluster, geneOI,annotationToUse = "orig.cluster", assayToUse = "SCT",
                                 pal = "Paired", n.neighbor = 30){
  scHolography.sc <- scHolography.obj$scHolography.sc
  scHolography.sc[[annotationToUse]][[1]] <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]),
                                                           function(x) paste0(strsplit(x, split = " ")[[1]], collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]] <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]),
                                                           function(x) paste0(strsplit(x, split = "/")[[1]], collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]] <- factor(scHolography.sc[[annotationToUse]][[1]],
                                                    levels = stringr::str_sort(unique(scHolography.sc[[annotationToUse]][[1]]),
                                                                               numeric = T))
  query.cluster <- unlist(lapply(as.character(query.cluster),
                                 function(x) paste0(strsplit(x, split = " ")[[1]], collapse = "_")))
  query.cluster <- unlist(lapply(as.character(query.cluster),
                                 function(x) paste0(strsplit(x, split = "/")[[1]], collapse = "_")))
  ref.cluster <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,
                                                                                      split = " ")[[1]], collapse = "_")))
  ref.cluster <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,
                                                                                      split = "/")[[1]], collapse = "_")))
  if (sum(is.na(as.numeric(as.character(levels(scHolography.sc[["orig.cluster"]][[1]]))))) ==
      0) {
    scHolography.sc[[annotationToUse]][[1]] <- factor(paste("c",
                                                            scHolography.sc[[annotationToUse]][[1]], sep = ""),
                                                      levels = paste("c", levels(scHolography.sc[[annotationToUse]][[1]]),
                                                                     sep = ""))
    query.cluster <- paste("c", query.cluster, sep = "")
    ref.cluster <- paste("c", ref.cluster, sep = "")
  }
  query.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]] %in%
                               query.cluster)
  ref.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]] %in%
                             ref.cluster)
  graph <- igraph::graph_from_adjacency_matrix(scHolography.obj$adj.mtx,
                                               mode = "undirected")
  dist <- igraph::distances(graph, mode = "out")
  clus.dist <- dist[query.cluster.ind, ref.cluster.ind]
  query.to.ref.dis <- (colMeans(apply(clus.dist, 1, sort)[1:n.neighbor,
  ]))

  names(query.to.ref.dis) <- colnames(scHolography.sc)[query.cluster.ind]


  plot.ls<-lapply(geneOI,function(x){
    dat <-data.frame(dist=query.to.ref.dis,expr=as.numeric(scHolography.sc@assays[[assayToUse]]@data[x,query.cluster.ind]))
    #dat$expr <- dat$expr/max(dat$expr)
    suppressWarnings(p<-ggplot(dat, aes(x=dist, y=expr)) + geom_point(alpha=.5) +stat_smooth(method="glm", se=T, method.args = list(family="poisson")) +theme_classic()+ggtitle(as.character(x))
    )
    p
  })
  names(plot.ls)<-geneOI
  plot.ls

}


#' Spatial Dynamics Heatmap
#' @import Seurat
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import igraph
#' @import viridis
#' @import colorspace
#' @import stringr
#' @export

spatialDynamicsFeaturePlot<-function (scHolography.obj, query.cluster, ref.cluster,
                                      geneOI = NULL, annotationToUse = "orig.cluster", assayToUse = "RNA",
                                      heatmapGp = NULL, pal = "Paired", n.neighbor = 30,
                                      heatmap.pal = "viridis")
{
  scHolography.sc <- scHolography.obj$scHolography.sc
  scHolography.sc[[annotationToUse]][[1]] <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]),
                                                           function(x) paste0(strsplit(x, split = " ")[[1]], collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]] <- unlist(lapply(as.character(scHolography.sc[[annotationToUse]][[1]]),
                                                           function(x) paste0(strsplit(x, split = "/")[[1]], collapse = "_")))
  scHolography.sc[[annotationToUse]][[1]] <- factor(scHolography.sc[[annotationToUse]][[1]],
                                                    levels = stringr::str_sort(unique(scHolography.sc[[annotationToUse]][[1]]),
                                                                               numeric = T))
  query.cluster <- unlist(lapply(as.character(query.cluster),
                                 function(x) paste0(strsplit(x, split = " ")[[1]], collapse = "_")))
  query.cluster <- unlist(lapply(as.character(query.cluster),
                                 function(x) paste0(strsplit(x, split = "/")[[1]], collapse = "_")))
  ref.cluster <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,
                                                                                      split = " ")[[1]], collapse = "_")))
  ref.cluster <- unlist(lapply(as.character(ref.cluster), function(x) paste0(strsplit(x,
                                                                                      split = "/")[[1]], collapse = "_")))
  if (heatmap.pal == "rdbu") {
    heat.col <- rev(RColorBrewer::brewer.pal(n = 25, name = "RdBu"))
  }else if (heatmap.pal == "magma") {
    heat.col <- viridis::viridis(25, option = "A")
  }else if (heatmap.pal == "rdbu_1") {
    heat.col <- colorspace::diverge_hsv(25)
  }else {
    heat.col <- viridis::viridis(25)
  }
  if (sum(is.na(as.numeric(as.character(levels(scHolography.sc[[annotationToUse]][[1]]))))) == 0) {
    scHolography.sc[[annotationToUse]][[1]] <- factor(paste("c", scHolography.sc[[annotationToUse]][[1]], sep = ""),
                                                      levels = paste("c", levels(scHolography.sc[[annotationToUse]][[1]]), sep = ""))
    query.cluster <- paste("c", query.cluster, sep = "")
    ref.cluster <- paste("c", ref.cluster, sep = "")
  }
  query.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]] %in%
                               query.cluster)
  ref.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]] %in%
                             ref.cluster)
  graph <- igraph::graph_from_adjacency_matrix(scHolography.obj$adj.mtx,
                                               mode = "undirected")
  dist <- igraph::distances(graph, mode = "out")
  clus.dist <- dist[query.cluster.ind, ref.cluster.ind]
  query.to.ref.dis <- (colMeans(apply(clus.dist, 1, sort)[1:n.neighbor, ]))

  nL <- length(unique(query.to.ref.dis))
  layer.seq <- rep(NA, length(query.cluster.ind))
  for (i in 1:(nL)) {
    l.ind <- which(query.to.ref.dis %in% (sort(unique(query.to.ref.dis))[i]) == T)
    layer.seq[l.ind] <- paste(paste(query.cluster,  collapse = "_"), i, sep = "_")
  }

  query.cluster.sub <- subset(scHolography.sc, cells = c(query.cluster.ind,
                                                         ref.cluster.ind))
  query.cluster.sub$cal.dist <- (c(query.to.ref.dis, rep(NA,
                                                         length(ref.cluster.ind))))
  query.cluster.sub[[paste(paste(query.cluster, collapse = "_"),
                           paste(ref.cluster, collapse = "_"), sep = "To")]] <- factor(c(layer.seq,
                                                                                         as.character(scHolography.sc[[annotationToUse]][[1]][ref.cluster.ind])),
                                                                                       levels = c(stringr::str_sort(unique(layer.seq), numeric = T),
                                                                                                  ref.cluster))
  if(length(intersect(query.cluster,ref.cluster))>0){
    names.ls <-colnames(query.cluster.sub)
    ind.rg <-((length(layer.seq)+1):ncol(query.cluster.sub))
    names.ls[ind.rg]<-paste(names.ls[ind.rg],"_ref",sep = "")
    query.cluster.sub <- Seurat::RenameCells(query.cluster.sub,new.names =names.ls)
  }

  query.cluster.sub <- Seurat::SetIdent(query.cluster.sub,
                                        value = paste(paste(query.cluster, collapse = "_"), paste(ref.cluster,
                                                                                                  collapse = "_"), sep = "To"), sep = "To")
  if (is.null(geneOI)) {
    show("Need the Gene of Interest list if ploting by ordered distance")
  }
  else {
    obj.sub <- subset(query.cluster.sub, cells = which((query.cluster.sub@active.ident %in%
                                                          ref.cluster) == F))
    if(assayToUse=="RNA"){
      Seurat::DefaultAssay(obj.sub) <- "RNA"
      obj.sub <-Seurat::NormalizeData(obj.sub)
      obj.sub <- Seurat::ScaleData(obj.sub)
      obj.sub <- Seurat::FindVariableFeatures(obj.sub,nfeatures = 5000)
    }

    p1 <- DoMultiBarHeatmap(obj.sub, features = geneOI,
                            assay = assayToUse, draw.lines = F, size = 3, disPlot = T, additional.group.by = annotationToUse,
                            label = F, palette = pal) + Seurat::NoLegend() +
      ggplot2::scale_fill_gradientn(colors = heat.col) +
      ggplot2::theme(legend.title = ggplot2::element_blank())
    fea <- Seurat::VariableFeatures(obj.sub)
    obj.sub[[annotationToUse]][[1]] <- droplevels(obj.sub[[annotationToUse]][[1]])
    labal.col <- colorRampPalette(brewer.pal(12, pal))(length(levels(scHolography.sc[[annotationToUse]][[1]])))
    p <- Seurat::DoHeatmap(obj.sub, features = fea[1], label=FALSE,
                           group.by = annotationToUse, group.colors = labal.col[unlist(lapply(levels(obj.sub[[annotationToUse]][[1]]),
                                                                                              function(x) {
                                                                                                which(levels(scHolography.sc[[annotationToUse]][[1]]) %in%
                                                                                                        x)
                                                                                              }))]) + ggplot2::scale_fill_gradientn(colors = colorRampPalette(brewer.pal(12, pal))(100), name = "Distance", labels = c("Proximal", "Distal"), n.breaks = 2)
    grid::grid.newpage()
    legend <- cowplot::get_legend(p)
    (cowplot::plot_grid(p1, NULL, legend, rel_widths = c(1, -0.1, 1), align = "hv", nrow = 1))

  }

}


#' Find Spatial Neighborhood
#' @import Seurat
#' @import dplyr
#' @import RColorBrewer
#' @import ggplot2
#' @import d3heatmap
#' @import textmineR
#' @import factoextra
#' @import RColorBrewer
#' @import stringr
#' @export

findSpatialNeighborhood <-function(scHolography.obj,annotationToUse,query.cluster,nclus=NULL, pal="Set3"){

  docs <- lapply(which(scHolography.obj$scHolography.sc[[annotationToUse]][[1]]%in%query.cluster),function(x){as.vector(scHolography.obj$scHolography.sc[[annotationToUse]][[1]][which(scHolography.obj$adj.mtx[x,]>0)])})
  vocab <- unique(unlist(docs))
  ls <- lapply(vocab, function(x){
    unlist(lapply(docs, function(y){
      sum(y%in%x)
    }))
  })
  dtm <- matrix(unlist(ls), ncol = length(vocab),nrow = length(docs),byrow = F)
  colnames(dtm) <- vocab

  tf_mat <- TermDocFreq(dtm)
  # TF-IDF and cosine similarity
  tf_mat$idf[which(tf_mat$idf==0) ] <- 10^-18
  tfidf <- t(dtm[ , tf_mat$term ]) * tf_mat$idf
  tfidf <- t(tfidf)
  zero.sum <- which(rowSums(tfidf)==0)
  nonzero.sum <- which(rowSums(tfidf)!=0)
  tfidf.sub <- tfidf[nonzero.sum,]
  csim <- tfidf.sub / sqrt(rowSums(tfidf.sub * tfidf.sub))
  csim <- csim %*% t(csim)
  cdist <- as.dist(1 - csim)
  silhouette <- fviz_nbclust((as.matrix(cdist)), kmeans, method = "silhouette")+labs(subtitle = "Silhouette method")
  if(is.null(nclus)){
    nclus <- as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)])
  }
  heatmap <- d3heatmap(as.matrix(cdist),k_col = nclus )
  heatmap.order <- as.numeric(heatmap$x$matrix$rows)
  node.ls <- unlist(heatmap$x$cols)
  entry.names <- unlist(lapply(names(node.ls),function(x){
    temp <-unlist(strsplit(x,split = "[.]"))
    temp[length(temp)]
  } ))
  label <-which(entry.names%in%"label")
  col <- label+1
  label.c <- node.ls[label]
  col.c <- node.ls[col]

  label.c <- as.numeric(label.c)
  names(label.c) <- col.c
  label.c <- sort(label.c)

  sub.color <-colorRampPalette(brewer.pal(12, pal))(length(unique(names(label.c))))
  clustering <- sub.color[as.numeric(as.factor(names(label.c)))]
  heat.matrix <- heatmap$x$params$x
  heat.matrix.reorder <- heat.matrix[heatmap.order,heatmap.order]
  clustering.reorder <- clustering[heatmap.order]
  heatmap.out <- heatmap(heat.matrix.reorder, Colv = NA, Rowv = NA, RowSideColors=clustering.reorder, ColSideColors=clustering.reorder, scale="none",revC=T,labRow=F,labCol=F)


  sub.clus <- rep("NA",ncol(scHolography.obj$scHolography.sc))
  sub.clus <- as.character(scHolography.obj$scHolography.sc[[annotationToUse]][[1]])
  sub.clus[which(scHolography.obj$scHolography.sc[[annotationToUse]][[1]]%in%query.cluster)] <- clustering
  names(sub.clus) <-colnames(scHolography.obj$scHolography.sc)
  scHolography.obj$scHolography.sc[["sub.clus"]] <- sub.clus


  annotationToUse2 <- "sub.clus"
  clus <- scHolography.obj$scHolography.sc[[annotationToUse]][[1]]
  uni.clus <- stringr::str_sort(unique(clustering),numeric = T)

  sig.class.ls <- lapply(uni.clus, function(ind.clus){

    q.cluster <- list(c(ind.clus),uni.clus[-which(uni.clus%in%c(ind.clus))])
    matrix.ls <-vector("list",2)

    for (i in 1:2) {
      ind <- which(scHolography.obj$scHolography.sc[[annotationToUse2]][[1]]%in%q.cluster[[i]])

      matrix.this<- matrix(unlist(lapply(stringr::str_sort(unique(clus),numeric = T), function(x){
        this.clus <- which(clus==x)
        rowSums(scHolography.obj$adj.mtx[,this.clus])
      })),ncol=length(stringr::str_sort(unique(clus),numeric = T)),byrow=F)
      colnames(matrix.this) <- stringr::str_sort(unique(clus),numeric = T)
      matrix.ls[[i]] <- matrix.this[ind,]

    }

    pval.ls <- lapply(stringr::str_sort(unique(clus),numeric = T),function(x){
      wilcox.test(matrix.ls[[1]][,x],matrix.ls[[2]][,x],alternative="greater")$p.value
    })
    names(pval.ls) <- stringr::str_sort(unique(clus),numeric = T)
    unlist(pval.ls)[names(which(unlist(pval.ls)<.05))]
  })

  names(sig.class.ls) <- uni.clus

  list(scHolography.obj=scHolography.obj,summary=sig.class.ls,cluster=clustering,heatmap=heatmap.out,cdist=cdist,order=heatmap.order)
}



















#' Find Space Driver Gene
#' @import dplyr
#' @import pracma
#' @import stats
#' @import enrichR

findDriverGene <- function(scHolography.obj,query.cluster,ref.cluster,k1=1.96,k2=1.96,annotationToUse="orig.cluster",n.neighbor=30, assayToUse="SCT",bandwidth=NULL){
  scHolography.sc<-scHolography.obj$scHolography.sc
  query.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%query.cluster)
  ref.cluster.ind <- which(scHolography.sc[[annotationToUse]][[1]]%in%ref.cluster)
  graph <- igraph::graph_from_adjacency_matrix(scHolography.obj$adj.mtx,mode = "undirected")
  dist <- igraph::distances(graph, mode="out")
  clus.dist <- dist[query.cluster.ind,ref.cluster.ind]
  query.to.ref.dis <- (colMeans(apply(clus.dist, 1, sort)[1:n.neighbor,]))

  if(is.null(bandwidth)){
    default.dist <- hist(query.to.ref.dis,plot = F)
    bandwidth <- default.dist$breaks[2]-default.dist$breaks[1]
  }
  sub.query.sc.obj<-subset(scHolography.sc, cells = query.cluster.ind)
  p<-hist(query.to.ref.dis,breaks =seq(from=0,to=max(ceiling(query.to.ref.dis/bandwidth)*bandwidth),by=bandwidth),plot = F)
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
  # setEnrichrSite("Enrichr") # Human genes
  # dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
  # enriched <- enrichr(rownames(regression)[which(regression$`z value`<(-k1))], dbs)
  # enriched[["GO_Biological_Process_2021"]]
  # show(plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = "Proximal Gene Ontology Analysis")
  # )
  # enriched <- enrichr(rownames(regression)[which(regression$`z value`>(k2))], dbs)
  # enriched[["GO_Biological_Process_2021"]]
  # show(plotEnrich(enriched[[3]], showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", title = "Distal Gene Ontology Analysis")
  # )
  list(regression=regression,avg.expr.sub=avg.expr.sub,count=count,bandwidth=bandwidth)

}

#' Visualization of Space Driver Genes
#' @import ggplot2

driverGenePlot <- function(findDG.obj,gene){
  avg.expr.sub <- findDG.obj[[2]]
  count <- findDG.obj[[3]]
  bandwidth <- findDG.obj[[4]]
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

