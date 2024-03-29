
#' Data Normalization
#' @param  seuratObj Seurat Object
#' @param  ndims Number of PCs used
#' @import Seurat
seuratSCT<-function(seuratObj,ndims ){
  if("Spatial"%in%names(seuratObj@assays)){
    assayUsed <- "Spatial"
  }else{
    assayUsed <- "RNA"
  }
  seuratObj <- Seurat::SCTransform(seuratObj, assay = assayUsed, verbose = FALSE)
  seuratObj <- Seurat::RunPCA(seuratObj, assay = "SCT", verbose = FALSE)
  seuratObj <- Seurat::RunUMAP(seuratObj,dims=1:ndims,verbose = FALSE)
  seuratObj <- Seurat::FindNeighbors(seuratObj, reduction = "pca", dims = 1:ndims,verbose = FALSE)
  seuratObj <- Seurat::FindClusters(seuratObj,resolution = 0.5 ,verbose = FALSE)
  cat("Data Normalized \n")
  seuratObj
}



#' Data Integration
#' @export
#' @import Seurat
#' @param  low.res.sp Seurat Object of ST (low resolution) data
#' @param  high.res.sp Seurat Object of SC (high resolution) data
#' @param  stProcessed Is low.res.sp processed? Default is False
#' @param  scProcessed Is low.res.sp processed? Default is False
#' @param  nFeatureToUse Number of features to use for downstream analysis (e.g., PCA). Default is 3000
#' @param  whichReference Use which dataset as the reference for Seurat integration. Default using SC data (whichReference=2). Set whichReference=1 to use ST data as the integration reference
#' @param  nPCtoUse Number of PCs used for analysis. Default is 32
#'
dataAlign<-function(low.res.sp,high.res.sp, stProcessed=F, scProcessed=F, nFeatureToUse=3000, whichReference=2,nPCtoUse=32,future.size=4000){
  if(stProcessed==0){low.res.sp<-seuratSCT(low.res.sp,ndims = nPCtoUse)}
  if(scProcessed==0){high.res.sp<-seuratSCT(high.res.sp,ndims = nPCtoUse)}
  options(future.globals.maxSize = future.size * 1024^2)
  low.res.sp$orig.cluster<-low.res.sp$seurat_clusters
  high.res.sp$orig.cluster<-high.res.sp$seurat_clusters
  low.res.sp[["type.assay"]]<-"sp"
  high.res.sp[["type.assay"]]<-"sc"
  integrate.list<-list(low.res.sp,high.res.sp)
  comp.umap <- merge(x=low.res.sp@reductions$umap,y=high.res.sp@reductions$umap)
  rm(low.res.sp)
  rm(high.res.sp)
  features <- Seurat::SelectIntegrationFeatures(object.list = integrate.list,nfeatures = nFeatureToUse, verbose = FALSE)
  integrate.list <- Seurat::PrepSCTIntegration(object.list = integrate.list, anchor.features = features, verbose = FALSE)
  sp.anchors <- Seurat::FindIntegrationAnchors(object.list = integrate.list, anchor.features = features,reference = whichReference, normalization.method = "SCT", verbose = FALSE)
  rm(integrate.list)
  sp.integrated <- Seurat::IntegrateData(anchorset = sp.anchors, normalization.method = "SCT", verbose = FALSE)
  sp.integrated <- Seurat::ScaleData(sp.integrated, verbose = FALSE)
  sp.integrated <- Seurat::RunPCA(sp.integrated, npcs = nPCtoUse, verbose = FALSE)
  sp.integrated <- Seurat::RunUMAP(sp.integrated, reduction = "pca", dims = 1:nPCtoUse, verbose = FALSE)
  sp.integrated <- Seurat::FindNeighbors(sp.integrated, reduction = "pca", dims = 1:nPCtoUse, verbose = FALSE)
  sp.integrated <- Seurat::FindClusters(sp.integrated, verbose = FALSE)
  sp.integrated$orig.cluster <-as.factor(as.numeric(sp.integrated$orig.cluster))
  sp.integrated@reductions[["ind.umap"]] <- comp.umap
  sp.integrated
}


#' Generate Data for Processing
#' @import SeuratObject
#' @import Seurat
#' @import gmodels
#'
getData <- function (sp.integrated, nPCtoOut, is.fov = F, fov)
{
  this.fov <- fov
  if (length(this.fov) > 1) {
    this.fov <- this.fov[1]
  }
  if (is.fov == F) {
    if("imagerow"%in%colnames(sp.integrated@images[[1]]@coordinates)){
      sp.coord <- sp.integrated@images[[1]]@coordinates[, c("imagerow","imagecol")]
    }else if("x"%in%colnames(sp.integrated@images[[1]]@coordinates)){
      sp.coord <- sp.integrated@images[[1]]@coordinates[, c("x","y")]
    }

  }else {
    sp.coord <- ( SeuratObject::GetTissueCoordinates(sp.integrated[[as.character(this.fov)]][["centroids"]]))[, c("x", "y")]
  }
  cat("Get Data \n")
  pixel.dist <- dist(sp.coord)
  pixel.dist <- as.matrix(pixel.dist)
  spatial.pca <-gmodels::fast.prcomp(pixel.dist, scale = TRUE)
  spatial.pca.32 <- spatial.pca$x[, 1:nPCtoOut]
  spatial.pca.32.submin <- t(t(spatial.pca.32) - apply(spatial.pca.32,MARGIN = 2, min))
  spatial.pca.32.submin.divmax <- t(t(spatial.pca.32.submin)/apply(spatial.pca.32.submin,
                                                                   MARGIN = 2, max))
  train_X <- as.matrix(sp.integrated@reductions$pca@cell.embeddings[which(sp.integrated$type.assay == "sp"), ])
  if (is.fov) {
    train_X <- train_X[which(names(which(sp.integrated$type.assay == "sp"))%in% GetTissueCoordinates(sp.integrated[[as.character(this.fov)]][["centroids"]])[, "cell"]), ]
  }
  test_X <- as.matrix(sp.integrated@reductions$pca@cell.embeddings[which(sp.integrated$type.assay == "sc"), ])
  train_Y <- spatial.pca.32.submin.divmax
  list(train_X, test_X, train_Y)
}




#' Build Neural Network Model
#' @import dplyr
#' @import keras
build_model <- function(nPCtoUse, nPCtoOut) {
  model <- keras::keras_model_sequential()
  model %>%
    keras::layer_dense(units = 32, input_shape = c(nPCtoUse)) %>%
    keras::layer_dropout(0.2) %>%
    keras::layer_activation_leaky_relu(alpha = 0.2) %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dense(units = 32) %>%
    keras::layer_dropout(0.2) %>%
    keras::layer_activation_leaky_relu(alpha = 0.2) %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dense(units = 8) %>%
    keras::layer_dropout(0.2) %>%
    keras::layer_activation_leaky_relu(alpha = 0.2) %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dense(units = 32) %>%
    keras::layer_dropout(0.2) %>%
    keras::layer_activation_leaky_relu(alpha = 0.2) %>%
    keras::layer_batch_normalization() %>%
    keras::layer_dense(units = nPCtoOut) %>%
    keras::layer_activation_relu()

  model %>% keras::compile(
    loss = "mse",
    optimizer = keras::keras$optimizers$legacy$Adam(learning_rate = 0.001,decay = 0.00005),
    metrics = list("mean_absolute_error")
  )
  model
}


#' Apply the NN Model
#' @import keras
#' @import tensorflow
applyModel <- function(x_train, x_test, train_y, nRepeat ,seed,n.PCtoUse,n.PCtoOut){
  corr.tensor <- array(NA, dim = c(nRepeat, nrow(x_train)+nrow(x_test), nrow(x_train)+nrow(x_test)))
  early_stop <- keras::callback_early_stopping(monitor = "val_loss", patience = 20)
  tensorflow::set_random_seed(seed)

  for (i in 1:nRepeat) {
    cat(paste(i,"th repeat run \n",sep = ""))
    y_train <- as.matrix(train_y)
    model <- build_model(nPCtoUse = n.PCtoUse,nPCtoOut = n.PCtoOut)
    history <- model %>% keras::fit(
      x = x_train,
      y = y_train,
      epochs = 500,
      validation_split = 0.2,
      verbose = 0,
      callbacks = list(early_stop)
    )
    show(plot(history))
    train_predictions <- model %>% predict(x_train)
    test_predictions <- model %>% predict(x_test)
    train_test_comb <- rbind(train_predictions,test_predictions)
    cormat <- as.matrix(dist(train_test_comb))
    corr.tensor[i,,]<-cormat

  }




  cat("Complete of NN step \n")
  list(corr.tensor=corr.tensor,test_predictions=test_predictions)
}



#' Stable Matching Neighbor Analysis
#' @import matchingR
#' @export
#' @param  est.array Distance matrix between individual cells/pixels
#' @param  nslot Maximum number of neighbors to find for stable matching
#' @param  is3D If est.array is a 3D tensor. Default is True
interpretDecomp<-function(est.array, nslot, is3D=T,isFilter){
  if(is3D){

    var.vec<-unlist(lapply(1:dim(est.array)[2], function(x){
      sum(unlist(lapply(1:dim(est.array)[2],function(ent){var(est.array[,x,ent])})))
    }))


    filter.index <- NULL
    if(isFilter){
      filter.index <- which(var.vec<min(boxplot.stats(var.vec)$out))
      est.array <- est.array[, filter.index, filter.index]

    }


    flatten.array <-apply(est.array,c(1),function(x){
      cell2cell.dia <- (x+t(x))/2
      cell2cell.dia.rowM<-mean(cell2cell.dia)
      cell2cell.dia.rowSD<-sqrt(var(as.vector(cell2cell.dia)))
      cell2cell.dia.norm<-cell2cell.dia-cell2cell.dia.rowM
      cell2cell.dia.norm<-cell2cell.dia.norm/cell2cell.dia.rowSD
      cell2cell.dia.norm
    })

    cell2cell.dia.norm<- matrix(rowMeans(flatten.array),nrow = ncol(est.array),byrow = F)

  }else{
    cell2cell<-est.array

    cell2cell.dia <- (cell2cell+t(cell2cell))/2
    cell2cell.dia.rowM<-mean(cell2cell.dia)
    cell2cell.dia.rowSD<-sqrt(var(as.vector(cell2cell.dia)))
    cell2cell.dia.norm<-cell2cell.dia-cell2cell.dia.rowM
    cell2cell.dia.norm<-cell2cell.dia.norm/cell2cell.dia.rowSD

  }
  cell2cell.ut<- -cell2cell.dia.norm
  diag(cell2cell.ut)<-(-99999)

  assign.matrix<-matrix(nrow = nrow(cell2cell.ut),ncol = nslot)


  for (k in 1:nslot) {
    cell2cell_results <- matchingR::galeShapley.collegeAdmissions(collegeUtils =cell2cell.ut, studentUtils =cell2cell.ut,slot=1, studentOptimal = T)
    assign.matrix[,k] <-cell2cell_results$matched.colleges[,1]
    for (i in 1:nrow(cell2cell.ut)) {
      cell2cell.ut[i,cell2cell_results$matched.colleges[i,1]]<-(-99999)
      cell2cell.ut[cell2cell_results$matched.colleges[i,1],i]<-(-99999)
    }
  }

  adj.mtx <- matrix(0, ncol = ncol(cell2cell.ut), nrow = nrow(cell2cell.ut))
  for (i in 1:nrow(cell2cell.ut)) {
    adj.mtx[i,assign.matrix[i,]]<-1
  }
  rm(cell2cell.ut)
  if(is3D){
    cat("Complete of intepretation step \n")
    list(adj.mtx,var.vec,filter.index)
  }else{
    cat("Complete of intepretation step \n")
    list(adj.mtx,(-cell2cell.dia.norm))
  }
}



#' Single cell 3D reconstruction
#' @export
#' @import igraph
#' @param  sp.integrated Seurat Object of integrated SC and ST data
#' @param  n.repeat Number of neural network models to train. Default is 30
#' @param  n.slot Maximum number of neighbors to find for stable matching. Default is 30
#' @param  n.pcUse Number of PCs used for NN input. Default is 32
#' @param  n.pcOut Number of PCs used for NN output. Default is 32
#' @param  vSeed Random seed. Default is 60611
#'
trainHolography<-function (sp.integrated, fov = NULL, which.fov = NULL, n.repeat = 30,
                           n.slot = 30, n.pcUse = 32, n.pcOut = 32, vSeed = 60611,filter=T)
{
  dat.list <- getData(sp.integrated, nPCtoOut = n.pcOut, is.fov = (1 %in% which.fov), fov = fov)
  train_y <- dat.list[[3]]
  test_x <- dat.list[[2]]
  train_x <- dat.list[[1]]
  x_train <- as.matrix(train_x)
  x_test <- as.matrix(test_x)
  rm(dat.list)
  rm(train_x)
  rm(test_x)
  apply.model.rslt <- applyModel(x_train, x_test, train_y, seed = vSeed, n.PCtoUse = n.pcUse, n.PCtoOut = n.pcOut, nRepeat = n.repeat)
  cp_decomp <- apply.model.rslt[[1]]
  indS <- nrow(x_train) + 1
  indE <- nrow(x_train) + nrow(x_test)
  cp_decomp_sc <- cp_decomp[, indS:indE, indS:indE]
  interpretDecomp.rslt <- interpretDecomp(cp_decomp_sc, nslot = n.slot,isFilter = filter)
  adj.mtx <- interpretDecomp.rslt[[1]]
  var.vec <- interpretDecomp.rslt[[2]]
  filter.index <- interpretDecomp.rslt[[3]]
  g <- igraph::graph.adjacency(adj.mtx, mode = "undirected")
  set.seed(vSeed)
  l <- igraph::layout_with_fr(g)
  l_3d <- igraph::layout_with_fr(g, dim = 3)
  if (is.null(which.fov)) {
    sc.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay ==  "sc"))
    if(is.null(filter.index)==F){
      sc.in.int <- subset(sc.in.int, cells = filter.index)
    }
    sp.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay == "sp"))
  }else if(length(which.fov)==2){
    fov.ls.1 <- sp.integrated[[as.character(fov[1])]]
    fov.ls.2 <- sp.integrated[[as.character(fov[2])]]
    sp.integrated[[as.character(fov[1])]] <- NULL
    sp.integrated[[as.character(fov[2])]] <- NULL
    sc.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay ==  "sc"))
    sp.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay == "sp"))
    sp.in.int[[as.character(fov[1])]] <- fov.ls.1
    sc.in.int[[as.character(fov[2])]] <- fov.ls.2
    if(is.null(filter.index)==F){
      sc.in.int <- subset(sc.in.int, cells = filter.index)
    }
  }else if (which.fov == 2) {
    sc.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay ==  "sc"))
    sp.integrated[[as.character(fov)]] <- NULL
    sp.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay == "sp"))
    if(is.null(filter.index)==F){
      sc.in.int <- subset(sc.in.int, cells = filter.index)
    }
  }else if (which.fov == 1) {
    sp.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay == "sp"))
    sp.integrated[[as.character(fov)]] <- NULL
    sc.in.int <- subset(sp.integrated, cells = which(sp.integrated$type.assay == "sc"))
    if(is.null(filter.index)==F){
      sc.in.int <- subset(sc.in.int, cells = filter.index)
    }
  }
  sc.in.int[["x_sp"]] <- l[, 1]
  sc.in.int[["y_sp"]] <- l[, 2]
  sc.in.int[["x3d_sp"]] <- l_3d[, 1]
  sc.in.int[["y3d_sp"]] <- l_3d[, 2]
  sc.in.int[["z3d_sp"]] <- l_3d[, 3]

  if(is.null(filter.index)){
    sc.in.int[["motility"]] <- (var.vec)
  }else{
    sc.in.int[["motility"]] <- (var.vec)[filter.index]
    }

  list(scHolography.sc = sc.in.int, scHolography.sp = sp.in.int,  adj.mtx = adj.mtx, est.array = colMeans(cp_decomp),raw.prediction = apply.model.rslt[[2]],filter.index=filter.index)
}



#' Find scHolography Inferred Distance between Cell Clusters
#' @export
#' @import igraph
#' @import stringr
#' @import Seurat
#' @param  scHolography.obj scHolography object list
#' @param  query.cluster A vector of query identity types
#' @param  ref.cluster A vector of reference identity types
#' @param  annotationToUse Which annotation to call identities from. Default is orig.cluster
#' @param  n.neighbor Number of nearest cells to use to define distance. Default is 30

findDistance<- function (scHolography.obj, query.cluster, ref.cluster, annotationToUse = "orig.cluster",n.neighbor = 30)
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
  query.to.ref.dis
}

#' Reconstruct scHolography graph from distance matrix
#' @export
#' @import igraph
#' @import Seurat
#' @param  high.res.sp Seurat Object of SC (high resolution) data
#' @param  low.res.sp Seurat Object of ST (low resolution) data
#' @param  dist.mat Estimated spatial distance matrix between cells
distToscHolography <- function(high.res.sp,low.res.sp,dist.mat){
  dist.mat.ls <-interpretDecomp(dist.mat,is3D = F,nslot = 30)
  adj.mtx <- dist.mat.ls[[1]]
  g <- igraph::graph.adjacency(adj.mtx, mode="undirected")
  set.seed(60611)
  l_3d <- igraph::layout_with_fr(g, dim=3)
  inf.obj <- list(scHolography.sc=high.res.sp, scHolography.sp=low.res.sp,adj.mtx=adj.mtx,est.array=dist.mat)
  inf.obj$scHolography.sc$x3d_sp<-l_3d[,1]
  inf.obj$scHolography.sc$y3d_sp<-l_3d[,2]
  inf.obj$scHolography.sc$z3d_sp<-l_3d[,3]
  inf.obj
}

