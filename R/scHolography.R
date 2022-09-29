
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
  print("Data Normalized")
  seuratObj
}



#' Data Integration
#' @export
#' @import Seurat
#'
dataAlign<-function(low.res.sp,high.res.sp, stProcessed=F, scProcessed=F, nFeatureToUse=3000, whichReference=2,nPCtoUse=32){
  if(stProcessed==0){low.res.sp<-seuratSCT(low.res.sp,ndims = nPCtoUse)}
  if(scProcessed==0){high.res.sp<-seuratSCT(high.res.sp,ndims = nPCtoUse)}
  options(future.globals.maxSize = 4000 * 1024^2)
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

### Generate Data for Processing
getData <- function(sp.integrated,nPCtoOut){
  sp.in.int<-subset(sp.integrated,cells=which(sp.integrated$type.assay=="sp"))
  rna.in.int<-subset(sp.integrated,cells=which(sp.integrated$type.assay=="sc"))
  sp.coord<-sp.in.int@images[[1]]@coordinates[,c('row', 'col')]
  pixel.dist<-dist(sp.coord)
  pixel.dist<-as.matrix(pixel.dist)
  spatial.pca <- prcomp(pixel.dist, scale = TRUE)
  spatial.pca.32 <- spatial.pca$x[,1:nPCtoOut]
  spatial.pca.32.submin <- t(t(spatial.pca.32)-apply(spatial.pca.32,MARGIN = 2,min))
  spatial.pca.32.submin.divmax <- t(t(spatial.pca.32.submin)/apply(spatial.pca.32.submin,MARGIN = 2,max))

  train_X <- (as.matrix(sp.in.int@reductions$pca@cell.embeddings))
  test_X <- (as.matrix(rna.in.int@reductions$pca@cell.embeddings))
  train_Y <- spatial.pca.32.submin.divmax
  list(train_X,test_X,train_Y)
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
    optimizer = keras::optimizer_adam(learning_rate = 0.001,decay = 0.00005),
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
    show(paste(i,"th repeat run",sep = ""))
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

  show("Complete of NN step")
  corr.tensor
}



#' Decomposition Analysis
#' @import matchingR
interpretDecomp<-function(est.array, nslot, is3D=T){
  if(is3D){
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
    var.vec<-unlist(lapply(1:dim(est.array)[2], function(x){
      sum(unlist(lapply(1:dim(est.array)[2],function(ent){var(est.array[,x,ent])})))

    }))
    show("Complete of intepretation step")
    list(adj.mtx,var.vec)
  }else{
    show("Complete of intepretation step")
    list(adj.mtx,(-cell2cell.dia.norm))
  }
}



#' Single cell 3D reconstruction
#' @export
#' @import igraph
trainHolography<-function(sp.integrated,vSeed=60611,n.repeat=50,n.slot=30,n.pcUse=32,n.pcOut=32){
  dat.list<-getData(sp.integrated,nPCtoOut = n.pcOut)
  train_y<-dat.list[[3]]
  test_x<-dat.list[[2]]
  train_x<-dat.list[[1]]
  x_train <- as.matrix(train_x)
  x_test <- as.matrix(test_x)
  rm(dat.list)
  rm(train_x)
  rm(test_x)
  cp_decomp<-applyModel(x_train,x_test,train_y,seed = vSeed,n.PCtoUse = n.pcUse,n.PCtoOut = n.pcOut,nRepeat = n.repeat)
  indS<-nrow(x_train)+1
  indE<-nrow(x_train)+nrow(x_test)
  cp_decomp_sc<-cp_decomp[,indS:indE,indS:indE]
  interpretDecomp.rslt <- interpretDecomp(cp_decomp_sc,nslot = n.slot)
  adj.mtx <- interpretDecomp.rslt[[1]]
  var.vec <- interpretDecomp.rslt[[2]]

  g <- igraph::graph.adjacency(adj.mtx, mode="undirected")
  set.seed(vSeed)
  l <- igraph::layout_with_fr(g)
  l_3d <- igraph::layout_with_fr(g, dim=3)

  sc.in.int <- subset(sp.integrated,cells=which(sp.integrated$type.assay=="sc"))
  sp.in.int<-subset(sp.integrated,cells=which(sp.integrated$type.assay=="sp"))
  sc.in.int[["x_sp"]]<-l[,1]
  sc.in.int[["y_sp"]]<-l[,2]
  sc.in.int[["x3d_sp"]]<-l_3d[,1]
  sc.in.int[["y3d_sp"]]<-l_3d[,2]
  sc.in.int[["z3d_sp"]]<-l_3d[,3]
  sc.in.int[["motility"]]<-log(var.vec)
  list(scHolography.sc=sc.in.int, scHolography.sp=sp.in.int,adj.mtx=adj.mtx,est.array=colMeans(cp_decomp))
}



