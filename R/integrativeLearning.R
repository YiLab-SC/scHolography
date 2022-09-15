
#' scHolography Joint Learning from Multiple ST references
#' @export
#' @import matchingR
#' @import igraph

scHolographySpatialIntegration<-function(spatial.int.list,nslot=30,seed=60611){

  sc.dist.sum<-matrix(0, nrow=nrow(spatial.int.list[[1]]$adj.mtx),ncol =  ncol(spatial.int.list[[1]]$adj.mtx))
  motility<-0

  for (i in 1:length(spatial.int.list)) {
    indS<-1+ncol(spatial.int.list[[i]]$scHolography.sp)
    indE<-ncol(spatial.int.list[[i]]$scHolography.sp)+ncol(spatial.int.list[[i]]$scHolography.sc)
    adj.tensor<-spatial.int.list[[i]]$est.array[indS:indE,indS:indE]
    adj.tensor.dia <- (adj.tensor+t(adj.tensor))/2
    adj.tensor.dia.rowM<-mean(adj.tensor.dia)
    adj.tensor.dia.rowSD<-sqrt(var(as.vector(adj.tensor.dia)))
    adj.tensor.dia.norm<-adj.tensor.dia-adj.tensor.dia.rowM
    adj.tensor.dia.norm<-adj.tensor.dia.norm/adj.tensor.dia.rowSD
    sc.dist.sum<-sc.dist.sum+adj.tensor.dia.norm
    motility <- motility +exp(spatial.int.list[[i]]$scHolography.sc$motility)
  }
  adj.tensor.ut<- -sc.dist.sum
  diag(adj.tensor.ut)<-0
  assign.matrix<-matrix(nrow = nrow(adj.tensor.ut),ncol = nslot)
  for (k in 1:nslot) {
    adj.tensor_results <- matchingR::galeShapley.collegeAdmissions(collegeUtils =adj.tensor.ut, studentUtils =adj.tensor.ut,slot=1, studentOptimal = T)
    assign.matrix[,k] <-adj.tensor_results$matched.colleges[,1]
    for (i in 1:nrow(adj.tensor.ut)) {
      adj.tensor.ut[i,adj.tensor_results$matched.colleges[i,1]]<-0
      adj.tensor.ut[adj.tensor_results$matched.colleges[i,1],i]<-0
    }
  }

  adj.mtx <- matrix(0, ncol = ncol(adj.tensor), nrow = nrow(adj.tensor))
  for (i in 1:nrow(adj.tensor)) {
    adj.mtx[i,assign.matrix[i,]]<-1
  }

  g <- igraph::graph.adjacency(adj.mtx, mode="undirected")
  set.seed(seed)
  l <- igraph::layout_with_fr(g)
  l_3d <- igraph::layout_with_fr(g, dim=3)

  sc.in.int <- spatial.int.list[[1]]$scHolography.sc
  sc.in.int[["x_sp"]]<-l[,1]
  sc.in.int[["y_sp"]]<-l[,2]
  sc.in.int[["x3d_sp"]]<-l_3d[,1]
  sc.in.int[["y3d_sp"]]<-l_3d[,2]
  sc.in.int[["z3d_sp"]]<-l_3d[,3]

  sc.in.int[["motility"]]<-log(motility/length(spatial.int.list))

  sp.in.int.ls<-lapply(spatial.int.list, function(x){
    x$scHolography.sp
  })

  spatial.int.list[[1]]$scHolography.sc<-sc.in.int
  spatial.int.list[[1]]$adj.mtx<-adj.mtx
  spatial.int.list[[1]]$est.array<-sc.dist.sum
  spatial.int.list[[1]]$scHolography.sp<-sp.in.int.ls
  spatial.int.list[[1]]
}
