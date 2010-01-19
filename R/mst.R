
mst<-function(X,rebuild=sqrt(nrow(X))/4){

  redo<-as.integer((1:rebuild)*nrow(X)/(rebuild+1))
  rval<-.Call("call_primq",X, redo)
  rval[[1]]<-rval[[1]][-1]+1
  rval[[2]]<-rval[[2]][-1]+1
  rval[[3]]<-rval[[3]][-1]
  names(rval)<-c("from","to","dist")
  rval$n<-NROW(X)-1
  rval$p<-NCOL(X)
  class(rval)<-"mst"
  rval
}
mst_restart<-function(X,rebuild=sqrt(nrow(X))/4,threshold=Inf,start=NULL){

  redo<-as.integer((1:rebuild)*nrow(X)/(rebuild+1))
  if(is.null(start)) start<-as.integer(sample(nrow(X),1))-1
  rval<-.Call("call_primq_restart",X, redo,threshold, as.integer(start))
  n<-min(which(rval[[1]][-1]==-1))
  if (n==Inf) n<-nrow(X)
  rval[[1]]<-rval[[1]][2:n]+1
  rval[[2]]<-rval[[2]][2:n]+1
  rval[[3]]<-rval[[3]][2:n]
  names(rval)<-c("from","to","dist")
  rval$n<-n
  rval$p<-NCOL(X)
  class(rval)<-c("mst")
  rval
}

print.mst<-function(x,...){
  cat("Minimum spanning tree on", x$n,"points in",x$p,"dimensions.\n")
  invisible(x)
}
