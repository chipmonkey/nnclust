nnfind <-
function(from,to){
  nfrom<-as.integer(NROW(from))
  np<-as.integer(NCOL(from))
  from<-as.double(from)
  if (missing(to) || is.null(to)){
    within<-TRUE
  }else{
    within<-FALSE
    nto<-as.integer(NROW(to))
    to<-as.double(to)
  }
  if (within)
    rval<-.C("within_neighbours",from,count=nfrom,np,neighbour=integer(nfrom),dist=double(nfrom))
  else
    rval<-.C("between_neighbours",from,count=nfrom,to,nto,np,neighbour=integer(nto),dist=double(nto))
  rval$neighbour<-rval$neighbour+1 ## C code has zero-based arrays
  class(rval)<-"neighbours"
  rval$call<-sys.call()
  rval$dim<-np
  rval[c("count","neighbour","dist")]
}

print.neighbours<-function(x,...){
  cat("Nearest-neighbours in",x$dim," dimensions:")
  print(x$call)
  invisible(x)
}

summary.neighbours<-function(object,...){
  rval<-list(nn=object)
  rval$distsummary<-summary(object$dist)
  class(rval)<-"summary.neighbours"
  rval
}

print.summary.neighbours<-function(x,...){
  print(x$nn)
  cat("Distances\n")
  print(x$distsummary)
  invisible(x)
}
