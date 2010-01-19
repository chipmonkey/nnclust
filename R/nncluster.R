
##clustering by restarted mst

#library(nnclust)

nncluster<-function(x, threshold, fill=0.95, maxclust=20, give.up=500,verbose=FALSE,start=NULL){
  forest<-list()
  i<-1
  n<-nrow(x)
  m<-0
  if (is.null(start)){
    nn<-nnfind(x)
    start<-which.min(nn$dist)
  }
  rows<-1:n
  repeat({
    if (length(threshold)<i)
      thresh<-threshold[length(threshold)]
    else
      thresh<-threshold[i]
    tree<-list(mst=mst_restart(x, start=start-1, threshold=thresh))
    tree$x<-x[c(start,tree$mst$to),,drop=FALSE]
    tree$rows<-rows[c(start,tree$mst$to)]
    rows<-rows[-c(start,tree$mst$to)]
    x<-x[-c(start,tree$mst$to),,drop=FALSE]
    m<-m+nrow(tree$x)
    forest[[i]]<-tree
    nn<-nnfind(x)
    if (!any(nn$dist<thresh)) break;
    if (m/n > fill) break;
    if (sum(nn$dist<thresh)<give.up) break;
    if (i==maxclust) break;
    if (verbose)
       print(c(m, nrow(x), sum(nn$dist<thresh)))
    i<-i+1
    start<-which.min(nn$dist)
  })
  forest[[i+1]]<-list(x=x,rows=rows)
  attr(forest,"threshold")<-threshold
  class(forest)<-"nncluster"
  forest		     
}

trimCluster<-function(nnclust, size=10){
  n<-length(nnclust)
  sizes<-sapply(nnclust[-n], function(t) t$mst$n)
  small<-sizes<size
  if (all(small)) stop(paste("There are no clusters larger than ",size))
  if (any(small)){
    tmp<-nnclust[which(small)]
    smallx<-do.call(rbind,lapply(tmp, function(t) t$x))
    smallrows<-do.call(c,lapply(tmp, function(t) t$rows))
    nnclust[[n]]$x<-rbind(nnclust[[n]]$x,smallx)
    nnclust[[n]]$rows<-c(nnclust[[n]]$rows,smallrows)
    nnclust<-nnclust[-which(small)]
    class(nnclust)<-"nncluster"
  }
  nnclust
}

print.nncluster<-function(x,...){
  cat("MST-based clustering in ",NCOL(x[[1]]$x)," dimensions\n")
  cat("Cluster sizes:",sapply(x[-length(x)], function(t) t$mst$n),"\n")
  cat(" and ",NROW(x[[length(x)]]$x)," outliers\n")
  d<-attr(x,"threshold")
  if (length(d)==1){
    cat("Threshold =",d,"\n")
  }else{
    cat("Thresholds")
    print(d)
  }
  invisible(x)
}

clusterMember<-function(nnclust,outlier=TRUE){
  allrows<-lapply(nnclust, function(t) t$rows)
  n<-max(sapply(allrows,NROW))
  index<-integer(n)
  for(i in 1:length(allrows)){
    index[allrows[[i]]]<-i
  }
  if (!outlier) index[index==i]<-NA
  index
}


nearestCluster<-function(nnclust, threshold=Inf,outlier=FALSE){
  incluster<-clusterMember(nnclust, outlier=FALSE)
  m<-length(nnclust)
  indata<-do.call(rbind, lapply(nnclust[-m], function(cluster) cluster$x))
  inrows<-do.call(c, lapply(nnclust[-m], function(cluster) cluster$rows))
  dists<-nnfind(indata,nnclust[[m]]$x)
  dists$neighbour[dists$dist>threshold]<-NA
  nearest<- incluster[inrows[dists$neighbour]]
  incluster[nnclust[[m]]$rows]<-nearest
  if (outlier)
    incluster[is.na(incluster)]<-max(incluster,na.rm=TRUE)+1
  incluster
}
