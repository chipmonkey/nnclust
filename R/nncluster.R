
##clustering by restarted mst

#library(nnclust)

nncluster<-function(x, threshold, fill=0.95, maxclust=20, give.up=500,verbose=FALSE){
  forest<-list()
  i<-1
  n<-nrow(x)
  m<-0
  start<-1
  rows<-1:n
  repeat({
    tree<-list(mst=mst_restart(x, start=start-1, threshold=threshold))
    tree$x<-x[c(start,tree$mst$to),,drop=FALSE]
    tree$rows<-rows[c(start,tree$mst$to)]
    rows<-rows[-c(start,tree$mst$to)]
    x<-x[-c(start,tree$mst$to),,drop=FALSE]
    m<-m+nrow(tree$x)
    forest[[i]]<-tree
    nn<-nnfind(x)
    if (!any(nn$dist<threshold)) break;
    if (m/n > fill) break;
    if (sum(nn$dist<threshold)<give.up) break;
    if (i==maxclust) break;
    if (verbose)
       print(c(m, nrow(x), sum(nn$dist<threshold)))
    i<-i+1
    start<-which.min(nn$dist)
  })
  forest[[i+1]]<-list(x=x,rows=rows)
  class(forest)<-"nncluster"
  forest		     
}

trimCluster<-function(nnclust, size=10){
  n<-length(nnclust)
  sizes<-sapply(nnclust[-n],function(t) t$mst$n)
  small<-sizes<size
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

nnSampleCluster<-function(x, threshold, samplesize=5000, fill=0.95, maxclust=20, give.up=20, nsamples=10,verbose=FALSE){
  if (nrow(x) < 2*samplesize+1) stop("samplesize is too large")
  if (samplesize<1 || samplesize<maxclust || samplesize<give.up) stop("samplesize is too small")

  meandist<-Inf
  best<-list()
  for(i in 1:nsamples){
    keep <- sample(nrow(x), samplesize,replace=FALSE)
    cluster <- nncluster(x[keep,], threshold, fill=fill, maxclust=maxclust, give.up=give.up,verbose=FALSE)
    link <- nnfind(x[keep,],x[-keep,])
    distances<-c(link$dist, do.call(c,lapply(cluster, function(tree) tree$dist)))
    if (mean(distances)<meandist) {
      meandist<-mean(distances)
      best<-list(keep,cluster,link)
    }
  }
  find<-nnfind(x[-keep,])
  boost<- sample(nrow(x)-samplesize, samplesize, prob=exp(link$dist-find$dist))
  cluster2 <- nncluster(x[-keep,][boost,], threshold, fill=fill, maxclust=maxclust, give.up=give.up,verbose=FALSE)
  if (nrow(x)>2*samplesize+3)
    link2<-nnfind(x[-keep,][boost,],x[-keep,][-boost,])
  else
    link2<-NULL
  clusters<-list(keep, cluster,link, boost,cluster2,link2)
  clusters
}
