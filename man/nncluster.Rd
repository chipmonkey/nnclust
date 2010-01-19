\name{nncluster}
\Rdversion{1.1}
\alias{nncluster}
\alias{nearestCluster}
\alias{trimCluster}
\alias{clusterMember}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast clustering with restarted minimum spanning tree.
}
\description{
Uses Prim's algorithm to build a minimum spanning tree for each
cluster, stopping when the nearest-neighbour distance rises above a
specified threshold.  Returns a set of clusters and a set of 'outliers'
not in any cluster. \code{trimCluster} tidies up the output by removing
small clusters, \code{clusterMember} returns cluster membership for the
original data points.

}
\usage{
nncluster(x, threshold, fill = 0.95, maxclust = 20, give.up = 500,verbose=FALSE,start=NULL)
trimCluster(nnclust, size=10)
clusterMember(nnclust, outlier=TRUE)
nearestCluster(nnclust, threshold=Inf,outlier=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data matrix
}
  \item{threshold}{
Threshold for stopping the tree building within a cluster. The tree
stops when the squared euclidean distance to the closest point to the
tree is greater than this. If \code{threshold} is a vector, the elements 
will be used in succession, with the last element repeated as necessary.
}
  \item{fill}{
Stop when the clusters make up this fraction of the data.
}
  \item{maxclust}{
Stop at this many clusters
}
  \item{give.up}{
Stop when fewer than this many pairs have nearest-neighbour distance
less than \code{threshold}.
}
\item{verbose}{Print some cluster summaries before restarting?}
\item{nnclust}{
  An object of class \code{nncluster}, returned by \code{nncluster}
}
\item{size}{Clusters smaller than this are added to the 'outlier' set}
\item{outlier}{If \code{FALSE}, use \code{NA} for the cluster identifier
for outliers}
\item{start}{integer index to start the minimum spanning tree at this observation}
}
\details{
   Works best for well-separated clusters in up to 8
   dimensions, and sample sizes up to hundreds of thousands.
   
If you want a complete minimum spanning tree, run \code{mst} on the
outlier set and then use \code{nnfind} to find the shortest links
connecting the clusters. When there are well-separated clusters this
will be faster than running \code{mst} once on the whole data set.

\code{clusterMember} returns a vector of integers indicating cluster membership. Outliers are treated as a separate cluster if \code{outlier} is \code{TRUE}, otherwise they code as \code{NA}.  \code{nearestCluster} assigns outliers at distance less than \code{threshold} from a cluster to the cluster whose nearest member is closest.

\code{trimCluster} returns a new \code{nncluster} object with small clusters converted to outliers. There must be at least one cluster larger than \code{size}.
}
\value{
A list of class \code{nncluster}. Each element but the last
describes a cluster, with components \code{mst} containing the tree,
\code{x} containing the data, and \code{rows} containing row numbers in
the initial data set.

The last element describes the unclustered outliers and has no
\code{mst} component.

}

\author{
Thomas Lumley
}
\note{
  The performance of this algorithm depends critically on the
  performance of the nearest-neighbour finder, and can decay
  catastrophically if too many uninformative variables are added.

  The performance can also be poor if the data are close to being ordered on some of the variables.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mst}}, \code{\link{nnfind}}
}
\examples{
x<-scale(faithful)
a<-nncluster(x, threshold=0.1, give.up=0, fill=1)
a
id<-clusterMember(a)
plot(faithful, col=id, pch=19)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

