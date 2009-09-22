\name{mst}
\Rdversion{1.1}
\alias{mst}
\alias{mst_restart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum spanning trees
}
\description{
Minimum spanning tree by Prim's algorithm, using a fast
nearest-neighbour search to find candidate points. The variant
\code{mst_restart} stops when the link length is greater than a
specified threshold.
}
\usage{
mst(X, rebuild = sqrt(nrow(X))/4)
mst_restart(X,rebuild=sqrt(nrow(X))/4,threshold=Inf,start=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Matrix of data points
}
  \item{rebuild}{
How often to rebuild the k-d tree for nearest neighbour lookup.
}
\item{threshold}{Return when the next link is longer than this}
\item{start}{Start here rather than at the first row of \code{X}}
}

\value{
\item{from,to}{Row numbers of \code{X} for each link}
\item{dist}{Squared length of link}
\item{n}{For \code{mst_restart} only, the number of links added}
}
\seealso{
\code{\link{nncluster}} for a wrapper for \code{mst_restart} to do clustering.
}
\author{
Thomas Lumley
}
\references{
Bentley JL, Friedman JH. "Fast algorithms for constructing minimal spanning trees in coordinate spaces" IEEE Transactions on Computers C27(2) Feb 1978

}

\examples{
x<-scale(faithful)
a<-mst(x)
plot(faithful)
segments(faithful[a$from,1], faithful[a$from,2], faithful[a$to,1],
faithful[a$to,2],col="blue")


## restarting 
plot(a$dist)
a<-mst_restart(x,threshold=0.04,start=0)
plot(faithful)
segments(faithful[a$from,1], faithful[a$from,2], faithful[a$to,1],
faithful[a$to,2],col="red")

b<-mst_restart(x[-a$to,],threshold=0.04)
ff<-faithful[-a$to,]
segments(ff[b$from,1], ff[b$from,2], ff[b$to,1],ff[b$to,2],col="purple")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
