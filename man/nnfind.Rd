\name{nnfind}
\Rdversion{1.1}
\alias{nnfind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nearest-neighbour distances}
\description{
Find the nearest neighbours of points in one data set from another data set. Useful for Mallows-type distance metrics.
}
\usage{
nnfind(from, to)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{
  A matrix giving the first sample, with rows specifying points and columns specifying dimensions.
}
  \item{to}{
  Optional matrix with the same number of columns as \code{from} giving the second sample.
}
}
\details{
If \code{to} is specified, for each point in \code{to} find the nearest neighbour in \code{from}.  If \code{to} is not specified, for each point in \code{from} find the nearest distinct neighbour in \code{from}.

The algorithm builds a k-d tree on \code{from} and drops points down the tree to find the nearest neighbour.  This is much more efficient than a brute-force search as long as the dimension is low enough.  For a million points in each data set and five dimensions, the code does only about 0.03\% of the possible distance computations.

}
\value{
 \item{count}{Number of distance computations performed }
 \item{neighbour}{Vector of row numbers in \code{from} containing the nearest neighbour of each point in \code{to}}
 \item{dist}{Vector of distances to the nearest neighbour for each point in \code{to}}
}

\examples{
data(faithful)
nn<-nnfind(as.matrix(faithful))
plot(faithful)
segments(faithful[,1], faithful[,2], faithful[nn$neighbour,1], faithful[nn$neighbour,2], col="blue")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

