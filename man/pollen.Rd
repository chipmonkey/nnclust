\name{pollen}
\Rdversion{1.1}
\alias{pollen}
\docType{data}
\title{
Artificial "pollen-grain" data from 1986 Data Expo.
}
\description{
Artificial data described as measurements of pollen grains. Created for the Data Expo at the 1986 Joint Statistical Meetings by David Coleman of RCA Labs.
}
\usage{data(pollen)}
\format{
  A data frame with 3848 observations on the following 6 variables.
  \describe{
    \item{\code{ridge}}{a numeric vector}
    \item{\code{nub}}{a numeric vector}
    \item{\code{crack}}{a numeric vector}
    \item{\code{weight}}{a numeric vector}
    \item{\code{density}}{a numeric vector}
    \item{\code{id}}{a numeric vector}
  }
}

\source{
\url{http://lib.stat.cmu.edu/data-expo/}
}

\examples{
data(pollen)
scaled <- scale(as.matrix(pollen[,1:5]))

## clearly at least two populations, based on neighbour distances
d <- nnfind(scaled)
plot(density(d$dist))
abline(v=0.08)

## cluster to extract two populations
nn <- nncluster(scaled, threshold=c(0.08,1), give.up=0)
nn

## tight cluster is the hidden message
plot(ridge~nub, data=pollen, subset=clusterMember(nn)==1)

## data set is a hollow ellipsoid with a treat in the middle.
coplot(density~ridge|crack*nub, data=pollen,pch=".",
    col=clusterMember(nn),cex=2, n=5)

}
\keyword{datasets}
