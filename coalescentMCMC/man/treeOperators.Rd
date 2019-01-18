\name{treeOperators}
\alias{treeOperators}
\alias{NeighborhoodRearrangement}
\alias{TipInterchange}
\title{Trees Operators for Running MCMC}
\description{
  These functions provide tools for tree rearrangement to be used as
  operators in a MCMC run.
}
\usage{
NeighborhoodRearrangement(phy, n, THETA, brtimes)
TipInterchange(phy, n)
}
\arguments{
  \item{phy}{a tree of class \code{"phylo"}.}
  \item{n}{the number of tips in \code{phy}.}
  \item{THETA}{The estimate of \eqn{\Theta}{THETA} for \code{phy} at the
    node 'target'.}
  \item{brtimes}{the branching times of \code{phy}.}
}
\details{
  \code{NeighborhoodRearrangement} performs a rearrangement as described
  by Kuhner et al. (1995).

  \code{TipInterchange} interchanges two tips under the condition that
  they are not sisters.

  \code{EdgeLengthJittering} alters the branch lengths by adding a
  random value from a uniform distribution defined by
  \code{range(phy$edge.length)} (the ultrametric nature of the tree is
  conserved).
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Kuhner, M. K., Yamato, J. and Felsenstein, J. (1995) Estimating
  effective population size and mutation rate from sequence data using
  Metropolis-Hastings sampling. \emph{Genetics}, \bold{140}, 1421--1430.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{coalescentMCMC}}, \code{\link{dcoal}}
}
\examples{
tr <- rcoal(10)
ts <- NeighborhoodRearrangement(tr, 10, 1, branching.times(tr))
layout(matrix(1:2, 2))
plot(tr); plot(ts)
layout(1)
}
\keyword{manip}
