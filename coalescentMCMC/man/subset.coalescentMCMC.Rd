\name{subset.coalescentMCMC}
\alias{subset.coalescentMCMC}
\title{Subset MCMC Output}
\description{
  This function helps to manipulate outputs from
  \code{\link{coalescentMCMC}}, particularly it sets the attributes
  correctly (unlike if you would use \code{[}).
}
\usage{
\method{subset}{coalescentMCMC}(x, burnin = 1000, thinning = 10, end = NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{"coalescentMCMC"}.}
  \item{burnin}{an integer: the number of generations to drop.}
  \item{thinning}{an integer: the sampling frequency.}
  \item{end}{an integer: the generations after \code{end} are removed.}
  \item{\dots}{(unused)}
}
\details{
  This function has more or less the same effect than the standard
  indexing operator (\code{[}). The main difference is that the
  attributes attached to the output from \code{\link{coalescentMCMC}}
  are not lost.

  The options \code{end} may help to focus on one part of the MCMC (see
  examples).
}
\value{an object of class \code{"coalescentMCMC"}.}
\author{Emmanuel Paradis}
\note{
  The default values of \code{burnin} and \code{thinning} are only
  indicative: it is recommended to use functions in the package
  \pkg{coda} to help find appropriate values (see examples).
}
\seealso{
  \code{\link[coda]{acfplot}}, \code{\link[coda]{effectiveSize}}
}
\examples{
\dontrun{
data(woodmouse)
res <- coalescentMCMC(woodmouse, 1e6, moves = c(1, 3)) # ~ 1 hr
plot(res) # surely hard to read
plot(subset(res, end = 1e3)) # plot only the first 1000 generations

acfplot(res)
acfplot(subset(res, 1e4, 100))
}
}
\keyword{manip}
