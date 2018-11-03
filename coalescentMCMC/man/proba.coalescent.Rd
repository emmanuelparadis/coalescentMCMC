\name{proba.coalescent}
\alias{proba.coalescent}
\title{Probability of Coalescence}
\description{
  This function calculates the probability that two lineages coalesce
  out of a sample of size \eqn{n} in a population of size \eqn{N} after
  \eqn{t} generations.
}
\usage{
proba.coalescent(t, N = 1e4, n = 2, exact = TRUE)
}
\arguments{
  \item{t}{a vector of generations (rounded to integers if needed).}
  \item{N}{the size of population (10,000 by default).}
  \item{n}{the sample size (2 by default).}
  \item{exact}{a logical value specifying whether exact calculation
    should be done or an approximation (Hudson, 1991, eq. 3).}
}
\references{
  Hudson, R. R. (1991). Gene genealogies and the coalescent
  process. \emph{Oxford Surveys in Evolutionary Biology}, \bold{7},
  1--44.
}
\author{Emmanuel Paradis}
\examples{
proba.coalescent(1:10)
## the approximate formula doesn't work well when n is
## not small compared to N:
proba.coalescent(1, 100, 50, exact = FALSE)
proba.coalescent(1, 100, 50)
}
\keyword{datagen}
