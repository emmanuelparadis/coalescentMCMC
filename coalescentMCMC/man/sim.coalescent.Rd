\name{sim.coalescent}
\alias{sim.coalescent}
\title{Coalescent Simulation and Visualisation}
\description{
  This is pedagogic function to show what is the coalescent in a simple
  population model with discrete generations and asexual reproduction.
}
\usage{
sim.coalescent(n = 5, TIME = 50, growth.rate = NULL, N.0 = 50, N.final = 20,
               col.lin = "grey", col.coal = "blue", pch = NULL, ...)
}
\arguments{
  \item{n}{the sample size.}
  \item{TIME}{the number of generations.}
  \item{growth.rate}{the growth rate of the population.}
  \item{N.0}{the initial size of the population.}
  \item{N.final}{the final size of the population (i.e., at present).}
  \item{col.lin}{the colour used to show links of ancestry in the
    population.}
  \item{col.coal}{the colour used to show the coalescent of the \eqn{n}
    individuals.}
  \item{pch}{the symbol used to show individuals (none by default).}
  \item{\dots}{further arguments passed to \code{points} if \code{pch}
    is used.}
}
\details{
  The simulation works along the following steps. The number of
  individuals at each generation is calculated. For each individual, a
  (unique) parent is randomly chosen at the previous generation. All
  individuals are then plotted and the ancestry lines are shown; the
  individuals are eventually ordered to avoid line-crossings. A sample
  of \eqn{n} individuals are randomly chosen from the last generation,
  and their shared ancestry is shown with thicker lines.

  The first (oldest) generation is at the bottom, and the final
  (present) one is at the top of the plot.

  The population size at each generation is determined from the four
  arguments: \code{TIME}, \code{growth.rate}, \code{N.0}, and
  \code{N.final}. At least three of them must be given by the user. If
  \code{TIME} is not given, its value is calculated with
  \code{log(N.final/N.0) / growth.rate}.

  This code was used to make the figures in Emerson et al. (2001).
}
\references{
  Emerson, B., Paradis, E. and Thebaud C. (2001). Revealing the
  demographic histories of species using DNA sequences. \emph{Trends in
    Ecology and Evolution}, \bold{16}, 707--716.
}
\author{Emmanuel Paradis}
\examples{
sim.coalescent()
sim.coalescent(N.0 = 20) # constant population size
}
\keyword{datagen}
