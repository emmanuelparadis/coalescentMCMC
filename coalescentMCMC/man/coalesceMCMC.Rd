\name{coalescentMCMC}
\alias{coalescentMCMC}
\alias{logLik.coalescentMCMC}
\alias{AIC.coalescentMCMC}
\alias{BIC.coalescentMCMC}
\alias{anova.coalescentMCMC}
\title{Run MCMC for Coalescent Trees}
\description{
  These are the main function of the package to run a Markov chain Monte
  Carlo (MCMC) to generate a set of trees which is returned with their
  likelihoods, the coalescent likelihoods and the respective parameter(s).

  The \code{logLik} method returns the average log-likelihood of the
  coalescent model. \code{AIC}, \code{BIC}, and \code{anova} use this
  average log-likelihood.
}
\usage{
coalescentMCMC(x, ntrees = 3000, model = "constant", tree0 = NULL,
               printevery = 100, degree = 1, nknots = 0,
               knot.times = NULL, moves = 1:6)
\method{logLik}{coalescentMCMC}(object, ...)
\method{AIC}{coalescentMCMC}(object, ..., k = 2)
\method{BIC}{coalescentMCMC}(object, ...)
\method{anova}{coalescentMCMC}(object, ...)
}
\arguments{
  \item{x}{a set of DNA sequences, typically an object of class
    \code{"DNAbin"} or \code{"phyDat"}.}
  \item{ntrees}{the number of trees to output.}
  \item{tree0}{the initial tree of the chain; by default, a UPGMA
    tree with a JC69 distance is generated.}
  \item{model}{the coalescent model to be used for resampling. By
    default, a constant-THETA is used.}
  \item{printevery}{an integer specifying the frequency at which to
    print the numbers of trees proposed and accepted; set to 0 to cancel
    all printings.}
  \item{degree, nknots, knot.times}{parameters used if \code{model =
    "splines"}.}
  \item{moves}{the tree moves used by the MCMC (see details).}
  \item{\dots}{options passed to other methods.}
  \item{object}{an bject of class \code{"coalescentMCMC"}.}
  \item{k}{the coefficient used to calculate the AIC (see
    \code{\link[stats]{AIC}}).}
}
\details{
  Six tree moves are programmed and one is chosen randomly at each step
  of the MCMC. The steps are: (1) NeighborhoodRearrangement (Kuhner et
  al., 1995), (2) ScalingMove, (3) branchSwapping, (4) subtreeExchange,
  (5) NodeAgeMove, and (6) randomWalkThetaMu (all five from Drummond et
  al., 2002). In practice, it appears that in many situations
  \code{moves = c(1, 3)} is a good selection resulting in around 50\% acceptance rate.
}
\value{
  \code{coalescentMCMC} returns an object of class
  \code{c("coalescentMCMC", "coda")} with the log-likelihood and the
  parameters of each tree.

  \code{logLik}, \code{AIC} and \code{BIC} return a numeric vector.

  \code{anova} return an object of class \code{"anova"}.
}
\references{
  Drummond, A. J., Nicholls, G. K., Rodrigo, A. G. and Solomon,
  W. (2002) Estimating mutation parameters, population history and
  genealogy simultaneously from temporally spaced sequence
  data. \emph{Genetics}, \bold{161}, 1307--1320.

  Hastings, W. K. (1970) Monte Carlo sampling methods using Markov
  chains and their applications. \emph{Biometrika}, \bold{57}, 97--109.

  Kuhner, M. K., Yamato, J. and Felsenstein, J. (1995) Estimating
  effective population size and mutation rate from sequence data using
  Metropolis-Hastings sampling. \emph{Genetics}, \bold{140}, 1421--1430.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{getMCMCtrees}}, \code{\link{dcoal}}, \code{\link{treeOperators}}
}
\examples{
\dontrun{
data(woodmouse)
out <- coalescentMCMC(woodmouse)
plot(out)
getMCMCtrees() # returns 3000 trees
}}
\keyword{models}
