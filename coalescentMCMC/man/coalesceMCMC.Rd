\name{coalescentMCMC}
\alias{coalescentMCMC}
\alias{getMCMCtrees}
\alias{saveMCMCtrees}
\alias{cleanMCMCtrees}
\alias{getLastTree}
\alias{getMCMCstats}
\title{Run MCMC for Coalescent Trees}
\description{
  This function runs a Markov chain Monte Carlo (MCMC) algorithm to
  generate a set of trees which is returned with their likelihoods.

  \code{getMCMCtrees} extracts the trees from previous MCMC runs.

  \code{saveMCMCtrees} saves the lists of trees from previous runs on
  the user's hard disk.

  \code{cleanMCMCtrees} deletes the lists of trees from previous runs
  (the files possibly on the hard disk are not changed).

  \code{getLastTree} extracts the last tree from a list of trees (object
  of class \code{"multiPhylo"}).

  \code{getMCMCstats} returns the summary data for the different chains
  run during a session.
}
\usage{
coalescentMCMC(x, ntrees = 3000, burnin = 1000, frequency = 1,
               tree0 = NULL, model = NULL, printevery = 100)
getMCMCtrees(chain = NULL)
saveMCMCtrees(destdir = ".", format = "RDS", ...)
cleanMCMCtrees()
getLastTree(X)
getMCMCstats()
}
\arguments{
  \item{x}{a set of DNA sequences, typically an object of class
    \code{"DNAbin"} or \code{"phyDat"}.}
  \item{ntrees}{the number of trees to output.}
  \item{burnin}{the number of trees to discard as ``burn-in''.}
  \item{frequency}{the frequency at which trees are sampled.}
  \item{tree0}{the initial tree of the chain; by default, a UPGMA
    tree with a JC69 distance is generated.}
  \item{model}{the coalescent model to be used for resampling. By
    default, a constant-THETA is used.}
  \item{printevery}{an integer specifying the frequency at which to
    print the numbers of trees proposed and accepted; set to 0 to cancel
    all printings.}
  \item{chain}{an integer giving which lists of trees to extract}
  \item{destdir}{a character string giving the location where to save
    the files; by default, this is the current working directory.}
  \item{format}{the format of the tree files. Three choices are
    possible (cae-insensitive): \code{"RDS"}, \code{"Newick"},
    \code{"NEXUS"}, or any unambiguous abbreviation of these.}
  \item{\dots}{options passed to the function used to write the tree
    files (see below).}
  \item{X}{an bject of class \code{"multiPhylo"}.}
}
\details{
  A simple MCMC algorithm is programmed using at each step the
  ``neighborhood rearrangement'' operation (Kuhner et al., 1995) and
  Hastings's ratio for acceptance/rejection of the proposed tree.

  The number of generations of the chain is determined by: `frequency'
  times `ntrees' plus `burnin'. Only the `ntrees' trees are output
  whereas all the log-likelihood values are output.

  The list of trees is returned in a specific environment and can be
  extracted with \code{getMCMCtrees}.

  \code{saveMCMCtrees} saves the files with, by default, the RDS format
  using \code{\link[base]{saveRDS}}. If \code{format = "Newick"},
  \code{\link[ape]{write.tree}} is used.; if \code{format = "NEXUS"},
  \code{\link[ape]{write.nexus}} is used. Options can be passed to any
  of these functions with \code{\dots}.

  \code{getLastTree(X)} is a short-cut to \code{X[[length(X)]]}.
}
\value{
  \code{coalescentMCMC} returns an object of class \code{"coda"} with
  the log-likelihood and the parameters of each tree.

  \code{getLastTree} returns an object of class \code{"phylo"}.

  \code{getMCMCstats} returns a data frame.
}
\references{
  Hastings, W. K. (1970) Monte Carlo sampling methods using Markov
  chains and their applications. \emph{Biometrika}, \bold{57}, 97--109.

  Kuhner, M. K., Yamato, J. and Felsenstein, J. (1995) Estimating
  effective population size and mutation rate from sequence data using
  Metropolis-Hastings sampling. \emph{Genetics}, \bold{140}, 1421--1430.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{dcoal}}, \code{\link{treeOperators}}
}
\examples{
\dontrun{
data(woodmouse)
out <- coalescentMCMC(woodmouse)
plot(out)
getMCMCtrees() # returns 3000 trees
}}
\keyword{models}
