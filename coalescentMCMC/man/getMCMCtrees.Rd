\name{getMCMCtrees}
\alias{getMCMCtrees}
\alias{saveMCMCtrees}
\alias{cleanMCMCtrees}
\alias{getLastTree}
\alias{getMCMCstats}
\title{Managing List of Trees From MCMC}
\description{
  These functions help to manage trees output from MCMCs.

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
getMCMCtrees(chain = NULL)
saveMCMCtrees(destdir = ".", format = "RDS", ...)
cleanMCMCtrees()
getLastTree(X)
getMCMCstats()
}
\arguments{
  \item{chain}{an integer giving which lists of trees to extract}
  \item{destdir}{a character string giving the location where to save
    the files; by default, this is the current working directory.}
  \item{format}{the format of the tree files. Three choices are
    possible (cae-insensitive): \code{"RDS"}, \code{"Newick"},
    \code{"NEXUS"}, or any unambiguous abbreviation of these.}
  \item{\dots}{options passed to the function used to write the tree
    files (see below) or passed to other methods.}
  \item{X}{an bject of class \code{"multiPhylo"}.}
}
\details{
  The list of trees is returned in a specific environment and can be
  extracted with \code{getMCMCtrees}.

  \code{saveMCMCtrees} saves the files with, by default, the RDS format
  using \code{\link[base]{saveRDS}}. If \code{format = "Newick"},
  \code{\link[ape]{write.tree}} is used.; if \code{format = "NEXUS"},
  \code{\link[ape]{write.nexus}} is used. Options can be passed to any
  of these functions with \code{\dots}.

  \code{getLastTree(X)} is a short-cut to \code{X[[length(X)]]}.

  Most functions from the package \pkg{coda} can also be used to
  analyse the MCMC outputs.
}
\value{
  \code{getLastTree} returns an object of class \code{"phylo"}.

  \code{getMCMCstats} returns a data frame.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{coalescentMCMC}}, \code{\link{treeOperators}},
  \code{\link{subset.coalescentMCMC}}
}
\keyword{manip}
