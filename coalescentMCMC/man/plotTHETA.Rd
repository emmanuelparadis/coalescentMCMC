\name{plotTHETA}
\alias{plotTHETA}
\title{Plot THETA From coalescentMCMC Output}
\description{
  This function plots the values of \eqn{\Theta}{THETA} predicted from
  model fitted with \code{\link{coalescentMCMC}}.
}
\usage{
plotTHETA(x, phy, add = FALSE, rightwards = TRUE, col = "blue",
          transparency = 50/length(phy), xlab = "Time",
          ylab = expression(Theta), ylim = NULL, x.scale = 1,
          y.scale = 1, show.present = TRUE, ...)
}
\arguments{
  \item{x}{an output from.}
  \item{phy}{an object of class \code{"phylo"} or  \code{"multiPhylo"}.}
  \item{add}{a logical value: whether to add the curves on an existing
    plot?}
  \item{rightwards}{a logical value: whether to draw the time
    (horizontal) axis rightwards (the default).}
  \item{col}{the base colour of the curves.}
  \item{transparency}{the degree of transparency of the curves.}
  \item{xlab, ylab}{character strings giving the labels for the x- and
    y-axes.}
  \item{ylim}{a numeric vector with two values giving the limits of the
    y-axis (useful if several curves are drawn).}
  \item{x.scale}{the scaling factor for the x-axis; typically, the
    inverse of the mutation rate.}
  \item{y.scale}{the scaling factor for the y-axis; typically, the
    inverse of twice the mutation rate.}
  \item{show.present}{a logical value: whether to indicate the present
    in italics (near the zero coordinate on the x-axis).}
  \item{\dots}{further arguments passed to \code{plot}.}
}
\value{NULL}
\author{Emmanuel Paradis}
\seealso{\code{\link{coalescentMCMC}}}
\keyword{hplot}
