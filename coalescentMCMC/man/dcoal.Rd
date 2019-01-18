\name{dcoal}
\encoding{utf8}
\alias{dcoal}
\alias{dcoal.step}
\alias{dcoal.linear}
\alias{dcoal.time}
\alias{dcoal.time2}
\title{Density Functions of Some Time-Dependent Coalescent Models}
\description{
  These functions compute the (log-)likelihood values for various
  coalescent models, including the constant-\eqn{\Theta}{THETA} model
  and various time-dependent models.
}
\usage{
dcoal(bt, theta, log = FALSE)
dcoal.step(bt, theta0, theta1, tau, log = FALSE)
dcoal.linear(bt, theta0, thetaT, log = FALSE)
dcoal.time(bt, theta0, rho, log = FALSE)
dcoal.time2(bt, theta0, rho1, rho2, tau, log = FALSE)
}
\arguments{
  \item{bt}{a vector of coalescent times (typically from
    \code{\link[ape]{branching.times}}).}
  \item{theta}{population parameter THETA.}
  \item{log}{a logical value specifying whether the probabilities should
    be returned log-transformed.}
  \item{theta0, theta1, thetaT}{THETA parameter for the time-dependent
    models.}
  \item{tau}{breakpoint in time when the parameters change.}
  \item{rho, rho1, rho2}{population (exponential) growth rates.}
}
\details{
  The models are detailed in a vignette: \code{vignette("CoalescentModels")}.
}
\value{
  a numeric vector with (log-)likelihood values.
}
\references{
  Griffiths, R. C. and Tavaré, S. (1994) Sampling theory for neutral
  alleles in a varying environment. \emph{Philosophical Transactions of
    the Royal Society of London. Series B. Biological Sciences},
  \bold{344}, 403--410.

  Kuhner, M. K., Yamato, J. and Felsenstein, J. (1998) Maximum
  likelihood estimation of population growth rates based on the
  coalescent. \emph{Genetics}, \bold{149}, 429--434.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{coalescentMCMC}}, \code{\link[ape]{branching.times}}
}
\keyword{models}
