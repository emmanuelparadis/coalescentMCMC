## dcoal.R (2019-01-18)

##   pdf of Various Time-Dependent Coalescent Models

## Copyright 2012-2019 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

dcoal.step <- function(bt, theta0, theta1, tau, log = FALSE)
{
    ## if theta0 = theta1 use dcoal():
    if (theta0 == theta1) return(dcoal(bt, theta0, log))
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(bt))
    n <- length(x)
    f <- F <- numeric(n)
    s <- x <= tau
    f[s] <- theta0 # theta(t) for x <= tau
    f[!s] <- theta1 # theta(t) for x > tau
    ## the primitive of '1/f' is 'x/theta0' for x <= tau...
    F[s] <- x[s]/theta0
    ## ... see vignette for x > tau:
    F[!s] <- tau/theta0 + (x[!s] - tau)/theta1
    Ncomb <- choose(n:2, 2)
    p <- sum(log(Ncomb) - log(f[-1]) - Ncomb * (F[-1] - F[-n]))
    if (!log) p <- exp(p)
    p
}

dcoal.linear <- function(bt, theta0, thetaT, log = FALSE)
{
    if (theta0 == thetaT) return(dcoal(bt, theta0, log))
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(bt))
    TMRCA <- x[length(x)]
    kappa <- (thetaT - theta0)/TMRCA
    f <- theta0 + x*kappa # = theta(t)
    ## don't need to compute the inverse of theta(t)
    n <- length(f)
    K <- n:2
    Ncomb <- K*(K - 1)/2 # Ncomb <- choose(n:2, 2)
    ## the primitive of 1/f is log(f)/f'
    lnf <- log(f)
    p <- sum(log(Ncomb) - lnf[-1] - Ncomb * (lnf[-1] - lnf[-n])/kappa)
    if (!log) p <- exp(p)
    p
}

dcoal.time2 <- function(bt, theta0, rho1, rho2, tau, log = FALSE)
{
    ## if rho1 = rho2 = 0 use dcoal():
    if (rho1 == rho2) {
        if (!rho1) return(dcoal(bt, theta0, log))
        else return(dcoal.time(bt, theta0, rho1, log))
    }
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(bt))
    ## \theta(t) = \theta_0 e^{\rho_1 t} \quad t \le \tau
    ## \theta(t) = \theta_0 e^{\rho_1 \tau} e^{\rho_2 t} \quad t \gt \tau
    ##           = \theta_0 e^{\rho_1 \tau + \rho_2 t}
    n <- length(x)
    f <- F <- numeric(n)
    s <- x <= tau
    f[s] <- exp(-rho1*x[s])/theta0 # inverse of theta(t) for x <= tau
    f[!s] <- exp(-rho2*x[!s] - (rho1 - rho2)*tau)/theta0 # inverse of theta(t) for x > tau
    ## the primitive of 'f' is '-(f - 1/theta0)/rho1' for x <= tau...
    F[s] <- -(f[s] - 1/theta0)/rho1
    ## ... see paper for x > tau:
    A <- exp(-rho1*tau)/theta0
    F[!s] <- -(A - 1/theta0)/rho1 - (f[!s] - A)/rho2
    Ncomb <- choose(n:2, 2)
    p <- sum(log(Ncomb) + log(f[-1]) - Ncomb * (F[-1] - F[-n]))
    if (!log) p <- exp(p)
    p
}

dcoal.time <- function(bt, theta0, rho, log = FALSE)
{
    ## would return NaN if rho = 0; use dcoal() instead:
    if (!rho) return(dcoal(bt, theta0, log))
    ## branching times from the most recent to the oldest:
    x <- c(0, sort(bt))
    ## \theta(t) = \theta_0 e^{\rho t}
    f <- exp(-rho * x)/theta0 # inverse of theta(t)
    n <- length(f)
    k <- n:2
    Ncomb <- k * (k - 1)/2 # choose(n:2, 2)
    ## the primitive of 'f' is '-f/rho'
    p <- sum(log(Ncomb) + log(f[-1]) + Ncomb * (f[-1] - f[-n])/rho)
    ## gr.theta0 <- -(n - 1)/theta0 - sum(Ncomb * (f[-1] - f[-n])/(theta0 * rho))
    ## gr.rho <- sum(-x[-1] + Ncomb * (-(f[-1] - f[-n])/rho^2 + (x[-1]*f[-1] - x[-n]*f[-n])/rho))
    ## if (log) attr(p, "gradient") <- c(gr.theta0, gr.rho)
    ## else p <- exp(p)
    if (!log) p <- exp(p)
    p
}

dcoal <- function(bt, theta, log = FALSE)
{
    x <- c(0, sort(bt)) # coalescent times sorted
    n <- length(x)
    k <- n:2
    Ncomb <- k * (k - 1)/2
    p <- sum(log(Ncomb)) - (n - 1) * log(theta) - sum((x[-1] - x[-n]) * Ncomb)/theta
    if (!log) p <- exp(p)
    p
}

dcoal.bsplines <- function(bt, beta, knots = NULL, degree = 3, log = FALSE)
{
    ## slightly faster than lchoose(n, 2) and equally vectorized:
    f <- function(n) log(n) + log(n - 1) - 0.6931471805599452862268
    ## coalescent times from the most recent to the oldest:
    x <- c(0, sort(bt))
    names(x) <- NULL
    n <- length(x)
    bsx <- bs(x, knots = knots, degree = degree, intercept = TRUE)
    THETA <- drop(bsx %*% beta)
    if (any(is.na(THETA) | !is.finite(THETA) | THETA < 0))
        return(if (log) -1e100 else 0)
    invTHETA <- 1/THETA
    ## simple trapeze integrals:
    INT <- (x[-1] - x[-n]) * (invTHETA[-1] + invTHETA[-n])/2
    res <- sum(f(n:2) - log(THETA[-1]) - 0.5*(n:2)*((n - 1):1) * INT)
    if (!log) res <- exp(res)
    res
}
