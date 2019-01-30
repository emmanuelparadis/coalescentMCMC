## plotTHETA.R (2019-01-28)

##   Plot THETA From coalescentMCMC Output

## Copyright 2019 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

plotTHETA <- function(x, phy, add = FALSE, rightwards = TRUE, col = "blue",
                      transparency = 50/length(phy), xlab = "Time",
                      ylab = expression(Theta), ylim = NULL, x.scale = 1,
                      y.scale = 1, show.present = TRUE, ...)
{
    if (!inherits(x, "mcmc")) stop('x should be of class "mcmc"')
    if (!inherits(phy, "phylo")) phy <- c(phy)
    else if (!inherits(phy, "multiPhylo"))
        stop('phy should be of class "multiPhylo" or "phylo"')

    colrgb <- as.vector(col2rgb(col))/255
    co <- rgb(colrgb[1], colrgb[2], colrgb[3], transparency)
    Ntrees <- length(phy)

    model <- attr(x, "model")
    if (model == "bsplines") {
        degree <- attr(x, "degree")
        nknots <- attr(x, "nknots")
    }

    tmrca <- sapply(phy, function(x) .branching_times(x)[1])
    xl <- c(0, max(tmrca))
    para.nms <- colnames(x)[-1]

    if (model == "constant") {
        if (is.null(ylim)) ylim <- range(x[, 3])
    }
    if (model %in% c("linear", "step")) {
        if (is.null(ylim)) ylim <- range(x[, 3:4])
    }
    if (model == "time") {
        t <- seq(0, xl[2], length.out = 100)
        Y <- matrix(NA_real_, 100, nrow(x))
        for (i in 1:nrow(x)) Y[, i] <- x[i, 3]*exp(x[i, 4]*t)
        if (is.null(ylim)) ylim <- range(Y)
    }
    if (model == "bsplines") {
        t <- seq(0, xl[2], length.out = 100)
        if (nknots) {
            knots <- x[, 3:(nknots + 1)]
            beta <- x[, -(1:(nknots + 1))]
        } else {
            knots <- NULL
            beta <- x[, -(1:2)]
            bst <- bs(t, knots = NULL, degree = degree, intercept = TRUE)
        }
        Y <- matrix(NA_real_, 100, nrow(x))
        for (i in 1:nrow(x)) {
            if (nknots)
                bst <- bs(t, knots = knots[i, ], degree = degree, intercept = TRUE)
            Y[, i] <- drop(beta[i, ] %*% t(bst))
        }
        if (is.null(ylim)) ylim <- range(Y)
    }

    if (rightwards) xl <- rev(-xl)

    if (x.scale != 1) {
        xl <- x.scale * xl
        tmrca <- x.scale * tmrca
    }
    if (y.scale != 1) ylim <- y.scale * ylim

    if (!add) {
        plot(NA, type = "n", xlim = xl, ylim = ylim, xlab = xlab, ylab = ylab,
             xaxt = "n", ...)
        if (rightwards) {
            pxl <- pretty(xl)
            axis(1, at = pxl, labels = -pxl)
        } else axis(1)
    }
    if (model == "constant")
        for (i in 1:nrow(x)) lines(xl, y.scale * rep(x[i, 3], 2), col = co)

    if (model == "linear") {
        if (y.scale != 1) x[, 3:4] <- y.scale * x[, 3:4]
        if (rightwards) {
            for (i in 1:nrow(x))
                lines(c(-tmrca[i], 0), x[i, 4:3], col = co)
        } else {
            for (i in 1:nrow(x))
                lines(c(0, tmrca[i]), x[i, 3:4], col = co)
        }
    }

    if (model == "step") {
        if (x.scale != 1) x[, 5] <- x.scale * x[, 5]
        if (y.scale != 1) x[, 3:4] <- y.scale * x[, 3:4]
        if (rightwards) {
            for (i in 1:nrow(x))
                lines(c(-tmrca[i], -x[i, 5], -x[i, 5], 0),
                      c(x[i, c(4, 4, 3, 3)]), col = co)
        } else {
            for (i in 1:nrow(x))
                lines(c(0, x[i, 5], x[i, 5], tmrca[i]),
                      c(x[i, c(3, 3, 4, 4)]), col = co)
        }
    }

    if (model == "time" || model == "bsplines") {
        if (x.scale != 1) t <- x.scale * t
        if (y.scale != 1) Y <- y.scale * Y
        if (rightwards) t <- -t
        for (i in 1:nrow(x)) lines(t, Y[, i], col = co)
    }

    if (show.present && !add) mtext("Present", 1, 2, at = 0, font = 3)
}
