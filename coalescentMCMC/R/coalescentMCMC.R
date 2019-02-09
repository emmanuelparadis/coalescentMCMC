## coalescentMCMC.R (2019-02-01)

##   Run MCMC for Coalescent Trees

## Copyright 2012-2019 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

coalescentMCMC <- function(x, ntrees = 3000, model = "constant",
                           tree0 = NULL, printevery = 100, degree = 1,
                           nknots = 0, knot.times = NULL, moves = 1:6)
{
    MODELS <- c("constant", "time", "step", "linear", "bsplines")
    model <- match.arg(model, MODELS)
    if (is.na(model)) stop(paste("model must one of:", MODELS))
    if (packageVersion("phangorn") >= "1.99.5") {
        edQt <- phangorn::edQt
        lli <- phangorn::lli
        pml.fit <- phangorn::pml.fit
        pml.free <- phangorn::pml.free
        pml.init <- phangorn::pml.init
    }

    on.exit({
        if (packageVersion("phangorn") >= "1.99.5") pml.free()
        if (i < ntrees) {
            one2i <- seq_len(i)
            TREES <- TREES[one2i]
            LL <- LL[one2i]
            params <- params[one2i, ]
            LLmod <- LLmod[one2i]
            warning(paste("MCMC interrupted after", i, "generations"))
        }
        ## compress the list of trees:
        attr(TREES, "TipLabel") <- TREES[[1L]]$tip.label
        lapply(TREES, function(x) x$tip.label <- NULL)
        class(TREES) <- "multiPhylo"

        suffix <- 1
        list.trees <- .get.list.trees()
        if (l <- length(list.trees))
            suffix <- 1 + as.numeric(sub("TREES_", "", list.trees[l]))
        suffix <- sprintf("%03d", suffix)
        assign(paste("TREES", suffix, sep = "_"), TREES,
               envir = .coalescentMCMCenv)

        i <- i - 1L

        MCMCstats <- get("MCMCstats", envir = .coalescentMCMCenv)
        MCMCstats[[suffix]] <- c(i, j)
        assign("MCMCstats", MCMCstats, envir = .coalescentMCMCenv)

        ## LLmod stores the half deviances, hence the minus sign:
        LL <- cbind(LL, -LLmod, params)
        colnames(LL) <- c("logLik.tree", "logLik.coal", para.nms)
        LL <- mcmc(LL, start = 1, end = i)
        class(LL) <- c("coalescentMCMC", "mcmc")
        attr(LL, "model") <- model
        attr(LL, "call") <- match.call()
        attr(LL, "nobs") <- n
        if (model == "bsplines") {
            attr(LL, "degree") <- degree
            attr(LL, "nknots") <- nknots
        }
        return(LL)
    })

    verbose <- as.logical(printevery)

    if (is.null(tree0)) {
        d <- dist.dna(x, "JC69")
        tree0 <- upgma(d, "average")
    }

    X <- if (inherits(x, "phyDat")) x else phyDat(x)
    n <- length(tree0$tip.label)
    nodeMax <- 2*n - 1

    ## log-lik of the tree given the genetic data
    getlogLik <- function(phy, X) pml(phy, X)$logLik

    if (packageVersion("phangorn") >= "1.99.5") {
        INV <- Matrix(lli(X, tree0), sparse = TRUE)
        ll.0 <- numeric(attr(X, "nr"))
        ## by Klaus (ensures that tip labels of tree and data have same order):
        X <- subset(X, tree0$tip.label)
        ##
        bf <- rep(0.25, 4)
        eig <- edQt()
        pml.init(X)
        getlogLik <- function(phy, X) {
            phy <- reorder(phy, "postorder")
            pml.fit(phy, X, bf = bf, eig = eig, INV = INV, ll.0 = ll.0)
        }
    }

    TREES <- vector("list", ntrees)
    LL <- LLmod <- numeric(ntrees)
    TREES[[1L]] <- tree0
    lnL0 <- getlogLik(tree0, X)
    LL[1L] <- lnL0

    ## For each model, calculate:
    ## np: number of parameters of the coalescent model
    ## para.nms: names of these parameters (output in the "mcmc" object)
    ## lo: lower bounds of these parameters (except for "constant" model)
    ##     for nlminb()
    ## up: id. for the upper bounds
    ## getparams: a function that returns the coalescent parameters and
    ##            the log-likelihood given the coalescent times (bt)

    switch(model, constant = {
        np <- 1L
        para.nms <- "theta"
        ## quantities to calculate THETA:
        n2two <- n:2
        Ncomb <- n2two * (n2two - 1)/2 # choose(two2n, 2)
        getparams <- function(bt) {
            x <- diff(c(0, sort(bt)))
            par <- sum(x * Ncomb)/(n - 1)
            obj <- sum(log(Ncomb)) - (n - 1) * log(par) - sum(x * Ncomb)/par
            list(par = par, objective = -obj)
        }
    }, time = {
        np <- 2L
        para.nms <- c("theta0", "rho")
        getparams <- function(bt) {
            halfdev <- function(p) {
                if (any(is.nan(p))) return(1e100)
                -dcoal.time(bt, p[1], p[2], log = TRUE)
            }
            lo <- c(1e-8, -1e6)
            up <- c(100, 1000)
            out <- nlminb(c(0.02, 0.1), halfdev, lower = lo, upper = up)
            out[c("par", "objective")]
        }
    }, step = {
        np <- 3L
        para.nms <- c("theta0", "theta1", "tau")
        getparams <- function(bt) {
            halfdev <- function(p) {
                if (any(p <= 0) || any(is.nan(p))) return(1e100)
                -dcoal.step(bt, p[1], p[2], p[3], log = TRUE)
            }
            out <- nlminb(c(0.02, 0.02, bt[1]/2), halfdev)
            out[c("par", "objective")]
        }
    }, linear = {
        np <- 2L
        para.nms <- c("theta0", "thetaT")
        getparams <- function(bt) {
            halfdev <- function(p) {
                if (any(p <= 0) || any(is.nan(p))) return(1e100)
                -dcoal.linear(bt, p[1], p[2], log = TRUE)
            }
            out <- nlminb(c(0.02, 0.02), halfdev)
            out[c("par", "objective")]
        }
    }, bsplines = {
        free.knot.times <- TRUE
        if (!is.null(knot.times)) {
            nknots <- length(knot.times)
            free.knot.times <- FALSE
        }
        np <- nknots + degree + 1L
        if (free.knot.times) np <- np + nknots
        para.nms <- paste0("b", 1:np)
        getparams <- function(bt) {
            up <- rep(1e3, np)
            lo <- -up
            if (nknots) {
                up[1:nknots] <- bt
                lo[1:nknots] <- 0
                ip <- c(1:nknots * bt/(nknots + 1), runif(np))
            } else {
                ip <- runif(np)
            }
            halfdev <- function(p) {
                if (any(p <= 0) || any(is.nan(p))) return(1e100)
                if (nknots) {
                    knots <- p[1:nknots]
                    beta <- p[-(1:nknots)]
                } else {
                    knots <- NULL
                    beta <- p
                }
                -dcoal.bsplines(bt, beta, knots = knots, degree = degree, log = TRUE)
            }
            out <- nlminb(ip, halfdev, lower = lo, upper = up)
            out[c("par", "objective")]
        }
    })

    params <- matrix(0, ntrees, np)

    i <- 2L # number of sampled trees = number of generations
    j <- 0L # number of accepted moves

    if (verbose) {
        cat("Running the Markov chain:\n")
        cat("  Number of generations to run:", ntrees, "\n")
        cat("Generation    Nb of accepted moves\n")
    }

    bt0 <- .branching_times(tree0)
    fitcoal <- getparams(bt0)
    params[1L, ] <- para0 <- fitcoal$par
    LLmod[1L] <- lnLcoal0 <- fitcoal$objective

    getTHETAct <- function(n, bt) {
        x <- diff(c(0, sort(bt)))
        n2two <- n:2
        Ncomb <- n2two * (n2two - 1)/2
        sum(x * Ncomb)/(n - 1)
    }

    while (i <= ntrees) {
        if (verbose) if (! i %% printevery)
            cat("\r  ", i, "                ", j, "           ")

        move <- sample(moves, 1L)
        if (move == 1)
            THETA <- ifelse(model == "constant", para0, getTHETAct(n, bt0))
        tr.b <- switch(move,
                       NeighborhoodRearrangement(tree0, n, THETA, bt0),
                       TipInterchange(tree0, n),
                       ScalingMove(tree0),
                       branchSwapping(tree0, n, bt0),
                       subtreeExchange(tree0, n, bt0),
                       NodeAgeMove(tree0, n, bt0))

        lnL.b <- getlogLik(tr.b, X)
        ## calculate coalescent parameters for the proposed tree:
        bt <- .branching_times(tr.b)
        fitcoal <- getparams(bt)
        para <- fitcoal$par
        lnLcoal.b <- fitcoal$objective
        if (is.na(lnL.b) || is.na(lnLcoal.b)) {
            ACCEPT <- FALSE
        } else {
            ## R <- lnL.b + lnLcoal.b - lnL0 - lnLcoal0
            R <- lnL.b - lnL0
            ACCEPT <- if (R >= 0) TRUE else rbinom(1, 1, exp(R))
        }
        if (ACCEPT) {
            j <- j + 1L
            lnL0 <- lnL.b
            lnLcoal0 <- lnLcoal.b
            tree0 <- tr.b
            para0 <- para
            bt0 <- bt
        }
        TREES[[i]] <- tree0
        LL[i] <- lnL0
        params[i, ] <- para0
        LLmod[i] <- lnLcoal0
        i <- i + 1L
    }
    if (verbose) cat("\nDone.\n")
}

.get.list.trees <- function()
    ls(envir = .coalescentMCMCenv, pattern = "^TREES_")

getMCMCtrees <- function(chain = NULL)
{
    list.trees <- .get.list.trees()
    l <- length(list.trees)
    if (is.null(chain)) {
        if (!l) return(NULL)
        if (l == 1)
            return(get(list.trees, envir = .coalescentMCMCenv))
        ## l > 1:
        cat("Several lists of MCMC trees are stored:\n\n")
        for (i in 1:l) cat(i, ":", list.trees[i], "\n")
        cat("\nReturn which number? ")
        chain <- as.numeric(readLines(n = 1))
    } else {
        if (!l) {
            warning("no list of MCMC trees stored")
            return(NULL)
        }
        if (l < chain) {
            warning("no enough lists of MCMC trees stored")
            return(NULL)
        }
    }
    get(paste("TREES", sprintf("%03d", chain), sep = "_"),
        envir = .coalescentMCMCenv)
}

saveMCMCtrees <- function(destdir = ".", format = "RDS", ...)
{
    format <- match.arg(toupper(format), c("RDS", "NEWICK", "NEXUS"))
    switch(format, RDS = {
        FUN <- saveRDS
        suffix <- ".rds"
    }, NEWICK = {
        FUN <- write.tree
        suffix <- ".tre"
    }, NEXUS = {
        FUN <- write.nexus
        suffix <- ".nex"
    })
    list.trees <- .get.list.trees()
    l <- length(list.trees)
    if (!l) warning("no list of trees to save") else {
        for (i in 1:l) {
            f <- list.trees[i]
            outfile <- paste(destdir, "/", f, suffix, sep = "")
            FUN(get(f, envir = .coalescentMCMCenv), outfile, ...)
        }
    }
}

cleanMCMCtrees <- function()
    rm(list = .get.list.trees(), envir = .coalescentMCMCenv)

getLastTree <- function(X) X[[length(X)]]

getMCMCstats <- function()
{
    cat("MCMC chain summaries (chains as columns):\n\n")
    get("MCMCstats", envir = .coalescentMCMCenv)
}

logLik.coalescentMCMC <- function(object, ...)
{
    model <- attr(object, "model")
    if (model == "bsplines") {
        degree <- attr(object, "degree")
        knots <- attr(object, "knots")
    }

    if (model == "bsplines") {
        beta <- vector("list", nrow(object))
        for (i in 1:nrow(object)) beta[[i]] <- object[i, -1]
    }
    para.nms <- colnames(object)[-(1:2)]
    res <- object[, "logLik.coal"]
    attributes(res) <- NULL
    ##if (useWeights) {
    ##    w <- object[, "logLik.tree"]
    ##    attributes(w) <- NULL
    ##    w <- w - min(w) + 1e-8 # all w > 0
    ##    w <- length(res)*w/sum(w) # sum(w) = length(res)
    ##    res <- w * res
    ##}
    attr(res, "nobs") <- attr(object, "nobs")
    attr(res, "df") <- length(para.nms) # number of parameters
    class(res) <- c("logLik")
    res
}

AIC.coalescentMCMC <- function(object, ..., k = 2)
{
    mod <- c(list(object), list(...))
    ll <- lapply(mod, logLik.coalescentMCMC)
    df <-  sapply(ll, attr, "df")
    -2 * sapply(ll, mean) + k * df
}

BIC.coalescentMCMC <- function(object, ...)
{
    mod <- c(list(object), list(...))
    ll <- lapply(mod, logLik.coalescentMCMC)
    df <-  sapply(ll, attr, "df")
    nobs <-  sapply(ll, attr, "nobs")
    -2 * sapply(ll, mean) + df * nobs
}

anova.coalescentMCMC <- function(object, ...)
{
    mod <- c(list(object), list(...))
    ll <- lapply(mod, logLik.coalescentMCMC)
    df <- sapply(ll, attr, "df")
    ll <- sapply(ll, mean)
    dev <- c(NA, 2 * diff(ll)) # LRT's
    ddf <- c(NA, diff(df))
    res <- data.frame(ll, df, ddf, dev, pchisq(dev, ddf, lower.tail = FALSE))
    dimnames(res) <- list(1:length(mod),
                          c("Log lik.", "Resid. df", "df", "Chisq", "Pr(>|Chi|)"))
    structure(res, heading = "Likelihood Ratio Test Table",
              class = c("anova", "data.frame"))
}

subset.coalescentMCMC <- function(x, burnin = 1000, thinning = 10, end = NULL, ...)
{
    oc <- oldClass(x)
    class(x) <- NULL
    attr.x <- attributes(x)
    mcpar <- attr.x$mcpar
    attr.x$old.call <- attr.x$call
    attr.x$call <-  match.call()

    if (mcpar[2] < burnin)
        stop("'burnin' longer than the number of generations")
    from <- burnin + 1
    n <- nrow(x)
    if (!is.null(end)) {
        if (end < from) stop("argument 'end' too small")
        if (end > n) {
            warning("argument 'end' greater than the number of generations: it was ignored")
            end <- n
        }
    } else  end <- n
    x <- x[from:end, , drop = FALSE]
    if (thinning > 1) {
        i <- seq(thinning, nrow(x), thinning)
        x <- x[i, , drop = FALSE]
        ## mcpar[c(1, 3)] <- thinning
    }

    mcpar[2] <- nrow(x)
    attr.x$mcpar <- mcpar
    attr.x$dim <- dim(x)
    attributes(x) <- attr.x
    class(x) <- oc
    x
}
