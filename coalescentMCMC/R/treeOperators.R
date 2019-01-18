## treeOperators.R (2019-01-17)

##   Trees Operators for Running MCMC

## Copyright 2012-2019 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

.branching_times <- function(phy)
{
    phy <- reorder(phy)
    .Call(branchingTimesCall, phy$edge,  phy$edge.length)
}

NeighborhoodRearrangement <- function(phy, n, THETA, brtimes)
{
    e1 <- phy$edge[, 1L] # local copy
    e2 <- phy$edge[, 2L]

    ## sample a node excluding the root:
    target <- sample.int(n - 2L, 1L) + n + 1L

### i1, i2, and i3 are indices of some edges
### target, anc, and sister are indices of some nodes

    i1 <- which(e2 == target)
    anc <- e1[i1] # the ancestor of 'target'
    i2 <- which(e1 == target) # the 2 edges where 'target' is basal
    i3 <- which(e1 == anc) # this includes i1, so:
    i3 <- i3[i3 != i1]
    sister <- e2[i3] # the sister-node of 'target'
    sel <- sample.int(2L, 1L)
    i2.move <- i2[sel]
    i2.stay <- i2[-sel]
    phy$edge[i3, 2L] <- child2move <- e2[i2.move]
    child2stay <- e2[i2.stay]
    phy$edge[i2.move, 2L] <- sister

    bt.anc <- brtimes[anc - n]
    bt.child2move <- if (child2move > n) brtimes[child2move - n] else 0
    bt.sister <- if (sister > n) brtimes[sister - n] else 0
    bt.child2stay <- if (child2stay > n) brtimes[child2stay - n] else 0

    ## now adjust branch lengths:
    ## adjust the branch length that was subtending 'sister':
    phy$edge.length[i3] <- bt.anc - bt.child2move

    ## random age for 'target' between the ones of 'sister' and 'anc':
    ## newage <- runif(1, max(bt.sister, bt.child2stay), bt.anc)

    ## random coalescent time:
    ## (gives a slightly better acceptance rate than the previous one)
    agemin <- max(bt.sister, bt.child2stay)
    pmax <- 1 - exp(-THETA * (bt.anc - agemin))
    p <- runif(1, 0, pmax)
    newage <- agemin - log(1 - p) / THETA

    phy$edge.length[i1] <- bt.anc - newage
    phy$edge.length[i2.move] <- newage - bt.sister
    ## adjust the branch length below the child that has not been moved:
    phy$edge.length[i2.stay] <- newage - bt.child2stay

    attr(phy, "order") <- NULL
    phy <- reorder(phy)
    newNb <- integer(2*n - 1L)
    newNb[n + 1L] <- n + 1L
    sndcol <- phy$edge[, 2L] > n
    phy$edge[sndcol, 2L] <- newNb[phy$edge[sndcol, 2L]] <- 1:(n - 2L) + n + 1L
    phy$edge[, 1L] <- newNb[phy$edge[, 1L]]
    phy
}

TipInterchange <- function(phy, n)
{
    e1 <- phy$edge[, 1L]
    e2 <- phy$edge[, 2L]
    repeat {
        ij <- sample.int(n, size = 2L)
        k <- match(ij, e2)
        ## check that the two tips in 'k' are not sisters
        if (e1[k[1]] != e1[k[2]]) break
    }
    phy$edge[k, 2L] <- ij[2:1]
    phy
}

###EdgeLengthJittering <- function(phy)
###### all edge lengths are added a random value on U[-MIN, MAX]
###### (the ultrametric nature of the tree is kept)
###{
###    el <- phy$edge.length
###    MIN <- min(el)
###    MAX <- max(el)
###    phy$edge.length <- el + runif(1, -MIN, MAX)
###    phy
###}

## moves from Drummond et al. (2002, Genetics):

ScalingMove <- function(phy)
{
    ## values decreased compared to Drummond et al.'s
    phy$edge.length <- runif(1, 0.95, 1.05) * phy$edge.length
    phy
}

branchSwapping <- function(phy, n, brtimes)
{
    e1 <- phy$edge[, 1L]
    e2 <- phy$edge[, 2L]
    Nedge <- length(e1)
    ROOT <- n + 1L
    ## sample directly the edges:
    x <- sample.int(Nedge, 2L)
    wi <- x[1L]
    wj <- x[2L]
    i <- e2[wi]
    j <- e2[wj]
    ip <- e1[wi]
    if (ip == ROOT) return(phy)
    jp <- e1[wj]
    if (ip == j || ip == jp || i == jp) return(phy)
    if (jp > n && i > n && brtimes[jp - n] < brtimes[i - n]) return(phy)
    wip <- which(e2 == ip)
    ipp <- e1[wip]
    tmp <- which(e1 == ip)
    witild <- tmp[tmp != wi]
    itild <- e2[witild]
    ## adjust edges
    phy$edge[witild, 2L] <- j
    phy$edge[wj, 2L] <- ip
    phy$edge[wip, 2L] <- itild
    ##if (j != ROOT) { # j cannot be the root if phy is rooted
    ## adjust branch lengths
    phy$edge.length[wip] <- phy$edge.length[wip] + phy$edge.length[witild]
    bti <- if (i > n) brtimes[i - n] else 0
    btj <- if (j > n) brtimes[j - n] else 0
    btjp <- brtimes[jp - n]
    new.age.ip <- runif(1, max(bti, btj), btjp)
    phy$edge.length[wj] <- btjp - new.age.ip
    phy$edge.length[wi] <- new.age.ip - bti
    phy$edge.length[witild] <- new.age.ip - btj
    ##}
    reorder(phy, "postorder")
}

subtreeExchange <- function(phy, n, brtimes)
{
    phybak <- phy
    e1 <- phy$edge[, 1L]
    e2 <- phy$edge[, 2L]
    Nedge <- length(e1)
    ROOT <- n + 1L
    ## sample directly an edge:
    wi <- sample.int(Nedge, 1L)
    i <- e2[wi]
    if (i == ROOT) return(phy)
    ip <- e1[wi]
    if (ip == ROOT) return(phy)
    wip <- which(e2 == ip)
    jp <- e1[wip]
    tmp <- which(e1 == jp)
    wj <- tmp[tmp != wip]
    j <- e2[wj]
    if (j > n) # no need to test the following condition if j is a tip
        if (brtimes[ip - n] < brtimes[j - n]) return(phy)
    phy$edge[wi, 2L] <- j
    phy$edge[wj, 2L] <- i
    ## no need to adjust edge lengths if both i and j are tips
    if (i > n || j > n) {
        a <- phy$edge.length[wi]
        b <- phy$edge.length[wip]
        c <- phy$edge.length[wj]
        ab <- a + b
        phy$edge.length[wj] <- ab
        phy$edge.length[wi] <- c * a/ab
        phy$edge.length[wip] <- c * b/ab
        tmp <- which(e1 == ip)
        wib <- tmp[tmp != wi]
        phy$edge.length[wib] <- phy$edge.length[wib] - phy$edge.length[wip] + b
    }
    if (any(phy$edge.length <= 0)) phybak else reorder(phy, "postorder")
}

NodeAgeMove <- function(phy, n, brtimes)
{
    phybak <- phy
    e1 <- phy$edge[, 1L]
    e2 <- phy$edge[, 2L]
    Nedge <- length(e1)
    ROOT <- n + 1L
    ## sample a node:
    i <- sample.int(n - 1L, 1L) + n
    wi <- which(e1 == i)
    j <- e2[wi[1]]
    k <- e2[wi[2]]
    ti <- brtimes[i - n]
    tj <- if (j > n) brtimes[j - n] else 0
    tk <- if (k > n) brtimes[k - n] else 0
    mx <- max(tj, tk)
    if (i == ROOT) {
        delta <- runif(1, 0.8333333, 1.2)
        newti <- mx + delta * (ti - delta * mx)
        dif <- newti - ti
    } else {
        wip <- which(e2 == i)
        ip <- e1[wip]
        newti <- runif(1, mx, brtimes[ip - n])
        dif <- newti - ti
        phy$edge.length[wip] <- phy$edge.length[wip] - dif
    }
    phy$edge.length[wi] <- phy$edge.length[wi] + dif
    if (any(phy$edge.length <= 0)) phybak else phy
}

randomWalkThetaMu <-
    function(theta, mu, range.theta = c(0, Inf), range.mu = c(0, Inf))
{
    wtheta <- 1e-3
    wmu <- 1e-3
    newtheta <- theta + runif(1, -wtheta, wtheta)
    newmu <- mu + runif(1, -wmu, wmu)
    if (newtheta < range.theta[1] || newtheta > range.theta[2]) newtheta <- theta
    if (newmu < range.mu[1] || newmu > range.mu[2]) newmu <- mu
    c(theta, mu)
}
