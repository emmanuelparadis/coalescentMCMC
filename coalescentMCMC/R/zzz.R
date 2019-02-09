## zzz.R (2019-02-01)

##   Library Loading

## Copyright 2013-2019 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

.coalescentMCMCenv <- new.env()

assign("MCMCstats",
       data.frame(row.names = c("Nb of generations", "Nb of accepted moves")),
       envir = .coalescentMCMCenv)
