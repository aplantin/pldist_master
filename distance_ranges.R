#setwd("~/Documents/Research/Methods_indiv/2018_LUniFrac/Code/pldist/pldist/R/")
setwd("~/Documents/Research/pldist/pldist/R/")
source("./pltransform.R")
source("./pl_braycurtis.R")
source("./pl_jaccard.R")
source("./pl_kulczynski.R")
source("./pl_gower.R")
source("./paired_unifrac.R")
source("./longit_unifrac.R")
source("./pldist.R")
library(ape)


#### Tests for paired data 

gen.data.paired <- function(seed, nsubj, notus, nzero, maxct = 2500) {
  set.seed(seed) 
  ncells = nsubj * notus * 2 
  toy.otus <- matrix(sample(c(sample(1:maxct, (ncells - nzero)), rep(0, nzero))), nrow = (2*nsubj), ncol = notus)
  while (any(c(apply(toy.otus, 1, FUN = function(x) all(x == 0)), 
               apply(toy.otus, 2, FUN = function(x) all(x == 0))))) {
    toy.otus <- matrix(sample(c(sample(1:maxct, (ncells - nzero), replace = TRUE), rep(0, nzero))), nrow = (2*nsubj), ncol = notus)
  }
  toy.props <- counts2props(toy.otus)
  toy.meta <- data.frame(subjID = rep(paste("subj", 1:nsubj, sep = ""), each = 2), 
                         sampID = paste(rep(paste("subj", 1:nsubj, sep = ""), each = 2), rep(c("a", "b"), nsubj), sep = ""), 
                         time  = rep(c(1,2), nsubj), stringsAsFactors = FALSE)
  rownames(toy.otus) = rownames(toy.props) = toy.meta$sampID
  colnames(toy.otus) = paste("otu", 1:notus, sep = "")
  return(list(otus = toy.otus, metadata = toy.meta))
}

gen.tree <- function(seed, notus) {
  set.seed(seed)
  sim.tree = rtree(n=notus)
  sim.tree$tip.label <- paste("otu", sample(1:notus), sep = "")
  return(sim.tree)
}

do.many.paired <- function(bigB, nsubj, notus, nzero, maxct) {
  out <- vector("list", 10)
  for (i in 1:10) { out[[i]] <- vector() }
  
  for (i in 1:bigB) {
    dat <- gen.data.paired(seed = i, nsubj = nsubj, notus = notus, nzero = nzero, maxct = maxct)
    tree <- gen.tree(i, notus)
    bray.pq <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[1]] <- c(out[[1]], bray.pq[upper.tri(bray.pq)])
    bray.pb <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = TRUE, method = "bray", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[2]] <- c(out[[2]], bray.pb[upper.tri(bray.pb)])
    gow.pq <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = FALSE, method = "gow", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[3]] <- c(out[[3]], gow.pq[upper.tri(gow.pq)]) 
    gow.pb <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = TRUE, method = "gow", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[4]] <- c(out[[4]], gow.pb[upper.tri(gow.pb)]) 
    jac.pq <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = FALSE, method = "jac", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[5]] <- c(out[[5]], jac.pq[upper.tri(jac.pq)]) 
    jac.pb <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = TRUE, method = "jac", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[6]] <- c(out[[6]], jac.pb[upper.tri(jac.pb)]) 
    kul.pq <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = FALSE, method = "kul", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[7]] <- c(out[[7]], kul.pq[upper.tri(kul.pq)]) 
    kul.pb <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = TRUE, method = "kul", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[8]] <- c(out[[8]], kul.pb[upper.tri(kul.pb)]) 
    unifracs <- pldist(dat$otus, dat$metadata, paired = TRUE, binary = NULL, method = "unifrac", tree = tree, gam = c(0.5))$D
    uf.pq <- unifracs[,,"d_0.5"] 
    uf.pb <- unifracs[,,"d_UW"]
    out[[9]] <- c(out[[9]], uf.pq[upper.tri(uf.pq)])
    out[[10]] <- c(out[[10]], uf.pb[upper.tri(uf.pb)])
  }
  names(out) <- c("bray_quant", "bray_bin", "gow_quant", "gow_bin", "jac_quant", "jac_bin", "kul_quant", "kul_bin", "gen_unifrac", "unweighted_unifrac")
  return(out)
}

simres <- do.many.paired(10, nsubj = 10, notus = 20, nzero = 100, maxct = 2500)  
lapply(simres, range)


#### Tests for balanced longitudinal data

gen.data.longit <- function(seed, nsubj, ntimes = 4, notus, propzero = 0.5, maxct = 2500) {
  set.seed(seed)
  ncells = nsubj * notus * ntimes 
  nzero = floor(ncells*propzero)
  toy.otus <- matrix(0, nrow = nsubj*ntimes, ncol = notus)
  while (any(c(apply(toy.otus, 1, FUN = function(x) all(x == 0)), 
               apply(toy.otus, 2, FUN = function(x) all(x == 0))))) {
    toy.otus <- matrix(sample(c(sample(1:maxct, (ncells - nzero), replace = TRUE), rep(0, nzero))), nrow = nsubj*ntimes, ncol = notus)
  }
  toy.props <- counts2props(toy.otus) 
  toy.meta <- data.frame(subjID = rep(paste("subj", 1:nsubj, sep = ""), each = ntimes), 
                         sampID = paste(rep(paste("subj", 1:nsubj, sep = ""), each = ntimes), 
                                        rep(letters[1:ntimes], nsubj), sep = ""), 
                         time  = rep(1:ntimes, nsubj), stringsAsFactors = FALSE)
  rownames(toy.otus) = toy.meta$sampID
  colnames(toy.otus) = paste("otu", 1:notus, sep = "")
  return(list(otus = toy.otus, metadata = toy.meta))
}

do.many.longit <- function(bigB, nsubj, ntimes, notus, propzero, maxct) {
  out <- vector("list", 10)
  for (i in 1:10) { out[[i]] <- vector() }
  
  for (i in 1:bigB) {
    dat <- gen.data.longit(seed = i, nsubj = nsubj, ntimes = ntimes, notus = notus, propzero = propzero, maxct = maxct)
    tree <- gen.tree(i, notus)
    bray.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[1]] <- c(out[[1]], bray.pq[upper.tri(bray.pq)])
    bray.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "bray", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[2]] <- c(out[[2]], bray.pb[upper.tri(bray.pb)])
    gow.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "gow", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[3]] <- c(out[[3]], gow.pq[upper.tri(gow.pq)]) 
    gow.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "gow", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[4]] <- c(out[[4]], gow.pb[upper.tri(gow.pb)]) 
    jac.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "jac", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[5]] <- c(out[[5]], jac.pq[upper.tri(jac.pq)]) 
    jac.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "jac", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[6]] <- c(out[[6]], jac.pb[upper.tri(jac.pb)]) 
    kul.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "kul", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[7]] <- c(out[[7]], kul.pq[upper.tri(kul.pq)]) 
    kul.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "kul", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[8]] <- c(out[[8]], kul.pb[upper.tri(kul.pb)]) 
    unifracs <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = NULL, method = "unifrac", tree = tree, gam = c(0.5))$D
    uf.pq <- unifracs[,,"d_0.5"] 
    uf.pb <- unifracs[,,"d_UW"]
    out[[9]] <- c(out[[9]], uf.pq[upper.tri(uf.pq)])
    out[[10]] <- c(out[[10]], uf.pb[upper.tri(uf.pb)])
  }
  names(out) <- c("bray_quant", "bray_bin", "gow_quant", "gow_bin", "jac_quant", "jac_bin", "kul_quant", "kul_bin", "gen_unifrac", "unweighted_unifrac")
  return(out)
}

simres <- do.many.longit(1000, nsubj = 10, ntimes = 4, notus = 20, propzero = 0.8, maxct = 2500)
lapply(simres, range)





#### Tests for unbalanced longitudinal data

gen.data.unbal.longit <- function(seed, nsubj, maxtimes = 4, maxdiff = 5, notus, propzero = 0.5, maxct = 2500) {
  set.seed(seed)
  ntimes <- sample(2:maxtimes, nsubj, replace = TRUE)
  ncells = sum(ntimes) * notus
  nzero = floor(ncells*propzero)
  toy.otus <- matrix(0, nrow = sum(ntimes), ncol = notus)
  while (any(c(apply(toy.otus, 1, FUN = function(x) all(x == 0)), 
               apply(toy.otus, 2, FUN = function(x) all(x == 0))))) {
    toy.otus <- matrix(sample(c(sample(1:maxct, (ncells - nzero), replace = TRUE), rep(0, nzero))), nrow = sum(ntimes), ncol = notus)
  }
  toy.props <- counts2props(toy.otus) 
  toy.meta <- data.frame(subjID = unlist(sapply(1:nsubj, FUN = function(i) rep(paste("subj", i, sep = ""), ntimes[i]), simplify = FALSE)), 
                         sampID = paste(unlist(sapply(1:nsubj, FUN = function(i) rep(paste("subj", i, sep = ""), ntimes[i]), simplify = FALSE)), 
                                        unlist(sapply(1:nsubj, FUN = function(i) letters[1:ntimes[i]], simplify = FALSE)), sep = ""), 
                         time  = unlist(sapply(1:nsubj, FUN = function(i) cumsum(c(1, sample(1:maxdiff, ntimes[i]-1, replace = TRUE))), simplify = FALSE)), 
                         stringsAsFactors = FALSE)
  rownames(toy.otus) = toy.meta$sampID
  colnames(toy.otus) = paste("otu", 1:notus, sep = "")
  return(list(otus = toy.otus, metadata = toy.meta))
}

do.many.unbal.longit <- function(bigB, nsubj, maxtimes, maxdiff, notus, propzero, maxct) {
  out <- vector("list", 10)
  for (i in 1:10) { out[[i]] <- vector() }
  
  for (j in 1:bigB) {
    dat <- gen.data.unbal.longit(seed = j, nsubj = nsubj, maxtimes = maxtimes, maxdiff = maxdiff, notus = notus, propzero = propzero, maxct = maxct)
    tree <- gen.tree(seed = j, notus)
    bray.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[1]] <- c(out[[1]], bray.pq[upper.tri(bray.pq)])
    bray.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "bray", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[2]] <- c(out[[2]], bray.pb[upper.tri(bray.pb)])
    gow.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "gow", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[3]] <- c(out[[3]], gow.pq[upper.tri(gow.pq)]) 
    gow.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "gow", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[4]] <- c(out[[4]], gow.pb[upper.tri(gow.pb)]) 
    jac.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "jac", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[5]] <- c(out[[5]], jac.pq[upper.tri(jac.pq)]) 
    jac.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "jac", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[6]] <- c(out[[6]], jac.pb[upper.tri(jac.pb)]) 
    kul.pq <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = FALSE, method = "kul", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[7]] <- c(out[[7]], kul.pq[upper.tri(kul.pq)]) 
    kul.pb <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = TRUE, method = "kul", tree = NULL, gam = c(0, 0.5, 1))$D
    out[[8]] <- c(out[[8]], kul.pb[upper.tri(kul.pb)]) 
    unifracs <- pldist(dat$otus, dat$metadata, paired = FALSE, binary = NULL, method = "unifrac", tree = tree, gam = c(0.5))$D
    uf.pq <- unifracs[,,"d_0.5"] 
    uf.pb <- unifracs[,,"d_UW"]
    out[[9]] <- c(out[[9]], uf.pq[upper.tri(uf.pq)])
    out[[10]] <- c(out[[10]], uf.pb[upper.tri(uf.pb)])
  }
  names(out) <- c("bray_quant", "bray_bin", "gow_quant", "gow_bin", "jac_quant", "jac_bin", "kul_quant", "kul_bin", "gen_unifrac", "unweighted_unifrac")
  return(out)
}

simres <- do.many.unbal.longit(1000, nsubj = 10, maxtimes = 4, maxdiff = 5, notus = 20, propzero = 0.5, maxct = 2500)
lapply(simres, range)


