setwd("~/Documents/Research/pldist")
source("./pltransform.R")

## Test paired transformations
set.seed(1)
toy.otus <- matrix(sample(c(1:9, 0, 0, 0)), nrow = 4, ncol = 3)
toy.props <- counts2props(toy.otus)
toy.meta <- data.frame(subjID = c("subj1", "subj1", "subj2", "subj2"), 
                       sampID = c("samp1a", "samp1b", "samp2b", "samp2a"), 
                       group  = c(1, 2, 2, 1), stringsAsFactors = FALSE); toy.meta 
rownames(toy.otus) = rownames(toy.props) = toy.meta$sampID


toy.props 
pl.transform(toy.props, toy.meta, paired = FALSE)
pl.transform(toy.otus, toy.meta)


## Test longitudinal (un)balanced transformations 
set.seed(1)
toy.otus <- matrix(sample(c(1:18, rep(0,6))), nrow = 8, ncol = 3)
toy.meta <- data.frame(subjID = rep(paste("subj", 1:2, sep = ""), each = 4), 
                       sampID = paste(rep(paste("subj", 1:2, sep = ""), each = 4), 
                                      rep(c("a","b","c","d"), 2), sep = ""), 
                       group  = c(1, 3, 5, 7, 1, 2, 3, 8), stringsAsFactors = FALSE); toy.meta 
rownames(toy.otus) = toy.meta$sampID
colnames(toy.otus) = paste("otu", 1:3, sep = "")

pl.transform(toy.otus, toy.meta, paired = FALSE)


# Manual subj 2 
toy.props <- counts2props(toy.otus)
t1t2 <- (toy.props[6,] - toy.props[5,])/(toy.props[6,] + toy.props[5,])
t1t2[2] = 0
t2t3 <- (toy.props[7,] - toy.props[6,])/(toy.props[7,] + toy.props[6,])
t3t4 <- (toy.props[8,] - toy.props[7,])/(toy.props[8,] + toy.props[7,])

(abs(t1t2) + abs(t2t3) + abs(t3t4)/5)/8  ## Checks out! 

