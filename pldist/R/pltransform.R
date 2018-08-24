## Input: 
##    Matrix of OTU counts or proportions 
##        - Will be transformed to proportions if it's not already
##        - Row names must be sample identifiers (matching metadata) 
##        - Columns are OTUs 
##    Metadata (subject ID, sample ID, group/time) 
## 
## Output: 
##    Matrices (n x p) of transformed data (both quantitative and binary) 
##        - Row names are subject identifiers 
##    Note of which transformations were used (paired, balanced, or unbalanced)
##        - With a warning if unbalanced 
## 

### Paired transformation 
tsf_paired <- function(otus, metadata, out.data) {
  out.binary = out.quant = out.avgprop = out.data 
  for (i in 1:nrow(out.data)) {
    t1.idx <- which(metadata$subjID == rownames(out.data)[i] & metadata$time == 1)
    t2.idx <- which(metadata$subjID == rownames(out.data)[i] & metadata$time == 2)
    out.binary[i, ] <- 0.5 * (as.numeric(otus[t2.idx,] > 0) - as.numeric(otus[t1.idx,] > 0)) 
    nonz <- which(otus[t2.idx,] != 0 | otus[t1.idx,] != 0) 
    out.quant[i, nonz] <- 0.5 * (otus[t2.idx, nonz] - otus[t1.idx, nonz]) / (otus[t2.idx, nonz] + otus[t1.idx, nonz]) 
    out.avgprop[i, ] <- 0.5 * (otus[t2.idx, ] + otus[t1.idx, ])
  }
  return(list(dat.binary = out.binary, dat.quant = out.quant, avg.prop = out.avgprop))   
} 

### Longitudinal transformation
tsf_long <- function(otus, metadata, out.data) {
  out.binary = out.quant = out.data 
  out.avgprop = out.data 
  
  for (i in 1:nrow(out.data)) {
    ## Prep subject 
    subj.idx <- which(metadata$subjID == rownames(out.data)[i])
    subj.otu <- otus[subj.idx, ]
    subj.times <- metadata$time[subj.idx] 
    
    ord <- order(metadata$time[subj.idx])
    subj.otu <- subj.otu[ord, ]
    subj.times <- subj.times[ord]
    qi <- nrow(subj.otu)
    
    ## Calculate both 
    dk.uw <- rep(0, ncol(otus))
    dk.g <- rep(0, ncol(otus))
    cumprop <- subj.otu[1,] 
    for (j in 1:(qi-1)) {
      dk.uw = dk.uw + (1/(subj.times[j+1] - subj.times[j])) * abs(as.numeric(subj.otu[(j+1), ] > 0) - as.numeric(subj.otu[j, ] > 0))
      nonz <- which(subj.otu[(j+1), ] != 0 | subj.otu[j, ] != 0)
      dk.g[nonz] = dk.g[nonz] + (1/(subj.times[j+1] - subj.times[j])) * 
        abs((subj.otu[(j+1), nonz] - subj.otu[j, nonz])/(subj.otu[(j+1), nonz] + subj.otu[j, nonz]))
      cumprop = cumprop + subj.otu[(j+1), ]
    }
    dk.uw = dk.uw/(2*qi)
    dk.g = dk.g/(2*qi)
    cumprop = cumprop/qi 
    
    ## Fill row 
    out.binary[i, ] <- dk.uw 
    out.quant[i, ] <- dk.g 
    out.avgprop[i, ] <- cumprop 
  }
  return(list(dat.binary = out.binary, dat.quant = out.quant, avg.prop = out.avgprop)) 
}

### Convert counts to proportions 
counts2props <- function(x) {
  return(t(apply(x, 1, FUN = function(y) y/sum(y))))
}

### Full transformation function 
pl.transform <- function(otus, metadata, paired) {
  ## Check that otus are proportions 
  if (all(apply(otus, 1, sum) != 1)) {
    otus <- counts2props(otus) 
  }
 
  ## prepare output 
  n <- length(unique(metadata$subjID))
  out.data <- matrix(0, nrow = n, ncol = ncol(otus))
  rownames(out.data) <- unique(metadata$subjID)
  colnames(out.data) <- colnames(otus) 
  
  ## calculate appropriate transformations 
  if (paired) {
    res <- tsf_paired(otus, metadata, out.data)
  } else {
    res <- tsf_long(otus, metadata, out.data)
    if (length(unique(table(metadata$time))) != 1) { balanced = FALSE } else { balanced = TRUE } 
  }
  if (paired) { type = "paired" 
  } else if (balanced) { type = "balanced longitudinal"
  } else { 
    type = "unbalanced longitudinal (WARNING: this transformation is not recommended for strongly unbalanced designs!)" 
    warning("WARNING: this transformation is not recommended for strongly unbalanced designs!")
    }
   
  ## return 
  return(list(tsf.data = res, type = type))
}


