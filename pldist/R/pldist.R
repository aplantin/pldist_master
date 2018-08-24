setwd("~/Documents/Research/pldist")
source("./pltransform.R")
source("./pl_braycurtis.R")
source("./pl_jaccard.R")
source("./pl_kulczynski.R")
source("./pl_gower.R")
source("./paired_unifrac.R")
source("./longit_unifrac.R")


pldist <- function(otus, metadata, paired = FALSE, binary = FALSE, method, tree = NULL, gam = c(0, 0.5, 1)) {
  ## Find desired method 
  method.opts = c("braycurtis", "jaccard", "kulczynski", "gower", "unifrac")
  this.method = pmatch(trimws(tolower(method)), method.opts, nomatch = NA)
  if (is.na(this.method)) {
    stop("Method does not match any expected methods. Please see list of options in documentation.")
  } 
  method = method.opts[this.method]
  
  ## Check and prepare all input 
  if (nrow(otus) != nrow(metadata)) {
    stop("Number of rows of metadata and OTUs should match") } 
  
  if (any(apply(otus, 1, FUN = function(x) all(x == 0)))) {
    stop("At least one subject has uniformly zero OTU counts. Please exclude.") }
  
  if (any(apply(otus, 2, FUN = function(x) all(x == 0)))) {
    warning("Some OTUs have count zero for all subjects and are being excluded.")
    otus <- otus[, which(apply(otus, 2, FUN = function(x) !all(x == 0)))]
  }
  
  if (!all(colnames(metadata) == c("subjID", "sampID", "time"))) {
    stop("Please format metadata with columns \"subjID\", \"sampID\", \"time\" in that order") }
  
  if (is.null(rownames(otus)) | !all(rownames(otus) == metadata$sampID)) {
    stop("Please ensure rownames of OTU matrix exactly match sample IDs in metadata") }
  
  if (!all(apply(otus, 1, sum) == 1))  otus <- counts2props(otus) 
  
  if (paired) { 
    if (length(unique(metadata$time)) > 2) {
      stop("Paired dissimilarities were requested, but >2 unique time points/groups were provided.")
    } else if (length(unique(metadata$time)) < 2) {
      stop("Paired dissimilarities were requested, but <2 unique time points/groups were provided.")
    } 
    
    persubj <- aggregate(metadata$time, by = list(metadata$subjID), FUN = function(x) length(x))$x
    if (length(unique(persubj)) != 1) {
      stop("Paired dissimilarities were requested, but some groups/subjects do not have 2 observations. \n
           Please check for missing or miscoded data and exclude any unpaired observations.")
    }
    metadata$time = as.numeric(as.factor(metadata$time))
  } else {      
    metadata$time = as.numeric(metadata$time) 
  }
  
  ## Calculate distances/dissimilarities 
  
  if (method.opts[this.method] != "unifrac") {
    ## Calculate transformed data and apply distance (all except UniFrac) 
    tsf.res <- pl.transform(otus = otus, metadata = metadata, paired = paired)
    D <- switch(method, 
                braycurtis = braycurtis(tsf.res$tsf.data, binary = binary), 
                jaccard = jaccard(tsf.res$tsf.data, paired = paired, binary = binary), 
                kulczynski = kulczynski(tsf.res$tsf.data, paired = paired, binary = binary), 
                gower = gower(tsf.res$tsf.data, binary = binary)
                ) 
  } else {
    ## Calculate paired/longitudinal UniFrac dissimilarities 
    if (is.null(tree)) stop("Tree is required for UniFrac family metrics.")
    if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!") 
    if (paired) {
      D <- PUniFrac(otu.tab = otus, tree = tree, gam = gam, metadata = metadata)
    } else {
      D <- LUniFrac(otu.tab = otus, tree = tree, gam = gam, metadata = metadata)
    }
  } 
  
  if (paired) {
    if (method != "unifrac") {
      if (binary) {
        type = paste("Method: ", method, "; Paired, Binary") 
      } else {
        type = paste("Method: ", method, "; Paired, Quantitative") 
      }
    } else {
      type = paste("UniFrac family; Paired") 
    } 
  } else {
    if (method != "unifrac") {
      if (binary) {
        type = paste("Method: ", method, "; Longitudinal, Binary") 
      } else {
        type = paste("Method: ", method, "; Longitudinal, Quantitative") 
      }
    } else {
      type = paste("UniFrac family; Longitudinal") 
    } 
  }
  return(list(D = D, type = type)) 
}