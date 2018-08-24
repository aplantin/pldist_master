flexsign <- function(v1, v2) {
  return(as.numeric( (v1 >= 0 & v2 >= 0) | (v1 <= 0 & v2 <= 0) ))
}


#' Paired or longitudinal Jaccard distances 
#'
#' @param tsf.data Transformed OTU table and metadata (from function pl.transform)
#' @param paired Logical indicating whether paired analysis is desired 
#' @param binary Logical indicating whether to use the binary version of the distance 
#' @return Returns an n x n distance matrix. 
#'
#' @export
#' 
jaccard <- function(tsf.data, paired, binary) {
  if (binary) { dat = tsf.data$dat.binary
  } else { dat = tsf.data$dat.quant }
  
  n = nrow(dat); m = ncol(dat) 
  out.D <- matrix(0, n, n) 
  
  if (paired) {
    if (binary) {
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          num = sum(as.numeric(dat[i, ] == dat[j, ] & dat[i, ] != 0))
          denom = sum(as.numeric(dat[i, ] != 0)) + sum(as.numeric(dat[j, ] != 0))
          if (denom != 0) {
            out.D[i,j] = out.D[j,i] = 1 - num/denom
          } else {
            out.D[i,j] = out.D[j,i] = 0
          }
          
        }
      }
    } else { # paired, not binary 
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          idx = which(pmax(abs(dat[i,]), abs(dat[j,])) != 0)
          num = sum(pmin(abs(dat[i,idx]), abs(dat[j,idx])) * flexsign(dat[i,idx], dat[j,idx]))
          denom = sum(pmax(abs(dat[i,idx]), abs(dat[j,idx])))
          out.D[i, j] = out.D[j, i] = 1 - num/denom
        }
      }
    }
  } else {  # not paired 
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        idx = which(pmax(dat[i,], dat[j,]) != 0)
        out.D[i,j] = out.D[j,i] = 1 - sum(pmin(dat[i,idx],dat[j,idx])) / sum(pmax(dat[i,idx],dat[j,idx]))
      }
    }
  }
  return(out.D) 
}