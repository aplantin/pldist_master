gower <- function(tsf.data, binary) {
  if (binary) { dat = tsf.data$dat.binary
  } else { dat = tsf.data$dat.quant }
  
  n = nrow(dat); m = ncol(dat) 
  out.D <- matrix(0, n, n) 
  
  taxmax <- apply(dat, 2, max)
  taxmin <- apply(dat, 2, min)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      idx = which(taxmax != taxmin)
      out.D[i, j] = out.D[j, i] = sum( abs(dat[i,idx] - dat[j,idx]) / (taxmax[idx] - taxmin[idx]) )/m
    }
  }

  return(out.D) 
}