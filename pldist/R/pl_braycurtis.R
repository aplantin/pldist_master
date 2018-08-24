braycurtis <- function(tsf.data, binary) {
  if (binary) { dat = tsf.data$dat.binary
  } else { dat = tsf.data$dat.quant }
  
  n = nrow(dat); m = ncol(dat) 
  out.D <- matrix(0, n, n) 
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      out.D[i, j] = out.D[j, i] = 1/m * sum(abs(dat[i, ] - dat[j, ]))
    }
  }
  
  return(out.D) 
}