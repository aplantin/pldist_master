kulczynski <- function(tsf.data, paired, binary) {
  if (binary) { dat = tsf.data$dat.binary
  } else { dat = tsf.data$dat.quant }
  
  n = nrow(dat); m = ncol(dat) 
  out.D <- matrix(0, n, n) 
  
  if (paired) {
    if (binary) {
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          out.D[i, j] = out.D[j, i] = 1 - sum(dat[i, ] == dat[j, ] & dat[i, ] != 0)/m 
        }
      }
    } else {   # paired but not binary 
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          num = sum(pmin(abs(dat[i,]), abs(dat[j,])) * flexsign(dat[i, ], dat[j, ]))
          out.D[i, j] = out.D[j, i] = 1 - 2*num/m
        }
      }
    }
  } else {  # not paired (longitudinal) 
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        out.D[i,j] = out.D[j,i] = 1 - sum(pmin(dat[i,],dat[j,])) / m
      }
    }
  }
  return(out.D) 
}