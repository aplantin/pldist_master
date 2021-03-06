% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pl_gower.R
\name{gower}
\alias{gower}
\title{Paired or longitudinal Gower distances}
\usage{
gower(tsf.data, binary)
}
\arguments{
\item{tsf.data}{Transformed OTU table and metadata (from function pl.transform)}

\item{binary}{Logical indicating whether to use the binary version of the distance}
}
\value{
Returns an n x n distance matrix.
}
\description{
The distances are calculated as follows, where d_k^X is the within-subject 
    measure of change appropriate to the setting (paired/longitudinal and 
    quantitative/qualitative), as described in the full package documentation 
    and vignette. 
    
\eqn{D_{AB} = (1/m) \sum_k (|d_k^A - d_k^B|)/(\max d_k - \min d_k)}
}
