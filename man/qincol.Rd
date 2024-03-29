% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mincol.R
\name{qincol}
\alias{qincol}
\title{Convert a symmetric matrix into a vector}
\usage{
qincol(m)
}
\arguments{
\item{m}{A symmetric matrix or dimension \code{rxr}.}
}
\value{
A vector of length \code{r(r+1)/2}.
}
\description{
Calling \code{\link[=qincol]{qincol()}} on a symmetric matrix with dimension \code{rxr} converts
this matrix a vector of length \code{r(r+1)/2}. This function is used to
handle the MCMC sampling output from multivariate finite mixture models. To
save storage the symmetric variance-covariance matrices of multivariate
mixtures are stored vector form. If the covariance matrices are needed for
calculations the functions \code{\link[=qinmatr]{qinmatr()}} and \code{\link[=qinmatrmult]{qinmatrmult()}} helps to restore
these matrices from the storage vectors.
}
\examples{
# Define a vector.
q <- rnorm(n = 6, mean = 0.5, sd = 2)
# Generate the symmetric matrix.
mat <- qinmatr(q)
# Convert the matrix back into the vector.
qincol(mat)

}
\seealso{
\itemize{
\item \code{\link[=qinmatr]{qinmatr()}} for converting a single vector into a symmetric matrix
\item \code{\link[=qinmatrmult]{qinmatrmult()}} for converting multiple vectors into symmetric matrices
\item \code{\link[=qincolmult]{qincolmult()}} for converting multiple symmetric matrice into vectors
}
}
