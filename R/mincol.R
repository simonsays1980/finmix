#' Convert vector into matrix.
#' 
#' @description 
#' Calling [qinmatr()] on a vector of dimension `r(r+1)/2x1` 
#' converts the vector into a symmetric matrix of dimension `rxr`. This 
#' function is used to handle the MCMC sampling output from multivariate finite 
#' mixture models. To save storage the symmetric variance-covariance matrices 
#' of multivariate mixtures are stored vector form. If the covariance matrices 
#' are needed for calculations this function helps to restore these matrices 
#' from the storage vectors.
#' 
#' @param q A vector of dimension `r(r+1)/2x1`.
#' @return A symmetric matrix of dimension `rxr`.
#' @export
#' 
#' @examples 
#' # Define a vector.
#' q <- rnorm(n = 6, mean = 0.5, sd = 2)
#' # Generate the symmetric matrix.
#' qinmatr(q)
#' 
#' @seealso 
#' * [qinmatrmult()] 
#' * [qincol()]
#' * [qincolmult()]
"qinmatr" <- function(q) {
  if (length(dim(q)) > 0) {
    stop(paste("The argument 'q' has to be an array or matrix in column or 
               row format, i.e. one dimension is 1.",
      sep = ""
    ),
    call. = TRUE
    )
  }
  r <- -.5 + sqrt(.25 + 2 * length(q))
  tmp <- matrix(numeric(), nrow = r, ncol = r)
  tmp[upper.tri(tmp, diag = TRUE)] <- q
  tmp[lower.tri(tmp)] <- t(tmp[upper.tri(tmp)])
  return(tmp)
}

#' Convert array of vectors into array of matrices.
#' 
#' @description 
#' Calling [qinmatrmult()] on multiple vectors of dimension `r(r+1)/2x1` 
#' converts these vectors into an array of symmetric matrices of dimension 
#' `rxr`. This function is used to handle the MCMC sampling output from 
#' multivariate finite mixture models. To save storage the symmetric 
#' variance-covariance matrices of multivariate mixtures are stored vector 
#' form. If the covariance matrices are needed for calculations this function 
#' helps to restore these matrices from the storage vectors.
#' 
#' @param q A matrix or array of vectors of dimension `r(r+1)/2x1`. 
#' @return An array of symmetric matrices, all of dimension `rxr`.
#' @export
#' 
#' @examples 
#' # Convert a matrix of vectors
#' qinmatrmult(matrix(rnorm(36), nrow = 6))
#' 
#' @seealso 
#' * [qinmatr()] for converting a single vector into a symmetric matrix
#' * [qincol()] for converting a symmetric matrix into a vector
#' * [qincolmult()] for converting an array of symmetric matrices into vectors
"qinmatrmult" <- function(m) {
  r <- -.5 + sqrt(.25 + 2 * nrow(m))
  tmp.array <- array(numeric(), dim = c(r, r, ncol(m)))
  for (k in 1:ncol(m)) {
    tmp.array[, , k] <- qinmatr(m[, k])
  }
  return(tmp.array)
}

#' Convert a symmetric matrix into a vector
#' 
#' @description 
#' Calling [qincol()] on a symmetric matrix with dimension `rxr` converts 
#' this matrix a vector of length `r(r+1)/2`. This function is used to 
#' handle the MCMC sampling output from multivariate finite mixture models. To 
#' save storage the symmetric variance-covariance matrices of multivariate 
#' mixtures are stored vector form. If the covariance matrices are needed for 
#' calculations the functions [qinmatr()] and [qinmatrmult()] helps to restore 
#' these matrices from the storage vectors.
#' 
#' @param q A symmetric matrix or dimension `rxr`. 
#' @return A vector of length `r(r+1)/2`.
#' @export
#' 
#' @examples 
#' # Define a vector.
#' q <- rnorm(n = 6, mean = 0.5, sd = 2)
#' # Generate the symmetric matrix.
#' mat <- qinmatr(q)
#' # Convert the matrix back into the vector.
#' qincol(mat)
#' 
#' @seealso 
#' * [qinmatr()] for converting a single vector into a symmetric matrix
#' * [qinmatrmult()] for converting multiple vectors into symmetric matrices
#' * [qincolmult()] for converting multiple symmetric matrice into vectors
"qincol" <- function(m) {
  r <- ncol(m)
  index <- 0
  s <- r * (r + 1) / 2
  qcol <- vector("numeric", s)
  for (rr in 1:r) {
    qcol[(index + 1):(index + rr)] <- m[1:rr, rr]
    index <- index + rr
  }
  return(qcol)
}

#' Convert multiple symmetric matrices into vectors
#' 
#' @description 
#' Calling [qincolmult()] on an array of symmetric matrices all with dimension 
#' `rxr` converts these matrices into an array of vectors with length 
#' `r(r+1)/2`. This function is used to handle the MCMC sampling output from 
#' multivariate finite mixture models. To save storage the symmetric 
#' variance-covariance matrices of multivariate mixtures are stored vector 
#' form. If the covariance matrices are needed for calculations the functions 
#' [qinmatr()] and [qinmatrmult()] helps to restore these matrices from the 
#' storage vectors.
#' 
#' @param q A symmetric matrix or dimension `rxr`. 
#' @return A vector of length `r(r+1)/2`.
#' @export
#' 
#' @examples 
#' # Convert a matrix of vectors
#' matrices <- qinmatrmult(matrix(rnorm(36), nrow = 6))
#' # Convert these matrices back into vectors.
#' qincolmult(matrices) 
#' 
#' @seealso 
#' * [qinmatr()] for converting a single vector into a symmetric matrix
#' * [qinmatrmult()] for converting multiple vectors into symmetric matrices
#' * [qincol()] for converting a single symmetric matrix into a vector
"qincolmult" <- function(a) {
  r <- dim(a)[1]
  K <- dim(a)[3]
  s <- r * (r + 1) / 2
  tmp.mat <- matrix(numeric(), nrow = s, ncol = K)
  for (k in 1:K) {
    tmp.mat[, k] <- qincol(a[, , k])
  }
  return(tmp.mat)
}
