#Functions to generate simulated data.
#===============================================================================
#'Function to generate parameters given number of subjects n and rank parameter K.
#'@param n number of players/subjects.
#'@param K rank of parameter matrix.
#'
#'@return A list containing relevant parameters, including
#' \item{Pi}{Matrix of winning probabilities.}
#' \item{T_max}{Maximum number of comparisons between any pair.}
#'
#' @export
genpar.func <- function( n,  K) {
#T_max: maximum number of comparisons
 T_max = 5

 # Generate Theta: orthonormal matrix from iid normal distribution using QR decomposition
 Z <- matrix(rnorm(n*2*K), nrow = n, ncol = 2*K)
 qr_decomp <- qr(Z)
 Theta <- qr.Q(qr_decomp)

 #Generate J and the parameter matrix M
 J <- getJmat.func(K)*n
 M <- Theta %*% J %*%t(Theta)


 # Compute winning probabilities using the logit function
 Pi <- g.func(M)

  # Returning the output list.
  out <- list(
    "Pi" = Pi, "T_max"=T_max
  )
  return(out)
}

#===============================================================================
#' Function to generate example data given T_max, winning probabilites, n and sparse level.
#' @param  Pi Matrix of winning probabilities.
#' @param T_max Maximum number of comparisons between any pair.
#' @param n number of players/subjects.
#' @param sparse_lv different levels of sparsity. Sparse(1), Less sparse (2) and dense (3).
#'
#' @return  \item{Y}{Matrix of comparison outcome}
#'
#' @export
gendata.func <- function(T_max,Pi,  n, sparse_lv) {
#Set p_n and q_n, the min and max comparison rates relative to different sparse levels.
  sparse_lv_vec <- c(n^(-1)*log(n), n^(-1/2), 1/4)
  p_n <- sparse_lv_vec[sparse_lv]
  q_n <- 4*p_n

  #Compute number of distinct pairs.
  n_pairs <- n*(n-1)/2

  #generate comparison rate pijs for each pair and number of comparisons nijs
  pijs<- runif(n_pairs,p_n,q_n )
  nijs <- rbinom(n_pairs,size = T_max,prob =  pijs)

  #Compute Y, the matrix of comparison results.
  Y <- matrix(0, nrow=n,ncol=n)
  Y_upper <- rbinom(n_pairs,size = nijs,prob =Pi[upper.tri(Pi)])
  Y[upper.tri(Y)] <- Y_upper
  Y <- t(Y)
  Y[upper.tri(Y)] <- nijs - Y_upper
  Y <- t(Y)

  return(Y)
}


#============================================================================
#Function to generate the matrix J for a given rank K.
getJmat.func = function(K){
  J <- matrix(0, nrow=2* K, ncol= 2* K)
  for(s in 1: K){
    J[2*s-1, 2*s]<- 1
    J[2*s,2*s-1]<- -1
  }
  return(J)
}


