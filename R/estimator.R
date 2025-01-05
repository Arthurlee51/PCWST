#' @useDynLib PCWST, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#Script contataining main functions to implement the proposed estimator in PCWST using the nonmonotone spectral projected gradient algorithm. Parts of the codes are adapted form the original code available at https://mdav.ece.gatech.edu/software/ for the paper "1-Bit Matrix Completion" to Tailor for skew-symmetric matrix.

#' @title Function to decide scale parameter controlling the nuclear norm constraint.
#' @description Function to choose scale parameter based on training and validation data. Here the scale parameter is equivalent to Cn in the paper.
#'
#' @param scaleVals A numeric vector of candidate values for the tuning parameter.
#' @param Y_train A comparison matrix from the training set, used for estimating the model.
#' @param Y_valid A comparison matrix from the validation set, used for evaluating the likelihood.
#'
#' @return A list containing:
#' \item{out}{The model output corresponding to the optimal scale.}
#' \item{valid_like}{A numeric vector of validation likelihoods for each candidate scale.}
#' \item{scale}{The optimal scale value that maximizes the validation likelihood.}
#'
#'@examples
#'set.seed(123)
#' n=300
#' K=2
#' scaleVals =10^seq(-1, 1, length.out = 20)
#' parameters = genpar.func(n, K)
#' Y_train = gendata.func(T_max = parameters$T_max, Pi = parameters$Pi, n,  sparse_lv =1 )
#' Y_valid = gendata.func(T_max = parameters$T_max, Pi = parameters$Pi, n,  sparse_lv =1 )
#' result = pcwst_est_tune(Y_train,Y_valid, scaleVals)
#'
#' @export
pcwst_est_tune <- function(Y_train, Y_valid,scaleVals) {
  # Get the number of players (or subjects) from the training matrix
  n <- nrow(Y_train)

  # Initialize storage for outputs and validation likelihoods
  out_list  <- list()
  valid_like <- rep(0,length(scaleVals) )

  # Iterate over each candidate scale value
  for(i in 1:length(scaleVals)){
    scale <- scaleVals[i]
    out <- pcwst_est(Y_train,scale)
    out_list[[i]] <- out
    Mhat <- g.func(out$M_mat)
    valid_like[i]<-  rcpp_like(Y_valid, Mhat, n)
  }
  # Identify the scale that maximizes validation likelihood and get the corresponding output
  ind <- which( valid_like==max(valid_like))
  scale <- scaleVals[ind]
  out <- out_list[[ind]]

  # Return the results
  return(list(out=out, valid_like = valid_like, scale = scale))
}

#============================================================================
#' @title Main estimation function
#' @description Function to implement the proposed estimator based on given comparison matrix and scale parameter.
#'
#' @param Y a n x n matrix representing comparison outcomes.
#' @param scale Scale parameter corresponding to C_n in the paper.
#'
#' @return A list containing the result, including
#' \item{M_mat}{Matrix of estimated parameters.}
#' \item{iterations}{Number of iterations.}
#' \item{objective}{Value of objectvice function.}
#' \item{exit}{Exit status. 0 for successful convergence and 1 otherwise.}
#' @examples
#' # Generate example dataset
#' set.seed(123)
#' n=300
#' K=2
#' scale = 2*K
#' parameters <- genpar.func(n, K)
#' Y = gendata.func(T_max = parameters$T_max, Pi = parameters$Pi, n,  sparse_lv =1 )
#' result <-pcwst_est(Y,scale)
#' @export
pcwst_est <- function(Y,scale) {
  # Get the number of players (or subjects)
  n <- ncol(Y)

  #Compute nij_mat, matrix storing the number of comparison between pairs.
  nij_mat <- Y + t(Y)

  # Vectorize nij_mat and Y to prepare inputs for the solver
   n_vec <- nij_mat[upper.tri(nij_mat)]
   y_vec <- Y[upper.tri(Y)]

   # Initialize m0: the starting point for optimization
   m0 <- rep(0, length( n_vec ))

   #Perform estimation using the solver.
  out <- spg_solver(m0, scale,n_vec, y_vec,n,options = list())

  # Return the results
  return(out)
}

#============================================================================
# Define the Spectral Projected-Gradient (SPG) Solver in R, tailored for the proposed method
spg_solver <- function(m0, scale,n_vec, y_vec,n,options = list()) {
  # Default options (if not provided)
  options <- modifyList(list(
    verbosity = 2,        # Verbosity level (0: No output, 1: Major, 2: Detailed)
    iterations = 1000,    # Maximum number of iterations
    optTol = 1e-9,        # Optimality tolerance in z
    stepMin = 1e-16,      # Minimum spectral step
    stepMax = 1e5         # Maximum spectral step
  ), options)

  #  Indicator of convergence results
  EXIT_CONVERGED <- 0
  EXIT_ITERATIONS <- 1#Does not converge after iteration.s

  exit <- 0

  # Initialize variables based on options
  maxIts <- options$iterations
  optTol <- options$optTol
  stepMin <- options$stepMin
  stepMax <- options$stepMax

  # Initialize tracking variables
  iter <- 0
  dxNorm <- 0
  nLineTot <- 0
  condition = TRUE#condition to track whether the loop is continued.
  gStep <- 1  # Initial gradient step size. Equivalent to gamma_k in paper.

  # Initialize solver parameters
  tau <- scale*n # Nuclear norm constraint parameter
  m <- proj_nucnorm(m0,  n, tau) # Project initial point
  pos_ind <- which(n_vec>0)# Indices with non-zero comparisons

  # Compute initial gradient and objective value
  object <- rcpp_grad_val_mat(m[pos_ind], n_vec[pos_ind], y_vec[pos_ind])
  grad <- object$grad
  obj_val <- object$objectiveValue

  # Non-monotone strategy (tracking last function values)
  fBest <- obj_val
  mBest <- m

  # Compute projected gradient direction
  newm <- m
  newm[pos_ind]<- newm[pos_ind] -  gStep*grad
  dm <- proj_nucnorm(newm  ,  n, tau) -m;

  #Compute variables for later convergence checking.
  gtd <- sum(grad * dm[pos_ind]);#Inner product of grad and dm[pos_ind]
  dmNorm <- max(abs(dm))#infinity norm for convergence checking.

  # Main optimization loop
  while (iter < maxIts && condition ) {
    iter <- iter + 1

    # Logging
    if (options$verbosity >= 1) {
      cat(sprintf("Iter %d | Objective: %.6f | dmNorm: %.6e\n", iter, obj_val, dmNorm))
    }

    #save mold, obj_valOld and grad_old in case convergence cannot be reached in the first type of line search.
    mOld <- m
    obj_valOld <- obj_val
    if (!is.null(grad)) {
      grad_old <- grad
    }

    # Projected Gradient Step and Linesearch
    if (!is.null(grad)) {
      spgLineResult <- spgLine(m, dm, gtd, obj_valOld, obj_val, n_vec, y_vec, pos_ind)
      fNew <- spgLineResult$fNew
      mNew <- spgLineResult$mNew
      nLine <- spgLineResult$iter
      lnErr <- (spgLineResult$stat == EXIT_ITERATIONS)
    } else {
      # If gradient is NULL, skip linesearch
      fNew <- fOld
      xNew <- xOld
      nLine <- 0
      lnErr <- FALSE
    }

    #Update values
    obj_val <- fNew
    m <- mNew
    nLineTot <- nLineTot + nLine

    # Handle Line Search Errors by retrying with Curvilinear Linesearch
    if (lnErr) {
      m <- mOld
      obj_val <- mOld

      spgLineCurvyResult <- spgLineCurvy(m, gStep * grad, obj_valOld, tau,n_vec, y_vec, pos_ind,n )
      fNew <- spgLineCurvyResult$fNew
      mNew <- spgLineCurvyResult$mNew
      nLine <- spgLineCurvyResult$iter
      stepG <- spgLineCurvyResult$step
      lnErr <- (spgLineCurvyResult$stat == EXIT_ITERATIONS)

      obj_val <- fNew
      m <- mNew
      nLineTot <- nLineTot + nLine
    }

    # If still error after curvilinear linesearch, break the loop and return result from previous iteration.
    if (lnErr) {
      m <- mOld
      obj_val <- obj_valOld
      if (fBest > obj_val) {
        fBest <- obj_val
        mBest <- m
      }

      #stat <- EXIT_LINE_ERROR
      print("Warning: Does not converge after line search")
      exit <- EXIT_ITERATIONS
      break
      # stat <- EXIT_LINE_ERROR
    }

    # ---------------------------
    # Update Gradient and Step Scaling
    # ---------------------------
    #if (stat == 0) {
    obj <-rcpp_grad_val_mat(m[pos_ind], n_vec[pos_ind], y_vec[pos_ind])
    grad <- obj$grad
    obj_val <- obj$objectiveValue

    s <- m - mOld
    y <- grad - grad_old
    sts <- sum(s^2)
    sty <- sum(s[pos_ind] * y)
    if (sty <= 0) {
      gStep <- stepMax
    } else {
      gStep <- min(stepMax, max(stepMin, sts / sty))
    }

    #Update dm and dmNorm. Make sure step size is great enough.
    newm <- m
    newm[pos_ind]<- newm[pos_ind] -  max(gStep,1)*grad
    dm <- proj_nucnorm(newm  ,  n, tau) -m;
    dmNorm <- max(abs(dm))#infinity norm for convergence checking.

    #Check whether condition 1 or 2 is satisfied for convergence. #Only continue when both are true.
    mNorm <- sqrt(sum(m^2))
    condition1 <- (max(abs(m - mOld)) / max(1, mNorm) > 1e-5)
    condition2 <- (dmNorm > optTol * max(1, mNorm))
    condition <- condition1 & condition2

    # Update best solution found
    if (fBest > obj_val) {
      fBest <- obj_val
      mBest <- m
    }
  }

  # Recover the matrix form of M for output:
  M_mat_prep <- matrix(0, nrow = n,ncol =n)
  M_mat_prep[upper.tri(M_mat_prep)] <- mBest
  M_mat <- M_mat_prep - t(M_mat_prep)

  list(M_mat = M_mat, iterations = iter, objective = fBest, exit = exit)
}

#============================================================================
# Nonmonotone Line Search
spgLine <- function(m, d, gtd, f, fMax, n_vec, y_vec, pos_ind) {
  #initialisation
  # Exit status constants
  EXIT_CONVERGED <- 0
  EXIT_ITERATIONS <- 1

  # Line search parameters
  maxIts_line <- 20 # Maximum iterations
  step <- 1 # Initial step size (alpha in the paper)
  iter_line <- 0 # Line search iteration counter
  gamma <- 1e-4

  while (TRUE) {
    # Evaluate trial point and function value
    mNew <- m + step * d
    fNew <- rcpp_val(mNew[pos_ind],  n_vec[pos_ind],  y_vec[pos_ind])

    # Check exit conditions
    if (fNew < fMax + gamma * step * gtd) {
      stat_line <- EXIT_CONVERGED
      break
    } else if (iter_line >= maxIts_line) {
      stat_line <- EXIT_ITERATIONS
      break
    }

    # Safeguarded quadratic interpolation
    denominator <- 2 * (fNew - f - step * gtd)
    if (denominator == 0) {
      tmp <- step / 2
    } else {
      tmp <- (-gtd * step^2) / denominator
      if (tmp < 0.1 || tmp > 0.9 * step || is.na(tmp)) {
        tmp <- step / 2
      }
    }

    step <- tmp # Update step size
    iter_line <- iter_line + 1
  }

  # Return updated values and status
  return(list(fNew = fNew, mNew = mNew, iter = iter_line, stat = stat_line))
}

#============================================================================
# Curvilinear Line Search
spgLineCurvy <- function(m, grad, f, tau,n_vec, y_vec, pos_ind,n ) {
  #initialisation
  # Exit status constants
  EXIT_CONVERGED <- 0
  EXIT_ITERATIONS <- 1

  # Line search parameters
  maxIts_curvy <- 20 # Maximum iterations
  step <- 1 # Initial step size
  alpha <- 1 # Scaling factor for the gradient
  iter_curvy <- 0 # Iteration counter
  gamma <- 1e-4

  while (TRUE) {
    # Evaluate trial point and function value
    mTrial <- m
    mTrial[pos_ind] <- m[pos_ind] - step * alpha * grad
    mNew <-proj_nucnorm(mTrial,  n, tau)
    fNew <- rcpp_val(mNew[pos_ind],  n_vec[pos_ind],  y_vec[pos_ind])
    s <- mNew - m
    gts <- alpha * sum(grad * s[pos_ind])

    # Normalized inner product for descent check
    gtsBar <- gts / (sqrt(sum(grad^2)) * sqrt(sum(s^2)))

    # Check exit conditions
    if (gtsBar < -1e-6) {
      if (fNew < f + gamma * step * gts) {
        stat_curvy <- EXIT_CONVERGED
        break
      }
    }

    if (iter_curvy >= maxIts_curvy) {
      stat_curvy <- EXIT_ITERATIONS
      break
    }

    # Safeguarded step reduction
    step <- step / 2
    iter_curvy <- iter_curvy + 1

    # Adjust scaling factor if step size becomes too small
    sNorm <- max(abs(s))
    if (sNorm <= 1e-6 * max(sqrt(sum(m^2)), 1)) {
      alpha <- alpha / 10
    }
  }

  # Return updated values and status
  return(list(fNew = fNew, mNew = mNew, iter = iter_curvy, step = step, stat = stat_curvy))
}


#============================================================================
# Projection of matrix Z onto the nuclear norm ball with given tau.
proj_nucnorm <- function(M_vec,  n, tau) {
  # Reshape M_vec into n x n matrix
  M_matrix_prep <- matrix(0, nrow = n,ncol =n)
  M_matrix_prep[upper.tri(M_matrix_prep)] <- M_vec
  M_matrix <- M_matrix_prep - t(M_matrix_prep)

  # Compute the Singular Value Decomposition (SVD)
  svd_result <- svd(M_matrix)
  U <- svd_result$u
  S <- svd_result$d
  V <- svd_result$v

  # Project the singular values onto the nuclear norm ball
  s2 <- nuclear_norm_projection(S, tau)

  # Reconstruct the projected matrix
  P <- U %*% diag(s2) %*% t(V)
  P_vec <- P[upper.tri(P)]  # Convert matrix back to vector form

  # Return result
  return(P_vec)
}

#============================================================================
# Helper function for projecting singular values onto nuclear norm ball.
nuclear_norm_projection <- function(singular_values, tau) {
  odd = FALSE#Default;even
  # Check if the nuclear norm (sum of singular values) exceed the tau.
  if (sum(singular_values) <= tau) {
    return(singular_values)  # No need to project if already within tau
  }

  # Handle odd number of singular values by appending a zero
  if (length(singular_values) %% 2 != 0) {
    singular_values <- c(singular_values, 0)
    odd = TRUE
  }

  # Group into pairs
  m <- length(singular_values) / 2
  paired_sigma <- singular_values[seq(1, length(singular_values), by = 2)]

  # Compute cumulative sums for pairs
  S_paired <- cumsum(2 * paired_sigma)

  # Initialize l_star
  l_star <- 1

  # Determine the threshold (lambda) for projection
  for (l in 1:m) {
    lambda_candidate <- (S_paired[l] - tau) / (2 * l)
    if(lambda_candidate <0){
      next
    }
    sigma_next <- ifelse(l < m, paired_sigma[l + 1], 0)

    if (lambda_candidate > sigma_next) {#make sure the terms are 0 after l
      l_star <- l
      break
    }

    # If reached the end without satisfying the condition
    if (l == m) {
     print("warning: Cannnot find lambda")
    }
  }

  # Compute lambda
  lambda <- lambda_candidate
  if(odd){#Remove last input if odd = TRUE.
    singular_values <- singular_values[-length(singular_values)]
  }
  return(pmax(singular_values - lambda, 0))
}

#================================================================================
#g.func: logit function
g.func = function(x){
  # Compute the logit function
  out=1/(1 + exp(-x))

  # Safeguards for extreme values
  out[x> 700]=1
  out[x< -800]=0
  return(out)
}




