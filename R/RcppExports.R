# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpp_grad_val_mat <- function(M_vec, n_vec, w_vec) {
    .Call(`_PCWST_rcpp_grad_val_mat`, M_vec, n_vec, w_vec)
}

rcpp_val <- function(M_vec, n_vec, w_vec) {
    .Call(`_PCWST_rcpp_val`, M_vec, n_vec, w_vec)
}

rcpp_like <- function(W, Mhat, n) {
    .Call(`_PCWST_rcpp_like`, W, Mhat, n)
}

