//Script containing support functions written in Rcpp.
#include <Rcpp.h>
using namespace Rcpp;

//Function to compute gradient and value of objective function given input vector M_vec.
// [[Rcpp::export]]
List rcpp_grad_val_mat(NumericVector M_vec, NumericVector n_vec, NumericVector w_vec) {
  // Initialize variables
  int M_vec_len  = M_vec.length();
  NumericVector grad(M_vec_len);
  double objectiveValue = 0.0;
  
  // Loop over all elements of M_vec
  for (int j = 0; j < M_vec_len; j++) {
     // Compute the logistic function (gmj) 
      double gmj = 1 / (1 + exp(-M_vec[j]));
      
      // Compute (-)gradient
      grad[j] = w_vec[j] - n_vec[j] * gmj;
      
      //update (-)objective value
      objectiveValue += -(n_vec[j]) *log(1 + exp(-M_vec[j]))  - (n_vec[j] - w_vec[j] ) *M_vec[j];
  }
  // Return the gradient and objective value 
  return List::create(Named("grad") = -grad,Named("objectiveValue") = -objectiveValue);
}

//Function to compute objective value only.
// [[Rcpp::export]]
double rcpp_val(NumericVector M_vec, NumericVector n_vec, NumericVector w_vec) {
  // Initialize variables
  int M_vec_len  = M_vec.length();
  double objectiveValue = 0.0;
  
  // Loop over all elements of M_vec
  for (int j = 0; j < M_vec_len; j++) {
    //update (-)objective value
    objectiveValue += -(n_vec[j]) *log(1 + exp(-M_vec[j]))  - (n_vec[j] - w_vec[j] ) *M_vec[j];
  }
 
 //Return objective value 
  return  -objectiveValue;
}



//function to compute likelihood based on estimated probabilities
// [[Rcpp::export]]
double rcpp_like(NumericMatrix W, NumericMatrix Mhat, int n) {
  // Initialize variables
  double objectiveValue = 0.0;
  
  // Iterate over unique pairs
  for (int i = 0; i < (n - 1); i++) {
    for (int j = i + 1; j < n; j++) {
      // Get the observed counts of wins for i over j (wij) and j over i (wji)
      double wij = W(i, j);
      double wji = W(j, i);
      
      if (!(wij == 0 && wji == 0)) {  // Skip pairs where no matches were observed
        // Update likelihood value 
        objectiveValue += wij * log(Mhat(i,j)) + wji * log(Mhat(j,i));
      }
    }
  }
  
  // Return the computed likelihood value
  return objectiveValue;
}