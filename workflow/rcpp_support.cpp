#include <Rcpp.h>

using namespace Rcpp;

//function to compute test accuracy
// [[Rcpp::export]]
double rcpp_acc(NumericMatrix W, NumericMatrix Mhat, int n) {
  double objectiveValue = 0.0;
  int count = 0;
  // Iterate over all pairs (i, j)
  for (int i = 0; i < (n - 1); i++) {
    for (int j = i + 1; j < n; j++) {
      // Get number of times i won against j and j won against i
      double wij = W(i, j);
      double wji = W(j, i);
      
      if (!(wij == 0 && wji == 0)) {  // Proceed only if both are not zero
        count = count + wij + wji;
        // Compute factor component: parts that exclude intercepts
        // Update likelihood value 
        if(Mhat(i,j) >= 0.5){
          objectiveValue = objectiveValue + wij; 
        }else{
          objectiveValue = objectiveValue + wji; 
        }
      }
    }
  }
  
  // Normalize the output
  objectiveValue /= count;
  
  return objectiveValue;
}

//function to compute likelihood based on estimated probabilities
// [[Rcpp::export]]
double rcpp_like(NumericMatrix W, NumericMatrix Mhat, int n) {
  // Initialize variables
  double objectiveValue = 0.0;
  int count = 0;
  // Iterate over unique pairs
  for (int i = 0; i < (n - 1); i++) {
    for (int j = i + 1; j < n; j++) {
      // Get the observed counts of wins for i over j (wij) and j over i (wji)
      double wij = W(i, j);
      double wji = W(j, i);
      if (!(wij == 0 && wji == 0)) {  // Skip pairs where no matches were observed
        // Update likelihood value 
        count = count + wij + wji ;
        objectiveValue += wij * log(Mhat(i,j)) + wji * log(Mhat(j,i));
      }
    }
  }
  
  // Return the computed likelihood value
  return objectiveValue/count;
}
