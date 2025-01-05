# PCWST
Pairwise comparisons without stochastic transitivity, as proposed by Lee et al. (2025+)

To install this package in R, run the following commands:  

```R
library(devtools) 
install_github("Arthurlee51/PCWST")
```

## Usage 

pcwst_est(): Performs parameter estimation given a pairwise comparison matrix.

pcwst_est_tune(): Decide the scale parameter (Cn) using training and validation datasets.


## Brief description
The PCWST package implements the estimator developed by Lee et al. (2025+), which analyses models pairwise comparison data without the traditional stochastic transitivity assumption. The function pcwst_est() performs parameter estimation, given the matrix of comparison outcome Y, and the scale parameter (equivalent to Cn in the paper) that controls the nuclear norm of the parameter matrix. When the scale parameter is unknown and needs to be estimated, we can apply the function pcwst_est_tune() to choose scale from a range of possible values specified by the user, given traingin data and validation data.

```R
?pcwst_est
?pcwst_est_tune
```

## Reference 
Sze Ming Lee and Yunxiao Chen. **Pairwise Comparisons without Stochastic Transitivity: Model, Theory and Applications**. 2025+. (manuscript)


