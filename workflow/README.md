# Workflow

This directory contains the code and workflow descriptions for reproducing the numerical results presented in the article.

---

## Simulation Results
- The main simulation script is `simulation.R`, which generates the estimated parameters, runtime, and Frobenius loss for each setting.
- The outputs from this script are saved as `.Rdata` files.

After generating the `.Rdata` files, the following scripts can be used to produce the figures reported in the manuscript:

- **Figure 1**  
  - Run `figure1.R`.
  - 
- **Figure 2**  
  - Run `figure2.R`. 


## Real Data Analysis

## Real Data Analysis

The real data analysis for both the ATP and HotS datasets can be reproduced by running `Data_analysis.R`. This script uses the processed datasets `ATPdata.Rdata` and `HotS.Rdata`.

The raw datasets are publicly available from the following sources:

- ATP tennis match data: http://www.tennis-data.co.uk  
- StarCraft II match data: https://www.kaggle.com/datasets/alimbekovkz/starcraft-ii-matches-history
