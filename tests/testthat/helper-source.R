# Load the packages and source the repository's functions so that the tests
# exercise exactly the code used by the analysis scripts.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source(test_path("..", "..", "R", "Functions.R"))
source(test_path("..", "..", "R", "Functions_cSTM_time_indep.R"))

# Base-case parameters of the Sick-Sicker model (same as in
# analysis/cSTM_time_indep.R)
l_params_base <- list(
  r_HD = 0.002, r_HS1 = 0.15, r_S1H = 0.5, r_S1S2 = 0.105,
  hr_S1 = 3, hr_S2 = 10, hr_S1S2_trtB = 0.6,
  c_H = 2000, c_S1 = 4000, c_S2 = 15000, c_D = 0,
  c_trtA = 12000, c_trtB = 13000,
  u_H = 1, u_S1 = 0.75, u_S2 = 0.5, u_D = 0, u_trtA = 0.95,
  n_age_init = 25, n_age_max = 100,
  d_c = 0.03, d_e = 0.03,
  cycle_length = 1
)
