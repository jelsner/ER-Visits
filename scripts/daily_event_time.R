# Current entry point: confirmed 0--5 suppression, observed counts >=6.
# Default central mean is 3.8 visits, with sensitivity from zero to five.
message("Using corrected suppression models for Figures 1 and S1.")
source("scripts/regenerate_suppression_figures.R")
source("scripts/plot_suppression_figures.R")
source("scripts/write_suppression_supplement.R")
