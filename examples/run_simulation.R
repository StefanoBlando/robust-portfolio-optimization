# ========================================================================
# Run Simulation Study for Robust Portfolio Optimization
# ========================================================================

# Load required modules
source("src/estimation.R")
source("src/portfolio.R")
source("src/utils.R")
source("analysis/simulation.R")

# -------------------------------------------------------------------------
# Set Simulation Parameters
# -------------------------------------------------------------------------

# Basic simulation settings
n <- 504                  # Number of observations (~ 2 years of daily data)
p <- 100                  # Number of assets (dimension)
replications <- 10        # Number of simulation replications
contamination_levels <- c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15)  # Contamination percentages to test
heavy_tailed <- TRUE      # Use t-distribution instead of normal
launch_dashboard <- FALSE # Set to TRUE to launch interactive dashboard

# Covariance methods to test
cov_methods <- c(
  "Sample",          # Standard sample covariance (baseline)
  "MCD",             # Minimum Covariance Determinant
  "HybridROBPCA",    # Hybrid MCD + ROBPCA method
  "HybridRobGlasso", # Hybrid MCD + ROBPCA + GLasso regularization
  "Tyler",           # Tyler M-estimator
  "AdaptiveLW",      # Adaptive Ledoit-Wolf shrinkage
  "HFBRE"            # Hybrid Factor-Based Robust Estimator (novel method)
)

# Portfolio optimization methods to test
opt_methods <- c(
  "EqualWeight",     # Equal weight (1/N) portfolio
  "MinVar",          # Minimum variance portfolio
  "MaxSharpe",       # Maximum Sharpe ratio portfolio
  "MinCVaR",         # Minimum Conditional Value-at-Risk portfolio
  "DE",              # Differential Evolution for minimum variance
  "HRP"              # Hierarchical Risk Parity
)

# -------------------------------------------------------------------------
# Run Simulation
# -------------------------------------------------------------------------

cat("Starting robust portfolio optimization simulation study...\n")
cat(paste("Testing", length(cov_methods), "covariance methods and", 
          length(opt_methods), "optimization methods\n"))
cat(paste("Using", n, "observations,", p, "assets, and", replications, "replications\n"))
cat(paste("Contamination levels:", paste(contamination_levels, collapse = ", "), "\n"))
cat(paste("Distribution:", ifelse(heavy_tailed, "Heavy-tailed (t)", "Normal"), "\n"))

# Main simulation
results <- run_robust_portfolio_analysis(
  n = n,
  p = p, 
  replications = replications,
  contamination_levels = contamination_levels,
  compare_contamination = TRUE,  # Compare balanced vs unbalanced contamination
  heavy_tailed = heavy_tailed,
  launch_dashboard = launch_dashboard,
  cov_methods = cov_methods,
  opt_methods = opt_methods
)

# -------------------------------------------------------------------------
# Save Results and Generate Report
# -------------------------------------------------------------------------

# Save results and generate report
output_file <- paste0("sim_study_results_", format(Sys.Date(), "%Y%m%d"))
save_output <- save_sim_study_results(results, filename = output_file)

cat("\nSimulation study complete.\n")
cat("Results saved to:", save_output$saved_file, "\n")
cat("Report saved to:", save_output$report_file, "\n")

# -------------------------------------------------------------------------
# (Optional) Launch Dashboard
# -------------------------------------------------------------------------

# If dashboard wasn't launched during simulation, you can launch it afterward with:
if(!launch_dashboard) {
  cat("\nTo launch the interactive dashboard, run:\n")
  cat("create_portfolio_dashboard(results)\n")
}
