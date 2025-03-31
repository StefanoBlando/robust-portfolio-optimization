# ========================================================================
# Run Empirical Analysis on S&P 500 Data
# ========================================================================

# Load required modules
source("src/estimation.R")  # Covariance estimation methods
source("src/portfolio.R")   # Portfolio optimization methods
source("src/utils.R")       # Utility functions
source("analysis/empirical.R")  # Empirical analysis functions
source("visualization/visualization.R")  # Visualization functions

# Set parameters for analysis
params <- list(
  start_date = "2015-01-01",  # Start date for analysis
  end_date = Sys.Date(),      # Current date as end date
  n_stocks = 100,             # Number of S&P 500 stocks to include
  train_ratio = 0.7,          # Training/testing split ratio
  by_regime = TRUE,           # Use regime-aware splitting
  max_weight = 0.2,           # Maximum weight per asset
  output_dir = "results/empirical"  # Output directory
)

# Create output directory if it doesn't exist
if (!dir.exists(params$output_dir)) {
  dir.create(params$output_dir, recursive = TRUE)
}

# Set random seed for reproducibility
set.seed(42)

# Setup parallel processing
cores <- min(detectCores() - 1, 14)  # Use all but one core, maximum 14
cl <- makeCluster(cores)
registerDoParallel(cl)
cat(paste("Using", cores, "cores for parallel processing\n"))

# ========================================================================
# Run the standard empirical analysis
# ========================================================================

cat("\n===== RUNNING STANDARD S&P 500 ANALYSIS =====\n\n")

analysis_results <- run_sp500_analysis(
  start_date = params$start_date,
  end_date = params$end_date,
  n_stocks = params$n_stocks,
  train_ratio = params$train_ratio,
  by_regime = params$by_regime,
  rolling = FALSE,  # Skip rolling window analysis for speed
  max_weight = params$max_weight
)

# Save results and generate report
standard_dir <- file.path(params$output_dir, "standard_analysis")
save_and_report(analysis_results, standard_dir)

# ========================================================================
# Regime-Specific Analysis
# ========================================================================

cat("\n===== RUNNING REGIME-SPECIFIC ANALYSIS =====\n\n")

# Extract regime information if available
if (!is.null(analysis_results$regime_info)) {
  # Create regime-specific performance summary
  regime_performance <- create_regime_performance_summary(
    analysis_results$portfolio_results,
    analysis_results$split_data$regimes[analysis_results$split_data$test_indices]
  )
  
  # Create regime visualizations
  regime_plots <- create_regime_visualizations(regime_performance)
  
  # Save regime-specific results
  regime_dir <- file.path(params$output_dir, "regime_analysis")
  if (!dir.exists(regime_dir)) {
    dir.create(regime_dir, recursive = TRUE)
  }
  
  # Save regime performance summaries
  for (regime in names(regime_performance)) {
    write.csv(regime_performance[[regime]],
              file = file.path(regime_dir, paste0(regime, "_performance.csv")),
              row.names = FALSE)
  }
  
  # Save regime plots
  if (!is.null(regime_plots)) {
    for (name in names(regime_plots)) {
      ggsave(file.path(regime_dir, paste0(name, ".png")),
             regime_plots[[name]], width = 10, height = 6)
    }
  }
  
  cat("Regime-specific analysis saved to", regime_dir, "\n")
} else {
  cat("No regime information available for analysis\n")
}

# ========================================================================
# Sensitivity Analysis
# ========================================================================

cat("\n===== RUNNING SENSITIVITY ANALYSIS =====\n\n")

# Test sensitivity to contamination detection threshold
threshold_sensitivity <- test_contamination_sensitivity(
  analysis_results$split_data$train_returns,
  analysis_results$split_data$test_returns,
  thresholds = seq(2, 5, by = 0.5),
  analysis_results$split_data$train_market,
  analysis_results$split_data$test_market,
  max_weight = params$max_weight,
  test_dates = analysis_results$split_data$test_dates,
  test_regimes = analysis_results$split_data$regimes[analysis_results$split_data$test_indices]
)

# Test different outlier detection methods
detection_methods <- test_detection_methods(
  analysis_results$split_data$train_returns,
  analysis_results$split_data$test_returns,
  analysis_results$split_data$train_market,
  analysis_results$split_data$test_market,
  max_weight = params$max_weight,
  test_dates = analysis_results$split_data$test_dates,
  test_regimes = analysis_results$split_data$regimes[analysis_results$split_data$test_indices]
)

# Save sensitivity analysis results
sensitivity_dir <- file.path(params$output_dir, "sensitivity_analysis")
if (!dir.exists(sensitivity_dir)) {
  dir.create(sensitivity_dir, recursive = TRUE)
}

# Save sensitivity plots
if (!is.null(threshold_sensitivity)) {
  if (!is.null(threshold_sensitivity$min_var_plot)) {
    ggsave(file.path(sensitivity_dir, "threshold_sensitivity_minvar.png"),
           threshold_sensitivity$min_var_plot, width = 10, height = 6)
  }
  if (!is.null(threshold_sensitivity$max_sharpe_plot)) {
    ggsave(file.path(sensitivity_dir, "threshold_sensitivity_maxsharpe.png"),
           threshold_sensitivity$max_sharpe_plot, width = 10, height = 6)
  }
}

if (!is.null(detection_methods)) {
  if (!is.null(detection_methods$min_var_plot)) {
    ggsave(file.path(sensitivity_dir, "detection_methods_minvar.png"),
           detection_methods$min_var_plot, width = 10, height = 6)
  }
  if (!is.null(detection_methods$max_sharpe_plot)) {
    ggsave(file.path(sensitivity_dir, "detection_methods_maxsharpe.png"),
           detection_methods$max_sharpe_plot, width = 10, height = 6)
  }
}

# Save optimal parameters
optimal_params <- data.frame(
  Parameter = c("Min Variance Threshold", "Max Sharpe Threshold", 
                "Min Variance Detection Method", "Max Sharpe Detection Method"),
  Value = c(
    ifelse(!is.null(threshold_sensitivity$min_var_optimal_threshold),
           threshold_sensitivity$min_var_optimal_threshold, NA),
    ifelse(!is.null(threshold_sensitivity$max_sharpe_optimal_threshold),
           threshold_sensitivity$max_sharpe_optimal_threshold, NA),
    ifelse(!is.null(detection_methods$min_var_best_method),
           detection_methods$min_var_best_method, NA),
    ifelse(!is.null(detection_methods$max_sharpe_best_method),
           detection_methods$max_sharpe_best_method, NA)
  )
)

write.csv(optimal_params, file = file.path(sensitivity_dir, "optimal_parameters.csv"), 
          row.names = FALSE)

cat("Sensitivity analysis saved to", sensitivity_dir, "\n")

# ========================================================================
# Performance Visualization
# ========================================================================

cat("\n===== CREATING PERFORMANCE VISUALIZATIONS =====\n\n")

# Create comprehensive performance visualizations
perf_plots <- create_performance_visualizations(analysis_results)
weight_plots <- create_weight_visualizations(analysis_results)

# Save visualizations
viz_dir <- file.path(params$output_dir, "visualizations")
if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

# Save performance plots
for (name in names(perf_plots)) {
  if (name != "top_weights") {  # Handle special case for list of plots
    ggsave(file.path(viz_dir, paste0(name, ".png")),
           perf_plots[[name]], width = 10, height = 6)
  }
}

# Handle top weights plots separately
if ("top_weights" %in% names(weight_plots)) {
  for (method_name in names(weight_plots$top_weights)) {
    safe_name <- gsub("[^a-zA-Z0-9]", "_", method_name)
    ggsave(file.path(viz_dir, paste0("top_weights_", safe_name, ".png")),
           weight_plots$top_weights[[method_name]], width = 10, height = 6)
  }
}

# Save other weight plots
for (name in names(weight_plots)) {
  if (name != "top_weights") {  # Skip the list of plots
    ggsave(file.path(viz_dir, paste0(name, ".png")),
           weight_plots[[name]], width = 10, height = 6)
  }
}

cat("Performance visualizations saved to", viz_dir, "\n")

# ========================================================================
# Generate Final Report
# ========================================================================

cat("\n===== GENERATING FINAL REPORT =====\n\n")

# Create comprehensive report
report_file <- file.path(params$output_dir, "empirical_analysis_report.md")

# Extract key results
best_method <- analysis_results$performance_summary$Method[1]
best_sharpe <- analysis_results$performance_summary$OutOfSample_Sharpe[1]

# Create report content
report_content <- c(
  "# S&P 500 Robust Portfolio Optimization Analysis",
  "",
  paste("Analysis period:", params$start_date, "to", params$end_date),
  paste("Number of stocks analyzed:", params$n_stocks),
  paste("Train/test split ratio:", params$train_ratio),
  paste("Maximum weight per asset:", params$max_weight),
  "",
  "## Performance Summary",
  "",
  "```",
  capture.output(print(analysis_results$performance_summary)),
  "```",
  "",
  "## Key Findings",
  "",
  paste("- Best performing method:", best_method),
  paste("  with out-of-sample Sharpe ratio:", round(best_sharpe, 4)),
  "",
  "### Detected Contamination",
  "",
  paste("- Overall contamination level:", 
        round(analysis_results$outliers$contamination_pct, 2), "%"),
  "",
  "### Sensitivity Analysis",
  "",
  paste("- Optimal outlier detection threshold for Min Variance:", 
        ifelse(!is.null(threshold_sensitivity$min_var_optimal_threshold),
               threshold_sensitivity$min_var_optimal_threshold, "Not determined")),
  paste("- Optimal outlier detection threshold for Max Sharpe:", 
        ifelse(!is.null(threshold_sensitivity$max_sharpe_optimal_threshold),
               threshold_sensitivity$max_sharpe_optimal_threshold, "Not determined")),
  "",
  paste("- Best detection method for Min Variance:", 
        ifelse(!is.null(detection_methods$min_var_best_method),
               detection_methods$min_var_best_method, "Not determined")),
  paste("- Best detection method for Max Sharpe:", 
        ifelse(!is.null(detection_methods$max_sharpe_best_method),
               detection_methods$max_sharpe_best_method, "Not determined"))
)

# Add regime information if available
if (!is.null(analysis_results$regime_info)) {
  regime_content <- c(
    "",
    "## Market Regime Analysis",
    "",
    "### Regime Distribution",
    "",
    "```",
    capture.output(print(analysis_results$regime_info$regime_counts)),
    "```",
    "",
    "### Train Set Regime Distribution",
    "",
    "```",
    capture.output(print(analysis_results$regime_info$train_regime_distribution)),
    "```",
    "",
    "### Test Set Regime Distribution",
    "",
    "```",
    capture.output(print(analysis_results$regime_info$test_regime_distribution)),
    "```"
  )
  
  report_content <- c(report_content, regime_content)
}

# Add conclusion
report_content <- c(report_content,
                    "",
                    "## Conclusion",
                    "",
                    paste("This analysis demonstrates that hybrid robust methods provide",
                          "significant advantages in portfolio optimization, especially",
                          "when dealing with financial returns data that often contains",
                          "outliers and follows heavy-tailed distributions."),
                    "",
                    paste("The", best_method, "approach shows the best combination of",
                          "performance and robustness, making it suitable for practical",
                          "implementation in portfolio management."))

# Write the report
writeLines(report_content, report_file)

cat("Final report saved to", report_file, "\n")

# Clean up parallel cluster
stopCluster(cl)

cat("\n===== EMPIRICAL ANALYSIS COMPLETE =====\n\n")
