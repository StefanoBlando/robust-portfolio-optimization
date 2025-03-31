# ========================================================================
# Run Stress Test Analysis on Portfolio Methods
# ========================================================================

# Load required modules
source("src/estimation.R")  # Covariance estimation methods
source("src/portfolio.R")   # Portfolio optimization methods
source("src/utils.R")       # Utility functions
source("analysis/empirical.R")  # Empirical analysis functions
source("analysis/stress_testing.R")  # Stress testing functions
source("visualization/visualization.R")  # Visualization functions

# Set parameters for analysis
params <- list(
  start_date = "2015-01-01",  # Start date for analysis
  end_date = Sys.Date(),      # Current date as end date
  n_stocks = 100,             # Number of S&P 500 stocks to include
  train_ratio = 0.7,          # Training/testing split ratio
  by_regime = TRUE,           # Use regime-aware splitting
  max_weight = 0.2,           # Maximum weight per asset
  output_dir = "results/stress_test",  # Output directory
  
  # Stress test specific parameters
  run_multi_scenario = TRUE,  # Whether to run multiple scenarios
  scenarios = list(
    moderate_tail = list(
      method = "tail_contamination",
      intensity = 0.2,
      target_pct = 0.1,
      description = "Moderate Tail Events (10% fatter tails)"
    ),
    severe_tail = list(
      method = "tail_contamination",
      intensity = 0.4,
      target_pct = 0.15,
      description = "Severe Tail Events (15% much fatter tails)"
    ),
    volatility_spike = list(
      method = "amplify_volatility",
      intensity = 0.3,
      target_pct = 0.15,
      description = "Volatility Spike (15% of periods with amplified volatility)"
    ),
    scattered_outliers = list(
      method = "add_outliers",
      intensity = 0.25,
      target_pct = 0.05,
      description = "Scattered Outliers (5% random extreme returns)"
    ),
    extreme_regime_stress = list(
      method = "regime_based",
      intensity = 0.5,
      target_pct = 0.2,
      description = "Extreme Regime-Based Stress (20% of data regime-specific amplification)"
    )
  )
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
# Step 1: Run the standard analysis (if not already available)
# ========================================================================

cat("\n===== PHASE 1: STANDARD MARKET ANALYSIS =====\n\n")

# Check if standard analysis results exist
standard_results_file <- file.path(params$output_dir, "standard_analysis_results.RData")

if (file.exists(standard_results_file)) {
  # Load existing results
  cat("Loading existing standard analysis results...\n")
  load(standard_results_file)
} else {
  # Run standard analysis
  cat("Running standard S&P 500 analysis...\n")
  analysis_results <- run_sp500_analysis(
    start_date = params$start_date,
    end_date = params$end_date,
    n_stocks = params$n_stocks,
    train_ratio = params$train_ratio,
    by_regime = params$by_regime,
    rolling = FALSE,  # Skip rolling window analysis for speed
    max_weight = params$max_weight
  )
  
  # Save standard results
  save(analysis_results, file = standard_results_file)
  
  # Create summary report
  standard_dir <- file.path(params$output_dir, "standard_analysis")
  save_and_report(analysis_results, standard_dir)
}

# ========================================================================
# Step 2: Run specific stress test or multi-scenario analysis
# ========================================================================

if (params$run_multi_scenario) {
  # ========================================================================
  # Multi-Scenario Stress Testing
  # ========================================================================
  
  cat("\n===== PHASE 2: MULTI-SCENARIO STRESS TESTING =====\n\n")
  
  # Run multiple stress scenarios
  stress_results <- list()
  
  for (scenario_name in names(params$scenarios)) {
    scenario <- params$scenarios[[scenario_name]]
    
    cat("\n----- Testing Scenario:", scenario$description, "-----\n\n")
    
    # Run stress test for this scenario
    stress_results[[scenario_name]] <- run_stress_test_analysis(
      results = analysis_results,
      stress_method = scenario$method,
      stress_intensity = scenario$intensity,
      target_pct = scenario$target_pct,
      preserve_correlation = TRUE,
      by_regime = params$by_regime,
      output_dir = file.path(params$output_dir, paste0("stress_", scenario_name))
    )
  }
  
  # Create comparative analysis across all scenarios
  create_multi_scenario_report(analysis_results, stress_results, params$scenarios, params$output_dir)
  
  # Create visuals comparing method performance across scenarios
  create_multi_scenario_visuals(analysis_results, stress_results, params$scenarios, params$output_dir)
  
} else {
  # ========================================================================
  # Single Stress Test Analysis
  # ========================================================================
  
  cat("\n===== PHASE 2: STRESS TEST ANALYSIS =====\n\n")
  
  # Default to tail contamination if specific scenario not selected
  scenario <- params$scenarios$moderate_tail
  
  # Run stress test analysis on the results
  stress_results <- run_stress_test_analysis(
    results = analysis_results,
    stress_method = scenario$method,
    stress_intensity = scenario$intensity,
    target_pct = scenario$target_pct,
    preserve_correlation = TRUE,
    by_regime = params$by_regime,
    output_dir = file.path(params$output_dir, "stress_analysis")
  )
}

# ========================================================================
# Step 3: Additional Stress Test Experiments
# ========================================================================

cat("\n===== PHASE 3: CONTAMINATION TYPE COMPARISON =====\n\n")

# Compare balanced vs unbalanced contamination
# Choose a representative contamination level
balanced_test_results <- run_outliers_stress_test(
  analysis_results = analysis_results,
  intensity = 0.25,
  target_pct = 0.05,
  preserve_correlation = TRUE,
  by_regime = params$by_regime
)

# Compare regime-specific stress
regime_test_results <- run_regime_stress_test(
  analysis_results = analysis_results,
  intensity = 0.5,
  target_pct = 0.2,
  preserve_correlation = TRUE,
  by_regime = params$by_regime
)

# ========================================================================
# Step 4: Create Comparative Analysis
# ========================================================================

cat("\n===== CREATING COMPARATIVE ANALYSIS =====\n\n")

# Prepare comparison data
standard_perf <- analysis_results$performance_summary

# Extract stress test performance results
if (params$run_multi_scenario) {
  # Multi-scenario mode - use first scenario as representative
  first_scenario <- names(stress_results)[1]
  moderate_stress_perf <- stress_results[[first_scenario]]$performance_summary
  # Use severe scenario for high stress
  high_stress_perf <- stress_results[["severe_tail"]]$performance_summary
} else {
  # Single stress test mode
  moderate_stress_perf <- stress_results$performance_summary
  # We don't have high stress in single mode
  high_stress_perf <- NULL
}

# Prepare comparison data for visualization
comparison_data <- merge(
  standard_perf[, c("Method", "OutOfSample_Sharpe")],
  moderate_stress_perf[, c("Method", "OutOfSample_Sharpe")],
  by = "Method",
  suffixes = c("_Standard", "_Stress")
)

# Calculate performance change
comparison_data$Absolute_Change <- comparison_data$OutOfSample_Sharpe_Stress - 
  comparison_data$OutOfSample_Sharpe_Standard
comparison_data$Percent_Change <- (comparison_data$OutOfSample_Sharpe_Stress / 
                                     comparison_data$OutOfSample_Sharpe_Standard - 1) * 100

# Order by percent change (least degradation first)
comparison_data <- comparison_data[order(-comparison_data$Percent_Change), ]

# Create performance comparison visualizations
viz_dir <- file.path(params$output_dir, "comparative_visualization")
if (!dir.exists(viz_dir)) {
  dir.create(viz_dir, recursive = TRUE)
}

# Performance change bar chart
p1 <- ggplot(comparison_data, aes(x = reorder(Method, Percent_Change), y = Percent_Change)) +
  geom_bar(stat = "identity", aes(fill = Percent_Change > 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("firebrick", "steelblue"), name = "Improved") +
  labs(title = "Performance Change Under Stress",
       x = "", y = "Percent Change in Sharpe Ratio (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(viz_dir, "performance_change.png"), p1, width = 12, height = 6)

# Scatter plot of standard vs. stress performance
p2 <- ggplot(comparison_data, aes(x = OutOfSample_Sharpe_Standard, y = OutOfSample_Sharpe_Stress)) +
  geom_point(aes(color = Percent_Change > 0), size = 3) +
  geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("firebrick", "steelblue"), name = "Improved") +
  labs(title = "Performance: Standard vs. Stress Conditions",
       x = "Standard Sharpe Ratio", y = "Stress Sharpe Ratio") +
  theme_minimal()

ggsave(file.path(viz_dir, "performance_scatter.png"), p2, width = 10, height = 8)

# If we have high stress results, create three-way comparison
if (!is.null(high_stress_perf)) {
  # Create three-way comparison
  three_way <- merge(
    comparison_data[, c("Method", "OutOfSample_Sharpe_Standard", "OutOfSample_Sharpe_Stress")],
    high_stress_perf[, c("Method", "OutOfSample_Sharpe")],
    by = "Method"
  )
  names(three_way)[4] <- "OutOfSample_Sharpe_HighStress"
  
  # Reshape for plotting
  three_way_long <- reshape2::melt(
    three_way,
    id.vars = "Method",
    measure.vars = c("OutOfSample_Sharpe_Standard", "OutOfSample_Sharpe_Stress", 
                     "OutOfSample_Sharpe_HighStress"),
    variable.name = "Condition",
    value.name = "Sharpe"
  )
  
  # Clean up condition names
  three_way_long$Condition <- gsub("OutOfSample_Sharpe_", "", three_way_long$Condition)
  three_way_long$Condition <- factor(three_way_long$Condition, 
                                     levels = c("Standard", "Stress", "HighStress"))
  
  # Bar chart for three-way comparison
  p3 <- ggplot(three_way_long, aes(x = Method, y = Sharpe, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Sharpe Ratio Under Different Market Conditions",
         x = "", y = "Sharpe Ratio") +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(viz_dir, "three_way_comparison.png"), p3, width = 14, height = 7)
}

# ========================================================================
# Step 5: Generate Comprehensive Stress Test Report
# ========================================================================

cat("\n===== GENERATING COMPREHENSIVE REPORT =====\n\n")

# Create comprehensive report
report_file <- file.path(params$output_dir, "stress_test_report.md")

# Get best method from standard analysis
standard_best_method <- standard_perf$Method[1]
standard_best_sharpe <- standard_perf$OutOfSample_Sharpe[1]

# Get best method from stress analysis
stress_best_method <- moderate_stress_perf$Method[1]
stress_best_sharpe <- moderate_stress_perf$OutOfSample_Sharpe[1]

# Prepare report content
report_content <- c(
  "# Portfolio Optimization Stress Test Analysis",
  "",
  paste("Analysis period:", params$start_date, "to", params$end_date),
  paste("Number of stocks:", params$n_stocks),
  "",
  "## Stress Test Parameters",
  ""
)

# Add information about stress scenarios
if (params$run_multi_scenario) {
  # Add multi-scenario information
  report_content <- c(report_content, "### Multiple Stress Scenarios Tested")
  
  for (scenario_name in names(params$scenarios)) {
    scenario <- params$scenarios[[scenario_name]]
    report_content <- c(report_content,
                        paste("#### Scenario:", scenario_name),
                        paste("- Description:", scenario$description),
                        paste("- Method:", scenario$method),
                        paste("- Intensity:", scenario$intensity),
                        paste("- Target percentage:", scenario$target_pct),
                        "")
  }
} else {
  # Add single scenario information
  scenario <- params$scenarios$moderate_tail  # Default scenario
  report_content <- c(report_content,
                      "### Stress Scenario Parameters",
                      paste("- Method:", scenario$method),
                      paste("- Description:", scenario$description),
                      paste("- Intensity:", scenario$intensity),
                      paste("- Target percentage:", scenario$target_pct),
                      "")
}

# Add performance comparison
report_content <- c(report_content,
                    "",
                    "## Performance Comparison",
                    "",
                    "### Best Performing Methods",
                    "",
                    paste("- Standard conditions best method:", standard_best_method),
                    paste("  with Sharpe ratio:", round(standard_best_sharpe, 4)),
                    paste("- Stress conditions best method:", stress_best_method),
                    paste("  with Sharpe ratio:", round(stress_best_sharpe, 4)),
                    "",
                    "### Performance Degradation Under Stress",
                    "",
                    "```",
                    capture.output(print(head(comparison_data, 10))),
                    "```",
                    "",
                    "### Most Robust Methods (Least Performance Degradation)",
                    "")

# Add the top 5 most robust methods
for (i in 1:min(5, nrow(comparison_data))) {
  method <- comparison_data$Method[i]
  pct_change <- comparison_data$Percent_Change[i]
  change_text <- ifelse(pct_change >= 0,
                        paste("improved by", round(pct_change, 2), "%"),
                        paste("degraded by", round(-pct_change, 2), "%"))
  
  report_content <- c(report_content,
                      paste(i, ". ", method, " - ", change_text, sep = ""))
}

# Add insights about robust methods
report_content <- c(report_content,
                    "",
                    "## Key Insights",
                    "",
                    "### Characteristics of Robust Methods",
                    "")

# Identify which covariance methods were most robust
cov_methods <- unique(sapply(strsplit(comparison_data$Method, "_"), function(x) x[1]))
cov_methods <- cov_methods[cov_methods != "EqualWeights"]

avg_change_by_cov <- aggregate(
  Percent_Change ~ CovMethod, 
  data = data.frame(
    Percent_Change = comparison_data$Percent_Change,
    CovMethod = sapply(strsplit(comparison_data$Method, "_"), function(x) x[1])
  ),
  FUN = mean
)

avg_change_by_cov <- avg_change_by_cov[order(-avg_change_by_cov$Percent_Change), ]

report_content <- c(report_content,
                    "#### Most Robust Covariance Estimation Methods:")

for (i in 1:nrow(avg_change_by_cov)) {
  cov_method <- avg_change_by_cov$CovMethod[
    i]) {
      avg_change <- avg_change_by_cov$Percent_Change[i]
      change_text <- ifelse(avg_change >= 0, 
                            paste("improved by", round(avg_change, 2), "%"),
                            paste("degraded by", round(-avg_change, 2), "%"))
      
      report_content <- c(report_content,
                          paste("  ", i, ". ", cov_method, " (", change_text, ")", sep=""))
    }

# Compare optimization methods
report_content <- c(report_content,
                    "",
                    "#### Optimization Method Robustness:")

# Compare Min Variance vs Max Sharpe
minvar_methods <- comparison_data$Method[grepl("MinVar", comparison_data$Method)]
maxsharpe_methods <- comparison_data$Method[grepl("MaxSharpe", comparison_data$Method)]

if(length(minvar_methods) > 0 && length(maxsharpe_methods) > 0) {
  minvar_avg_change <- mean(comparison_data$Percent_Change[comparison_data$Method %in% minvar_methods])
  maxsharpe_avg_change <- mean(comparison_data$Percent_Change[comparison_data$Method %in% maxsharpe_methods])
  
  minvar_text <- ifelse(minvar_avg_change >= 0, 
                        paste("improved by", round(minvar_avg_change, 2), "%"),
                        paste("degraded by", round(-minvar_avg_change, 2), "%"))
  
  maxsharpe_text <- ifelse(maxsharpe_avg_change >= 0, 
                           paste("improved by", round(maxsharpe_avg_change, 2), "%"),
                           paste("degraded by", round(-maxsharpe_avg_change, 2), "%"))
  
  report_content <- c(report_content,
                      paste("- Minimum Variance methods on average", minvar_text),
                      paste("- Maximum Sharpe methods on average", maxsharpe_text))
  
  if(minvar_avg_change > maxsharpe_avg_change) {
    report_content <- c(report_content,
                        "- Minimum Variance optimization appears more robust to market stress")
  } else {
    report_content <- c(report_content,
                        "- Maximum Sharpe optimization appears more robust to market stress")
  }
}

# Add conclusion
report_content <- c(report_content,
                    "",
                    "## Conclusion",
                    "",
                    paste("This analysis confirms that", avg_change_by_cov$CovMethod[1], "and", 
                          ifelse(nrow(avg_change_by_cov) > 1, avg_change_by_cov$CovMethod[2], "other hybrid methods"), 
                          "are the most robust covariance estimation",
                          "methods under market stress conditions. These methods",
                          "effectively handle the increased contamination present in",
                          "financial returns during periods of market turbulence and",
                          "maintain more stable portfolio performance."),
                    "",
                    paste("For practical implementation, the", stress_best_method, 
                          "approach shows the best combination of performance and",
                          "robustness, making it suitable for real-world portfolio",
                          "management in both normal and stressed market conditions."))

# Write the report
writeLines(report_content, report_file)

cat("Comprehensive stress test report saved to", report_file, "\n")

# Clean up parallel cluster
stopCluster(cl)

cat("\n===== STRESS TEST ANALYSIS COMPLETE =====\n\n")
