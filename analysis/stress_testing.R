# ========================================================================
# Stress Testing Framework for Robust Portfolio Optimization
# ========================================================================

# Add controlled contamination to financial data to simulate market stress
add_market_stress <- function(returns, dates, 
                              method = c("amplify_volatility", "add_outliers", "regime_based", "tail_contamination"),
                              intensity = 0.2,  # Controls overall stress level
                              target_pct = 0.1, # Percentage of data to contaminate
                              preserve_correlation = TRUE, # Whether to preserve correlation structure
                              seed = 42) {
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Basic validation
  method <- match.arg(method)
  n <- nrow(returns)
  p <- ncol(returns)
  
  # Initialize contaminated returns
  stressed_returns <- returns
  
  # Create mask to track which values are contaminated
  contamination_mask <- matrix(1, nrow = n, ncol = p)
  
  cat("Applying", method, "stress method with intensity", intensity, "\n")
  
  if(method == "amplify_volatility") {
    # This method increases the volatility of the most volatile periods
    
    # Calculate row volatility
    row_vol <- apply(returns, 1, function(x) sd(x, na.rm = TRUE))
    
    # Identify most volatile periods (top target_pct percentile)
    vol_threshold <- quantile(row_vol, 1 - target_pct)
    high_vol_rows <- which(row_vol > vol_threshold)
    
    # Amplify returns in those periods
    for(i in high_vol_rows) {
      stress_factor <- 1 + intensity * (row_vol[i] / max(row_vol))
      stressed_returns[i, ] <- returns[i, ] * stress_factor
      contamination_mask[i, ] <- 0  # Mark as contaminated
    }
    
    cat("Amplified volatility in", length(high_vol_rows), "periods (", 
        round(length(high_vol_rows)/n*100, 1), "% of data)\n")
    
  } else if(method == "add_outliers") {
    # This method adds synthetic outliers at random positions
    
    # Calculate each asset's volatility
    asset_vol <- apply(returns, 2, sd, na.rm = TRUE)
    
    # Determine number of outliers to add
    n_outliers <- round(n * p * target_pct)
    
    # Generate random positions for outliers
    outlier_rows <- sample(1:n, n_outliers, replace = TRUE)
    outlier_cols <- sample(1:p, n_outliers, replace = TRUE)
    
    # Add outliers
    for(i in 1:n_outliers) {
      row <- outlier_rows[i]
      col <- outlier_cols[i]
      
      # Generate outlier magnitude (scaled by asset volatility and intensity)
      direction <- sample(c(-1, 1), 1)
      magnitude <- asset_vol[col] * (3 + intensity * 7) * direction
      
      stressed_returns[row, col] <- returns[row, col] + magnitude
      contamination_mask[row, col] <- 0  # Mark as contaminated
    }
    
    cat("Added", n_outliers, "outliers (", 
        round(n_outliers/(n*p)*100, 1), "% of data points)\n")
    
  } else if(method == "regime_based") {
    # This method applies different stress patterns based on market regimes
    
    # Determine simple market regimes based on market-wide returns
    market_returns <- rowMeans(returns)
    market_vol <- sd(market_returns) 
    
    # Define regimes
    up_market <- market_returns > market_vol/2
    down_market <- market_returns < -market_vol/2
    sideways_market <- !up_market & !down_market
    
    # Number of days to contaminate
    n_days <- round(n * target_pct)
    
    # Select days to stress with higher probability in down markets
    down_prob <- 0.6   # 60% probability for down markets
    up_prob <- 0.1     # 10% probability for up markets
    side_prob <- 0.3   # 30% probability for sideways markets
    
    # Create sampling weights
    weights <- rep(side_prob, n)
    weights[up_market] <- up_prob
    weights[down_market] <- down_prob
    weights <- weights / sum(weights) * n
    
    # Sample days
    stress_days <- sample(1:n, n_days, prob = weights)
    
    # Apply different stress patterns based on the regime
    for(day in stress_days) {
      if(down_market[day]) {
        # In down markets, amplify negative returns
        neg_assets <- returns[day, ] < 0
        if(any(neg_assets)) {
          stressed_returns[day, neg_assets] <- returns[day, neg_assets] * (1 + intensity * 2)
          contamination_mask[day, neg_assets] <- 0
        }
      } else if(up_market[day]) {
        # In up markets, introduce some reversals
        reversal_factor <- -intensity * 1.5
        stressed_returns[day, ] <- returns[day, ] * (1 + reversal_factor)
        contamination_mask[day, ] <- 0
      } else {
        # In sideways markets, increase volatility
        vol_increase <- runif(p, 1, 1 + intensity * 3)
        stressed_returns[day, ] <- returns[day, ] * vol_increase
        contamination_mask[day, ] <- 0
      }
    }
    
    cat("Applied regime-based stress to", length(stress_days), "days (", 
        round(length(stress_days)/n*100, 1), "% of data)\n")
    
  } else if(method == "tail_contamination") {
    # This method focuses on making fat tails fatter - contaminating the extremes
    
    # Calculate extreme returns for each asset (e.g., 5% most extreme)
    extreme_threshold <- quantile(abs(returns), 1 - target_pct/2)
    extreme_points <- abs(returns) > extreme_threshold
    
    # Make extremes more extreme
    extreme_scale <- 1 + intensity * 2
    
    # Apply contamination
    stressed_returns[extreme_points] <- returns[extreme_points] * extreme_scale
    contamination_mask[extreme_points] <- 0
    
    cat("Applied tail contamination to", sum(extreme_points), "data points (", 
        round(sum(extreme_points)/(n*p)*100, 1), "% of data)\n")
    
    # If requested, preserve correlation structure while making tails fatter
    if(preserve_correlation) {
      # Calculate original correlation
      orig_cor <- cor(returns, use = "pairwise.complete.obs")
      
      # Calculate correlation after contamination
      stress_cor <- cor(stressed_returns, use = "pairwise.complete.obs")
      
      # If correlation structure changed significantly, adjust
      if(norm(orig_cor - stress_cor, "F") > 0.1 * norm(orig_cor, "F")) {
        cat("Adjusting contamination to better preserve correlation structure\n")
        
        # Use factor model to adjust while preserving fatter tails
        # This is a simplified adjustment that blends the original and stressed returns
        blend_factor <- 0.7  # Preserve 70% of original correlation structure
        
        # Perform SVD on original correlation
        svd_result <- svd(orig_cor)
        
        # Reconstruct correlation with first few factors
        k <- min(10, ncol(returns)/2, nrow(returns)/5)
        reconstructed_cor <- svd_result$u[, 1:k] %*% diag(svd_result$d[1:k]) %*% t(svd_result$v[, 1:k])
        
        # Create adjusted returns
        adjusted_returns <- matrix(0, nrow = n, ncol = p)
        for(j in 1:p) {
          # Blend the stressed returns with a projection that preserves correlation
          contaminated <- which(contamination_mask[, j] == 0)
          if(length(contaminated) > 0) {
            # Only adjust contaminated points
            adjusted_returns[, j] <- stressed_returns[, j]
            
            # Calculate adjustment factor
            target_cors <- reconstructed_cor[j, ]
            adjustment <- rep(0, n)
            
            for(i in contaminated) {
              # Create adjustment based on correlation targets
              weighted_adj <- 0
              for(l in 1:p) {
                if(l != j) {
                  weighted_adj <- weighted_adj + target_cors[l] * returns[i, l] / p
                }
              }
              adjustment[i] <- weighted_adj
            }
            
            # Apply adjustment
            adjusted_returns[contaminated, j] <- blend_factor * stressed_returns[contaminated, j] + 
              (1 - blend_factor) * adjustment[contaminated]
          } else {
            adjusted_returns[, j] <- stressed_returns[, j]
          }
        }
        
        stressed_returns <- adjusted_returns
      }
    }
  }
  
  # Calculate overall contamination percentage
  contamination_pct <- 100 * (1 - sum(contamination_mask) / (n * p))
  cat("Total contamination applied:", round(contamination_pct, 2), "%\n")
  
  # Return stressed data and contamination information
  list(
    original_returns = returns,
    stressed_returns = stressed_returns,
    contamination_mask = contamination_mask,
    contamination_pct = contamination_pct,
    stress_method = method,
    stress_intensity = intensity,
    dates = dates
  )
}

# Run analysis with simulated market stress
run_stress_test_analysis <- function(results, 
                                     stress_method = "tail_contamination", 
                                     stress_intensity = 0.3,
                                     target_pct = 0.15,
                                     preserve_correlation = TRUE,
                                     by_regime = TRUE,
                                     output_dir = NULL) {
  
  # Extract original data from previous results
  original_data <- results$data
  
  # Apply simulated stress to the data
  stressed_data <- add_market_stress(
    returns = original_data$returns,
    dates = original_data$dates,
    method = stress_method,
    intensity = stress_intensity,
    target_pct = target_pct,
    preserve_correlation = preserve_correlation
  )
  
  # Create a copy of the original data with stressed returns
  stressed_data_full <- original_data
  stressed_data_full$returns <- stressed_data$stressed_returns
  stressed_data_full$returns_xts <- xts(stressed_data$stressed_returns, order.by = original_data$dates)
  
  # Identify market regimes (same as original analysis)
  regimes <- identify_market_regimes(stressed_data_full$dates)
  
  # Split into training and testing sets
  cat("Splitting stressed data into training and testing sets...\n")
  
  if (by_regime) {
    # Use regime-aware split
    split_data <- split_train_test_by_regime(
      stressed_data_full$returns, 
      stressed_data_full$dates, 
      stressed_data_full$market_returns, 
      train_ratio = 0.7, 
      regimes = regimes
    )
  } else {
    # Use original chronological split
    split_data <- split_train_test(
      stressed_data_full$returns, 
      stressed_data_full$dates, 
      stressed_data_full$market_returns, 
      train_ratio = 0.7
    )
    split_data$regimes <- regimes  # Add regimes to split data for later use
  }
  
  # Detect outliers in stressed training data
  cat("Detecting outliers in stressed training data...\n")
  outliers <- detect_outliers(split_data$train_returns, method = "MAD", threshold = 3)
  
  # Compare with true contamination (for validation)
  train_mask <- stressed_data$contamination_mask[split_data$train_indices, ]
  train_contamination_pct <- 100 * (1 - sum(train_mask) / length(train_mask))
  
  cat("True contamination in training data:", round(train_contamination_pct, 2), "%\n")
  cat("Detected contamination in training data:", round(outliers$contamination_pct, 2), "%\n")
  
  # Create and evaluate portfolios on stressed data
  cat("Creating and evaluating portfolios on stressed data...\n")
  portfolio_results <- create_and_evaluate_portfolios(
    split_data$train_returns, 
    split_data$test_returns,
    split_data$train_market,
    split_data$test_market,
    max_weight = 0.2,
    test_dates = split_data$test_dates,
    test_regimes = split_data$regimes[split_data$test_indices]
  )
  
  # Create regime-specific performance summary and visualizations
  if (by_regime) {
    cat("Creating regime-specific performance analysis for stressed data...\n")
    regime_performance <- create_regime_performance_summary(
      portfolio_results, 
      split_data$regimes[split_data$test_indices]
    )
    
    regime_plots <- create_regime_visualizations(regime_performance)
  } else {
    regime_performance <- NULL
    regime_plots <- NULL
  }
  
  # Create performance summary
  performance_summary <- data.frame(
    Method = names(portfolio_results),
    InSample_Return = sapply(portfolio_results, function(x) x$performance$train$ann_return),
    InSample_Volatility = sapply(portfolio_results, function(x) x$performance$train$ann_volatility),
    InSample_Sharpe = sapply(portfolio_results, function(x) x$performance$train$sharpe),
    OutOfSample_Return = sapply(portfolio_results, function(x) x$performance$test$ann_return),
    OutOfSample_Volatility = sapply(portfolio_results, function(x) x$performance$test$ann_volatility),
    OutOfSample_Sharpe = sapply(portfolio_results, function(x) x$performance$test$sharpe),
    MaxDrawdown = sapply(portfolio_results, function(x) x$performance$test$max_drawdown),
    Sortino = sapply(portfolio_results, function(x) x$performance$test$sortino),
    Beta = sapply(portfolio_results, function(x) x$performance$test$beta),
    Condition_Number = sapply(portfolio_results, function(x) {
      if("condition_number" %in% names(x)) x$condition_number else NA
    })
  )
  
  # Sort by out-of-sample Sharpe ratio
  performance_summary <- performance_summary[order(-performance_summary$OutOfSample_Sharpe), ]
  
  # Print performance summary
  cat("Stressed Data Performance Summary:\n")
  print(performance_summary)
  
  # Package results
  stress_results <- list(
    data = stressed_data_full,
    stress_info = stressed_data[c("stress_method", "stress_intensity", "contamination_pct")],
    split_data = split_data,
    outliers = outliers,
    portfolio_results = portfolio_results,
    performance_summary = performance_summary,
    regime_performance = regime_performance,
    regime_plots = regime_plots
  )
  
  # Set default output directory if not provided
  if(is.null(output_dir)) {
    output_dir <- paste0("sp500_stress_analysis_", stress_method, "_", format(Sys.Date(), "%Y%m%d"))
  }
  
  # Save results and create visualizations
  save_stress_results(stress_results, results, output_dir)
  
  return(stress_results)
}

# Function to save stress test results and generate comparative report
save_stress_results <- function(stress_results, original_results, output_dir) {
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save full results as RData
  save(stress_results, file = file.path(output_dir, "sp500_stress_analysis_results.RData"))
  
  # Save performance summary as CSV
  write.csv(stress_results$performance_summary, 
            file = file.path(output_dir, "stress_performance_summary.csv"),
            row.names = FALSE)
  
  # Create visualizations
  perf_plots <- create_performance_visualizations(stress_results)
  weight_plots <- create_weight_visualizations(stress_results)
  
  # Add a special function to create comparative visualizations
  comparative_plots <- create_comparative_visualizations(stress_results, original_results)
  
  # Save plots
  for(name in names(perf_plots)) {
    if(name != "top_weights") {  # Handle special case for list of plots
      ggsave(file.path(output_dir, paste0("stress_", name, ".png")), 
             perf_plots[[name]], width = 10, height = 6)
    }
  }
  
  # Save comparative plots
  for(name in names(comparative_plots)) {
    ggsave(file.path(output_dir, paste0("comparative_", name, ".png")), 
           comparative_plots[[name]], width = 10, height = 6)
  }
  
  # Generate comparative report
  report_file <- file.path(output_dir, "stress_analysis_report.md")
  
  # Get best method based on out-of-sample Sharpe in stress test
  stress_best_method <- stress_results$performance_summary$Method[1]
  stress_best_sharpe <- stress_results$performance_summary$OutOfSample_Sharpe[1]
  
  # Get best method from original analysis
  original_best_method <- original_results$performance_summary$Method[1]
  original_best_sharpe <- original_results$performance_summary$OutOfSample_Sharpe[1]
  
  # Prepare content
  report_content <- c(
    "# S&P 500 Portfolio Optimization Under Market Stress",
    "",
    "## Market Stress Simulation Parameters",
    "",
    paste("- Stress method:", stress_results$stress_info$stress_method),
    paste("- Stress intensity:", stress_results$stress_info$stress_intensity),
    paste("- Contamination level:", round(stress_results$stress_info$contamination_pct, 2), "%"),
    "",
    "## Performance Summary Under Stress",
    "",
    "```",
    capture.output(print(stress_results$performance_summary)),
    "```",
    "",
    "## Comparative Analysis",
    "",
    "### Best Performing Methods",
    "",
    paste("- Original data best method:", original_best_method),
    paste("  with out-of-sample Sharpe ratio:", round(original_best_sharpe, 4)),
    paste("- Stressed data best method:", stress_best_method),
    paste("  with out-of-sample Sharpe ratio:", round(stress_best_sharpe, 4)),
    "",
    "### Performance Degradation Under Stress",
    ""
  )
  
  # Calculate and add method-by-method performance degradation
  common_methods <- intersect(stress_results$performance_summary$Method, 
                              original_results$performance_summary$Method)
  
  degradation_data <- data.frame(
    Method = common_methods,
    Original_Sharpe = sapply(common_methods, function(m) {
      original_results$performance_summary$OutOfSample_Sharpe[
        original_results$performance_summary$Method == m]
    }),
    Stress_Sharpe = sapply(common_methods, function(m) {
      stress_results$performance_summary$OutOfSample_Sharpe[
        stress_results$performance_summary$Method == m]
    })
  )
  
  degradation_data$Absolute_Change <- degradation_data$Stress_Sharpe - degradation_data$Original_Sharpe
  degradation_data$Percent_Change <- (degradation_data$Stress_Sharpe / degradation_data$Original_Sharpe - 1) * 100
  
  # Sort by percent change (least degradation/most improvement first)
  degradation_data <- degradation_data[order(-degradation_data$Percent_Change), ]
  
  report_content <- c(report_content,
                      "```",
                      capture.output(print(degradation_data)),
                      "```",
                      "",
                      "### Most Robust Methods (Least Performance Degradation)",
                      "")
  
  # Add the top 5 most robust methods
  for(i in 1:min(5, nrow(degradation_data))) {
    method <- degradation_data$Method[i]
    pct_change <- degradation_data$Percent_Change[i]
    change_text <- ifelse(pct_change >= 0, 
                          paste("improved by", round(pct_change, 2), "%"),
                          paste("degraded by", round(-pct_change, 2), "%"))
    
    report_content <- c(report_content,
                        paste(i, ". ", method, " - ", change_text, sep=""))
  }
  
  # Add insights about robust methods
  report_content <- c(report_content,
                      "",
                      "## Key Insights",
                      "",
                      "### Characteristics of Robust Methods",
                      "")
  
  # Identify which covariance methods were most robust
  cov_methods <- unique(sapply(strsplit(degradation_data$Method, "_"), function(x) x[1]))
  cov_methods <- cov_methods[cov_methods != "EqualWeights"]
  
  avg_degradation_by_cov <- sapply(cov_methods, function(cov_method) {
    methods <- degradation_data$Method[grepl(paste0("^", cov_method), degradation_data$Method)]
    mean(degradation_data$Percent_Change[degradation_data$Method %in% methods])
  })
  
  ranked_cov_methods <- names(sort(avg_degradation_by_cov, decreasing = TRUE))
  
  report_content <- c(report_content,
                      "#### Most Robust Covariance Estimation Methods:")
  
  for(i in 1:length(ranked_cov_methods)) {
    avg_change <- avg_degradation_by_cov[ranked_cov_methods[i]]
    change_text <- ifelse(avg_change >= 0, 
                          paste("improved by", round(avg_change, 2), "%"),
                          paste("degraded by", round(-avg_change, 2), "%"))
    
    report_content <- c(report_content,
                        paste("  ", i, ". ", ranked_cov_methods[i], " (", change_text, ")", sep=""))
  }
  
  # Compare optimization methods
  report_content <- c(report_content,
                      "",
                      "#### Optimization Method Robustness:")
  
  # Compare Min Variance vs Max Sharpe
  minvar_methods <- degradation_data$Method[grepl("MinVar", degradation_data$Method)]
  maxsharpe_methods <- degradation_data$Method[grepl("MaxSharpe", degradation_data$Method)]
  
  if(length(minvar_methods) > 0 && length(maxsharpe_methods) > 0) {
    minvar_avg_change <- mean(degradation_data$Percent_Change[degradation_data$Method %in% minvar_methods])
    maxsharpe_avg_change <- mean(degradation_data$Percent_Change[degradation_data$Method %in% maxsharpe_methods])
    
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
                      paste("This analysis confirms that", ranked_cov_methods[1], "and", 
                            ranked_cov_methods[2], "are the most robust covariance estimation",
                            "methods under market stress conditions. These methods",
                            "effectively handle the increased contamination present in",
                            "financial returns during periods of market turbulence and",
                            "maintain more stable portfolio performance."),
                      "",
                      paste("For practical implementation, the", stress_best_method, 
                            "approach shows the best combination of performance and",
                            "robustness, making it suitable for real-world portfolio",
                            "management in both normal and stressed market conditions.")
  )
  
  # Write the report
  writeLines(report_content, report_file)
  
  cat("Stress analysis results and report saved to", output_dir, "\n")
}

# Create comparative visualizations between normal and stressed results
create_comparative_visualizations <- function(stress_results, original_results) {
  library(ggplot2)
  library(reshape2)
  
  plots <- list()
  
  # Get common methods between both analyses
  common_methods <- intersect(stress_results$performance_summary$Method, 
                              original_results$performance_summary$Method)
  
  # Prepare data for comparison
  comparison_data <- data.frame(
    Method = common_methods,
    Original_Sharpe = sapply(common_methods, function(m) {
      original_results$performance_summary$OutOfSample_Sharpe[
        original_results$performance_summary$Method == m]
    }),
    Stress_Sharpe = sapply(common_methods, function(m) {
      stress_results$performance_summary$OutOfSample_Sharpe[
        stress_results$performance_summary$Method == m]
    })
  )
  
  comparison_data$Change <- comparison_data$Stress_Sharpe - comparison_data$Original_Sharpe
  comparison_data$Percent_Change <- (comparison_data$Stress_Sharpe / comparison_data$Original_Sharpe - 1) * 100
  
  # Extract method types for coloring
  comparison_data$OptimizationType <- ifelse(grepl("MaxSharpe", comparison_data$Method), "Max Sharpe", 
                                             ifelse(grepl("MinVar", comparison_data$Method), "Min Variance", "Equal Weight"))
  
  comparison_data$CovMethod <- sapply(strsplit(as.character(comparison_data$Method), "_"), function(x) x[1])
  
  # 1. Sharpe Ratio Comparison Bar Chart
  # Reshape for side-by-side comparison
  sharpe_comparison <- melt(comparison_data, 
                            id.vars = c("Method", "OptimizationType", "CovMethod"),
                            measure.vars = c("Original_Sharpe", "Stress_Sharpe"),
                            variable.name = "Condition", 
                            value.name = "Sharpe_Ratio")
  
  # Clean up condition names
  sharpe_comparison$Condition <- ifelse(sharpe_comparison$Condition == "Original_Sharpe",
                                        "Normal Market", "Stressed Market")
  
  # Order by original performance
  method_order <- comparison_data$Method[order(-comparison_data$Original_Sharpe)]
  sharpe_comparison$Method <- factor(sharpe_comparison$Method, levels = method_order)
  
  # Create bar chart
  p1 <- ggplot(sharpe_comparison, aes(x = Method, y = Sharpe_Ratio, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("Normal Market" = "steelblue", "Stressed Market" = "firebrick")) +
    labs(title = "Portfolio Performance: Normal vs. Stressed Market",
         x = "", y = "Out-of-Sample Sharpe Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plots$sharpe_comparison <- p1
  
  # 2. Performance Degradation
  # Order by percent change (least degradation first)
  comparison_data <- comparison_data[order(-comparison_data$Percent_Change), ]
  comparison_data$Method <- factor(comparison_data$Method, levels = comparison_data$Method)
  
  # Create performance change plot
  p2 <- ggplot(comparison_data, aes(x = Method, y = Percent_Change, fill = OptimizationType)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Performance Change Under Market Stress",
         x = "", y = "Percent Change in Sharpe Ratio (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plots$performance_change <- p2
  
  # 3. Covariance Method Comparison
  # Average change by covariance method
  cov_method_data <- data.frame()
  
  for(cov_method in unique(comparison_data$CovMethod)) {
    if(cov_method == "EqualWeights") next
    
    # MinVar methods
    minvar_methods <- comparison_data$Method[comparison_data$CovMethod == cov_method & 
                                               comparison_data$OptimizationType == "Min Variance"]
    
    if(length(minvar_methods) > 0) {
      avg_change <- mean(comparison_data$Percent_Change[comparison_data$Method %in% minvar_methods])
      cov_method_data <- rbind(cov_method_data, data.frame(
        CovMethod = cov_method,
        OptimizationType = "Min Variance",
        Percent_Change = avg_change
      ))
    }
    
    # MaxSharpe methods
    maxsharpe_methods <- comparison_data$Method[comparison_data$CovMethod == cov_method & 
                                                  comparison_data$OptimizationType == "Max Sharpe"]
    
    if(length(maxsharpe_methods) > 0) {
      avg_change <- mean(comparison_data$Percent_Change[comparison_data$Method %in% maxsharpe_methods])
      cov_method_data <- rbind(cov_method_data, data.frame(
        CovMethod = cov_method,
        OptimizationType = "Max Sharpe",
        Percent_Change = avg_change
      ))
    }
  }
  
  # Order covariance methods by average performance change
  avg_by_cov <- aggregate(Percent_Change ~ CovMethod, cov_method_data, mean)
  method_order <- avg_by_cov$CovMethod[order(-avg_by_cov$Percent_Change)]
  cov_method_data$CovMethod <- factor(cov_method_data$CovMethod, levels = method_order)
  
  # Create covariance method comparison plot
  p3 <- ggplot(cov_method_data, aes(x = CovMethod, y = Percent_Change, fill = OptimizationType)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Robustness of Covariance Methods Under Stress",
         x = "Covariance Method", y = "Percent Change in Sharpe Ratio (%)") +
    theme_minimal()
  
  plots$cov_method_comparison <- p3
  
  # 4. Scatter plot of normal vs. stressed performance
  p4 <- ggplot(comparison_data, aes(x = Original_Sharpe, y = Stress_Sharpe, 
                                    color = CovMethod, shape = OptimizationType)) +
    geom_point(size = 3) +
    geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(title = "Performance Comparison: Normal vs. Stressed Market",
         x = "Normal Market Sharpe Ratio", y = "Stressed Market Sharpe Ratio") +
    theme_minimal()
  
  plots$performance_scatter <- p4
  
  # 5. Risk-Return Trade-off Comparison
  # Get volatility and return data
  risk_return_data <- data.frame(
    Method = common_methods,
    
    Original_Return = sapply(common_methods, function(m) {
      original_results$performance_summary$OutOfSample_Return[
        original_results$performance_summary$Method == m]
    }),
    Original_Volatility = sapply(common_methods, function(m) {
      original_results$performance_summary$OutOfSample_Volatility[
        original_results$performance_summary$Method == m]
    }),
    
    Stress_Return = sapply(common_methods, function(m) {
      stress_results$performance_summary$OutOfSample_Return[
        stress_results$performance_summary$Method == m]
    }),
    Stress_Volatility = sapply(common_methods, function(m) {
      stress_results$performance_summary$OutOfSample_Volatility[
        stress_results$performance_summary$Method == m]
    })
  )
  
  # Add method info
  risk_return_data$OptimizationType <- ifelse(grepl("MaxSharpe", risk_return_data$Method), "Max Sharpe", 
                                              ifelse(grepl("MinVar", risk_return_data$Method), "Min Variance", "Equal Weight"))
  
  risk_return_data$CovMethod <- sapply(strsplit(as.character(risk_return_data$Method), "_"), function(x) x[1])
  
  # Reshape for plotting
  risk_return_long <- melt(risk_return_data, 
                           id.vars = c("Method", "OptimizationType", "CovMethod"),
                           measure.vars = c("Original_Return", "Original_Volatility", 
                                            "Stress_Return", "Stress_Volatility"))
  
  # Create new variables for condition and metric
  risk_return_long$Condition <- ifelse(grepl("Original", risk_return_long$variable), 
                                       "Normal Market", "Stressed Market")
  risk_return_long$Metric <- ifelse(grepl("Return", risk_return_long$variable), 
                                    "Return", "Volatility")
  
  # Reshape to get one row per method and condition
  risk_return_wide <- dcast(risk_return_long, 
                            Method + OptimizationType + CovMethod + Condition ~ Metric, 
                            value.var = "value")
  
  # Create risk-return scatter plot
  p5 <- ggplot(risk_return_wide, aes(x = Volatility, y = Return, 
                                     color = CovMethod, shape = Condition)) +
    geom_point(size = 3) +
    geom_path(aes(group = Method), linetype = "dashed", alpha = 0.5) +
    labs(title = "Risk-Return Profile: Normal vs. Stressed Market",
         x = "Annualized Volatility", y = "Annualized Return") +
    theme_minimal()
  
  plots$risk_return_comparison <- p5
  
  return(plots)
}

# Run multi-scenario stress testing
run_multi_scenario_analysis <- function(start_date = NULL,
                                        end_date = Sys.Date(),
                                        n_stocks = 100,
                                        by_regime = TRUE,
                                        output_dir = NULL) {
  
  # Set start date if not provided (default to ~10 years ago)
  if(is.null(start_date)) {
    start_date <- as.Date(end_date) - 365*10  # 10 years of data
  }
  
  # Set default output directory if not provided
  if(is.null(output_dir)) {
    output_dir <- paste0("multi_scenario_analysis_", format(Sys.Date(), "%Y%m%d"))
  }
  
  cat("Running multi-scenario robust portfolio analysis from", format(start_date, "%Y-%m-%d"), 
      "to", format(end_date, "%Y-%m-%d"), "\n")
  
  # Step 1: Run the standard analysis
  cat("\n===== PHASE 1: STANDARD MARKET ANALYSIS =====\n\n")
  standard_results <- run_sp500_analysis(
    start_date = start_date,
    end_date = end_date,
    n_stocks = n_stocks,
    train_ratio = 0.7,
    by_regime = by_regime,
    rolling = FALSE,
    max_weight = 0.2
  )
  
  # Create standard output directory
  standard_dir <- file.path(output_dir, "standard_analysis")
  save_and_report(standard_results, standard_dir)
  
  # Step 2: Define stress scenarios
  scenarios <- list(
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
  
  # Step 3: Run each stress scenario
  stress_results <- list()
  
  for(scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    
    cat("\n===== PHASE 2: STRESS SCENARIO -", scenario$description, "=====\n\n")
    
    stress_results[[scenario_name]] <- run_stress_test_analysis(
      results = standard_results,
      stress_method = scenario$method,
      stress_intensity = scenario$intensity,
      target_pct = scenario$target_pct,
      preserve_correlation = TRUE,
      by_regime = by_regime,
      output_dir = file.path(output_dir, paste0("stress_", scenario_name))
    )
  }
  
  # Step 4: Create comparative analysis across all scenarios
  create_multi_scenario_report(standard_results, stress_results, scenarios, output_dir)
  
  # Step 5: Create visuals comparing method performance across scenarios
  create_multi_scenario_visuals(standard_results, stress_results, scenarios, output_dir)
  
  # Clean up parallel cluster
  stopCluster(cl)
  
  return(list(
    standard_results = standard_results,
    stress_results = stress_results
  ))
}

# Create report comparing method performance across multiple stress scenarios
create_multi_scenario_report <- function(standard_results, stress_results, scenarios, output_dir) {
  report_file <- file.path(output_dir, "multi_scenario_analysis_report.md")
  
  # Extract methods
  methods <- standard_results$performance_summary$Method
  
  # Create matrix of Sharpe ratios across scenarios
  sharpe_matrix <- matrix(NA, nrow = length(methods), ncol = 1 + length(scenarios))
  colnames(sharpe_matrix) <- c("Standard", names(scenarios))
  rownames(sharpe_matrix) <- methods
  
  # Fill in standard results
  for(i in 1:length(methods)) {
    method <- methods[i]
    idx <- which(standard_results$performance_summary$Method == method)
    if(length(idx) > 0) {
      sharpe_matrix[i, 1] <- standard_results$performance_summary$OutOfSample_Sharpe[idx]
    }
  }
  
  # Fill in stress scenario results
  for(j in 1:length(scenarios)) {
    scenario_name <- names(scenarios)[j]
    scenario_results <- stress_results[[scenario_name]]
    
    for(i in 1:length(methods)) {
      method <- methods[i]
      idx <- which(scenario_results$performance_summary$Method == method)
      if(length(idx) > 0) {
        sharpe_matrix[i, j+1] <- scenario_results$performance_summary$OutOfSample_Sharpe[idx]
      }
    }
  }
  
  # Calculate average performance across all scenarios
  avg_performance <- rowMeans(sharpe_matrix, na.rm = TRUE)
  
  # Calculate robustness score (lower standard deviation = more robust)
  robustness_score <- apply(sharpe_matrix, 1, function(x) -sd(x, na.rm = TRUE))
  
  # Calculate combined score (weighing both performance and robustness)
  combined_score <- scale(avg_performance) + scale(robustness_score)
  
  # Combine into a ranking data frame
  ranking_df <- data.frame(
    Method = methods,
    Avg_Sharpe = avg_performance,
    Robustness = -robustness_score,  # Convert back to std dev for reporting
    Combined_Score = combined_score
  )
  
  # Sort by combined score
  ranking_df <- ranking_df[order(-ranking_df$Combined_Score), ]
  
  # Prepare report content
  report_content <- c(
    "# Multi-Scenario Robust Portfolio Analysis",
    "",
    "## Scenarios Analyzed",
    ""
  )
  
  # Add scenario descriptions
  for(scenario_name in names(scenarios)) {
    scenario <- scenarios[[scenario_name]]
    report_content <- c(report_content,
                        paste("### Scenario:", scenario_name),
                        paste("- Description:", scenario$description),
                        paste("- Method:", scenario$method),
                        paste("- Intensity:", scenario$intensity),
                        paste("- Target percentage:", scenario$target_pct),
                        "")
  }
  
  # Add performance summary
  report_content <- c(report_content,
                      "## Performance Across Scenarios",
                      "",
                      "```",
                      capture.output(print(sharpe_matrix)),
                      "```",
                      "",
                      "## Method Ranking (Combined Performance and Robustness)",
                      "",
                      "```",
                      capture.output(print(ranking_df)),
                      "```",
                      "",
                      "## Top Performing Methods",
                      "")
  
  # Add top 5 methods
  for(i in 1:min(5, nrow(ranking_df))) {
    method <- ranking_df$Method[i]
    avg_sharpe <- ranking_df$Avg_Sharpe[i]
    robustness <- ranking_df$Robustness[i]
    
    report_content <- c(report_content,
                        paste(i, ". ", method, sep=""),
                        paste("   - Average Sharpe:", round(avg_sharpe, 4)),
                        paste("   - Robustness (Std Dev):", round(robustness, 4)),
                        "")
  }
  
  # Add recommendations
  best_overall <- ranking_df$Method[1]
  
  # Find best by optimization type
  best_minvar <- ranking_df$Method[grep("MinVar", ranking_df$Method)][1]
  best_maxsharpe <- ranking_df$Method[grep("MaxSharpe", ranking_df$Method)][1]
  
  report_content <- c(report_content,
                      "## Recommendations",
                      "",
                      paste("### Best Overall Method: **", best_overall, "**"),
                      paste("This method provides the best balance of performance and robustness",
                            "across all tested market scenarios."),
                      "",
                      paste("### Best Minimum Variance Strategy: **", best_minvar, "**"),
                      paste("### Best Maximum Sharpe Strategy: **", best_maxsharpe, "**"),
                      "",
                      "## Conclusions",
                      "",
                      paste("The multi-scenario analysis confirms that robust covariance estimation",
                            "methods provide significant advantages in portfolio optimization,",
                            "especially under stressed market conditions. The", best_overall,
                            "approach consistently outperforms other methods across various",
                            "market stress scenarios."),
                      "",
                      paste("For practical implementation, portfolio managers should consider",
                            "their specific risk preferences and market outlook. In periods of",
                            "anticipated market stability, methods optimized for performance may",
                            "be preferred. However, during uncertain times or for risk-averse",
                            "investors, methods with higher robustness scores would be more appropriate."))
  
  # Write the report
  writeLines(report_content, report_file)
  
  cat("Multi-scenario analysis report saved to", report_file, "\n")
}

# Create visualizations comparing method performance across scenarios
create_multi_scenario_visuals <- function(standard_results, stress_results, scenarios, output_dir) {
  library(ggplot2)
  library(reshape2)
  
  # Create directory for visualizations
  viz_dir <- file.path(output_dir, "visualizations")
  if(!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE)
  }
  
  # Extract methods
  methods <- standard_results$performance_summary$Method
  
  # Create data frame of Sharpe ratios across scenarios
  sharpe_df <- data.frame(Method = methods)
  
  # Add standard results
  sharpe_df$Standard <- NA
  for(i in 1:length(methods)) {
    method <- methods[i]
    idx <- which(standard_results$performance_summary$Method == method)
    if(length(idx) > 0) {
      sharpe_df$Standard[i] <- standard_results$performance_summary$OutOfSample_Sharpe[idx]
    }
  }
  
  # Add stress scenario results
  for(scenario_name in names(scenarios)) {
    scenario_results <- stress_results[[scenario_name]]
    sharpe_df[[scenario_name]] <- NA
    
    for(i in 1:length(methods)) {
      method <- methods[i]
      idx <- which(scenario_results$performance_summary$Method == method)
      if(length(idx) > 0) {
        sharpe_df[[scenario_name]][i] <- scenario_results$performance_summary$OutOfSample_Sharpe[idx]
      }
    }
  }
  
  # Add method type for coloring
  sharpe_df$OptimizationType <- ifelse(grepl("MaxSharpe", sharpe_df$Method), "Max Sharpe", 
                                       ifelse(grepl("MinVar", sharpe_df$Method), "Min Variance", "Equal Weight"))
  
  sharpe_df$CovMethod <- sapply(strsplit(as.character(sharpe_df$Method), "_"), function(x) x[1])
  
  # Calculate average performance and stability metrics
  sharpe_df$Avg_Sharpe <- rowMeans(sharpe_df[, c("Standard", names(scenarios))], na.rm = TRUE)
  sharpe_df$Stability <- -apply(sharpe_df[, c("Standard", names(scenarios))], 1, sd, na.rm = TRUE)
  
  # 1. Heatmap of performance across scenarios
  # Reshape for heatmap
  heatmap_data <- melt(sharpe_df[, c("Method", "Standard", names(scenarios))],
                       id.vars = "Method",
                       variable.name = "Scenario",
                       value.name = "Sharpe")
  
  # Order by average performance
  method_order <- sharpe_df$Method[order(-sharpe_df$Avg_Sharpe)]
  heatmap_data$Method <- factor(heatmap_data$Method, levels = method_order)
  
  # Create heatmap
  p1 <- ggplot(heatmap_data, aes(x = Scenario, y = Method, fill = Sharpe)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = median(heatmap_data$Sharpe, na.rm = TRUE)) +
    labs(title = "Performance Heatmap Across Market Scenarios",
         x = "Scenario", y = "Method") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # Save heatmap
  ggsave(file.path(viz_dir, "performance_heatmap.png"), p1, width = 12, height = 10)
  
  # 2. Performance vs. Stability scatter plot
  p2 <- ggplot(sharpe_df, aes(x = Stability, y = Avg_Sharpe, 
                              color = CovMethod, shape = OptimizationType)) +
    geom_point(size = 3) +
    geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
    labs(title = "Performance vs. Stability Across Scenarios",
         x = "Stability (negative std dev, higher is more stable)", 
         y = "Average Sharpe Ratio") +
    theme_minimal()
  
  # Save scatter plot
  ggsave(file.path(viz_dir, "performance_vs_stability.png"), p2, width = 12, height = 8)
  
  # 3. Bar chart of average performance by covariance method
  cov_summary <- aggregate(Avg_Sharpe ~ CovMethod + OptimizationType, 
                           data = sharpe_df, FUN = mean)
  
  p3 <- ggplot(cov_summary, aes(x = CovMethod, y = Avg_Sharpe, fill = OptimizationType)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Average Performance by Covariance Method",
         x = "Covariance Method", y = "Average Sharpe Ratio") +
    theme_minimal()
  
  # Save bar chart
  ggsave(file.path(viz_dir, "avg_performance_by_method.png"), p3, width = 10, height = 6)
  
  # 4. Line chart showing performance across scenarios
  # Select top 5 methods for clarity
  top_methods <- head(sharpe_df$Method[order(-sharpe_df$Avg_Sharpe)], 5)
  line_data <- heatmap_data[heatmap_data$Method %in% top_methods, ]
  
  p4 <- ggplot(line_data, aes(x = Scenario, y = Sharpe, color = Method, group = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(title = "Top 5 Methods: Performance Across Scenarios",
         x = "Scenario", y = "Sharpe Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save line chart
  ggsave(file.path(viz_dir, "top_methods_comparison.png"), p4, width = 10, height = 6)
  
  cat("Multi-scenario visualizations saved to", viz_dir, "\n")
}
        OptimizationType = "Min
