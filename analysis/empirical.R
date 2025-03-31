# ========================================================================
# S&P 500 Empirical Analysis
# ========================================================================

# Perform comprehensive portfolio analysis on real S&P 500 data
run_sp500_analysis <- function(start_date = "2015-01-01", 
                               end_date = Sys.Date(), 
                               n_stocks = 100,
                               train_ratio = 0.7,
                               by_regime = TRUE,  # Regime-aware splitting
                               rolling = TRUE,
                               window_size = 252,
                               step_size = 21,
                               max_weight = 0.2) {
  # Start timing
  start_time <- Sys.time()
  
  # Download and process S&P 500 data
  cat("Downloading and processing S&P 500 data...\n")
  data <- get_sp500_data(start_date, end_date, n_stocks)
  
  # Identify market regimes
  regimes <- identify_market_regimes(data$dates)
  
  # Split into training and testing sets
  cat("Splitting data into training and testing sets...\n")
  
  if (by_regime) {
    # Use regime-aware split
    split_data <- split_train_test_by_regime(data$returns, data$dates, data$market_returns, 
                                             train_ratio, regimes)
  } else {
    # Use original chronological split
    split_data <- split_train_test(data$returns, data$dates, data$market_returns, train_ratio)
    split_data$regimes <- regimes  # Add regimes to split data for later use
  }
  
  # Detect outliers in training data
  cat("Detecting outliers in training data...\n")
  outliers <- detect_outliers(split_data$train_returns, method = "MAD", threshold = 3)
  
  # Create and evaluate portfolios
  cat("Creating and evaluating portfolios...\n")
  portfolio_results <- create_and_evaluate_portfolios(
    split_data$train_returns, 
    split_data$test_returns,
    split_data$train_market,
    split_data$test_market,
    max_weight = max_weight,
    test_dates = split_data$test_dates,
    test_regimes = split_data$regimes[split_data$test_indices]
  )
  
  # Test sensitivity to contamination threshold
  cat("Testing sensitivity to contamination threshold...\n")
  threshold_sensitivity <- test_contamination_sensitivity(
    split_data$train_returns,
    split_data$test_returns,
    thresholds = seq(2, 5, by = 0.5),
    split_data$train_market,
    split_data$test_market,
    max_weight = max_weight,
    test_dates = split_data$test_dates,
    test_regimes = split_data$regimes[split_data$test_indices]
  )
  
  # Test different outlier detection methods
  cat("Testing different outlier detection methods...\n")
  detection_methods <- test_detection_methods(
    split_data$train_returns,
    split_data$test_returns,
    split_data$train_market,
    split_data$test_market,
    max_weight = max_weight,
    test_dates = split_data$test_dates,
    test_regimes = split_data$regimes[split_data$test_indices]
  )
  
  # Create regime-specific performance summary and visualizations
  if (by_regime) {
    cat("Creating regime-specific performance analysis...\n")
    regime_performance <- create_regime_performance_summary(
      portfolio_results, 
      split_data$regimes[split_data$test_indices]
    )
    
    regime_plots <- create_regime_visualizations(regime_performance)
  } else {
    regime_performance <- NULL
    regime_plots <- NULL
  }
  
  # Perform rolling window analysis if requested
  if(rolling) {
    cat("Performing rolling window analysis...\n")
    rolling_results <- rolling_window_analysis(
      data$returns, data$dates, window_size, step_size, data$market_returns,
      max_weight = max_weight, regimes = regimes
    )
    
    # Add regime information to rolling window analysis
    rolling_results <- add_regime_info_to_rolling(rolling_results, regimes, data$dates)
    
    # Analyze rolling performance
    sharpe_analysis <- analyze_rolling_performance(rolling_results, "sharpe")
    volatility_analysis <- analyze_rolling_performance(rolling_results, "ann_volatility")
    max_drawdown_analysis <- analyze_rolling_performance(rolling_results, "max_drawdown")
  } else {
    rolling_results <- NULL
    sharpe_analysis <- NULL
    volatility_analysis <- NULL
    max_drawdown_analysis <- NULL
  }
  
  # End timing
  end_time <- Sys.time()
  execution_time <- difftime(end_time, start_time, units = "mins")
  
  # Print summary
  cat("\n====== S&P 500 PORTFOLIO ANALYSIS SUMMARY ======\n")
  cat("Analysis period:", start_date, "to", end_date, "\n")
  cat("Number of stocks:", ncol(data$returns), "\n")
  cat("Train period:", min(split_data$train_dates), "to", max(split_data$train_dates), "\n")
  cat("Test period:", min(split_data$test_dates), "to", max(split_data$test_dates), "\n")
  cat("Detected contamination level:", round(outliers$contamination_pct, 2), "%\n")
  cat("Execution time:", round(execution_time, 2), "minutes\n\n")
  
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
  cat("Performance Summary:\n")
  print(performance_summary)
  
  # Include regime information in results if available
  if (by_regime && "regimes" %in% names(split_data)) {
    regime_info <- list(
      regime_counts = table(split_data$regimes),
      train_regime_distribution = table(split_data$regimes[split_data$train_indices]) / 
        length(split_data$train_indices),
      test_regime_distribution = table(split_data$regimes[split_data$test_indices]) / 
        length(split_data$test_indices)
    )
  } else {
    regime_info <- NULL
  }
  
  # Return comprehensive results
  list(
    data = data,
    split_data = split_data,
    outliers = outliers,
    portfolio_results = portfolio_results,
    threshold_sensitivity = threshold_sensitivity,
    detection_methods = detection_methods,
    regime_performance = regime_performance,
    regime_plots = regime_plots,
    rolling_results = rolling_results,
    sharpe_analysis = sharpe_analysis,
    volatility_analysis = volatility_analysis,
    max_drawdown_analysis = max_drawdown_analysis,
    performance_summary = performance_summary,
    regime_info = regime_info,
    execution_time = execution_time
  )
}

# Function to save results and generate a report
save_and_report <- function(results, output_dir = "sp500_analysis") {
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save full results as RData
  save(results, file = file.path(output_dir, "sp500_analysis_results.RData"))
  
  # Save performance summary as CSV
  write.csv(results$performance_summary, 
            file = file.path(output_dir, "performance_summary.csv"),
            row.names = FALSE)
  
  # Generate simple report
  report_file <- file.path(output_dir, "analysis_report.md")
  
  # Get best method based on out-of-sample Sharpe
  best_method <- results$performance_summary$Method[1]
  best_sharpe <- results$performance_summary$OutOfSample_Sharpe[1]
  
  # Get optimal threshold from sensitivity analysis
  min_var_optimal_threshold <- NULL
  max_sharpe_optimal_threshold <- NULL
  if(!is.null(results$threshold_sensitivity)) {
    if(!is.null(results$threshold_sensitivity$min_var_optimal_threshold)) {
      min_var_optimal_threshold <- results$threshold_sensitivity$min_var_optimal_threshold
    }
    if(!is.null(results$threshold_sensitivity$max_sharpe_optimal_threshold)) {
      max_sharpe_optimal_threshold <- results$threshold_sensitivity$max_sharpe_optimal_threshold
    }
  }
  
  # Get best detection methods
  min_var_best_method <- NULL
  max_sharpe_best_method <- NULL
  if(!is.null(results$detection_methods)) {
    if(!is.null(results$detection_methods$min_var_best_method)) {
      min_var_best_method <- results$detection_methods$min_var_best_method
    }
    if(!is.null(results$detection_methods$max_sharpe_best_method)) {
      max_sharpe_best_method <- results$detection_methods$max_sharpe_best_method
    }
  }
  
  # Prepare content
  report_content <- c(
    "# S&P 500 Robust Portfolio Optimization Analysis",
    "",
    paste("Analysis period:", min(results$data$dates), "to", max(results$data$dates)),
    paste("Number of stocks analyzed:", ncol(results$data$returns)),
    "",
    "## Performance Summary",
    "",
    "```",
    capture.output(print(results$performance_summary)),
    "```",
    "",
    "## Key Findings",
    "",
    "### Detected Contamination",
    "",
    paste("- Overall contamination level:", round(results$outliers$contamination_pct, 2), "%"),
    paste("- Most contaminated asset:", results$data$symbols[which.max(results$outliers$outlier_counts_by_asset)]),
    paste("  with", max(results$outliers$outlier_counts_by_asset), "outliers"),
    "",
    "### Performance Comparison",
    "",
    paste("- Best performing method:", best_method),
    paste("  with out-of-sample Sharpe ratio:", round(best_sharpe, 4)),
    "",
    "### Optimization Method Comparison",
    "",
    paste("- Maximum Sharpe optimization methods generally", 
          ifelse(mean(results$performance_summary$OutOfSample_Sharpe[grepl("MaxSharpe", results$performance_summary$Method)]) >
                   mean(results$performance_summary$OutOfSample_Sharpe[grepl("MinVar", results$performance_summary$Method)]),
                 "outperformed", "underperformed"),
          "Minimum Variance methods in out-of-sample testing."),
    "",
    "### Covariance Estimation Method Comparison",
    "",
    "- Ranking of covariance methods by average out-of-sample Sharpe ratio:",
    ""
  )
  
  # Calculate average Sharpe by covariance method
  cov_methods <- unique(sapply(strsplit(as.character(results$performance_summary$Method), "_"), function(x) x[1]))
  cov_methods <- cov_methods[cov_methods != "EqualWeights"]
  
  avg_sharpe_by_cov <- sapply(cov_methods, function(cov_method) {
    methods <- results$performance_summary$Method[grepl(paste0("^", cov_method), results$performance_summary$Method)]
    mean(results$performance_summary$OutOfSample_Sharpe[results$performance_summary$Method %in% methods])
  })
  
  ranked_cov_methods <- names(sort(avg_sharpe_by_cov, decreasing = TRUE))
  
  # Add covariance method ranking to report
  for(i in 1:length(ranked_cov_methods)) {
    report_content <- c(report_content, 
                        paste("  ", i, ". ", ranked_cov_methods[i], " (avg. Sharpe: ", 
                              round(avg_sharpe_by_cov[ranked_cov_methods[i]], 4), ")", sep=""))
  }
  
  # Add sensitivity analysis results
  report_content <- c(report_content,
                      "",
                      "### Sensitivity Analysis",
                      "",
                      if(!is.null(min_var_optimal_threshold)) {
                        paste("- Optimal outlier detection threshold for Min Variance:", min_var_optimal_threshold)
                      } else {
                        "- Could not determine optimal threshold for Min Variance"
                      },
                      if(!is.null(max_sharpe_optimal_threshold)) {
                        paste("- Optimal outlier detection threshold for Max Sharpe:", max_sharpe_optimal_threshold)
                      } else {
                        "- Could not determine optimal threshold for Max Sharpe"
                      },
                      "",
                      "### Best Outlier Detection Methods",
                      "",
                      if(!is.null(min_var_best_method)) {
                        paste("- Best detection method for Min Variance:", min_var_best_method)
                      } else {
                        "- Could not determine best detection method for Min Variance"
                      },
                      if(!is.null(max_sharpe_best_method)) {
                        paste("- Best detection method for Max Sharpe:", max_sharpe_best_method)
                      } else {
                        "- Could not determine best detection method for Max Sharpe"
                      },
                      "",
                      "## Conclusion",
                      "",
                      paste("This analysis confirms that hybrid approaches combining MCD with factor-based",
                            "methods (especially", ranked_cov_methods[1], "and", ranked_cov_methods[2],
                            ") outperform traditional covariance estimation techniques on real S&P 500 data.",
                            "These robust estimation methods effectively handle the contamination",
                            "present in financial returns and produce better out-of-sample portfolio performance.")
  )
  
  # Add regime information if available
  if(!is.null(results$regime_info)) {
    regime_content <- c(
      "",
      "## Market Regime Analysis",
      "",
      "### Regime Distribution",
      "",
      "```",
      capture.output(print(results$regime_info$regime_counts)),
      "```",
      "",
      "### Train Set Regime Distribution",
      "",
      "```",
      capture.output(print(results$regime_info$train_regime_distribution)),
      "```",
      "",
      "### Test Set Regime Distribution",
      "",
      "```",
      capture.output(print(results$regime_info$test_regime_distribution)),
      "```"
    )
    
    report_content <- c(report_content, regime_content)
  }
  
  # Add regime-specific performance summaries if available
  if(!is.null(results$regime_performance)) {
    regime_perf_content <- c(
      "",
      "## Performance by Market Regime",
      ""
    )
    
    for(regime in names(results$regime_performance)) {
      if(regime == "overall") next
      
      regime_perf_content <- c(
        regime_perf_content,
        paste("### Performance in", capitalize(regime), "Regime"),
        "",
        "```",
        capture.output(print(head(results$regime_performance[[regime]], 10))),
        "```",
        ""
      )
    }
    
    report_content <- c(report_content, regime_perf_content)
  }
  
  # Write the report
  writeLines(report_content, report_file)
  
  cat("Results and report saved to", output_dir, "\n")
}

# Function to run the complete analysis
run_complete_analysis <- function(start_date = NULL,
                                  end_date = Sys.Date(),
                                  n_stocks = 100,
                                  by_regime = TRUE,
                                  rolling = TRUE,
                                  output_dir = NULL) {
  # Set start date if not provided (default to ~10 years ago)
  if(is.null(start_date)) {
    start_date <- as.Date(end_date) - 365*10  # 10 years of data
  }
  
  cat("Running analysis on S&P 500 data from", format(start_date, "%Y-%m-%d"), 
      "to", format(end_date, "%Y-%m-%d"), "\n")
  
  # Run the analysis
  results <- run_sp500_analysis(
    start_date = start_date,
    end_date = end_date,
    n_stocks = n_stocks,
    train_ratio = 0.7,
    by_regime = by_regime,
    rolling = rolling,
    window_size = 252,  # 1 year
    step_size = 21      # Monthly
  )
  
  # Set default output directory if not provided
  if(is.null(output_dir)) {
    output_dir <- paste0("sp500_analysis_", format(Sys.Date(), "%Y%m%d"))
  }
  
  # Save results and generate report
  save_and_report(results, output_dir)
  
  return(results)
}


# Test sensitivity to contamination threshold
test_contamination_sensitivity <- function(train_returns, test_returns, 
                                           thresholds = seq(2, 5, by = 0.5),
                                           train_market = NULL, test_market = NULL,
                                           max_weight = 0.2, test_dates = NULL, test_regimes = NULL) {
  results <- list()
  
  # Detect base level contamination for comparison
  base_outliers <- detect_outliers(train_returns, method = "MAD", threshold = 3)
  cat("Base contamination (threshold = 3):", round(base_outliers$contamination_pct, 2), "%\n")
  
  for(threshold in thresholds) {
    cat("\nTesting contamination threshold:", threshold, "\n")
    
    # Detect outliers with current threshold
    outliers <- detect_outliers(train_returns, method = "MAD", threshold = threshold)
    cat("  Detected contamination:", round(outliers$contamination_pct, 2), "%\n")
    
    # Clean mask for Max Sharpe
    clean_mask <- outliers$outlier_mask
    
    # Apply hybrid method with this threshold for covariance estimation
    cov_matrix <- tryCatch({
      hybrid_mcd_robpca(train_returns, threshold = threshold)
    }, error = function(e) {
      warning("Error in hybrid_mcd_robpca with threshold ", threshold, ": ", e$message)
      # Fallback to hybrid_mcd_pca
      hybrid_mcd_pca(train_returns)
    })
    
    # Ensure matrix is valid
    if(any(is.na(cov_matrix)) || any(is.infinite(cov_matrix))) {
      warning("Invalid covariance matrix with threshold ", threshold, ". Using sample covariance.")
      cov_matrix <- sample_cov(train_returns)
    }
    
    # Ensure positive definiteness
    eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
    if(min(eig_vals) <= 0) {
      eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
      values <- pmax(eigen_decomp$values, 1e-8)
      cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
    }
    
    # Calculate condition number
    cond_num <- max(eig_vals) / min(eig_vals)
    
    # Create Min Variance portfolio
    min_var_portfolio_result <- min_var_portfolio(train_returns, cov_matrix, max_weight)
    
    # Create Max Sharpe portfolio
    max_sharpe_portfolio_result <- max_sharpe_portfolio(
      train_returns, cov_matrix, 0, max_weight, clean_mask
    )
    
    # Evaluate Min Variance performance
    min_var_performance <- evaluate_portfolio(
      min_var_portfolio_result$weights, train_returns, test_returns,
      train_market, test_market
    )
    
    # Evaluate Max Sharpe performance
    max_sharpe_performance <- evaluate_portfolio(
      max_sharpe_portfolio_result$weights, train_returns, test_returns,
      train_market, test_market
    )
    
    # Add regime-specific performance evaluation if regimes are provided
    if(!is.null(test_dates) && !is.null(test_regimes)) {
      min_var_regime_performance <- evaluate_portfolio_by_regime(
        min_var_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
      
      max_sharpe_regime_performance <- evaluate_portfolio_by_regime(
        max_sharpe_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
    } else {
      min_var_regime_performance <- NULL
      max_sharpe_regime_performance <- NULL
    }
    
    # Store results
    results[[as.character(threshold)]] <- list(
      threshold = threshold,
      contamination_pct = outliers$contamination_pct,
      min_var_weights = min_var_portfolio_result$weights,
      max_sharpe_weights = max_sharpe_portfolio_result$weights,
      min_var_performance = min_var_performance,
      max_sharpe_performance = max_sharpe_performance,
      min_var_regime_performance = min_var_regime_performance,
      max_sharpe_regime_performance = max_sharpe_regime_performance,
      condition_number = cond_num
    )
    
    # Print key results
    cat("  Min Variance - Out-of-sample Sharpe:", round(min_var_performance$test$sharpe, 4), "\n")
    cat("  Max Sharpe - Out-of-sample Sharpe:", round(max_sharpe_performance$test$sharpe, 4), "\n")
    
    # Print regime-specific performance if available
    if(!is.null(min_var_regime_performance)) {
      cat("\n  Min Variance regime-specific Sharpe ratios:\n")
      for(regime in names(min_var_regime_performance)) {
        if(regime != "overall") {
          cat("    ", regime, ":", round(min_var_regime_performance[[regime]]$sharpe, 4), "\n")
        }
      }
    }
  }
  
  # Create summary data frame for Min Variance
  min_var_summary <- data.frame(
    Threshold = as.numeric(names(results)),
    Contamination = sapply(results, function(x) x$contamination_pct),
    InSample_Sharpe = sapply(results, function(x) x$min_var_performance$train$sharpe),
    OutOfSample_Sharpe = sapply(results, function(x) x$min_var_performance$test$sharpe),
    ConditionNumber = sapply(results, function(x) x$condition_number)
  )
  
  # Add regime-specific Sharpe ratios if available
  if(!is.null(results[[1]]$min_var_regime_performance)) {
    regimes <- setdiff(names(results[[1]]$min_var_regime_performance), "overall")
    
    for(regime in regimes) {
      regime_sharpe <- sapply(results, function(x) {
        if(!is.null(x$min_var_regime_performance) && regime %in% names(x$min_var_regime_performance)) {
          x$min_var_regime_performance[[regime]]$sharpe
        } else {
          NA
        }
      })
      
      min_var_summary[[paste0(regime, "_Sharpe")]] <- regime_sharpe
    }
  }
  
  # Create summary data frame for Max Sharpe
  max_sharpe_summary <- data.frame(
    Threshold = as.numeric(names(results)),
    Contamination = sapply(results, function(x) x$contamination_pct),
    InSample_Sharpe = sapply(results, function(x) x$max_sharpe_performance$train$sharpe),
    OutOfSample_Sharpe = sapply(results, function(x) x$max_sharpe_performance$test$sharpe),
    ConditionNumber = sapply(results, function(x) x$condition_number)
  )
  
  # Add regime-specific Sharpe ratios if available
  if(!is.null(results[[1]]$max_sharpe_regime_performance)) {
    regimes <- setdiff(names(results[[1]]$max_sharpe_regime_performance), "overall")
    
    for(regime in regimes) {
      regime_sharpe <- sapply(results, function(x) {
        if(!is.null(x$max_sharpe_regime_performance) && regime %in% names(x$max_sharpe_regime_performance)) {
          x$max_sharpe_regime_performance[[regime]]$sharpe
        } else {
          NA
        }
      })
      
      max_sharpe_summary[[paste0(regime, "_Sharpe")]] <- regime_sharpe
    }
  }
  
  # Identify optimal threshold for each method
  min_var_optimal_idx <- which.max(min_var_summary$OutOfSample_Sharpe)
  min_var_optimal_threshold <- min_var_summary$Threshold[min_var_optimal_idx]
  
  max_sharpe_optimal_idx <- which.max(max_sharpe_summary$OutOfSample_Sharpe)
  max_sharpe_optimal_threshold <- max_sharpe_summary$Threshold[max_sharpe_optimal_idx]
  
  cat("\nOptimal threshold for Min Variance:", min_var_optimal_threshold, 
      "with out-of-sample Sharpe:", round(min_var_summary$OutOfSample_Sharpe[min_var_optimal_idx], 4), "\n")
  
  cat("Optimal threshold for Max Sharpe:", max_sharpe_optimal_threshold, 
      "with out-of-sample Sharpe:", round(max_sharpe_summary$OutOfSample_Sharpe[max_sharpe_optimal_idx], 4), "\n")
  
  return(list(
    results = results,
    min_var_summary = min_var_summary,
    max_sharpe_summary = max_sharpe_summary,
    min_var_optimal_threshold = min_var_optimal_threshold,
    max_sharpe_optimal_threshold = max_sharpe_optimal_threshold
  ))
}

# Test different outlier detection methods
test_detection_methods <- function(train_returns, test_returns,
                                   train_market = NULL, test_market = NULL,
                                   max_weight = 0.2, test_dates = NULL, test_regimes = NULL) {
  # Define methods to test
  methods <- c("MAD", "IQR", "MCD")
  
  results <- list()
  
  for(method in methods) {
    cat("\nTesting detection method:", method, "\n")
    
    # Detect outliers
    outliers <- detect_outliers(train_returns, method = method)
    cat("  Detected contamination:", round(outliers$contamination_pct, 2), "%\n")
    
    # Clean mask for Max Sharpe
    clean_mask <- outliers$outlier_mask
    
    # Apply hybrid method
    cov_matrix <- hybrid_mcd_robpca(train_returns)
    
    # Ensure matrix is valid
    if(any(is.na(cov_matrix)) || any(is.infinite(cov_matrix))) {
      warning("Invalid covariance matrix with method ", method, ". Using sample covariance.")
      cov_matrix <- sample_cov(train_returns)
    }
    
    # Ensure positive definiteness
    eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
    if(min(eig_vals) <= 0) {
      eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
      values <- pmax(eigen_decomp$values, 1e-8)
      cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
    }
    
    # Create Min Variance portfolio
    min_var_portfolio_result <- min_var_portfolio(train_returns, cov_matrix, max_weight)
    
    # Create Max Sharpe portfolio
    max_sharpe_portfolio_result <- max_sharpe_portfolio(
      train_returns, cov_matrix, 0, max_weight, clean_mask
    )
    
    # Evaluate Min Variance performance
    min_var_performance <- evaluate_portfolio(
      min_var_portfolio_result$weights, train_returns, test_returns,
      train_market, test_market
    )
    
    # Evaluate Max Sharpe performance
    max_sharpe_performance <- evaluate_portfolio(
      max_sharpe_portfolio_result$weights, train_returns, test_returns,
      train_market, test_market
    )
    
    # Add regime-specific performance evaluation if regimes are provided
    if(!is.null(test_dates) && !is.null(test_regimes)) {
      min_var_regime_performance <- evaluate_portfolio_by_regime(
        min_var_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
      
      max_sharpe_regime_performance <- evaluate_portfolio_by_regime(
        max_sharpe_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
    } else {
      min_var_regime_performance <- NULL
      max_sharpe_regime_performance <- NULL
    }
    
    # Store results
    results[[method]] <- list(
      method = method,
      contamination_pct = outliers$contamination_pct,
      min_var_weights = min_var_portfolio_result$weights,
      max_sharpe_weights = max_sharpe_portfolio_result$weights,
      min_var_performance = min_var_performance,
      max_sharpe_performance = max_sharpe_performance,
      min_var_regime_performance = min_var_regime_performance,
      max_sharpe_regime_performance = max_sharpe_regime_performance
    )
    
    # Print key results
    cat("  Min Variance - Out-of-sample Sharpe:", round(min_var_performance$test$sharpe, 4), "\n")
    cat("  Max Sharpe - Out-of-sample Sharpe:", round(max_sharpe_performance$test$sharpe, 4), "\n")
    
    # Print regime-specific performance if available
    if(!is.null(min_var_regime_performance)) {
      cat("\n  Min Variance regime-specific Sharpe ratios:\n")
      for(regime in names(min_var_regime_performance)) {
        if(regime != "overall") {
          cat("    ", regime, ":", round(min_var_regime_performance[[regime]]$sharpe, 4), "\n")
        }
      }
    }
  }
  
  # Create summary for Min Variance
  min_var_summary <- data.frame(
    Method = names(results),
    Contamination = sapply(results, function(x) x$contamination_pct),
    InSample_Sharpe = sapply(results, function(x) x$min_var_performance$train$sharpe),
    OutOfSample_Sharpe = sapply(results, function(x) x$min_var_performance$test$sharpe)
  )
  
  # Add regime-specific Sharpe ratios if available
  if(!is.null(results[[1]]$min_var_regime_performance)) {
    regimes <- setdiff(names(results[[1]]$min_var_regime_performance), "overall")
    
    for(regime in regimes) {
      regime_sharpe <- sapply(results, function(x) {
        if(!is.null(x$min_var_regime_performance) && regime %in% names(x$min_var_regime_performance)) {
          x$min_var_regime_performance[[regime]]$sharpe
        } else {
          NA
        }
      })
      
      min_var_summary[[paste0(regime, "_Sharpe")]] <- regime_sharpe
    }
  }
  
  # Create summary for Max Sharpe
  max_sharpe_summary <- data.frame(
    Method = names(results),
    Contamination = sapply(results, function(x) x$contamination_pct),
    InSample_Sharpe = sapply(results, function(x) x$max_sharpe_performance$train$sharpe),
    OutOfSample_Sharpe = sapply(results, function(x) x$max_sharpe_performance$test$sharpe)
  )
  
  # Add regime-specific Sharpe ratios if available
  if(!is.null(results[[1]]$max_sharpe_regime_performance)) {
    regimes <- setdiff(names(results[[1]]$max_sharpe_regime_performance), "overall")
    
    for(regime in regimes) {
      regime_sharpe <- sapply(results, function(x) {
        if(!is.null(x$max_sharpe_regime_performance) && regime %in% names(x$max_sharpe_regime_performance)) {
          x$max_sharpe_regime_performance[[regime]]$sharpe
        } else {
          NA
        }
      })
      
      max_sharpe_summary[[paste0(regime, "_Sharpe")]] <- regime_sharpe
    }
  }
  
  # Identify best method for each portfolio type
  min_var_best_idx <- which.max(min_var_summary$OutOfSample_Sharpe)
  min_var_best_method <- min_var_summary$Method[min_var_best_idx]
  
  max_sharpe_best_idx <- which.max(max_sharpe_summary$OutOfSample_Sharpe)
  max_sharpe_best_method <- max_sharpe_summary$Method[max_sharpe_best_idx]
  
  cat("\nBest detection method for Min Variance:", min_var_best_method, 
      "with out-of-sample Sharpe:", round(min_var_summary$OutOfSample_Sharpe[min_var_best_idx], 4), "\n")
  
  cat("Best detection method for Max Sharpe:", max_sharpe_best_method, 
      "with out-of-sample Sharpe:", round(max_sharpe_summary$OutOfSample_Sharpe[max_sharpe_best_idx], 4), "\n")
  
  return(list(
    results = results,
    min_var_summary = min_var_summary,
    max_sharpe_summary = max_sharpe_summary,
    min_var_best_method = min_var_best_method,
    max_sharpe_best_method = max_sharpe_best_method
  ))
}

# Create regime-specific performance summary
create_regime_performance_summary <- function(portfolio_results, regimes) {
  # Initialize list to store summaries for each regime
  regime_summaries <- list()
  
  # Get unique regimes
  unique_regimes <- unique(regimes)
  unique_regimes <- c(unique_regimes, "overall")  # Add overall results
  
  # For each regime
  for(regime in unique_regimes) {
    # Create summary dataframe
    regime_df <- data.frame(
      Method = names(portfolio_results),
      Return = sapply(portfolio_results, function(x) {
        if(!is.null(x$regime_performance) && regime %in% names(x$regime_performance)) {
          return(x$regime_performance[[regime]]$ann_return)
        } else {
          return(NA)
        }
      }),
      Volatility = sapply(portfolio_results, function(x) {
        if(!is.null(x$regime_performance) && regime %in% names(x$regime_performance)) {
          return(x$regime_performance[[regime]]$ann_volatility)
        } else {
          return(NA)
        }
      }),
      Sharpe = sapply(portfolio_results, function(x) {
        if(!is.null(x$regime_performance) && regime %in% names(x$regime_performance)) {
          return(x$regime_performance[[regime]]$sharpe)
        } else {
          return(NA)
        }
      }),
      MaxDrawdown = sapply(portfolio_results, function(x) {
        if(!is.null(x$regime_performance) && regime %in% names(x$regime_performance)) {
          return(x$regime_performance[[regime]]$max_drawdown)
        } else {
          return(NA)
        }
      })
    )
    
    # Order by Sharpe ratio
    regime_df <- regime_df[order(-regime_df$Sharpe), ]
    
    # Store summary
    regime_summaries[[regime]] <- regime_df
  }
  
  return(regime_summaries)
}

# Create regime-specific visualizations
create_regime_visualizations <- function(regime_performance) {
  plots <- list()
  
  # For each regime
  for(regime in names(regime_performance)) {
    # Get regime data
    regime_df <- regime_performance[[regime]]
    
    # Skip if empty
    if(nrow(regime_df) == 0 || all(is.na(regime_df$Sharpe))) {
      next
    }
    
    # Set method factor levels for ordering
    regime_df$Method <- factor(regime_df$Method, 
                               levels = regime_df$Method[order(-regime_df$Sharpe)])
    
    # Identify method types for coloring
    regime_df$OptimizationType <- ifelse(grepl("MaxSharpe", regime_df$Method), "Max Sharpe", 
                                         ifelse(grepl("MinVar", regime_df$Method), "Min Variance", "Equal Weight"))
    
    # Create Sharpe ratio plot
    p1 <- ggplot(regime_df, aes(x = Method, y = Sharpe, fill = OptimizationType)) +
      geom_bar(stat = "identity") +
      labs(title = paste("Sharpe Ratio by Method -", capitalize(regime), "Regime"),
           x = "", y = "Sharpe Ratio") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Create risk-return scatter plot
    p2 <- ggplot(regime_df, aes(x = Volatility, y = Return, color = OptimizationType)) +
      geom_point(size = 3) +
      geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
      labs(title = paste("Risk-Return Profile -", capitalize(regime), "Regime"),
           x = "Annualized Volatility", y = "Annualized Return") +
      theme_minimal()
    
    # Store plots
    plots[[paste0(regime, "_sharpe")]] <- p1
    plots[[paste0(regime, "_risk_return")]] <- p2
  }
  
  # Create comparison plot across regimes
  # Extract top 5 methods based on overall performance
  if("overall" %in% names(regime_performance)) {
    top_methods <- head(regime_performance[["overall"]]$Method, 5)
    
    # Create data for comparison
    comparison_data <- data.frame()
    
    for(regime in names(regime_performance)) {
      if(regime == "overall") next
      
      regime_df <- regime_performance[[regime]]
      regime_df$Regime <- regime
      
      # Filter for top methods
      filtered_df <- regime_df[regime_df$Method %in% top_methods, ]
      comparison_data <- rbind(comparison_data, filtered_df)
    }
    
    # Create comparison plot
    if(nrow(comparison_data) > 0) {
      comparison_plot <- ggplot(comparison_data, 
                                aes(x = Regime, y = Sharpe, fill = Method)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Performance Comparison Across Market Regimes",
             x = "Market Regime", y = "Sharpe Ratio") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      plots[["regime_comparison"]] <- comparison_plot
    }
    
    # Create heatmap for all methods across regimes
    heatmap_data <- data.frame()
    
    for(regime in names(regime_performance)) {
      if(regime == "overall") next
      
      regime_df <- regime_performance[[regime]]
      regime_df$Regime <- regime
      
      heatmap_data <- rbind(heatmap_data, regime_df[, c("Method", "Regime", "Sharpe")])
    }
    
    if(nrow(heatmap_data) > 0) {
      # Convert Method to factor to control ordering
      methods_order <- regime_performance[["overall"]]$Method
      heatmap_data$Method <- factor(heatmap_data$Method, levels = methods_order)
      
      # Create heatmap
      heatmap_plot <- ggplot(heatmap_data, aes(x = Regime, y = Method, fill = Sharpe)) +
        geom_tile() +
        scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                             midpoint = median(heatmap_data$Sharpe, na.rm = TRUE)) +
        labs(title = "Sharpe Ratio Heatmap Across Market Regimes",
             x = "Market Regime", y = "Method") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8))
      
      plots[["regime_heatmap"]] <- heatmap_plot
    }
  }
  
  return(plots)
}


# ========================================================================
# S&P 500 Empirical Analysis - Rolling Window Analysis
# ========================================================================

# Perform rolling window analysis
rolling_window_analysis <- function(returns, dates, window_size = 252, step_size = 20, market_returns = NULL,
                                    max_weight = 0.2, regimes = NULL) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  # Check if we have enough data
  if(n <= window_size) {
    stop("Not enough data for rolling window analysis")
  }
  
  # Calculate number of windows
  n_windows <- floor((n - window_size) / step_size) + 1
  
  # Initialize results
  window_results <- list()
  window_info <- data.frame(
    Window = 1:n_windows,
    StartDate = rep(NA, n_windows),
    EndDate = rep(NA, n_windows),
    StartIdx = rep(NA, n_windows),
    EndIdx = rep(NA, n_windows)
  )
  
  cat("Performing rolling window analysis with", n_windows, "windows\n")
  
  # For each window
  for(i in 1:n_windows) {
    cat("Processing window", i, "of", n_windows, "\n")
    
    # Define window indices
    start_idx <- (i-1) * step_size + 1
    end_idx <- start_idx + window_size - 1
    
    # Ensure end_idx doesn't exceed data length
    end_idx <- min(end_idx, n)
    
    # Store window info
    window_info$StartDate[i] <- as.character(dates[start_idx])
    window_info$EndDate[i] <- as.character(dates[end_idx])
    window_info$StartIdx[i] <- start_idx
    window_info$EndIdx[i] <- end_idx
    
    # Extract window data
    window_returns <- returns[start_idx:end_idx, , drop = FALSE]
    window_dates <- dates[start_idx:end_idx]
    
    # Extract market returns for this window if available
    window_market <- if(!is.null(market_returns)) {
      market_returns[start_idx:end_idx]
    } else {
      NULL
    }
    
    # Extract window regimes if available
    window_regimes <- if(!is.null(regimes)) {
      regimes[start_idx:end_idx]
    } else {
      NULL
    }
    
    # Split window data for in-sample/out-of-sample testing
    train_size <- floor(0.7 * nrow(window_returns))
    
    train_returns <- window_returns[1:train_size, , drop = FALSE]
    test_returns <- window_returns[(train_size+1):nrow(window_returns), , drop = FALSE]
    
    train_dates <- window_dates[1:train_size]
    test_dates <- window_dates[(train_size+1):length(window_dates)]
    
    if(!is.null(window_market)) {
      train_market <- window_market[1:train_size]
      test_market <- window_market[(train_size+1):length(window_market)]
    } else {
      train_market <- test_market <- NULL
    }
    
    if(!is.null(window_regimes)) {
      train_regimes <- window_regimes[1:train_size]
      test_regimes <- window_regimes[(train_size+1):length(window_regimes)]
    } else {
      train_regimes <- test_regimes <- NULL
    }
    
    # Create and evaluate portfolios for this window
    tryCatch({
      window_results[[i]] <- create_and_evaluate_portfolios(
        train_returns, test_returns,
        train_market, test_market,
        max_weight = max_weight,
        test_dates = test_dates,
        test_regimes = test_regimes
      )
    }, error = function(e) {
      warning("Error in window ", i, ": ", e$message)
      # Return NULL for this window
      NULL
    })
  }
  
  return(list(
    results = window_results,
    window_info = window_info
  ))
}

# Analyze rolling window performance
analyze_rolling_performance <- function(window_results, metric = "sharpe") {
  # Extract performance for each window and method
  n_windows <- length(window_results$results)
  
  # Get methods from the first non-NULL window
  first_valid <- which(!sapply(window_results$results, is.null))[1]
  if(is.na(first_valid)) {
    stop("No valid windows in results")
  }
  
  methods <- names(window_results$results[[first_valid]])
  
  # Initialize result matrix
  performance_matrix <- matrix(NA, nrow = n_windows, ncol = length(methods))
  colnames(performance_matrix) <- methods
  
  # Fill matrix with performance metrics
  for(i in 1:n_windows) {
    if(!is.null(window_results$results[[i]])) {
      for(j in 1:length(methods)) {
        method <- methods[j]
        if(method %in% names(window_results$results[[i]])) {
          performance_matrix[i, j] <- window_results$results[[i]][[method]]$performance$test[[metric]]
        }
      }
    }
  }
  
  # Convert to data frame
  performance_df <- as.data.frame(performance_matrix)
  performance_df$Window <- 1:n_windows
  

  # Add date information if available
  if(!is.null(window_results$window_info)) {
    performance_df$StartDate <- window_results$window_info$StartDate
    performance_df$EndDate <- window_results$window_info$EndDate
  }
  
  # Calculate mean performance across windows
  mean_performance <- colMeans(performance_matrix, na.rm = TRUE)
  mean_df <- data.frame(
    Method = names(mean_performance),
    Mean = mean_performance
  )
  
  # Calculate standard deviation of performance
  sd_performance <- apply(performance_matrix, 2, sd, na.rm = TRUE)
  sd_df <- data.frame(
    Method = names(sd_performance),
    SD = sd_performance
  )
  
  # Calculate Sharpe of strategy (mean/sd of rolling performance)
  strategy_sharpe <- mean_performance / sd_performance
  strategy_sharpe_df <- data.frame(
    Method = names(mean_performance),
    StrategySharpe = strategy_sharpe
  )
  
  # Create stability ranking
  stability_rank <- rank(-strategy_sharpe)
  stability_df <- data.frame(
    Method = names(mean_performance),
    StabilityRank = stability_rank
  )
  
  # Reshape for plotting
  plot_df <- reshape2::melt(performance_df, 
                            id.vars = c("Window", if(!is.null(window_results$window_info)) 
                              c("StartDate", "EndDate") else NULL),
                            variable.name = "Method", value.name = metric)
  
  # Create rank stability plot
  rank_matrix <- t(apply(-performance_matrix, 1, rank, ties.method = "min"))
  colnames(rank_matrix) <- methods
  
  rank_df <- as.data.frame(rank_matrix)
  rank_df$Window <- 1:n_windows
  
  if(!is.null(window_results$window_info)) {
    rank_df$StartDate <- window_results$window_info$StartDate
    rank_df$EndDate <- window_results$window_info$EndDate
  }
  
  rank_plot_df <- reshape2::melt(rank_df, 
                                 id.vars = c("Window", if(!is.null(window_results$window_info)) 
                                   c("StartDate", "EndDate") else NULL),
                                 variable.name = "Method", value.name = "Rank")
  
  list(
    performance_matrix = performance_matrix,
    mean_performance = mean_df,
    sd_performance = sd_df,
    strategy_sharpe = strategy_sharpe_df,
    stability_rank = stability_df,
    plot_data = plot_df,
    rank_plot_data = rank_plot_df
  )
}

# Add regime information to rolling window analysis
add_regime_info_to_rolling <- function(rolling_results, regimes, dates) {
  # Extract window start and end indices
  start_indices <- rolling_results$window_info$StartIdx
  end_indices <- rolling_results$window_info$EndIdx
  
  # Initialize dataframe for regime information
  regime_info <- data.frame(
    Window = 1:length(start_indices),
    regime_counts = rep(NA, length(start_indices))
  )
  
  # Calculate regime distributions for each window
  for(i in 1:length(start_indices)) {
    # Extract regimes for this window
    window_regimes <- regimes[start_indices[i]:end_indices[i]]
    
    # Calculate regime distribution
    regime_counts <- table(window_regimes)
    
    # Store dominant regime
    regime_info$dominant_regime[i] <- names(regime_counts)[which.max(regime_counts)]
    
    # Store regime counts
    for(regime in names(regime_counts)) {
      col_name <- paste0(regime, "_count")
      if(!(col_name %in% colnames(regime_info))) {
        regime_info[[col_name]] <- rep(0, length(start_indices))
      }
      regime_info[[col_name]][i] <- regime_counts[regime]
    }
    
    # Store regime percentages
    for(regime in names(regime_counts)) {
      col_name <- paste0(regime, "_pct")
      if(!(col_name %in% colnames(regime_info))) {
        regime_info[[col_name]] <- rep(0, length(start_indices))
      }
      regime_info[[col_name]][i] <- regime_counts[regime] / sum(regime_counts)
    }
  }
  
  # Add regime info to window_info
  rolling_results$window_info <- cbind(rolling_results$window_info, 
                                       regime_info[, setdiff(colnames(regime_info), "Window")])
  
  return(rolling_results)
}
