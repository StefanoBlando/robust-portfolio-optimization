# ========================================================================
# Statistical Validation Framework for Robust Portfolio Methods
# ========================================================================

run_statistical_tests <- function(standard_results, stress_results, 
                                  test_types = c("performance", "normality", "stability",
                                                 "robustness", "numerical", "regime"),
                                  bootstrap_replicates = 2000,
                                  confidence_level = 0.95,
                                  output_dir = "statistical_analysis") {
  
  cat("Running comprehensive statistical analysis of portfolio methods...\n")
  
  # Verify that required packages are installed and load them
  required_packages <- c("tseries", "boot", "car", "robustbase", "moments", "ggplot2", 
                         "reshape2", "dplyr", "corrplot")
  
  missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
  if(length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse=", "), "\n")
    install.packages(missing_packages)
  }
  
  # Load required packages
  invisible(sapply(required_packages, library, character.only = TRUE))
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize container for all analysis results
  analysis_results <- list()
  
  # ========================================================================
  # Extract and prepare data
  # ========================================================================
  
  cat("Extracting portfolio data...\n")
  
  # Extract method names from both analyses
  std_methods <- if(!is.null(standard_results$portfolio_results)) {
    names(standard_results$portfolio_results)
  } else {
    character(0)
  }
  
  stress_methods <- if(!is.null(stress_results$portfolio_results)) {
    names(stress_results$portfolio_results)
  } else {
    character(0)
  }
  
  all_methods <- unique(c(std_methods, stress_methods))
  common_methods <- intersect(std_methods, stress_methods)
  
  # Create data frame for methods with corresponding optimization and covariance types
  methods_info <- data.frame(
    Method = all_methods,
    stringsAsFactors = FALSE
  )
  
  # Extract optimization and covariance method types
  methods_info$OptimizationType <- sapply(methods_info$Method, function(m) {
    if(grepl("MaxSharpe", m)) {
      return("MaxSharpe")
    } else if(grepl("MinVar", m)) {
      return("MinVar")
    } else if(grepl("EqualWeights", m)) {
      return("EqualWeights")
    } else {
      return("Other")
    }
  })
  
  methods_info$CovMethod <- sapply(strsplit(methods_info$Method, "_"), function(x) {
    if(length(x) > 0) x[1] else "Unknown"
  })
  
  # Create data structure for portfolio returns
  portfolio_returns <- list(
    standard = list(),
    stress = list()
  )
  
  # Extract portfolio returns
  for(method in std_methods) {
    if(!is.null(standard_results$portfolio_results[[method]]$performance$test_returns)) {
      portfolio_returns$standard[[method]] <- 
        standard_results$portfolio_results[[method]]$performance$test_returns
    }
  }
  
  for(method in stress_methods) {
    if(!is.null(stress_results$portfolio_results[[method]]$performance$test_returns)) {
      portfolio_returns$stress[[method]] <- 
        stress_results$portfolio_results[[method]]$performance$test_returns
    }
  }
  
  # Create data structure for portfolio weights
  portfolio_weights <- list(
    standard = list(),
    stress = list()
  )
  
  # Extract portfolio weights
  for(method in std_methods) {
    if(!is.null(standard_results$portfolio_results[[method]]$weights)) {
      portfolio_weights$standard[[method]] <- 
        standard_results$portfolio_results[[method]]$weights
    }
  }
  
  for(method in stress_methods) {
    if(!is.null(stress_results$portfolio_results[[method]]$weights)) {
      portfolio_weights$stress[[method]] <- 
        stress_results$portfolio_results[[method]]$weights
    }
  }
  
  # Extract performance summaries
  perf_standard <- if(!is.null(standard_results$performance_summary)) {
    standard_results$performance_summary
  } else {
    NULL
  }
  
  perf_stress <- if(!is.null(stress_results$performance_summary)) {
    stress_results$performance_summary
  } else {
    NULL
  }
  
  # Extract regime information if available
  regime_info <- list(
    standard = if(!is.null(standard_results$regime_performance)) {
      standard_results$regime_performance
    } else {
      NULL
    },
    stress = if(!is.null(stress_results$regime_performance)) {
      stress_results$regime_performance
    } else {
      NULL
    }
  )
  
  # ========================================================================
  # 1. Distributional tests (normality, heavy tails)
  # ========================================================================
  
  if("normality" %in% test_types) {
    cat("Running distributional tests on portfolio returns...\n")
    
    # Initialize results
    normality_tests <- data.frame(
      Method = character(),
      JB_Statistic = numeric(),
      JB_PValue = numeric(),
      SW_Statistic = numeric(),
      SW_PValue = numeric(),
      AD_Statistic = numeric(),
      AD_PValue = numeric(),
      Skewness = numeric(),
      Kurtosis = numeric(),
      IsNormal = logical(),
      Condition = character(),
      stringsAsFactors = FALSE
    )
    
    # Function to run distributional tests
    run_dist_tests <- function(returns, method, condition) {
      # Check if we have enough data
      if(length(returns) >= 20 && !all(is.na(returns))) {
        # Calculate skewness and kurtosis
        skew <- tryCatch(moments::skewness(returns, na.rm = TRUE), error = function(e) NA)
        kurt <- tryCatch(moments::kurtosis(returns, na.rm = TRUE), error = function(e) NA)
        
        # Jarque-Bera test
        jb <- tryCatch({
          test <- tseries::jarque.bera.test(as.numeric(returns))
          list(statistic = as.numeric(test$statistic), p.value = test$p.value)
        }, error = function(e) {
          list(statistic = NA, p.value = NA)
        })
        
        # Shapiro-Wilk test (for smaller samples)
        sw <- tryCatch({
          if(length(returns) <= 5000) {
            test <- shapiro.test(as.numeric(returns))
            list(statistic = test$statistic, p.value = test$p.value)
          } else {
            list(statistic = NA, p.value = NA)
          }
        }, error = function(e) {
          list(statistic = NA, p.value = NA)
        })
        
        # Anderson-Darling test
        ad <- tryCatch({
          test <- nortest::ad.test(as.numeric(returns))
          list(statistic = test$statistic, p.value = test$p.value)
        }, error = function(e) {
          list(statistic = NA, p.value = NA)
        })
        
        # Determine if returns are normal (using multiple tests)
        is_normal <- TRUE
        if(!is.na(jb$p.value) && jb$p.value < 0.05) is_normal <- FALSE
        if(!is.na(sw$p.value) && sw$p.value < 0.05) is_normal <- FALSE
        if(!is.na(ad$p.value) && ad$p.value < 0.05) is_normal <- FALSE
        
        return(data.frame(
          Method = method,
          JB_Statistic = jb$statistic,
          JB_PValue = jb$p.value,
          SW_Statistic = sw$statistic,
          SW_PValue = sw$p.value,
          AD_Statistic = ad$statistic,
          AD_PValue = ad$p.value,
          Skewness = skew,
          Kurtosis = kurt,
          IsNormal = is_normal,
          Condition = condition,
          stringsAsFactors = FALSE
        ))
      } else {
        return(NULL)
      }
    }
    
    # Run tests for standard portfolio returns
    for(method in names(portfolio_returns$standard)) {
      result <- run_dist_tests(portfolio_returns$standard[[method]], method, "Standard")
      if(!is.null(result)) {
        normality_tests <- rbind(normality_tests, result)
      }
    }
    
    # Run tests for stressed portfolio returns
    for(method in names(portfolio_returns$stress)) {
      result <- run_dist_tests(portfolio_returns$stress[[method]], method, "Stress")
      if(!is.null(result)) {
        normality_tests <- rbind(normality_tests, result)
      }
    }
    
    # Add optimization and covariance method types
    if(nrow(normality_tests) > 0) {
      normality_tests$OptimizationType <- sapply(normality_tests$Method, function(m) {
        methods_info$OptimizationType[methods_info$Method == m]
      })
      
      normality_tests$CovMethod <- sapply(normality_tests$Method, function(m) {
        methods_info$CovMethod[methods_info$Method == m]
      })
      
      # Check for heavy tails
      normality_tests$HeavyTails <- normality_tests$Kurtosis > 3
      
      # Calculate statistics for both conditions
      normality_stats <- aggregate(
        cbind(IsNormal, HeavyTails) ~ Condition, 
        data = normality_tests, 
        FUN = function(x) mean(!x, na.rm = TRUE)
      )
      
      # Store results
      analysis_results$normality <- list(
        tests = normality_tests,
        stats = normality_stats
      )
    }
  }
  
  # ========================================================================
  # 2. Bootstrap confidence intervals for performance metrics
  # ========================================================================
  
  if("performance" %in% test_types) {
    cat("Bootstrapping performance metrics...\n")
    
    # Initialize results
    bootstrap_results <- data.frame(
      Method = character(),
      Condition = character(),
      Metric = character(),
      Original = numeric(),
      BootMean = numeric(),
      BootMedian = numeric(),
      BootSD = numeric(),
      CI_Lower = numeric(),
      CI_Upper = numeric(),
      RelPrecision = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Function to bootstrap a specific metric
    bootstrap_metric <- function(returns, method, condition, metric) {
      # Check if we have enough data
      if(length(returns) >= 10 && !all(is.na(returns))) {
        # Define function to calculate the metric
        calc_metric <- switch(metric,
                              "sharpe" = function(data) mean(data, na.rm = TRUE) / sd(data, na.rm = TRUE),
                              "sortino" = function(data) {
                                mean_ret <- mean(data, na.rm = TRUE)
                                downside <- data[data < 0]
                                if(length(downside) > 0) {
                                  downside_dev <- sd(downside, na.rm = TRUE)
                                  if(downside_dev > 0) return(mean_ret / downside_dev)
                                }
                                return(mean_ret / 0.0001)  # Avoid division by zero
                              },
                              "var95" = function(data) quantile(data, 0.05, na.rm = TRUE),
                              "cvar95" = function(data) {
                                var95 <- quantile(data, 0.05, na.rm = TRUE)
                                mean(data[data <= var95], na.rm = TRUE)
                              },
                              "maxdd" = function(data) {
                                cumret <- cumprod(1 + data)
                                peak <- cummax(cumret)
                                drawdown <- (cumret - peak) / peak
                                min(drawdown, na.rm = TRUE)
                              }
        )
        
        # Calculate original metric value
        original_value <- tryCatch(calc_metric(returns), error = function(e) NA)
        
        if(!is.na(original_value) && is.finite(original_value)) {
          # Bootstrap function
          boot_fn <- function(data, indices) {
            boot_sample <- data[indices]
            tryCatch(calc_metric(boot_sample), error = function(e) NA)
          }
          
          # Run bootstrap
          B <- bootstrap_replicates
          tryCatch({
            boot_results <- boot::boot(returns, boot_fn, R = B)
            
            # Check for enough valid results
            valid_boot <- !is.na(boot_results$t) & is.finite(boot_results$t)
            if(sum(valid_boot) < B * 0.8) {
              cat("Warning: Too many invalid bootstrap samples for", method, metric, "\n")
            }
            
            # Calculate bootstrap statistics
            boot_mean <- mean(boot_results$t[valid_boot], na.rm = TRUE)
            boot_median <- median(boot_results$t[valid_boot], na.rm = TRUE)
            boot_sd <- sd(boot_results$t[valid_boot], na.rm = TRUE)
            
            # Calculate percentile confidence intervals
            alpha <- (1 - confidence_level) / 2
            ci_lower <- quantile(boot_results$t[valid_boot], alpha, na.rm = TRUE)
            ci_upper <- quantile(boot_results$t[valid_boot], 1 - alpha, na.rm = TRUE)
            
            # Calculate relative precision (CI width / mean)
            rel_prec <- (ci_upper - ci_lower) / abs(boot_mean)
            if(boot_mean == 0) rel_prec <- NA
            
            return(data.frame(
              Method = method,
              Condition = condition,
              Metric = metric,
              Original = original_value,
              BootMean = boot_mean,
              BootMedian = boot_median,
              BootSD = boot_sd,
              CI_Lower = ci_lower,
              CI_Upper = ci_upper,
              RelPrecision = rel_prec,
              stringsAsFactors = FALSE
            ))
          }, error = function(e) {
            cat("Bootstrap error for", method, metric, ":", e$message, "\n")
            
            # Try manual bootstrap
            boot_samples <- numeric(B)
            for(i in 1:B) {
              idx <- sample(length(returns), replace = TRUE)
              boot_sample <- returns[idx]
              boot_samples[i] <- tryCatch(calc_metric(boot_sample), error = function(e) NA)
            }
            
            valid_boot <- !is.na(boot_samples) & is.finite(boot_samples)
            if(sum(valid_boot) < B * 0.8) {
              cat("Warning: Too many invalid manual bootstrap samples for", method, metric, "\n")
            }
            
            # Calculate bootstrap statistics
            boot_mean <- mean(boot_samples[valid_boot], na.rm = TRUE)
            boot_median <- median(boot_samples[valid_boot], na.rm = TRUE)
            boot_sd <- sd(boot_samples[valid_boot], na.rm = TRUE)
            
            # Calculate percentile confidence intervals
            alpha <- (1 - confidence_level) / 2
            ci_lower <- quantile(boot_samples[valid_boot], alpha, na.rm = TRUE)
            ci_upper <- quantile(boot_samples[valid_boot], 1 - alpha, na.rm = TRUE)
            
            # Calculate relative precision (CI width / mean)
            rel_prec <- (ci_upper - ci_lower) / abs(boot_mean)
            if(boot_mean == 0) rel_prec <- NA
            
            return(data.frame(
              Method = method,
              Condition = condition,
              Metric = metric,
              Original = original_value,
              BootMean = boot_mean,
              BootMedian = boot_median,
              BootSD = boot_sd,
              CI_Lower = ci_lower,
              CI_Upper = ci_upper,
              RelPrecision = rel_prec,
              stringsAsFactors = FALSE
            ))
          })
        } else {
          return(NULL)
        }
      } else {
        return(NULL)
      }
    }
    
    # Define metrics to bootstrap
    metrics <- c("sharpe", "sortino", "var95", "cvar95", "maxdd")
    
    # Bootstrap each metric for standard portfolio returns
    for(method in names(portfolio_returns$standard)) {
      returns <- portfolio_returns$standard[[method]]
      for(metric in metrics) {
        result <- bootstrap_metric(returns, method, "Standard", metric)
        if(!is.null(result)) {
          bootstrap_results <- rbind(bootstrap_results, result)
        }
      }
    }
    
    # Bootstrap each metric for stressed portfolio returns
    for(method in names(portfolio_returns$stress)) {
      returns <- portfolio_returns$stress[[method]]
      for(metric in metrics) {
        result <- bootstrap_metric(returns, method, "Stress", metric)
        if(!is.null(result)) {
          bootstrap_results <- rbind(bootstrap_results, result)
        }
      }
    }
    
    # Add optimization and covariance method types
    if(nrow(bootstrap_results) > 0) {
      bootstrap_results$OptimizationType <- sapply(bootstrap_results$Method, function(m) {
        methods_info$OptimizationType[methods_info$Method == m]
      })
      
      bootstrap_results$CovMethod <- sapply(bootstrap_results$Method, function(m) {
        methods_info$CovMethod[methods_info$Method == m]
      })
      
      # Store results
      analysis_results$bootstrap <- bootstrap_results
    }
  }
  
  # ========================================================================
  # 3. Hypothesis tests for method comparison
  # ========================================================================
  
  if("robustness" %in% test_types) {
    cat("Running hypothesis tests for method comparison...\n")
    
    # Create data frame with performance under both conditions
    comparison_data <- data.frame()
    
    if(!is.null(perf_standard) && !is.null(perf_stress)) {
      # Extract relevant columns from performance summary
      if("OutOfSample_Sharpe" %in% colnames(perf_standard) && 
         "OutOfSample_Sharpe" %in% colnames(perf_stress)) {
        
        std_sharpe <- perf_standard[, c("Method", "OutOfSample_Sharpe")]
        stress_sharpe <- perf_stress[, c("Method", "OutOfSample_Sharpe")]
        
        # Rename columns
        names(std_sharpe) <- c("Method", "Standard_Sharpe")
        names(stress_sharpe) <- c("Method", "Stress_Sharpe")
        
        # Merge data for common methods
        common_methods <- intersect(std_sharpe$Method, stress_sharpe$Method)
        
        if(length(common_methods) > 0) {
          # Filter for common methods
          std_common <- std_sharpe[std_sharpe$Method %in% common_methods, ]
          stress_common <- stress_sharpe[stress_sharpe$Method %in% common_methods, ]
          
          # Sort in the same order
          std_common <- std_common[match(common_methods, std_common$Method), ]
          stress_common <- stress_common[match(common_methods, stress_common$Method), ]
          
          # Calculate absolute and percentage changes
          comparison_data <- data.frame(
            Method = common_methods,
            Standard_Sharpe = std_common$Standard_Sharpe,
            Stress_Sharpe = stress_common$Stress_Sharpe,
            Absolute_Change = stress_common$Stress_Sharpe - std_common$Standard_Sharpe
          )
          
          # Calculate percent change with error handling
          comparison_data$Percent_Change <- ifelse(
            !is.na(comparison_data$Standard_Sharpe) & comparison_data$Standard_Sharpe != 0,
            (comparison_data$Stress_Sharpe / comparison_data$Standard_Sharpe - 1) * 100,
            NA
          )
          
          # Add optimization and covariance method types
          comparison_data$OptimizationType <- sapply(comparison_data$Method, function(m) {
            methods_info$OptimizationType[methods_info$Method == m]
          })
          
          comparison_data$CovMethod <- sapply(comparison_data$Method, function(m) {
            methods_info$CovMethod[methods_info$Method == m]
          })
          
          # Sort by percent change (less negative = more robust)
          comparison_data <- comparison_data[order(-comparison_data$Percent_Change), ]
        }
      }
    }
    
    # Run statistical tests on the comparison data
    if(nrow(comparison_data) > 0) {
      # Paired t-test for Sharpe differences
      paired_ttest <- NULL
      ttest_valid <- !is.na(comparison_data$Standard_Sharpe) & 
        !is.na(comparison_data$Stress_Sharpe) &
        is.finite(comparison_data$Standard_Sharpe) & 
        is.finite(comparison_data$Stress_Sharpe)
      
      if(sum(ttest_valid) >= 5) {  # Need at least 5 pairs for a reasonable test
        paired_ttest <- tryCatch({
          t.test(comparison_data$Standard_Sharpe[ttest_valid], 
                 comparison_data$Stress_Sharpe[ttest_valid], 
                 paired = TRUE)
        }, error = function(e) {
          cat("Error in paired t-test:", e$message, "\n")
          NULL
        })
      }
      
      # Run Wilcoxon signed-rank test (non-parametric alternative)
      wilcox_test <- NULL
      if(sum(ttest_valid) >= 5) {
        wilcox_test <- tryCatch({
          wilcox.test(comparison_data$Standard_Sharpe[ttest_valid], 
                      comparison_data$Stress_Sharpe[ttest_valid], 
                      paired = TRUE)
        }, error = function(e) {
          cat("Error in Wilcoxon test:", e$message, "\n")
          NULL
        })
      }
      
      # Compare performance by optimization type
      opt_type_comparison <- NULL
      if("OptimizationType" %in% colnames(comparison_data)) {
        opt_types <- unique(comparison_data$OptimizationType[!is.na(comparison_data$OptimizationType)])
        if(length(opt_types) > 1) {
          opt_type_comparison <- data.frame(
            OptimizationType = character(),
            Mean_Percent_Change = numeric(),
            SD_Percent_Change = numeric(),
            Count = numeric(),
            stringsAsFactors = FALSE
          )
          
          for(opt in opt_types) {
            subset_data <- comparison_data[comparison_data$OptimizationType == opt, ]
            if(nrow(subset_data) > 0) {
              opt_type_comparison <- rbind(opt_type_comparison, data.frame(
                OptimizationType = opt,
                Mean_Percent_Change = mean(subset_data$Percent_Change, na.rm = TRUE),
                SD_Percent_Change = sd(subset_data$Percent_Change, na.rm = TRUE),
                Count = nrow(subset_data)
              ))
            }
          }
          
          # Sort by mean percent change
          opt_type_comparison <- opt_type_comparison[order(-opt_type_comparison$Mean_Percent_Change), ]
          
          # Run ANOVA to test for significant differences
          anova_result <- NULL
          if(length(opt_types) > 1 && all(sapply(opt_types, function(opt) {
            sum(!is.na(comparison_data$Percent_Change[comparison_data$OptimizationType == opt])) >= 3
          }))) {
            anova_result <- tryCatch({
              aov(Percent_Change ~ OptimizationType, data = comparison_data)
            }, error = function(e) {
              cat("Error in ANOVA:", e$message, "\n")
              NULL
            })
          }
        }
      }
      
      # Compare performance by covariance method
      cov_method_comparison <- NULL
      if("CovMethod" %in% colnames(comparison_data)) {
        cov_methods <- unique(comparison_data$CovMethod[!is.na(comparison_data$CovMethod)])
        if(length(cov_methods) > 1) {
          cov_method_comparison <- data.frame(
            CovMethod = character(),
            Mean_Percent_Change = numeric(),
            SD_Percent_Change = numeric(),
            Count = numeric(),
            stringsAsFactors = FALSE
          )
          
          for(cov in cov_methods) {
            subset_data <- comparison_data[comparison_data$CovMethod == cov, ]
            if(nrow(subset_data) > 0) {
              cov_method_comparison <- rbind(cov_method_comparison, data.frame(
                CovMethod = cov,
                Mean_Percent_Change = mean(subset_data$Percent_Change, na.rm = TRUE),
                SD_Percent_Change = sd(subset_data$Percent_Change, na.rm = TRUE),
                Count = nrow(subset_data)
              ))
            }
          }
          
          # Sort by mean percent change
          cov_method_comparison <- cov_method_comparison[order(-cov_method_comparison$Mean_Percent_Change), ]
          
          # Run ANOVA to test for significant differences
          anova_result_cov <- NULL
          if(length(cov_methods) > 1 && all(sapply(cov_methods, function(cov) {
            sum(!is.na(comparison_data$Percent_Change[comparison_data$CovMethod == cov])) >= 3
          }))) {
            anova_result_cov <- tryCatch({
              aov(Percent_Change ~ CovMethod, data = comparison_data)
            }, error = function(e) {
              cat("Error in ANOVA:", e$message, "\n")
              NULL
            })
          }
        }
      }
      
      # Store results
      analysis_results$comparison <- list(
        data = comparison_data,
        paired_ttest = paired_ttest,
        wilcox_test = wilcox_test,
        opt_type_comparison = opt_type_comparison,
        cov_method_comparison = cov_method_comparison
      )
    }
  }
  
  # ========================================================================
  # 4. Weight stability analysis
  # ========================================================================
  
  if("stability" %in% test_types) {
    cat("Analyzing portfolio weight stability...\n")
    
    # Initialize weight stability dataframe
    weight_stability <- data.frame(
      Method = character(),
      Turnover = numeric(),
      MaxWeightChange = numeric(),
      WeightCorrelation = numeric(),
      WMAD = numeric(),  # Weighted Mean Absolute Deviation
      stringsAsFactors = FALSE
    )
    
    # Compare weights between standard and stress conditions
    for(method in common_methods) {
      if(method %in% names(portfolio_weights$standard) && 
         method %in% names(portfolio_weights$stress)) {
        
        std_weights <- portfolio_weights$standard[[method]]
        stress_weights <- portfolio_weights$stress[[method]]
        
        # Check if weights have the same length
        if(length(std_weights) == length(stress_weights)) {
          # Calculate turnover (sum of absolute changes divided by 2)
          turnover <- sum(abs(std_weights - stress_weights)) / 2
          
          # Calculate max weight change
          max_weight_change <- max(abs(std_weights - stress_weights))
          
          # Calculate correlation between weights
          weight_corr <- tryCatch({
            cor(std_weights, stress_weights)
          }, error = function(e) {
            cat("Error calculating weight correlation for", method, ":", e$message, "\n")
            NA
          })
          
          # Calculate weighted mean absolute deviation (weight changes weighted by original weights)
          wmad <- tryCatch({
            sum(abs(std_weights - stress_weights) * std_weights) / sum(std_weights)
          }, error = function(e) {
            cat("Error calculating WMAD for", method, ":", e$message, "\n")
            NA
          })
          
          # Add to dataframe
          weight_stability <- rbind(weight_stability, data.frame(
            Method = method,
            Turnover = turnover,
            MaxWeightChange = max_weight_change,
            WeightCorrelation = weight_corr,
            WMAD = wmad
          ))
        } else {
          cat("Weights have different lengths for method:", method, "\n")
        }
      }
    }
    
    # Add optimization and covariance method types
    if(nrow(weight_stability) > 0) {
      weight_stability$OptimizationType <- sapply(weight_stability$Method, function(m) {
        methods_info$OptimizationType[methods_info$Method == m]
      })
      
      weight_stability$CovMethod <- sapply(weight_stability$Method, function(m) {
        methods_info$CovMethod[methods_info$Method == m]
      })
      
      # Calculate stability score (higher = more stable)
      # Combine turnover, correlation, and max change into a single score
      weight_stability$StabilityScore <- with(weight_stability, {
        # Rescale metrics to 0-1 (1 = most stable)
        turnover_scaled <- 1 - (Turnover / max(Turnover, na.rm = TRUE))
        maxchange_scaled <- 1 - (MaxWeightChange / max(MaxWeightChange, na.rm = TRUE))
        correlation_scaled <- (WeightCorrelation + 1) / 2  # Scale from [-1,1] to [0,1]
        
        # Combine metrics (simple average)
        (turnover_scaled + maxchange_scaled + correlation_scaled) / 3
      })
      
      # Sort by stability score
      weight_stability <- weight_stability[order(-weight_stability$StabilityScore), ]
      
      # Calculate average stability by optimization type
      opt_stability <- NULL
      if(length(unique(weight_stability$OptimizationType)) > 1) {
        opt_stability <- aggregate(
          cbind(Turnover, StabilityScore) ~ OptimizationType, 
          data = weight_stability, 
          FUN = mean, 
          na.rm = TRUE
        )
        
        # Sort by stability score
        opt_stability <- opt_stability[order(-opt_stability$StabilityScore), ]
      }
      
      # Calculate average stability by covariance method
      cov_stability <- NULL
      if(length(unique(weight_stability$CovMethod)) > 1) {
        cov_stability <- aggregate(
          cbind(Turnover, StabilityScore) ~ CovMethod, 
          data = weight_stability, 
          FUN = mean, 
          na.rm = TRUE
        )
        
        # Sort by stability score
        cov_stability <- cov_stability[order(-cov_stability$StabilityScore), ]
      }
      
      # Store results
      analysis_results$weight_stability <- list(
        data = weight_stability,
        opt_stability = opt_stability,
        cov_stability = cov_stability
      )
    }
  }
  
  # ========================================================================
  # 5. Numerical stability analysis (condition numbers)
  # ========================================================================
  
  if("numerical" %in% test_types) {
    cat("Analyzing numerical stability (condition numbers)...\n")
    
    # Initialize condition numbers dataframe
    condition_numbers <- data.frame(
      Method = character(),
      ConditionNumber = numeric(),
      LogConditionNumber = numeric(),
      Condition = character(),
      stringsAsFactors = FALSE
    )
    
    # Extract condition numbers from standard results
    if(!is.null(standard_results$portfolio_results)) {
      for(method in names(standard_results$portfolio_results)) {
        cond_num <- NA
        
        # Try to extract directly if available
        if(!is.null(standard_results$portfolio_results[[method]]$condition_number)) {
          cond_num <- standard_results$portfolio_results[[method]]$condition_number
        } 
        # Try to calculate from covariance matrix
        else if(!is.null(standard_results$portfolio_results[[method]]$cov_matrix)) {
          cov_matrix <- standard_results$portfolio_results[[method]]$cov_matrix
          
          # Calculate condition number
          tryCatch({
            eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
            
            # Handle negative or zero eigenvalues
            if(min(eig_vals) <= 0) {
              cat("Warning: Non-positive eigenvalues for", method, "- adjusting\n")
              eig_vals <- eig_vals[eig_vals > 0]
              if(length(eig_vals) < 2) {
                cat("Too few positive eigenvalues for", method, "- skipping\n")
                next
              }
            }
            
            cond_num <- max(eig_vals) / min(eig_vals)
          }, error = function(e) {
            cat("Error in condition number calculation for", method, ":", e$message, "\n")
          })
        }
        
        if(!is.na(cond_num) && is.finite(cond_num)) {
          condition_numbers <- rbind(condition_numbers, data.frame(
            Method = method,
            ConditionNumber = cond_num,
            LogConditionNumber = log10(cond_num),
            Condition = "Standard",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Extract condition numbers from stress results
    if(!is.null(stress_results$portfolio_results)) {
      for(method in names(stress_results$portfolio_results)) {
        cond_num <- NA
        
        # Try to extract directly if available
        if(!is.null(stress_results$portfolio_results[[method]]$condition_number)) {
          cond_num <- stress_results$portfolio_results[[method]]$condition_number
        } 
        # Try to calculate from covariance matrix
        else if(!is.null(stress_results$portfolio_results[[method]]$cov_matrix)) {
          cov_matrix <- stress_results$portfolio_results[[method]]$cov_matrix
          
          # Calculate condition number
          tryCatch({
            eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
            
            # Handle negative or zero eigenvalues
            if(min(eig_vals) <= 0) {
              cat("Warning: Non-positive eigenvalues for", method, "- adjusting\n")
              eig_vals <- eig_vals[eig_vals > 0]
              if(length(eig_vals) < 2) {
                cat("Too few positive eigenvalues for", method, "- skipping\n")
                next
              }
            }
            
            cond_num <- max(eig_vals) / min(eig_vals)
          }, error = function(e) {
            cat("Error in condition number calculation for", method, ":", e$message, "\n")
          })
        }
        
        if(!is.na(cond_num) && is.finite(cond_num)) {
          condition_numbers <- rbind(condition_numbers, data.frame(
            Method = method,
            ConditionNumber = cond_num,
            LogConditionNumber = log10(cond_num),
            Condition = "Stress",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Add optimization and covariance method types
    if(nrow(condition_numbers) > 0) {
      condition_numbers$OptimizationType <- sapply(condition_numbers$Method, function(m) {
        methods_info$OptimizationType[methods_info$Method == m]
      })
      
      condition_numbers$CovMethod <- sapply(condition_numbers$Method, function(m) {
        methods_info$CovMethod[methods_info$Method == m]
      })
      
      # Calculate statistics by covariance method
      cov_condition <- NULL
      if(length(unique(condition_numbers$CovMethod)) > 1) {
        cov_condition <- aggregate(
          LogConditionNumber ~ CovMethod + Condition, 
          data = condition_numbers, 
          FUN = mean, 
          na.rm = TRUE
        )
        
        # Reshape for comparison
        cov_condition_wide <- reshape2::dcast(
          cov_condition, 
          CovMethod ~ Condition, 
          value.var = "LogConditionNumber"
        )
        
        # Calculate change in condition number
        if(all(c("Standard", "Stress") %in% colnames(cov_condition_wide))) {
          cov_condition_wide$Change <- cov_condition_wide$Stress - cov_condition_wide$Standard
          
          # Sort by Standard condition number
          cov_condition_wide <- cov_condition_wide[order(cov_condition_wide$Standard), ]
        }
        
        # Calculate mean by covariance method (across both conditions)
        cov_condition_mean <- aggregate(
          LogConditionNumber ~ CovMethod, 
          data = condition_numbers, 
          FUN = mean, 
          na.rm = TRUE
        )
        
        # Sort by log condition number
        cov_condition_mean <- cov_condition_mean[order(cov_condition_mean$LogConditionNumber), ]
      }
      
      # Store results
      analysis_results$condition_numbers <- list(
        data = condition_numbers,
        cov_condition = list(
          by_condition = cov_condition,
          wide = cov_condition_wide,
          mean = cov_condition_mean
        )
      )
    }
  }
  
  # ========================================================================
  # 6. Regime-specific performance analysis
  # ========================================================================
  
  if("regime" %in% test_types && !is.null(regime_info$standard) && !is.null(regime_info$stress)) {
    cat("Analyzing regime-specific performance...\n")
    
    # Extract regime names
    regime_names <- unique(c(
      names(regime_info$standard)[names(regime_info$standard) != "overall"],
      names(regime_info$stress)[names(regime_info$stress) != "overall"]
    ))
    
    # Create data frame for regime-specific performance comparison
    regime_comparison <- data.frame()
    
    for(regime in regime_names) {
      # Check if regime exists in both analyses
      if(regime %in% names(regime_info$standard) && regime %in% names(regime_info$stress
      )) {
        std_regime <- regime_info$standard[[regime]]
        stress_regime <- regime_info$stress[[regime]]
        
        # Extract Sharpe ratios
        if("sharpe" %in% colnames(std_regime) && "sharpe" %in% colnames(stress_regime)) {
          # Get common methods
          common_methods <- intersect(std_regime$Method, stress_regime$Method)
          
          if(length(common_methods) > 0) {
            # Filter for common methods
            std_common <- std_regime[std_regime$Method %in% common_methods, ]
            stress_common <- stress_regime[stress_regime$Method %in% common_methods, ]
            
            # Sort in the same order
            std_common <- std_common[match(common_methods, std_common$Method), ]
            stress_common <- stress_common[match(common_methods, stress_common$Method), ]
            
            # Calculate changes
            regime_data <- data.frame(
              Method = common_methods,
              Regime = regime,
              Standard_Sharpe = std_common$sharpe,
              Stress_Sharpe = stress_common$sharpe,
              Absolute_Change = stress_common$sharpe - std_common$sharpe
            )
            
            # Calculate percent change
            regime_data$Percent_Change <- ifelse(
              !is.na(regime_data$Standard_Sharpe) & regime_data$Standard_Sharpe != 0,
              (regime_data$Stress_Sharpe / regime_data$Standard_Sharpe - 1) * 100,
              NA
            )
            
            # Add to main data frame
            regime_comparison <- rbind(regime_comparison, regime_data)
          }
        }
      }
    }
    
    # Add optimization and covariance method types
    if(nrow(regime_comparison) > 0) {
      regime_comparison$OptimizationType <- sapply(regime_comparison$Method, function(m) {
        methods_info$OptimizationType[methods_info$Method == m]
      })
      
      regime_comparison$CovMethod <- sapply(regime_comparison$Method, function(m) {
        methods_info$CovMethod[methods_info$Method == m]
      })
      
      # Calculate average change by regime
      regime_summary <- aggregate(
        Percent_Change ~ Regime, 
        data = regime_comparison, 
        FUN = function(x) c(
          Mean = mean(x, na.rm = TRUE),
          Median = median(x, na.rm = TRUE),
          SD = sd(x, na.rm = TRUE),
          Count = sum(!is.na(x))
        )
      )
      
      # Format the regime summary
      regime_summary <- do.call(data.frame, c(
        list(Regime = regime_summary$Regime),
        regime_summary$Percent_Change
      ))
      
      # Calculate regime-specific robustness by covariance method
      regime_cov_summary <- NULL
      if("CovMethod" %in% colnames(regime_comparison)) {
        regime_cov_summary <- aggregate(
          Percent_Change ~ Regime + CovMethod, 
          data = regime_comparison, 
          FUN = mean, 
          na.rm = TRUE
        )
        
        # Reshape for comparison
        regime_cov_wide <- reshape2::dcast(
          regime_cov_summary, 
          CovMethod ~ Regime, 
          value.var = "Percent_Change"
        )
      }
      
      # Calculate regime-specific robustness by optimization type
      regime_opt_summary <- NULL
      if("OptimizationType" %in% colnames(regime_comparison)) {
        regime_opt_summary <- aggregate(
          Percent_Change ~ Regime + OptimizationType, 
          data = regime_comparison, 
          FUN = mean, 
          na.rm = TRUE
        )
        
        # Reshape for comparison
        regime_opt_wide <- reshape2::dcast(
          regime_opt_summary, 
          OptimizationType ~ Regime, 
          value.var = "Percent_Change"
        )
      }
      
      # Store results
      analysis_results$regime_analysis <- list(
        data = regime_comparison,
        summary = regime_summary,
        cov_summary = regime_cov_summary,
        cov_wide = regime_cov_wide,
        opt_summary = regime_opt_summary,
        opt_wide = regime_opt_wide
      )
    }
  }
  
  # ========================================================================
  # Generate comprehensive report
  # ========================================================================
  
  cat("Generating statistical analysis report...\n")
  
  # Create report file path
  report_file <- file.path(output_dir, "statistical_analysis_report.md")
  
  # Initialize report content
  report_content <- c(
    "# Statistical Analysis of Portfolio Optimization Results",
    "",
    paste("Report generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Executive Summary",
    ""
  )
  
  # Add executive summary
  exec_summary <- c(
    "This report provides a comprehensive statistical analysis of portfolio optimization methods, comparing their performance under normal and stressed market conditions."
  )
  
  # Add key findings
  key_findings <- c("### Key Findings:", "")
  
  # 1. Add normality test findings
  if("normality" %in% names(analysis_results)) {
    normality_results <- analysis_results$normality$tests
    if(nrow(normality_results) > 0) {
      non_normal_pct <- round(mean(!normality_results$IsNormal, na.rm = TRUE) * 100, 1)
      heavy_tailed_pct <- round(mean(normality_results$HeavyTails, na.rm = TRUE) * 100, 1)
      
      key_findings <- c(key_findings,
                        paste("- **Normality Tests:** ", non_normal_pct, "% of portfolio return distributions are non-normal."),
                        paste("- **Heavy Tails:** ", heavy_tailed_pct, "% of return distributions exhibit heavy tails (excess kurtosis).")
      )
      
      if(non_normal_pct > 50) {
        key_findings <- c(key_findings,
                          "  - This strongly supports the use of robust estimation methods that don't rely on normality assumptions."
        )
      }
    }
  }
  
  # 2. Add comparison findings
  if("comparison" %in% names(analysis_results)) {
    comparison_data <- analysis_results$comparison$data
    paired_ttest <- analysis_results$comparison$paired_ttest
    
    if(nrow(comparison_data) > 0) {
      # Most robust method
      most_robust <- comparison_data$Method[1]
      most_robust_pct <- round(comparison_data$Percent_Change[1], 2)
      
      key_findings <- c(key_findings,
                        paste("- **Most Robust Method:** ", most_robust, " with a performance change of ", 
                              ifelse(most_robust_pct >= 0, "+", ""), most_robust_pct, "% under stress.", sep = "")
      )
      
      # Statistical significance
      if(!is.null(paired_ttest)) {
        if(paired_ttest$p.value < 0.05) {
          key_findings <- c(key_findings,
                            paste("- **Statistical Significance:** The performance degradation under stress is statistically significant (p-value = ", 
                                  round(paired_ttest$p.value, 4), ").", sep = "")
          )
        } else {
          key_findings <- c(key_findings,
                            paste("- **Statistical Significance:** The performance changes under stress are not statistically significant (p-value = ", 
                                  round(paired_ttest$p.value, 4), ").", sep = "")
          )
        }
      }
      
      # Optimization method comparison
      if(!is.null(analysis_results$comparison$opt_type_comparison)) {
        opt_comparison <- analysis_results$comparison$opt_type_comparison
        if(nrow(opt_comparison) > 1) {
          best_opt <- opt_comparison$OptimizationType[1]
          worst_opt <- opt_comparison$OptimizationType[nrow(opt_comparison)]
          best_pct <- round(opt_comparison$Mean_Percent_Change[1], 2)
          worst_pct <- round(opt_comparison$Mean_Percent_Change[nrow(opt_comparison)], 2)
          
          key_findings <- c(key_findings,
                            paste("- **Optimization Methods:** ", best_opt, " strategies (", 
                                  ifelse(best_pct >= 0, "+", ""), best_pct, "%) are more robust than ", 
                                  worst_opt, " strategies (", ifelse(worst_pct >= 0, "+", ""), worst_pct, "%).", sep = "")
          )
        }
      }
      
      # Covariance method comparison
      if(!is.null(analysis_results$comparison$cov_method_comparison)) {
        cov_comparison <- analysis_results$comparison$cov_method_comparison
        if(nrow(cov_comparison) > 1) {
          best_cov <- cov_comparison$CovMethod[1]
          best_cov_pct <- round(cov_comparison$Mean_Percent_Change[1], 2)
          
          key_findings <- c(key_findings,
                            paste("- **Covariance Methods:** ", best_cov, " provides the most robust performance (", 
                                  ifelse(best_cov_pct >= 0, "+", ""), best_cov_pct, "%).", sep = "")
          )
        }
      }
    }
  }
  
  # 3. Add weight stability findings
  if("weight_stability" %in% names(analysis_results)) {
    weight_stability <- analysis_results$weight_stability$data
    if(nrow(weight_stability) > 0) {
      most_stable <- weight_stability$Method[1]
      most_stable_score <- round(weight_stability$StabilityScore[1], 2)
      avg_turnover <- round(mean(weight_stability$Turnover, na.rm = TRUE) * 100, 1)
      
      key_findings <- c(key_findings,
                        paste("- **Weight Stability:** ", most_stable, " shows the highest weight stability (score: ", 
                              most_stable_score, ").", sep = ""),
                        paste("- **Average Turnover:** Portfolio methods require an average of ", avg_turnover, 
                              "% turnover when adapting to stress conditions.", sep = "")
      )
      
      # Add optimization type comparison
      if(!is.null(analysis_results$weight_stability$opt_stability)) {
        opt_stability <- analysis_results$weight_stability$opt_stability
        if(nrow(opt_stability) > 1) {
          most_stable_opt <- opt_stability$OptimizationType[1]
          
          key_findings <- c(key_findings,
                            paste("- **Most Stable Optimization:** ", most_stable_opt, 
                                  " strategies show the highest weight stability.", sep = "")
          )
        }
      }
    }
  }
  
  # 4. Add condition number findings
  if("condition_numbers" %in% names(analysis_results)) {
    if(!is.null(analysis_results$condition_numbers$cov_condition$mean)) {
      cov_condition <- analysis_results$condition_numbers$cov_condition$mean
      if(nrow(cov_condition) > 0) {
        best_numerical <- cov_condition$CovMethod[1]
        worst_numerical <- cov_condition$CovMethod[nrow(cov_condition)]
        
        key_findings <- c(key_findings,
                          paste("- **Numerical Stability:** ", best_numerical, 
                                " provides the best numerical conditioning, while ", worst_numerical, 
                                " has the highest condition numbers.", sep = "")
        )
      }
    }
  }
  
  # 5. Add regime analysis findings
  if("regime_analysis" %in% names(analysis_results)) {
    if(!is.null(analysis_results$regime_analysis$summary)) {
      regime_summary <- analysis_results$regime_analysis$summary
      if(nrow(regime_summary) > 0) {
        most_affected_regime <- regime_summary$Regime[which.min(regime_summary$Mean)]
        most_affected_change <- round(min(regime_summary$Mean, na.rm = TRUE), 2)
        
        key_findings <- c(key_findings,
                          paste("- **Regime Analysis:** The ", most_affected_regime, 
                                " regime shows the largest performance impact under stress (", 
                                most_affected_change, "%).", sep = "")
        )
        
        # Add method recommendations for specific regimes
        if(!is.null(analysis_results$regime_analysis$data)) {
          regime_data <- analysis_results$regime_analysis$data
          for(regime in unique(regime_data$Regime)) {
            regime_subset <- regime_data[regime_data$Regime == regime, ]
            if(nrow(regime_subset) > 0) {
              best_method <- regime_subset$Method[which.max(regime_subset$Percent_Change)]
              best_change <- round(max(regime_subset$Percent_Change, na.rm = TRUE), 2)
              
              key_findings <- c(key_findings,
                                paste("  - For ", regime, " regime: ", best_method, 
                                      " is most robust (", ifelse(best_change >= 0, "+", ""), 
                                      best_change, "%).", sep = "")
              )
            }
          }
        }
      }
    }
  }
  
  # Combine executive summary and key findings
  exec_summary <- c(exec_summary, "", key_findings)
  
  # Add to report
  report_content <- c(report_content, exec_summary, "")
  
  # ---- Add detailed sections ----
  
  # 1. Normality Tests
  if("normality" %in% names(analysis_results)) {
    normality_section <- c(
      "## 1. Distributional Analysis of Returns",
      "",
      "This section examines the statistical properties of portfolio returns, testing for normality and heavy tails.",
      ""
    )
    
    normality_results <- analysis_results$normality$tests
    if(nrow(normality_results) > 0) {
      # Add summary statistics
      normality_stats <- analysis_results$normality$stats
      
      normality_section <- c(normality_section,
                             "### Summary of Normality Tests",
                             "",
                             "```",
                             capture.output(print(normality_stats)),
                             "```",
                             "",
                             paste("Under normal market conditions, ", 
                                   round(mean(!normality_results$IsNormal[normality_results$Condition == "Standard"], na.rm = TRUE) * 100, 1),
                                   "% of portfolio return distributions are non-normal.",
                                   sep = ""),
                             "",
                             paste("Under stressed market conditions, ", 
                                   round(mean(!normality_results$IsNormal[normality_results$Condition == "Stress"], na.rm = TRUE) * 100, 1),
                                   "% of portfolio return distributions are non-normal.",
                                   sep = ""),
                             "",
                             "### Skewness and Kurtosis Analysis",
                             "",
                             "```",
                             capture.output(print(aggregate(cbind(Skewness, Kurtosis) ~ Condition, data = normality_results, mean, na.rm = TRUE))),
                             "```",
                             ""
      )
      
      # Add interpretation
      if(mean(normality_results$Kurtosis, na.rm = TRUE) > 3) {
        normality_section <- c(normality_section,
                               "**Interpretation:** The portfolio returns exhibit excess kurtosis (heavy tails), which is typical of financial returns. This confirms the importance of using robust estimation methods that don't rely on normality assumptions.",
                               ""
        )
      } else {
        normality_section <- c(normality_section,
                               "**Interpretation:** The portfolio returns show moderate departures from normality. Standard estimation methods may be adequate for most conditions, but robust methods could provide additional stability under stress.",
                               ""
        )
      }
      
      # Add table of top 10 methods with highest non-normality
      non_normality_ranking <- normality_results[order(-normality_results$Kurtosis), ]
      if(nrow(non_normality_ranking) > 10) {
        non_normality_ranking <- non_normality_ranking[1:10, ]
      }
      
      normality_section <- c(normality_section,
                             "### Top Methods with Highest Non-Normality",
                             "",
                             "```",
                             capture.output(print(non_normality_ranking[, c("Method", "Condition", "Skewness", "Kurtosis", "IsNormal")])),
                             "```",
                             ""
      )
      
      report_content <- c(report_content, normality_section)
    }
  }
  
  # 2. Bootstrap Analysis
  if("bootstrap" %in% names(analysis_results)) {
    bootstrap_section <- c(
      "## 2. Bootstrap Analysis of Performance Metrics",
      "",
      paste("This section provides bootstrap confidence intervals (", confidence_level * 100, 
            "%) for key performance metrics.", sep = ""),
      ""
    )
    
    bootstrap_results <- analysis_results$bootstrap
    if(nrow(bootstrap_results) > 0) {
      # Create subsections for each metric
      metrics <- unique(bootstrap_results$Metric)
      
      for(metric in metrics) {
        metric_name <- switch(metric,
                              "sharpe" = "Sharpe Ratio",
                              "sortino" = "Sortino Ratio",
                              "var95" = "Value at Risk (95%)",
                              "cvar95" = "Conditional VaR (95%)",
                              "maxdd" = "Maximum Drawdown",
                              metric)
        
        metric_results <- bootstrap_results[bootstrap_results$Metric == metric, ]
        
        # Separate standard and stress conditions
        std_results <- metric_results[metric_results$Condition == "Standard", ]
        stress_results <- metric_results[metric_results$Condition == "Stress", ]
        
        bootstrap_section <- c(bootstrap_section,
                               paste("### ", metric_name, sep = ""),
                               "",
                               "#### Standard Market Conditions",
                               ""
        )
        
        if(nrow(std_results) > 0) {
          # Order by original value
          std_results <- std_results[order(-std_results$Original), ]
          
          # Take top 10 for table
          top_std <- std_results
          if(nrow(top_std) > 10) top_std <- top_std[1:10, ]
          
          bootstrap_section <- c(bootstrap_section,
                                 "```",
                                 capture.output(print(top_std[, c("Method", "Original", "BootMean", "CI_Lower", "CI_Upper", "RelPrecision")])),
                                 "```",
                                 "",
                                 paste("Best method: ", top_std$Method[1], " with ", metric_name, " of ", 
                                       round(top_std$Original[1], 4), " (", round(top_std$CI_Lower[1], 4), " - ", 
                                       round(top_std$CI_Upper[1], 4), ")", sep = ""),
                                 ""
          )
        } else {
          bootstrap_section <- c(bootstrap_section,
                                 "No results available for standard conditions.",
                                 ""
          )
        }
        
        bootstrap_section <- c(bootstrap_section,
                               "#### Stressed Market Conditions",
                               ""
        )
        
        if(nrow(stress_results) > 0) {
          # Order by original value
          stress_results <- stress_results[order(-stress_results$Original), ]
          
          # Take top 10 for table
          top_stress <- stress_results
          if(nrow(top_stress) > 10) top_stress <- top_stress[1:10, ]
          
          bootstrap_section <- c(bootstrap_section,
                                 "```",
                                 capture.output(print(top_stress[, c("Method", "Original", "BootMean", "CI_Lower", "CI_Upper", "RelPrecision")])),
                                 "```",
                                 "",
                                 paste("Best method: ", top_stress$Method[1], " with ", metric_name, " of ", 
                                       round(top_stress$Original[1], 4), " (", round(top_stress$CI_Lower[1], 4), " - ", 
                                       round(top_stress$CI_Upper[1], 4), ")", sep = ""),
                                 ""
          )
          
          # Check for overlap in top methods
          common_top <- intersect(std_results$Method[1:min(5, nrow(std_results))], 
                                  stress_results$Method[1:min(5, nrow(stress_results))])
          
          if(length(common_top) > 0) {
            bootstrap_section <- c(bootstrap_section,
                                   paste("Methods that perform well in both conditions: ", paste(common_top, collapse = ", "), sep = ""),
                                   ""
            )
          } else {
            bootstrap_section <- c(bootstrap_section,
                                   "No methods appear in the top 5 for both standard and stressed conditions.",
                                   ""
            )
          }
        } else {
          bootstrap_section <- c(bootstrap_section,
                                 "No results available for stressed conditions.",
                                 ""
          )
        }
      }
      
      report_content <- c(report_content, bootstrap_section)
    }
  }
  
  # 3. Robustness Comparison
  if("comparison" %in% names(analysis_results)) {
    comparison_section <- c(
      "## 3. Robustness Analysis (Standard vs. Stress)",
      "",
      "This section compares portfolio performance under normal and stressed market conditions, identifying the most robust methods.",
      ""
    )
    
    comparison_data <- analysis_results$comparison$data
    if(nrow(comparison_data) > 0) {
      # Add paired t-test results
      paired_ttest <- analysis_results$comparison$paired_ttest
      wilcox_test <- analysis_results$comparison$wilcox_test
      
      comparison_section <- c(comparison_section,
                              "### Statistical Tests for Performance Differences",
                              ""
      )
      if(!is.null(paired_ttest)) {
        comparison_section <- c(comparison_section,
                                "#### Paired t-test",
                                "",
                                "```",
                                capture.output(print(paired_ttest)),
                                "```",
                                "",
                                paste("**Interpretation:** The difference in Sharpe ratios between standard and stressed conditions is",
                                      ifelse(paired_ttest$p.value < 0.05, "", " not"), " statistically significant (p-value = ",
                                      round(paired_ttest$p.value, 4), ").", sep = ""),
                                ""
        )
      }
      
      if(!is.null(wilcox_test)) {
        comparison_section <- c(comparison_section,
                                "#### Wilcoxon Signed-Rank Test (non-parametric)",
                                "",
                                "```",
                                capture.output(print(wilcox_test)),
                                "```",
                                "",
                                paste("**Interpretation:** The non-parametric test ", 
                                      ifelse(wilcox_test$p.value < 0.05, "confirms", "does not confirm"), 
                                      " a significant difference between standard and stressed conditions (p-value = ",
                                      round(wilcox_test$p.value, 4), ").", sep = ""),
                                ""
        )
      }
      
      # Add overall robustness ranking
      comparison_section <- c(comparison_section,
                              "### Performance Robustness Ranking",
                              "",
                              "Methods ranked by percent change in Sharpe ratio under stress (higher = more robust):",
                              "",
                              "```",
                              capture.output(print(comparison_data[, c("Method", "Standard_Sharpe", "Stress_Sharpe", "Percent_Change")])),
                              "```",
                              ""
      )
      
      # Add optimization type comparison
      if(!is.null(analysis_results$comparison$opt_type_comparison)) {
        opt_comparison <- analysis_results$comparison$opt_type_comparison
        
        comparison_section <- c(comparison_section,
                                "### Optimization Method Comparison",
                                "",
                                "```",
                                capture.output(print(opt_comparison)),
                                "```",
                                "",
                                paste("**Interpretation:** ", opt_comparison$OptimizationType[1], " methods show the best robustness with an average change of ",
                                      ifelse(opt_comparison$Mean_Percent_Change[1] >= 0, "+", ""),
                                      round(opt_comparison$Mean_Percent_Change[1], 2), "% under stress.", sep = ""),
                                ""
        )
      }
      
      # Add covariance method comparison
      if(!is.null(analysis_results$comparison$cov_method_comparison)) {
        cov_comparison <- analysis_results$comparison$cov_method_comparison
        
        comparison_section <- c(comparison_section,
                                "### Covariance Method Comparison",
                                "",
                                "```",
                                capture.output(print(cov_comparison)),
                                "```",
                                "",
                                paste("**Interpretation:** ", cov_comparison$CovMethod[1], " covariance estimation shows the best robustness with an average change of ",
                                      ifelse(cov_comparison$Mean_Percent_Change[1] >= 0, "+", ""),
                                      round(cov_comparison$Mean_Percent_Change[1], 2), "% under stress.", sep = ""),
                                ""
        )
      }
      
      report_content <- c(report_content, comparison_section)
    }
  }
  
  # 4. Weight Stability Analysis
  if("weight_stability" %in% names(analysis_results)) {
    stability_section <- c(
      "## 4. Portfolio Weight Stability Analysis",
      "",
      "This section analyzes how portfolio weights change under stress, measuring turnover and identifying the most stable methods.",
      ""
    )
    
    weight_stability <- analysis_results$weight_stability$data
    if(nrow(weight_stability) > 0) {
      # Add overall stability ranking
      stability_section <- c(stability_section,
                             "### Weight Stability Ranking",
                             "",
                             "```",
                             capture.output(print(weight_stability[, c("Method", "Turnover", "MaxWeightChange", "WeightCorrelation", "StabilityScore")])),
                             "```",
                             "",
                             paste("**Most stable method:** ", weight_stability$Method[1], 
                                   " with turnover of ", round(weight_stability$Turnover[1] * 100, 2), 
                                   "% and stability score of ", round(weight_stability$StabilityScore[1], 2), sep = ""),
                             ""
      )
      
      # Add optimization type comparison
      if(!is.null(analysis_results$weight_stability$opt_stability)) {
        opt_stability <- analysis_results$weight_stability$opt_stability
        
        stability_section <- c(stability_section,
                               "### Stability by Optimization Type",
                               "",
                               "```",
                               capture.output(print(opt_stability)),
                               "```",
                               "",
                               paste("**Interpretation:** ", opt_stability$OptimizationType[1], 
                                     " methods show the highest weight stability with an average turnover of ", 
                                     round(opt_stability$Turnover[1] * 100, 2), "%.", sep = ""),
                               ""
        )
      }
      
      # Add covariance method comparison
      if(!is.null(analysis_results$weight_stability$cov_stability)) {
        cov_stability <- analysis_results$weight_stability$cov_stability
        
        stability_section <- c(stability_section,
                               "### Stability by Covariance Method",
                               "",
                               "```",
                               capture.output(print(cov_stability)),
                               "```",
                               "",
                               paste("**Interpretation:** ", cov_stability$CovMethod[1], 
                                     " covariance estimation provides the highest weight stability with an average turnover of ", 
                                     round(cov_stability$Turnover[1] * 100, 2), "%.", sep = ""),
                               ""
        )
      }
      
      report_content <- c(report_content, stability_section)
    }
  }
  
  # 5. Numerical Stability Analysis
  if("condition_numbers" %in% names(analysis_results)) {
    numerical_section <- c(
      "## 5. Numerical Stability Analysis",
      "",
      "This section examines the condition numbers of covariance matrices, indicating numerical stability and potential sensitivity to estimation errors.",
      ""
    )
    
    condition_numbers <- analysis_results$condition_numbers$data
    if(nrow(condition_numbers) > 0) {
      # Add condition number by covariance method
      if(!is.null(analysis_results$condition_numbers$cov_condition$mean)) {
        cov_condition <- analysis_results$condition_numbers$cov_condition$mean
        
        numerical_section <- c(numerical_section,
                               "### Condition Number by Covariance Method",
                               "",
                               "Lower values indicate better numerical conditioning (log scale):",
                               "",
                               "```",
                               capture.output(print(cov_condition)),
                               "```",
                               "",
                               paste("**Best numerically conditioned method:** ", cov_condition$CovMethod[1], 
                                     " with log10 condition number of ", round(cov_condition$LogConditionNumber[1], 2), sep = ""),
                               ""
        )
      }
      
      # Add change in condition number under stress
      if(!is.null(analysis_results$condition_numbers$cov_condition$wide)) {
        cond_wide <- analysis_results$condition_numbers$cov_condition$wide
        
        if(all(c("Standard", "Stress", "Change") %in% colnames(cond_wide))) {
          numerical_section <- c(numerical_section,
                                 "### Change in Condition Number Under Stress",
                                 "",
                                 "```",
                                 capture.output(print(cond_wide)),
                                 "```",
                                 ""
          )
          
          # Add interpretation
          if(any(cond_wide$Change > 0.5)) {
            numerical_section <- c(numerical_section,
                                   "**Interpretation:** Some covariance methods show significant degradation in numerical conditioning under stress conditions. This could lead to estimation errors and portfolio instability.",
                                   ""
            )
          } else {
            numerical_section <- c(numerical_section,
                                   "**Interpretation:** Covariance methods maintain relatively stable numerical conditioning under stress, which is favorable for reliable portfolio optimization.",
                                   ""
            )
          }
        }
      }
      
      report_content <- c(report_content, numerical_section)
    }
  }
  
  # 6. Regime Analysis
  if("regime_analysis" %in% names(analysis_results)) {
    regime_section <- c(
      "## 6. Market Regime Analysis",
      "",
      "This section examines how different market regimes affect the robustness of portfolio methods under stress.",
      ""
    )
    
    if(!is.null(analysis_results$regime_analysis$summary)) {
      regime_summary <- analysis_results$regime_analysis$summary
      
      regime_section <- c(regime_section,
                          "### Performance Change by Market Regime",
                          "",
                          "```",
                          capture.output(print(regime_summary)),
                          "```",
                          ""
      )
      
      # Find most and least affected regimes
      most_affected <- regime_summary$Regime[which.min(regime_summary$Mean)]
      least_affected <- regime_summary$Regime[which.max(regime_summary$Mean)]
      
      regime_section <- c(regime_section,
                          paste("- Most affected regime: ", most_affected, " with average change of ", 
                                round(min(regime_summary$Mean, na.rm = TRUE), 2), "%", sep = ""),
                          paste("- Least affected regime: ", least_affected, " with average change of ", 
                                round(max(regime_summary$Mean, na.rm = TRUE), 2), "%", sep = ""),
                          ""
      )
      
      # Add covariance method performance by regime
      if(!is.null(analysis_results$regime_analysis$cov_wide)) {
        cov_regime <- analysis_results$regime_analysis$cov_wide
        
        regime_section <- c(regime_section,
                            "### Covariance Method Performance by Regime",
                            "",
                            "```",
                            capture.output(print(cov_regime)),
                            "```",
                            ""
        )
        
        # Find best method for each regime
        best_methods <- character(0)
        for(regime in unique(analysis_results$regime_analysis$data$Regime)) {
          if(regime %in% colnames(cov_regime)) {
            best_method <- cov_regime[which.max(cov_regime[[regime]]), "CovMethod"]
            best_change <- max(cov_regime[[regime]], na.rm = TRUE)
            
            best_methods <- c(best_methods,
                              paste("- For ", regime, " regime: ", best_method, " (", 
                                    ifelse(best_change >= 0, "+", ""), round(best_change, 2), "%)", sep = "")
            )
          }
        }
        
        if(length(best_methods) > 0) {
          regime_section <- c(regime_section,
                              "#### Best Covariance Method by Regime:",
                              "",
                              best_methods,
                              ""
          )
        }
      }
      
      # Add optimization method performance by regime
      if(!is.null(analysis_results$regime_analysis$opt_wide)) {
        opt_regime <- analysis_results$regime_analysis$opt_wide
        
        regime_section <- c(regime_section,
                            "### Optimization Method Performance by Regime",
                            "",
                            "```",
                            capture.output(print(opt_regime)),
                            "```",
                            ""
        )
        
        # Find best optimization for each regime
        best_opts <- character(0)
        for(regime in unique(analysis_results$regime_analysis$data$Regime)) {
          if(regime %in% colnames(opt_regime)) {
            best_opt <- opt_regime[which.max(opt_regime[[regime]]), "OptimizationType"]
            best_change <- max(opt_regime[[regime]], na.rm = TRUE)
            
            best_opts <- c(best_opts,
                           paste("- For ", regime, " regime: ", best_opt, " (", 
                                 ifelse(best_change >= 0, "+", ""), round(best_change, 2), "%)", sep = "")
            )
          }
        }
        
        if(length(best_opts) > 0) {
          regime_section <- c(regime_section,
                              "#### Best Optimization Method by Regime:",
                              "",
                              best_opts,
                              ""
          )
        }
      }
      
      report_content <- c(report_content, regime_section)
    }
  }
  
  # Add final recommendations
  recommendations <- c(
    "## Recommendations",
    "",
    "Based on the comprehensive statistical analysis, the following recommendations are made:"
  )
  
  # Determine most robust method overall
  most_robust_method <- "Unknown"
  if("comparison" %in% names(analysis_results)) {
    if(nrow(analysis_results$comparison$data) > 0) {
      most_robust_method <- analysis_results$comparison$data$Method[1]
    }
  }
  
  # Determine most stable method
  most_stable_method <- "Unknown"
  if("weight_stability" %in% names(analysis_results)) {
    if(nrow(analysis_results$weight_stability$data) > 0) {
      most_stable_method <- analysis_results$weight_stability$data$Method[1]
    }
  }
  
  # Determine best numerical method
  best_numerical_method <- "Unknown"
  if("condition_numbers" %in% names(analysis_results)) {
    if(!is.null(analysis_results$condition_numbers$cov_condition$mean)) {
      if(nrow(analysis_results$condition_numbers$cov_condition$mean) > 0) {
        best_numerical_method <- analysis_results$condition_numbers$cov_condition$mean$CovMethod[1]
      }
    }
  }
  
  # Get best covariance method
  best_cov_method <- "Unknown"
  if("comparison" %in% names(analysis_results)) {
    if(!is.null(analysis_results$comparison$cov_method_comparison)) {
      if(nrow(analysis_results$comparison$cov_method_comparison) > 0) {
        best_cov_method <- analysis_results$comparison$cov_method_comparison$CovMethod[1]
      }
    }
  }
  
  # Get best optimization type
  best_opt_type <- "Unknown"
  if("comparison" %in% names(analysis_results)) {
    if(!is.null(analysis_results$comparison$opt_type_comparison)) {
      if(nrow(analysis_results$comparison$opt_type_comparison) > 0) {
        best_opt_type <- analysis_results$comparison$opt_type_comparison$OptimizationType[1]
      }
    }
  }
  
  # Add specific recommendations
  recommendations <- c(recommendations, "",
                       paste("1. **Best Overall Method:** ", most_robust_method, 
                             " - This method shows the highest robustness under stress conditions while maintaining strong performance.", sep = ""),
                       "",
                       paste("2. **Most Weight-Stable Method:** ", most_stable_method, 
                             " - This method requires the least portfolio rebalancing when adapting to changing market conditions.", sep = ""),
                       "",
                       paste("3. **Best Numerical Stability:** ", best_numerical_method, 
                             " covariance estimation provides the most reliable numerical properties, reducing sensitivity to estimation errors.", sep = ""),
                       "",
                       paste("4. **Recommended Covariance Estimation:** ", best_cov_method, 
                             " - This approach to covariance estimation consistently provides the best robustness across various optimization methods.", sep = ""),
                       "",
                       paste("5. **Recommended Optimization Approach:** ", best_opt_type, 
                             " optimization provides superior robustness compared to other approaches.", sep = "")
  )
  
  # Add regime-specific recommendations if available
  if("regime_analysis" %in% names(analysis_results)) {
    if(!is.null(analysis_results$regime_analysis$data)) {
      regime_recs <- c("", "6. **Regime-Specific Recommendations:**", "")
      
      regime_data <- analysis_results$regime_analysis$data
      for(regime in unique(regime_data$Regime)) {
        regime_subset <- regime_data[regime_data$Regime == regime, ]
        if(nrow(regime_subset) > 0) {
          best_method <- regime_subset$Method[which.max(regime_subset$Percent_Change)]
          
          regime_recs <- c(regime_recs,
                           paste("   - For ", regime, " regime: ", best_method, sep = "")
          )
        }
      }
      
      recommendations <- c(recommendations, regime_recs)
    }
  }
  
  report_content <- c(report_content, recommendations)
  
  # Add conclusion
  conclusion <- c(
    "",
    "## Conclusion",
    "",
    "This statistical analysis confirms that robust portfolio optimization methods provide significant advantages, especially under stressed market conditions. The analysis highlights the importance of:",
    "",
    "- Considering non-normality and heavy tails in return distributions",
    "- Using robust covariance estimation to improve numerical stability",
    "- Selecting optimization methods that maintain performance under stress",
    "- Adapting methods to specific market regimes when possible",
    "",
    paste("The ", most_robust_method, " approach demonstrates the best combination of performance, robustness, and stability, making it suitable for practical implementation in real-world portfolio management.", sep = "")
  )
  
  report_content <- c(report_content, conclusion)
  
  # Write the report
  writeLines(report_content, report_file)
  
  cat("Statistical analysis report generated at:", report_file, "\n")
  
  # ========================================================================
  # Create visualizations
  # ========================================================================
  
  cat("Creating visualizations...\n")
  
  # Create plots directory
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # 1. Normality tests visualization
  if("normality" %in% names(analysis_results)) {
    if(nrow(analysis_results$normality$tests) > 0) {
      normality_data <- analysis_results$normality$tests
      
      # Plot skewness vs kurtosis
      p1 <- ggplot(normality_data, aes(x = Skewness, y = Kurtosis, color = Condition, shape = CovMethod)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        labs(title = "Skewness vs. Kurtosis of Portfolio Returns",
             subtitle = "Normal distribution: Skewness = 0, Kurtosis = 3",
             x = "Skewness", y = "Kurtosis") +
        theme_minimal()
      
      ggsave(file.path(plots_dir, "skewness_kurtosis.png"), p1, width = 10, height = 8)
      
      # Plot normality test results
      normality_summary <- aggregate(IsNormal ~ Method + Condition, 
                                     data = normality_data, 
                                     FUN = function(x) as.numeric(!x))
      
      normality_wide <- reshape2::dcast(normality_summary, Method ~ Condition, value.var = "IsNormal")
      
      if(all(c("Standard", "Stress") %in% colnames(normality_wide))) {
        p2 <- ggplot(normality_wide, aes(x = Standard, y = Stress)) +
          geom_point(size = 3, alpha = 0.7) +
          geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
          labs(title = "Non-Normality Under Different Market Conditions",
               subtitle = "Higher values indicate stronger evidence against normality",
               x = "Non-Normality (Standard)", y = "Non-Normality (Stress)") +
          theme_minimal()
        
        ggsave(file.path(plots_dir, "normality_comparison.png"), p2, width = 10, height = 8)
      }
    }
  }
  
  # 2. Performance comparison visualization
  if("comparison" %in% names(analysis_results)) {
    if(nrow(analysis_results$comparison$data) > 0) {
      comparison_data <- analysis_results$comparison$data
      
      # Plot performance change
      p3 <- ggplot(comparison_data, aes(x = reorder(Method, Percent_Change), y = Percent_Change, fill = CovMethod)) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(title = "Portfolio Performance Change Under Stress",
             subtitle = "Percent change in Sharpe ratio (higher is better)",
             x = "", y = "Percent Change (%)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(plots_dir, "performance_change.png"), p3, width = 12, height = 8)
      
      # Plot standard vs. stress performance
      p4 <- ggplot(comparison_data, aes(x = Standard_Sharpe, y = Stress_Sharpe, color = CovMethod, shape = OptimizationType)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
        labs(title = "Performance Under Standard vs. Stressed Conditions",
             subtitle = "Points above the line perform better under stress than expected",
             x = "Sharpe Ratio (Standard)", y = "Sharpe Ratio (Stress)") +
        theme_minimal()
      
      ggsave(file.path(plots_dir, "performance_comparison.png"), p4, width = 10, height = 8)
      
      # Plot by covariance method
      if(!is.null(analysis_results$comparison$cov_method_comparison)) {
        cov_comparison <- analysis_results$comparison$cov_method_comparison
        
        p5 <- ggplot(cov_comparison, aes(x = reorder(CovMethod, Mean_Percent_Change), y = Mean_Percent_Change)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          geom_errorbar(aes(ymin = Mean_Percent_Change - SD_Percent_Change, 
                            ymax = Mean_Percent_Change + SD_Percent_Change), 
                        width = 0.2) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
          labs(title = "Performance Robustness by Covariance Method",
               subtitle = "Average percent change in Sharpe ratio under stress",
               x = "Covariance Method", y = "Mean Percent Change (%)") +
          theme_minimal()
        
        ggsave(file.path(plots_dir, "covariance_comparison.png"), p5, width = 10, height = 6)
      }
      
      # Plot by optimization type
      if(!is.null(analysis_results$comparison$opt_type_comparison)) {
        opt_comparison <- analysis_results$comparison$opt_type_comparison
        
        p6 <- ggplot(opt_comparison, aes(x = reorder(OptimizationType, Mean_Percent_Change), y = Mean_Percent_Change)) +
          geom_bar(stat = "identity", fill = "darkgreen") +
          geom_errorbar(aes(ymin = Mean_Percent_Change - SD_Percent_Change, 
                            ymax = Mean_Percent_Change + SD_Percent_Change), 
                        width = 0.2) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
          labs(title = "Performance Robustness by Optimization Method",
               subtitle = "Average percent change in Sharpe ratio under stress",
               x = "Optimization Method", y = "Mean Percent Change (%)") +
          theme_minimal()
        
        ggsave(file.path(plots_dir, "optimization_comparison.png"), p6, width = 10, height = 6)
      }
    }
  }
  
  # 3. Weight stability visualization
  if("weight_stability" %in% names(analysis_results)) {
    if(nrow(analysis_results$weight_stability$data) > 0) {
      weight_stability <- analysis_results$weight_stability$data
      
      # Plot turnover by method
      p7 <- ggplot(weight_stability, aes(x = reorder(Method, -Turnover), y = Turnover * 100, fill = CovMethod)) +
        geom_bar(stat = "identity") +
        labs(title = "Portfolio Weight Turnover Under Stress",
             subtitle = "Lower values indicate higher stability",
             x = "", y = "Turnover (%)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(plots_dir, "weight_turnover.png"), p7, width = 12, height = 8)
      
      # Plot weight correlation vs turnover
      p8 <- ggplot(weight_stability, aes(x = WeightCorrelation, y = Turnover * 100, color = CovMethod, shape = OptimizationType)) +
        geom_point(size = 3, alpha = 0.7) +
        geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
        labs(title = "Weight Correlation vs. Turnover",
             subtitle = "Upper left corner represents the most stable methods",
             x = "Weight Correlation", y = "Turnover (%)") +
        theme_minimal()
      
      ggsave(file.path(plots_dir, "stability_scatter.png"), p8, width = 10, height = 8)
    }
  }
  
  # 4. Condition number visualization
  if("condition_numbers" %in% names(analysis_results)) {
    if(nrow(analysis_results$condition_numbers$data) > 0) {
      condition_numbers <- analysis_results$condition_numbers$data
      
      # Plot condition numbers by method
      p9 <- ggplot(condition_numbers, aes(x = reorder(Method, -LogConditionNumber), 
                                          y = LogConditionNumber, fill = Condition)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Numerical Conditioning (Log10 Condition Number)",
             subtitle = "Lower values indicate better conditioning",
             x = "", y = "Log10 Condition Number") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggsave(file.path(plots_dir, "condition_numbers.png"), p9, width = 12, height = 8)
      
      # Plot by covariance method
      if(!is.null(analysis_results$condition_numbers$cov_condition$wide)) {
        cond_wide <- analysis_results$condition_numbers$cov_condition$wide
        
        if(all(c("Standard", "Stress") %in% colnames(cond_wide))) {
          p10 <- ggplot(cond_wide, aes(x = Standard, y = Stress)) +
            geom_point(size = 3, alpha = 0.7) +
            geom_text(aes(label = CovMethod), vjust = -1, size = 3, check_overlap = TRUE) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
            labs(title = "Condition Numbers: Standard vs. Stress",
                 subtitle = "Points below the line have better conditioning under stress",
                 x = "Log10 Condition Number (Standard)", 
                 y = "Log10 Condition Number (Stress)") +
            theme_minimal()
          
          ggsave(file.path(plots_dir, "condition_comparison.png"), p10, width = 10, height = 8)
        }
      }
    }
  }
  
  # 5. Regime analysis visualization
  if("regime_analysis" %in% names(analysis_results)) {
    if(!is.null(analysis_results$regime_analysis$data)) {
      regime_data <- analysis_results$regime_analysis$data
      
      # Plot by regime
      p11 <- ggplot(regime_data, aes(x = Regime, y = Percent_Change, fill = CovMethod)) +
        geom_boxplot() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(title = "Performance Change by Market Regime",
             subtitle = "Distribution of percent changes in Sharpe ratio",
             x = "Market Regime", y = "Percent Change (%)") +
        theme_minimal()
      
      ggsave(file.path(plots_dir, "regime_performance.png"), p11, width = 10, height = 8)
      
      # Create heatmap for regime-specific performance
      if(!is.null(analysis_results$regime_analysis$cov_wide)) {
        cov_wide <- analysis_results$regime_analysis$cov_wide
        
        # Extract regime columns
        regime_cols <- setdiff(colnames(cov_wide), "CovMethod")
        
        # Melt for heatmap
        cov_long <- reshape2::melt(cov_wide, id.vars = "CovMethod", 
                                   variable.name = "Regime", value.name = "Percent_Change")
        
        p12 <- ggplot(cov_long, aes(x = Regime, y = CovMethod, fill = Percent_Change)) +
          geom_tile() +
          scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
          labs(title = "Covariance Method Performance by Regime",
               subtitle = "Blue = better performance under stress",
               x = "Market Regime", y = "Covariance Method") +
          theme_minimal()
        
        ggsave(file.path(plots_dir, "regime_heatmap.png"), p12, width = 10, height = 8)
      }
    }
  }
  
  cat("Visualizations saved to", plots_dir, "\n")
  
  # Return all analysis results
  return(analysis_results)
}

# Function to create combined comparison visualizations for covariance methods
create_combined_visuals <- function(analysis_results, output_dir) {
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  
  # Check required data is available
  if(!all(c("comparison", "weight_stability", "condition_numbers") %in% names(analysis_results))) {
    cat("Missing required analysis results for combined visuals\n")
    return(NULL)
  }
  
  # Create directory for combined plots
  combined_dir <- file.path(output_dir, "combined_plots")
  if(!dir.exists(combined_dir)) {
    dir.create(combined_dir, recursive = TRUE)
  }
  
  # Prepare data for covariance methods comparison
  cov_methods <- unique(analysis_results$comparison$data$CovMethod)
  if(length(cov_methods) < 2) {
    cat("Not enough covariance methods for comparison\n")
    return(NULL)
  }
  
  # Extract metrics for each covariance method
  cov_metrics <- data.frame(
    CovMethod = cov_methods,
    stringsAsFactors = FALSE
  )
  
  # Add performance robustness
  if(!is.null(analysis_results$comparison$cov_method_comparison)) {
    perf_data <- analysis_results$comparison$cov_method_comparison
    
    cov_metrics$PerformanceChange <- sapply(cov_methods, function(cov) {
      idx <- which(perf_data$CovMethod == cov)
      if(length(idx) > 0) perf_data$Mean_Percent_Change[idx] else NA
    })
  }
  
  # Add weight stability
  if(!is.null(analysis_results$weight_stability$cov_stability)) {
    stab_data <- analysis_results$weight_stability$cov_stability
    
    cov_metrics$WeightStability <- sapply(cov_methods, function(cov) {
      idx <- which(stab_data$CovMethod == cov)
      if(length(idx) > 0) stab_data$StabilityScore[idx] else NA
    })
    
    cov_metrics$Turnover <- sapply(cov_methods, function(cov) {
      idx <- which(stab_data$CovMethod == cov)
      if(length(idx) > 0) stab_data$Turnover[idx] else NA
    })
  }
  
  # Add numerical conditioning
  if(!is.null(analysis_results$condition_numbers$cov_condition$mean)) {
    cond_data <- analysis_results$condition_numbers$cov_condition$mean
    
    cov_metrics$LogConditionNumber <- sapply(cov_methods, function(cov) {
      idx <- which(cond_data$CovMethod == cov)
      if(length(idx) > 0) cond_data$LogConditionNumber[idx] else NA
    })
  }
  
  # Create radar chart data (normalize all metrics to 0-1 scale)
  radar_data <- cov_metrics
  
  # Normalize performance change (higher = better)
  if("PerformanceChange" %in% colnames(radar_data)) {
    perf_range <- range(radar_data$PerformanceChange, na.rm = TRUE)
    if(diff(perf_range) > 0) {
      radar_data$NormPerformance <- (radar_data$PerformanceChange - perf_range[1]) / diff(perf_range)
    } else {
      radar_data$NormPerformance <- 0.5  # If all same value
    }
  }
  
  # Normalize weight stability (higher = better)
  if("WeightStability" %in% colnames(radar_data)) {
    stab_range <- range(radar_data$WeightStability, na.rm = TRUE)
    if(diff(stab_range) > 0) {
      radar_data$NormStability <- (radar_data$WeightStability - stab_range[1]) / diff(stab_range)
    } else {
      radar_data$NormStability <- 0.5  # If all same value
    }
  }
  
  # Normalize conditioning (lower = better, so invert)
  if("LogConditionNumber" %in% colnames(radar_data)) {
    cond_range <- range(radar_data$LogConditionNumber, na.rm = TRUE)
    if(diff(cond_range) > 0) {
      radar_data$NormConditioning <- 1 - (radar_data$LogConditionNumber - cond_range[1]) / diff(cond_range)
    } else {
      radar_data$NormConditioning <- 0.5  # If all same value
    }
  }
  
  # Calculate overall score (average of normalized metrics)
  norm_cols <- grep("^Norm", colnames(radar_data), value = TRUE)
  if(length(norm_cols) > 0) {
    radar_data$OverallScore <- rowMeans(radar_data[, norm_cols], na.rm = TRUE)
    
    # Sort by overall score
    radar_data <- radar_data[order(-radar_data$OverallScore), ]
  }
  
  # Create a bar chart of overall scores
  p1 <- ggplot(radar_data, aes(x = reorder(CovMethod, OverallScore), y = OverallScore)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Overall Performance of Covariance Methods",
         subtitle = "Combined score across performance, stability, and conditioning",
         x = "Covariance Method", y = "Overall Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(combined_dir, "overall_score.png"), p1, width = 10, height = 6)
  
  # Create a heatmap of normalized metrics
  heatmap_data <- reshape2::melt(radar_data, 
                                 id.vars = "CovMethod", 
                                 measure.vars = norm_cols,
                                 variable.name = "Metric", 
                                 value.name = "Score")
  
  # Clean up metric names
  heatmap_data$Metric <- gsub("Norm", "", heatmap_data$Metric)
  
  p2 <- ggplot(heatmap_data, aes(x = Metric, y = CovMethod, fill = Score)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = "Performance Metrics by Covariance Method",
         subtitle = "Higher values (darker blue) are better",
         x = "Metric", y = "Covariance Method") +
    theme_minimal()
  
  ggsave(file.path(combined_dir, "metrics_heatmap.png"), p2, width = 10, height = 8)
  
  # Create a bubble chart of performance vs. stability
  if(all(c("NormPerformance", "NormStability", "NormConditioning") %in% colnames(radar_data))) {
    p3 <- ggplot(radar_data, aes(x = NormPerformance, y = NormStability, size = NormConditioning)) +
      geom_point(alpha = 0.7) +
      geom_text(aes(label = CovMethod), vjust = -1.5, size = 4) +
      labs(title = "Multi-dimensional Comparison of Covariance Methods",
           subtitle = "Larger points indicate better numerical conditioning",
           x = "Normalized Performance", y = "Normalized Stability",
           size = "Numerical\nConditioning") +
      theme_minimal() +
      scale_size(range = c(2, 10))
    
    ggsave(file.path(combined_dir, "multidimensional_comparison.png"), p3, width = 10, height = 8)
  }
  
  # Write summary table
  summary_file <- file.path(combined_dir, "covariance_methods_summary.md")
  
  summary_content <- c(
    "# Covariance Method Performance Summary",
    "",
    "This table summarizes the performance of different covariance estimation methods across key metrics:",
    "",
    "| Covariance Method | Performance Change | Weight Stability | Turnover | Log Condition Number | Overall Score |",
    "|-------------------|-------------------|-----------------|----------|--------------------|--------------|"
  )
  
  for(i in 1:nrow(radar_data)) {
    row <- radar_data[i, ]
    summary_content <- c(summary_content,
                         paste("| ", row$CovMethod, " | ", 
                               ifelse(!is.na(row$PerformanceChange), 
                                      paste(ifelse(row$PerformanceChange >= 0, "+", ""), 
                                            round(row$PerformanceChange, 2), "%", sep = ""), "N/A"),
                               " | ", 
                               ifelse(!is.na(row$WeightStability), round(row$WeightStability, 2), "N/A"),
                               " | ", 
                               ifelse(!is.na(row$Turnover), paste(round(row$Turnover * 100, 1), "%", sep = ""), "N/A"),
                               " | ", 
                               ifelse(!is.na(row$LogConditionNumber), round(row$LogConditionNumber, 2), "N/A"),
                               " | ", 
                               round(row$OverallScore, 3),
                               " |", sep = "")
    )
  }
  
  # Add interpretation
  summary_content <- c(summary_content, "",
                       "## Interpretation",
                       "",
                       paste("- **Best Overall Method:** ", radar_data$CovMethod[1], " with a score of ", round(radar_data$OverallScore[1], 3), sep = ""),
                       paste("- **Most Robust to Stress:** ", radar_data$CovMethod[which.max(radar_data$NormPerformance)], 
                             " with performance change of ", 
                             ifelse(radar_data$PerformanceChange[which.max(radar_data$NormPerformance)] >= 0, "+", ""),
                             round(radar_data$PerformanceChange[which.max(radar_data$NormPerformance)], 2), "%", sep = ""),
                       paste("- **Most Weight-Stable:** ", radar_data$CovMethod[which.max(radar_data$NormStability)], 
                             " with stability score of ", round(radar_data$WeightStability[which.max(radar_data$NormStability)], 2), sep = ""),
                       paste("- **Best Numerical Conditioning:** ", radar_data$CovMethod[which.max(radar_data$NormConditioning)], 
                             " with log condition number of ", round(radar_data$LogConditionNumber[which.max(radar_data$NormConditioning)], 2), sep = "")
  )
  
  writeLines(summary_content, summary_file)
  
  cat("Combined visualizations and summary created in", combined_dir, "\n")
  
  # Return the combined metrics data
  return(radar_data)
}

# Example usage
create_example_code <- function(output_dir = "example_code") {
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create example code
  example_file <- file.path(output_dir, "example_usage.R")
  
  example_content <- c(
    "# Example Usage of Robust Portfolio Statistical Testing Framework",
    "",
    "# Load necessary libraries",
    "library(tseries)",
    "library(boot)",
    "library(car)",
    "library(robustbase)",
    "library(moments)",
    "library(ggplot2)",
    "library(reshape2)",
    "library(dplyr)",
    "",
    "# Load the statistical testing functions",
    "source('analysis/validation.R')",
    "",
    "# Option 1: Load existing analysis results from RData files",
    "load('sp500_analysis_results.RData')  # This loads your standard_results",
    "load('sp500_stress_analysis_results.RData')  # This loads your stress_results",
    "",
    "# Run the comprehensive statistical analysis",
    "test_results <- run_statistical_tests(",
    "  standard_results = analysis_results,  # Your standard market results",
    "  stress_results = stress_results,      # Your stress test results",
    "  test_types = c('performance', 'normality', 'stability', 'robustness', 'numerical', 'regime'),",
    "  bootstrap_replicates = 1000,",
    "  confidence_level = 0.95,",
    "  output_dir = 'statistical_analysis_results'",
    ")",
    "",
    "# Create combined visualizations",
    "combined_metrics <- create_combined_visuals(",
    "  analysis_results = test_results,",
    "  output_dir = 'statistical_analysis_results'",
    ")",
    "",
    "# Option 2: Run individual tests",
    "",
    "# Run only normality tests",
    "normality_results <- run_statistical_tests(",
    "  standard_results = analysis_results,",
    "  stress_results = stress_results,",
    "  test_types = 'normality',",
    "  output_dir = 'normality_tests'",
    ")",
    "",
    "# Run only weight stability analysis",
    "stability_results <- run_statistical_tests(",
    "  standard_results = analysis_results,",
    "  stress_results = stress_results,",
    "  test_types = 'stability',",
    "  output_dir = 'stability_tests'",
    ")",
    "",
    "# Run only performance comparison",
    "performance_results <- run_statistical_tests(",
    "  standard_results = analysis_results,",
    "  stress_results = stress_results,",
    "  test_types = c('performance', 'robustness'),",
    "  bootstrap_replicates = 2000,  # More replicates for higher precision",
    "  output_dir = 'performance_tests'",
    ")",
    "",
    "# How to work with the results",
    "",
    "# 1. Extract the most robust method",
    "if('comparison' %in% names(test_results)) {",
    "  most_robust_method <- test_results$comparison$data$Method[1]",
    "  cat('Most robust method:', most_robust_method, '\\n')",
    "}",
    "",
    "# 2. Check if performance degradation is statistically significant",
    "if('comparison' %in% names(test_results) && !is.null(test_results$comparison$paired_ttest)) {",
    "  p_value <- test_results$comparison$paired_ttest$p.value",
    "  is_significant <- p_value < 0.05",
    "  cat('Performance difference is', ifelse(is_significant, '', 'not'), 'statistically significant (p =', p_value, ')\\n')",
    "}",
    "",
    "# 3. Find the best covariance method for a specific regime",
    "if('regime_analysis' %in% names(test_results) && !is.null(test_results$regime_analysis$cov_wide)) {",
    "  regime_of_interest <- colnames(test_results$regime_analysis$cov_wide)[2]  # Skip the CovMethod column",
    "  best_cov_for_regime <- test_results$regime_analysis$cov_wide$CovMethod[",
    "    which.max(test_results$regime_analysis$cov_wide[[regime_of_interest]])",
    "  ]",
    "  cat('Best covariance method for', regime_of_interest, 'regime:', best_cov_for_regime, '\\n')",
    "}",
    "",
    "# 4. Compare non-normality across different portfolio methods",
    "if('normality' %in% names(test_results)) {",
    "  non_normal_pct <- aggregate(IsNormal ~ OptimizationType + Condition, ",
    "                             data = test_results$normality$tests,",
    "                             FUN = function(x) mean(!x))",
    "  print(non_normal_pct)",
    "}",
    "",
    "# 5. View weight stability rankings",
    "if('weight_stability' %in% names(test_results)) {",
    "  stability_ranking <- test_results$weight_stability$data[, c('Method', 'Turnover', 'StabilityScore')]",
    "  print(head(stability_ranking))",
    "}",
    "",
    "# 6. Extract bootstrap confidence intervals for Sharpe ratios",
    "if('bootstrap' %in% names(test_results)) {",
    "  sharpe_cis <- test_results$bootstrap[test_results$bootstrap$Metric == 'sharpe', ",
    "                                     c('Method', 'Condition', 'Original', 'CI_Lower', 'CI_Upper')]",
    "  print(head(sharpe_cis))",
    "}"
  )
  
  writeLines(example_content, example_file)
  
  cat("Example usage code created at", example_file, "\n")
  
  return(example_file)
}
