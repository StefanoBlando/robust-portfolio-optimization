# ========================================================================
# Portfolio Optimization Methods
# ========================================================================

# Equal weights portfolio (benchmark)
equal_weights_portfolio <- function(returns) {
  p <- ncol(returns)
  weights <- rep(1/p, p)
  
  cat("Creating equal weights portfolio with", p, "assets\n")
  
  return(list(
    weights = weights,
    method = "Equal-Weighted"
  ))
}

# Minimum variance portfolio using quadratic programming
min_var_portfolio <- function(returns, cov_matrix, max_weight = 0.2) {
  p <- ncol(returns)
  
  cat("Optimizing minimum variance portfolio with max weight =", max_weight, "\n")
  
  # Check if covariance matrix is valid
  if(any(is.na(cov_matrix)) || any(is.infinite(cov_matrix))) {
    warning("Invalid covariance matrix, using equal weights")
    return(equal_weights_portfolio(returns))
  }
  
  # Ensure positive definiteness
  eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  if(min(eig_vals) <= 0) {
    warning("Covariance matrix is not positive definite. Adding regularization.")
    # Add small value to diagonal
    cov_matrix <- cov_matrix + diag(1e-8, p)
  }
  
  # Try quadratic programming with constraints
  if(requireNamespace("quadprog", quietly = TRUE)) {
    tryCatch({
      # Set up quadratic programming problem
      Dmat <- 2 * cov_matrix
      dvec <- rep(0, p)
      
      # Constraints: 
      # 1. sum of weights = 1
      # 2. weights >= 0
      # 3. weights <= max_weight (diversification constraint)
      Amat <- cbind(rep(1, p), diag(p), -diag(p))
      bvec <- c(1, rep(0, p), rep(-max_weight, p))
      
      # Solve quadratic programming problem
      qp_solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
      weights <- qp_solution$solution
      
      # Ensure weights sum to 1 (handle numerical precision issues)
      weights <- weights / sum(weights)
      
      return(list(
        weights = weights,
        method = "Min-Variance (QP with constraints)"
      ))
    }, error = function(e) {
      warning("Constrained quadratic programming failed:", e$message)
      
      # Try again without max weight constraint
      tryCatch({
        # Simpler constraints: sum to 1 and non-negative
        Amat <- cbind(rep(1, p), diag(p))
        bvec <- c(1, rep(0, p))
        
        # Solve simpler QP
        qp_solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
        weights <- qp_solution$solution
        
        # Ensure weights sum to 1
        weights <- weights / sum(weights)
        
        # Apply weight caps as post-processing if needed
        if(max(weights) > max_weight) {
          exceeds <- weights > max_weight
          excess <- sum(weights[exceeds] - max_weight)
          weights[exceeds] <- max_weight
          weights[!exceeds] <- weights[!exceeds] + excess / sum(!exceeds)
          weights <- weights / sum(weights)  # Re-normalize
        }
        
        return(list(
          weights = weights,
          method = "Min-Variance (QP with post-processing)"
        ))
      }, error = function(e2) {
        warning("Basic quadratic programming failed:", e2$message)
        # Fallback to optimization
      })
    })
  }
  
  # Fallback to general optimization if QP fails
  result <- optim(
    par = rep(1/p, p),
    fn = function(w) {
      # Handle negative weights or weights that don't sum to 1
      w <- pmax(0, w)  # Ensure non-negative
      w <- w / sum(w)  # Normalize weights
      return(t(w) %*% cov_matrix %*% w)
    },
    method = "L-BFGS-B",
    lower = rep(0, p),
    upper = rep(max_weight * 1.2, p),  # Slightly higher max for optimization stability
    control = list(maxit = 1000)
  )
  
  # Check if optimization converged
  if(result$convergence != 0) {
    warning("L-BFGS-B optimization did not converge (code ", result$convergence, 
            "). ", result$message, ". Trying BFGS method.")
    
    # Try again with BFGS
    result <- optim(
      par = rep(1/p, p),
      fn = function(w) {
        w <- pmax(0, w)  # Ensure non-negative
        w <- w / sum(w)  # Normalize weights
        return(t(w) %*% cov_matrix %*% w)
      },
      method = "BFGS",
      control = list(maxit = 1000)
    )
    
    # If still no convergence, use equal weights
    if(result$convergence != 0) {
      warning("BFGS optimization also failed to converge. Using equal weights.")
      return(equal_weights_portfolio(returns))
    }
  }
  
  # Extract and normalize weights
  weights <- result$par
  weights <- pmax(0, weights)  # Ensure non-negative
  weights <- weights / sum(weights)  # Normalize
  
  # Apply weight caps
  if(max(weights) > max_weight) {
    exceeds <- weights > max_weight
    excess <- sum(weights[exceeds] - max_weight)
    weights[exceeds] <- max_weight
    if(sum(!exceeds) > 0) {
      weights[!exceeds] <- weights[!exceeds] + excess / sum(!exceeds)
    }
    weights <- weights / sum(weights)  # Re-normalize
  }
  
  return(list(
    weights = weights,
    method = "Min-Variance (BFGS with constraints)"
  ))
}

# Maximum Sharpe ratio portfolio using GA
max_sharpe_portfolio <- function(returns, cov_matrix, risk_free = 0, max_weight = 0.2, clean_mask = NULL) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  cat("Optimizing maximum Sharpe ratio portfolio using GA\n")
  
  # Validate inputs
  if(any(is.na(cov_matrix)) || any(is.infinite(cov_matrix))) {
    warning("Invalid covariance matrix, using equal weights")
    return(equal_weights_portfolio(returns))
  }
  
  # Ensure covariance matrix is positive definite
  eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
  if(min(eig_vals) <= 0) {
    warning("Covariance matrix is not positive definite. Adding regularization.")
    cov_matrix <- cov_matrix + diag(1e-8, p)
  }
  
  # Calculate mean returns with handling for contaminated data
  if(!is.null(clean_mask)) {
    # Use masked returns for mean estimation
    masked_returns <- returns * clean_mask
    
    # For contaminated values, replace with NA
    masked_returns[clean_mask == 0] <- NA
    
    # Compute robust means (ignoring NAs)
    mu <- apply(masked_returns, 2, mean, na.rm = TRUE)
    
    # Check for any columns with all NAs
    all_na_cols <- colSums(!is.na(masked_returns)) == 0
    if(any(all_na_cols)) {
      # For columns with all NAs, use overall mean
      overall_mean <- mean(returns, na.rm = TRUE)
      mu[all_na_cols] <- overall_mean
    }
  } else {
    # If no mask provided, use simple means
    mu <- colMeans(returns)
  }
  
  # Ensure GA package is available
  if(!requireNamespace("GA", quietly = TRUE)) {
    warning("GA package not available. Falling back to minimum variance portfolio.")
    return(min_var_portfolio(returns, cov_matrix, max_weight))
  }
  
  # Define Sharpe ratio function to maximize
  sharpe_fn <- function(w) {
    # Normalize weights to sum to 1
    w <- w / sum(w)
    
    # Calculate portfolio mean and standard deviation
    port_return <- sum(w * mu)
    port_sd <- sqrt(t(w) %*% cov_matrix %*% w)
    
    # Return Sharpe ratio (for maximization)
    if(port_sd <= 0 || is.na(port_sd)) {
      return(-999)  # Return very low value for invalid portfolios
    }
    return((port_return - risk_free) / port_sd)
  }
  
  # Run genetic algorithm optimization
  tryCatch({
    ga_result <- GA::ga(
      type = "real-valued",
      fitness = sharpe_fn,
      lower = rep(0, p),
      upper = rep(max_weight, p),
      maxiter = 100,
      run = 50,
      popSize = min(100, max(50, 5*p)),  # Adaptive population size
      pmutation = 0.2,                   # Higher mutation rate for diversity
      parallel = cores > 1,              # Use parallel if available
      seed = 42                          # For reproducibility
    )
    
    # Extract best solution
    weights <- ga_result@solution[1, ]
    
    # Ensure weights sum to 1
    weights <- weights / sum(weights)
    
    # Apply weight caps if needed
    if(max(weights) > max_weight) {
      exceeds <- weights > max_weight
      excess <- sum(weights[exceeds] - max_weight)
      weights[exceeds] <- max_weight
      if(sum(!exceeds) > 0) {
        weights[!exceeds] <- weights[!exceeds] + excess / sum(!exceeds)
      }
      weights <- weights / sum(weights)  # Re-normalize
    }
    
    # Calculate final Sharpe ratio for reporting
    port_return <- sum(weights * mu)
    port_sd <- sqrt(t(weights) %*% cov_matrix %*% weights)
    sharpe <- (port_return - risk_free) / port_sd
    
    cat("GA optimization complete. Portfolio Sharpe ratio:", round(sharpe, 4), "\n")
    
    return(list(
      weights = weights,
      method = "Max-Sharpe (GA)",
      sharpe = sharpe,
      port_return = port_return,
      port_sd = port_sd
    ))
  }, error = function(e) {
    warning("GA optimization failed:", e$message, ". Using min-variance portfolio instead.")
    return(min_var_portfolio(returns, cov_matrix, max_weight))
  })
}

# Minimum CVaR portfolio
min_cvar_portfolio <- function(returns, cov_matrix, clean_mask = NULL, alpha = 0.05, max_weight = 0.2) {
  p <- ncol(returns)
  
  # Function to calculate CVaR for a portfolio
  calc_cvar <- function(weights, rets, alpha = 0.05) {
    port_returns <- rets %*% weights
    var_alpha <- quantile(port_returns, alpha)
    cvar <- mean(port_returns[port_returns <= var_alpha])
    return(-cvar)  # Negative for minimization
  }
  
  # Prepare returns with clean mask if provided
  if (!is.null(clean_mask)) {
    effective_returns <- t(apply(returns, 1, function(row) {
      valid <- clean_mask[1,] == 1
      if (sum(valid) == 0) return(row)
      row
    }))
  } else {
    effective_returns <- returns
  }
  
  # Optimize
  result <- try({
    optim(
      par = rep(1/p, p),
      fn = function(w) {
        w <- w / sum(w)  # Normalize
        calc_cvar(w, effective_returns, alpha)
      },
      method = "L-BFGS-B",
      lower = rep(0, p),
      upper = rep(max_weight, p),
      control = list(maxit = 1000)
    )
  }, silent = TRUE)
  
  if (inherits(result, "try-error") || !result$convergence == 0) {
    warning("CVaR optimization failed, using equal weights")
    return(equal_weights_portfolio(returns))
  }
  
  weights <- result$par / sum(result$par)
  
  # Apply weight caps
  if(max(weights) > max_weight) {
    exceeds <- weights > max_weight
    excess <- sum(weights[exceeds] - max_weight)
    weights[exceeds] <- max_weight
    if(sum(!exceeds) > 0) {
      weights[!exceeds] <- weights[!exceeds] + excess / sum(!exceeds)
    }
    weights <- weights / sum(weights)  # Re-normalize
  }
  
  return(list(
    weights = weights,
    method = "MinCVaR"
  ))
}

# Minimum variance using Differential Evolution
de_portfolio <- function(returns, cov_matrix, max_weight = 0.2) {
  p <- ncol(returns)
  
  if (!requireNamespace("DEoptim", quietly = TRUE)) {
    warning("DEoptim package not available, using L-BFGS-B")
    return(min_var_portfolio(returns, cov_matrix, max_weight))
  }
  
  tryCatch({
    de_result <- DEoptim::DEoptim(
      fn = function(w) {
        w <- w / sum(w)  # Normalize
        as.numeric(t(w) %*% cov_matrix %*% w)
      },
      lower = rep(0, p),
      upper = rep(max_weight, p),
      control = DEoptim::DEoptim.control(
        NP = 10*p,
        itermax = 100,
        trace = FALSE
      )
    )
    
    weights <- de_result$optim$bestmem
    weights <- weights / sum(weights)
    
    # Apply weight caps
    if(max(weights) > max_weight) {
      exceeds <- weights > max_weight
      excess <- sum(weights[exceeds] - max_weight)
      weights[exceeds] <- max_weight
      if(sum(!exceeds) > 0) {
        weights[!exceeds] <- weights[!exceeds] + excess / sum(!exceeds)
      }
      weights <- weights / sum(weights)  # Re-normalize
    }
    
    return(list(
      weights = weights,
      method = "DE-MinVar"
    ))
  }, error = function(e) {
    warning("DE optimization failed, using L-BFGS-B")
    return(min_var_portfolio(returns, cov_matrix, max_weight))
  })
}

# ========================================================================
# Hierarchical Risk Parity Implementation
# ========================================================================

# Distance calculation between assets based on correlation
get_distance_matrix <- function(corr_matrix) {
  # Convert correlation to distance: sqrt(0.5*(1-correlation))
  dist_matrix <- sqrt(0.5 * (1 - corr_matrix))
  return(dist_matrix)
}

# Hierarchical clustering of assets based on distance
get_cluster_hierarchy <- function(dist_matrix) {
  # Perform hierarchical clustering using the distance matrix
  hc <- hclust(as.dist(dist_matrix), method = "complete")
  return(hc)
}

# Quasi-diagonalization of covariance matrix based on hierarchical clustering
get_quasi_diag <- function(hc) {
  # Reorder assets based on clustering
  asset_order <- hc$order
  return(asset_order)
}

# Recursive bisection for weight allocation
get_cluster_var <- function(cov_matrix, cluster_indices) {
  # Get subset of covariance matrix for the cluster
  cluster_cov <- cov_matrix[cluster_indices, cluster_indices, drop = FALSE]
  
  # Calculate cluster variance (portfolio variance with equal weights)
  weights <- rep(1/length(cluster_indices), length(cluster_indices))
  cluster_var <- t(weights) %*% cluster_cov %*% weights
  
  return(as.numeric(cluster_var))
}

# Recursive bisection to allocate weights based on inverse variance
get_recursive_bisection <- function(cov_matrix, sorted_indices, weights) {
  # Base case: if only one asset, return
  if(length(sorted_indices) <= 1) {
    return(weights)
  }
  
  # Perform recursive bisection
  # Split cluster in two
  n <- length(sorted_indices)
  mid <- ceiling(n/2)
  left_indices <- sorted_indices[1:mid]
  right_indices <- sorted_indices[(mid+1):n]
  
  # Get variance of clusters
  left_var <- get_cluster_var(cov_matrix, left_indices)
  right_var <- get_cluster_var(cov_matrix, right_indices)
  
  # Allocate weights inversely proportional to cluster variance
  alpha <- 1 - left_var/(left_var + right_var)
  
  # Update weights
  weights[left_indices] <- weights[left_indices] * alpha
  weights[right_indices] <- weights[right_indices] * (1-alpha)
  
  # Recurse
  if(length(left_indices) > 1)
    weights <- get_recursive_bisection(cov_matrix, left_indices, weights)
  if(length(right_indices) > 1)
    weights <- get_recursive_bisection(cov_matrix, right_indices, weights)
  
  return(weights)
}

# Main Hierarchical Risk Parity function
hrp_portfolio <- function(returns, cov_matrix, max_weight = 0.2) {
  p <- ncol(returns)
  
  # Calculate correlation matrix from covariance
  diag_inv <- diag(1/sqrt(diag(cov_matrix)))
  corr_matrix <- diag_inv %*% cov_matrix %*% diag_inv
  
  # Get distance matrix
  dist_matrix <- get_distance_matrix(corr_matrix)
  
  # Get clusters
  hc <- get_cluster_hierarchy(dist_matrix)
  
  # Get sorted indices
  sorted_indices <- get_quasi_diag(hc)
  
  # Initialize weights
  weights <- rep(1/p, p)
  
  # Get final weights using recursive bisection
  weights <- get_recursive_bisection(cov_matrix, sorted_indices, weights)
  
  # Apply weight caps if needed
  if(max(weights) > max_weight) {
    exceeds <- weights > max_weight
    excess <- sum(weights[exceeds] - max_weight)
    weights[exceeds] <- max_weight
    if(sum(!exceeds) > 0) {
      weights[!exceeds] <- weights[!exceeds] + excess / sum(!exceeds)
    }
    weights <- weights / sum(weights)  # Re-normalize
  }
  
  return(list(
    weights = weights,
    method = "HRP"
  ))
}

# ========================================================================
# Portfolio Weight Stability Analysis
# ========================================================================

# Calculate turnover between two sets of weights
calculate_turnover <- function(weights1, weights2) {
  # Total sum of absolute weight changes
  turnover <- sum(abs(weights1 - weights2)) / 2
  return(turnover)
}

# Calculate weight stability metrics
calculate_weight_stability <- function(all_weights) {
  # all_weights should be a matrix where each row represents a simulation
  # and each column represents an asset
  
  # Check if input is a matrix
  if(!is.matrix(all_weights)) {
    # Convert to matrix
    all_weights <- matrix(all_weights, ncol = length(all_weights))
  }
  
  # Average absolute deviation from mean weights
  mean_weights <- colMeans(all_weights)
  mad_weights <- mean(colMeans(abs(all_weights - matrix(mean_weights, 
                                                        nrow = nrow(all_weights), 
                                                        ncol = ncol(all_weights), 
                                                        byrow = TRUE))))
  
  # Average pairwise turnover
  n_sims <- nrow(all_weights)
  total_turnover <- 0
  count <- 0
  
  if(n_sims > 1) {
    for(i in 1:(n_sims-1)) {
      for(j in (i+1):n_sims) {
        total_turnover <- total_turnover + calculate_turnover(all_weights[i,], all_weights[j,])
        count <- count + 1
      }
    }
    avg_turnover <- total_turnover / count
  } else {
    avg_turnover <- 0
  }
  
  # Max absolute weight change for any asset across simulations
  max_weight_change <- max(apply(all_weights, 2, function(x) max(x) - min(x)))
  
  # Return stability metrics
  return(list(
    mad_weights = mad_weights,             # Mean absolute deviation
    avg_turnover = avg_turnover,           # Average pairwise turnover
    max_weight_change = max_weight_change  # Maximum weight change
  ))
}

# ========================================================================
# Portfolio Evaluation Functions
# ========================================================================

# Evaluate portfolio performance
evaluate_portfolio <- function(weights, train_returns, test_returns, 
                               train_market = NULL, test_market = NULL) {
  # Calculate portfolio returns
  train_port_returns <- calculate_portfolio_returns(weights, train_returns)
  test_port_returns <- calculate_portfolio_returns(weights, test_returns)
  
  # Calculate performance metrics
  train_metrics <- calculate_risk_measures(train_port_returns)
  test_metrics <- calculate_risk_measures(test_port_returns)
  
  # Calculate risk contribution
  risk_contrib <- tryCatch({
    # Get covariance of test returns
    test_cov <- cov(test_returns)
    
    # Calculate portfolio variance
    port_var <- as.numeric(t(weights) %*% test_cov %*% weights)
    
    if(port_var > 0) {
      # Calculate marginal contributions
      marginal_contrib <- (test_cov %*% weights) / sqrt(port_var)
      contrib <- weights * marginal_contrib
      
      # Normalize to percentage
      contrib / sum(contrib) * 100
    } else {
      rep(NA, length(weights))
    }
  }, error = function(e) {
    # Return NA if calculation fails
    rep(NA, length(weights))
  })
  
  # Return comprehensive evaluation
  list(
    train = train_metrics,
    test = test_metrics,
    train_returns = train_port_returns,
    test_returns = test_port_returns,
    risk_contribution = risk_contrib
  )
}

# Create and evaluate portfolios using different methods
create_and_evaluate_portfolios <- function(train_returns, test_returns, 
                                           train_market = NULL, test_market = NULL,
                                           max_weight = 0.2, risk_free = 0,
                                           test_dates = NULL, test_regimes = NULL) {
  # Define covariance estimation methods
  cov_methods <- list(
    Sample = function(r) sample_cov(r),
    MCD = function(r) mcd_cov(r),
    HybridMCD_PCA = function(r) hybrid_mcd_pca(r),
    HybridMCD_RobPCA = function(r) hybrid_mcd_robpca(r),
    HybridMCD_RobPCA_GLasso = function(r) hybrid_mcd_robpca_glasso(r),
    Tyler = function(r) tyler_m_cov(r),
    AdaptiveLW = function(r) adaptive_lw_cov(r),
    HFBRE = function(r) hybrid_factor_robust_cov(r)
  )
  
  # Initialize results
  results <- list()
  
  # Equal weights portfolio (benchmark)
  equal_weights <- equal_weights_portfolio(train_returns)
  results$EqualWeights <- list(
    weights = equal_weights$weights,
    method = equal_weights$method,
    performance = evaluate_portfolio(equal_weights$weights, train_returns, test_returns, 
                                     train_market, test_market)
  )
  
  # Add regime-specific performance evaluation for equal weights if regimes are provided
  if(!is.null(test_dates) && !is.null(test_regimes)) {
    results$EqualWeights$regime_performance <- evaluate_portfolio_by_regime(
      equal_weights$weights, test_returns, test_dates, test_regimes, test_market
    )
  }
  
  # For each covariance estimation method
  for(method_name in names(cov_methods)) {
    cat("\n---------------------------------------------\n")
    cat("Creating portfolio using", method_name, "method...\n")
    
    # Estimate covariance matrix
    cov_matrix <- cov_methods[[method_name]](train_returns)
    
    # Verify covariance matrix is valid
    if(any(is.na(cov_matrix)) || any(is.infinite(cov_matrix))) {
      warning("Method", method_name, "produced invalid covariance matrix. Skipping.")
      next
    }
    
    # Check positive definiteness
    eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
    if(min(eig_vals) <= 0) {
      warning("Method", method_name, "produced non-positive definite matrix. Fixing...")
      
      # Fix by ensuring minimum eigenvalue
      eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
      values <- pmax(eigen_decomp$values, 1e-8)
      cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
    }
    
    # Calculate and report condition number
    cond_num <- max(eig_vals) / min(eig_vals)
    cat("Condition number of", method_name, "covariance matrix:", format(cond_num, scientific = TRUE), "\n")
    
    # Detect outliers for Max Sharpe optimization (optional)
    outliers <- tryCatch({
      detect_outliers(train_returns, method = "MAD", threshold = 3)
    }, error = function(e) {
      warning("Outlier detection failed:", e$message)
      NULL
    })
    
    clean_mask <- if(!is.null(outliers)) outliers$outlier_mask else NULL
    
    # Min Variance portfolio
    min_var_portfolio_result <- min_var_portfolio(train_returns, cov_matrix, max_weight)
    
    # Max Sharpe portfolio
    max_sharpe_portfolio_result <- max_sharpe_portfolio(
      train_returns, cov_matrix, risk_free, max_weight, clean_mask
    )
    
    # Min CVaR portfolio
    min_cvar_portfolio_result <- min_cvar_portfolio(
      train_returns, cov_matrix, clean_mask, 0.05, max_weight
    )
    
    # HRP portfolio
    hrp_portfolio_result <- hrp_portfolio(train_returns, cov_matrix, max_weight)
    
    # Evaluate performance for all portfolios
    min_var_performance <- evaluate_portfolio(
      min_var_portfolio_result$weights, train_returns, test_returns, 
      train_market, test_market
    )
    
    max_sharpe_performance <- evaluate_portfolio(
      max_sharpe_portfolio_result$weights, train_returns, test_returns, 
      train_market, test_market
    )
    
    min_cvar_performance <- evaluate_portfolio(
      min_cvar_portfolio_result$weights, train_returns, test_returns,
      train_market, test_market
    )
    
    hrp_performance <- evaluate_portfolio(
      hrp_portfolio_result$weights, train_returns, test_returns,
      train_market, test_market
    )
    
    # Add regime-specific evaluation if regimes are provided
    if(!is.null(test_dates) && !is.null(test_regimes)) {
      # Add regime evaluation for each portfolio type
      min_var_regime_performance <- evaluate_portfolio_by_regime(
        min_var_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
      
      max_sharpe_regime_performance <- evaluate_portfolio_by_regime(
        max_sharpe_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
      
      min_cvar_regime_performance <- evaluate_portfolio_by_regime(
        min_cvar_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
      
      hrp_regime_performance <- evaluate_portfolio_by_regime(
        hrp_portfolio_result$weights, test_returns, test_dates, test_regimes, test_market
      )
    } else {
      min_var_regime_performance <- NULL
      max_sharpe_regime_performance <- NULL
      min_cvar_regime_performance <- NULL
      hrp_regime_performance <- NULL
    }
    
    # Store results for all portfolio types
    results[[paste0(method_name, "_MinVar")]] <- list(
      weights = min_var_portfolio_result$weights,
      method = paste(method_name, min_var_portfolio_result$method),
      cov_matrix = cov_matrix,
      performance = min_var_performance,
      condition_number = cond_num,
      regime_performance = min_var_regime_performance
    )
    
    results[[paste0(method_name, "_MaxSharpe")]] <- list(
      weights = max_sharpe_portfolio_result$weights,
      method = paste(method_name, max_sharpe_portfolio_result$method),
      cov_matrix = cov_matrix,
      performance = max_sharpe_performance,
      condition_number = cond_num,
      regime_performance = max_sharpe_regime_performance
    )
    
    results[[paste0(method_name, "_MinCVaR")]] <- list(
      weights = min_cvar_portfolio_result$weights,
      method = paste(method_name, min_cvar_portfolio_result$method),
      cov_matrix = cov_matrix,
      performance = min_cvar_performance,
      condition_number = cond_num,
      regime_performance = min_cvar_regime_performance
    )
    
    results[[paste0(method_name, "_HRP")]] <- list(
      weights = hrp_portfolio_result$weights,
      method = paste(method_name, hrp_portfolio_result$method),
      cov_matrix = cov_matrix,
      performance = hrp_performance,
      condition_number = cond_num,
      regime_performance = hrp_regime_performance
    )
    
    # Report key metrics
    cat("Min Variance - In-sample Sharpe:", round(min_var_performance$train$sharpe, 4), "\n")
    cat("Min Variance - Out-of-sample Sharpe:", round(min_var_performance$test$sharpe, 4), "\n")
    cat("Max Sharpe - In-sample Sharpe:", round(max_sharpe_performance$train$sharpe, 4), "\n")
    cat("Max Sharpe - Out-of-sample Sharpe:", round(max_sharpe_performance$test$sharpe, 4), "\n")
    cat("Min CVaR - In-sample Sharpe:", round(min_cvar_performance$train$sharpe, 4), "\n")
    cat("Min CVaR - Out-of-sample Sharpe:", round(min_cvar_performance$test$sharpe, 4), "\n")
    cat("HRP - In-sample Sharpe:", round(hrp_performance$train$sharpe, 4), "\n")
    cat("HRP - Out-of-sample Sharpe:", round(hrp_performance$test$sharpe, 4), "\n")
    
    # Print regime-specific performance if available
    if(!is.null(min_var_regime_performance)) {
      cat("\nMin Variance regime-specific Sharpe ratios:\n")
      for(regime in names(min_var_regime_performance)) {
        if(regime != "overall") {
          cat("  ", regime, ":", round(min_var_regime_performance[[regime]]$sharpe, 4), "\n")
        }
      }
    }
  }
  
  return(results)
}

# Evaluate portfolio performance by market regime
evaluate_portfolio_by_regime <- function(weights, returns, dates, regimes, market_returns = NULL) {
  # Initialize results list
  regime_results <- list()
  
  # Get unique regimes
  unique_regimes <- unique(regimes)
  
  # For each regime
  for(regime in unique_regimes) {
    # Get indices for this regime
    regime_idx <- which(regimes == regime)
    
    # Skip if not enough data
    if(length(regime_idx) < 10) {
      cat("Skipping regime", regime, "- insufficient data (", length(regime_idx), "observations)\n")
      next
    }
    
    # Extract regime-specific data
    regime_returns <- returns[regime_idx, , drop = FALSE]
    regime_dates <- dates[regime_idx]
    regime_market <- if(!is.null(market_returns)) market_returns[regime_idx] else NULL
    
    # Calculate portfolio returns for this regime
    port_returns <- calculate_portfolio_returns(weights, regime_returns)
    
    # Calculate performance metrics
    metrics <- calculate_risk_measures(port_returns)
    
    # Store results
    regime_results[[regime]] <- metrics
  }
  
  # Add overall metrics
  regime_results[["overall"]] <- calculate_risk_measures(
    calculate_portfolio_returns(weights, returns)
  )
  
  return(regime_results)
}
