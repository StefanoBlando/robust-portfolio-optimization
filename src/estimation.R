# ========================================================================
# Covariance Estimation Methods
# ========================================================================

# Sample covariance estimation with proper conditioning checks
sample_cov <- function(returns) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  cat("Computing sample covariance matrix (", n, "x", p, ")\n", sep="")
  
  tryCatch({
    # Compute standard sample covariance
    cov_mat <- cov(returns)
    
    # Check condition number
    eig_vals <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
    cond_num <- max(eig_vals) / min(eig_vals)
    
    if(cond_num > 1e6 || min(eig_vals) <= 0) {
      warning("Sample covariance matrix is ill-conditioned (cond = ", 
              format(cond_num, scientific = TRUE), 
              "). Adding minimal regularization.")
      
      # Add minimal shrinkage to ensure conditioning
      diag_mat <- diag(diag(cov_mat))
      shrinkage <- 0.05  # Small shrinkage to maintain most of the structure
      cov_mat <- (1 - shrinkage) * cov_mat + shrinkage * diag_mat
      
      # Verify improvement
      eig_vals_new <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
      cond_num_new <- max(eig_vals_new) / min(eig_vals_new)
      
      cat("Condition number improved from", format(cond_num, scientific = TRUE), 
          "to", format(cond_num_new, scientific = TRUE), "\n")
    }
    
    return(cov_mat)
  }, error = function(e) {
    warning("Error in sample covariance calculation:", e$message)
    # Return diagonal matrix as fallback but with proper scaling
    diag_var <- apply(returns, 2, var, na.rm = TRUE)
    # Ensure all variances are positive
    diag_var <- pmax(diag_var, 1e-8)
    return(diag(diag_var))
  })
}

# MCD robust covariance estimation with proper implementation
mcd_cov <- function(returns, alpha = NULL) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  # Determine appropriate alpha based on data dimensions
  if(is.null(alpha)) {
    # Use higher alpha for high-dimensional problems
    alpha <- min(0.875, max(0.5, 0.75 + 0.125 * (p/n)))
  }
  
  cat("Computing MCD covariance matrix with alpha =", alpha, "\n")
  
  # For very high-dimensional cases, MCD becomes unreliable
  if(n <= p * 1.2) {
    warning("Insufficient observations for reliable MCD (n ≤ 1.2p). Using robust diagonal covariance.")
    
    # Use robust diagonal covariance as a reliable alternative
    robust_var <- apply(returns, 2, function(x) {
      med_x <- median(x, na.rm = TRUE)
      mad_x <- median(abs(x - med_x), na.rm = TRUE)
      if(mad_x <= 0) mad_x <- sd(x, na.rm = TRUE) / 1.4826
      if(mad_x <= 0) mad_x <- 1e-8
      return((mad_x * 1.4826)^2)  # Convert MAD to variance
    })
    
    # Use robust pairwise correlations
    cors <- tryCatch({
      # Use spatial sign correlations for robustness
      robustbase::covSS(returns, cor = TRUE)$cov
    }, error = function(e) {
      # Fall back to Spearman correlations
      cor(returns, method = "spearman")
    })
    
    # Create covariance from correlations and robust variances
    cov_matrix <- diag(sqrt(robust_var)) %*% cors %*% diag(sqrt(robust_var))
    
    # Ensure positive definiteness
    eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
    values <- pmax(eigen_decomp$values, 1e-8)
    cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
    
    return(cov_matrix)
  }
  
  tryCatch({
    # Use deterministic MCD for reproducibility
    mcd_result <- rrcov::CovMcd(returns, alpha = alpha, nsamp = "deterministic")
    
    # Extract covariance matrix
    cov_matrix <- mcd_result@cov
    
    # Check condition number of MCD covariance
    eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
    cond_num <- max(eig_vals) / min(eig_vals)
    
    if(cond_num > 1e6 || min(eig_vals) <= 0) {
      warning("MCD covariance matrix is ill-conditioned (cond = ", 
              format(cond_num, scientific = TRUE), 
              "). Adding minimal regularization.")
      
      # Add small regularization to ensure positive definiteness
      diag_mat <- diag(diag(cov_matrix))
      shrinkage <- 0.05  # Small shrinkage
      cov_matrix <- (1 - shrinkage) * cov_matrix + shrinkage * diag_mat
      
      # Verify improvement in condition number
      eig_vals_new <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
      cond_num_new <- max(eig_vals_new) / min(eig_vals_new)
      
      cat("Condition number improved from", format(cond_num, scientific = TRUE), 
          "to", format(cond_num_new, scientific = TRUE), "\n")
    }
    
    return(cov_matrix)
  }, error = function(e) {
    warning("MCD estimation failed:", e$message)
    
    # Try alternative robust estimator
    tryCatch({
      # Try S-estimator which can be more stable
      s_est <- rrcov::CovSest(returns, method = "bisquare")
      return(s_est@cov)
    }, error = function(e2) {
      warning("S-estimator also failed:", e2$message, "Using robust diagonal covariance.")
      
      # Use diagonal robust covariance as ultimate fallback
      robust_var <- apply(returns, 2, function(x) {
        med_x <- median(x, na.rm = TRUE)
        mad_x <- median(abs(x - med_x), na.rm = TRUE)
        if(mad_x <= 0) mad_x <- sd(x, na.rm = TRUE) / 1.4826
        if(mad_x <= 0) mad_x <- 1e-8
        return((mad_x * 1.4826)^2)
      })
      
      return(diag(robust_var))
    })
  })
}

# Hybrid MCD + PCA covariance estimation
hybrid_mcd_pca <- function(returns, k = NULL, alpha = NULL) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  # Determine number of factors adaptively if not provided
  if(is.null(k)) {
    # Use parallel analysis to determine optimal k
    tryCatch({
      # Compute eigenvalues of original data
      cov_mat <- cov(returns, use = "pairwise.complete.obs")
      orig_eig <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
      
      # Generate random data with same dimensions
      n_simulations <- 10
      random_eigs <- matrix(0, n_simulations, p)
      
      for(i in 1:n_simulations) {
        random_data <- matrix(rnorm(n*p), nrow=n, ncol=p)
        random_cov <- cov(random_data)
        random_eigs[i,] <- eigen(random_cov, symmetric = TRUE, only.values = TRUE)$values
      }
      
      # Average the random eigenvalues
      mean_random_eigs <- colMeans(random_eigs)
      
      # Find where original eigenvalues > random eigenvalues
      k <- sum(orig_eig > mean_random_eigs)
      
      # Ensure k is reasonable
      k <- min(k, ceiling(p/4), floor(n/5))
      k <- max(k, 3)  # At least 3 factors
    }, error = function(e) {
      # Simple fallback
      k <- min(ceiling(p/10), floor(n/5), 20)
      k <- max(k, 3)
    })
    
    cat("Selected optimal number of factors k =", k, "\n")
  }
  
  # Determine appropriate alpha for MCD
  if(is.null(alpha)) {
    alpha <- min(0.875, max(0.5, 0.75 + 0.125 * (k/n)))
  }
  
  cat("Computing Hybrid MCD+PCA with k =", k, "and alpha =", alpha, "\n")
  
  tryCatch({
    # Step 1: PCA for dimensionality reduction with proper centering
    returns_centered <- scale(returns, center = TRUE, scale = FALSE)
    pca_result <- prcomp(returns_centered, center = FALSE, scale. = FALSE)
    
    # Extract factors and loadings
    factors <- pca_result$x[, 1:k, drop = FALSE]
    loadings <- pca_result$rotation[, 1:k, drop = FALSE]
    
    # Step 2: Apply MCD on factors
    if(n > 2*k) {  # Check if we have enough observations for MCD
      cat("Applying MCD to", k, "factors\n")
      mcd_result <- tryCatch({
        rrcov::CovMcd(factors, alpha = alpha)
      }, error = function(e) {
        warning("MCD on factors failed: ", e$message)
        # Fall back to standard covariance if MCD fails
        list(center = colMeans(factors), cov = cov(factors))
      })
      
      # Extract factor covariance matrix
      if("cov" %in% names(mcd_result)) {
        factor_cov <- mcd_result$cov
      } else {
        factor_cov <- mcd_result@cov
      }
    } else {
      # Not enough data for MCD - use robust correlations
      cat("Insufficient observations for MCD on factors, using robust correlation\n")
      robust_var <- apply(factors, 2, function(x) {
        mad_x <- median(abs(x - median(x)), na.rm = TRUE)
        return((mad_x * 1.4826)^2)
      })
      
      cors <- cor(factors, method = "spearman")  # Robust correlation
      factor_cov <- diag(sqrt(robust_var)) %*% cors %*% diag(sqrt(robust_var))
    }
    
    # Step 3: Calculate residuals and idiosyncratic variances robustly
    reconstructed <- factors %*% t(loadings)
    residuals <- returns_centered - reconstructed
    
    # Step 4: Detect outliers in residuals
    outlier_mask <- matrix(1, nrow = n, ncol = p)
    residual_threshold <- 3.0  # Fixed threshold for residuals
    
    for(j in 1:p) {
      res_j <- residuals[, j]
      med_j <- median(res_j, na.rm = TRUE)
      mad_j <- mad(res_j, na.rm = TRUE)
      if(mad_j <= 0) mad_j <- sd(res_j, na.rm = TRUE) / 1.4826
      if(mad_j <= 0) mad_j <- 1e-8
      
      # Flag outliers
      z_scores <- abs(res_j - med_j) / mad_j
      outliers <- z_scores > residual_threshold
      outlier_mask[outliers, j] <- 0
    }
    
    # Step 5: Compute robust residual variances with outlier handling
    D <- diag(p)
    for(j in 1:p) {
      # Use clean residuals for variance estimation
      res_clean <- residuals[outlier_mask[, j] == 1, j]
      if(length(res_clean) < 10) res_clean <- residuals[, j]  # Fall back if too few clean points
      
      # Robust variance estimate
      med_j <- median(res_clean, na.rm = TRUE)
      mad_j <- median(abs(res_clean - med_j), na.rm = TRUE)
      if(mad_j <= 0) mad_j <- sd(res_clean, na.rm = TRUE) / 1.4826
      if(mad_j <= 0) mad_j <- 1e-8
      
      D[j, j] <- (mad_j * 1.4826)^2
    }
    
    # Step 6: Reconstruct covariance matrix
    cov_matrix <- loadings %*% factor_cov %*% t(loadings) + D
    
    # Step 7: Ensure positive definiteness
    eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
    values <- pmax(eigen_decomp$values, 1e-8)  # Ensure minimum eigenvalue
    cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
    
    return(cov_matrix)
  }, error = function(e) {
    warning("Hybrid MCD+PCA failed:", e$message)
    # Return MCD covariance as fallback
    return(mcd_cov(returns, alpha))
  })
}

# Hybrid MCD + Robust PCA covariance estimation
hybrid_mcd_robpca <- function(returns, k = NULL, alpha = NULL, threshold = 3.0) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  # Determine number of factors adaptively if not provided
  if(is.null(k)) {
    tryCatch({
      # Try to use pcaPP for robust factor selection if available
      if(requireNamespace("pcaPP", quietly = TRUE)) {
        # Use robust scaling 
        scaled_returns <- scale(returns, 
                                center = apply(returns, 2, median, na.rm = TRUE),
                                scale = apply(returns, 2, mad, na.rm = TRUE))
        
        # Compute robust PCA with PCAgrid
        pca_result <- pcaPP::PCAgrid(scaled_returns, k = min(20, p/2, n/3), 
                                     method = "qn", center = FALSE)
        
        # Examine eigenvalue decay for elbow point
        evals <- pca_result$eigenvalues
        decay_rate <- -diff(evals) / evals[-length(evals)]
        
        # Find significant drop in eigenvalues
        k <- which(decay_rate > mean(decay_rate))[1]
        if(is.na(k)) k <- min(which(cumsum(evals) / sum(evals) >= 0.8))
      } else {
        # Detect outliers with fixed threshold
        outliers <- detect_outliers(returns, method = "MAD", threshold = threshold)
        
        # Apply winsorization
        cleaned_returns <- winsorize_returns(returns, outliers$outlier_mask, method = "shrinkage")
        
        # Run PCA on cleaned data
        pca_result <- prcomp(cleaned_returns, center = TRUE, scale. = FALSE)
        var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
        
        # Use enough components to explain 80% of variance
        k <- min(which(cumsum(var_explained) >= 0.8))
      }
      
      # Ensure k is reasonable
      k <- min(k, ceiling(p/4), floor(n/5))
      k <- max(k, 3)  # At least 3 factors
    }, error = function(e) {
      # Fallback k selection
      warning("Automatic factor selection failed:", e$message)
      k <- min(ceiling(p/10), floor(n/5), 15)
      k <- max(k, 3)
    })
    
    cat("Selected optimal number of factors k =", k, "\n")
  }
  
  # Determine appropriate alpha based on factor dimensionality
  if(is.null(alpha)) {
    alpha <- min(0.875, max(0.5, 0.75 + 0.125 * (k/n)))
  }
  
  cat("Computing Hybrid MCD+RobPCA with k =", k, "and alpha =", alpha, "\n")
  
  tryCatch({
    # Step 1: Perform preliminary outlier filtering
    outlier_result <- detect_outliers(returns, method = "MAD", threshold = threshold)
    clean_mask <- outlier_result$outlier_mask
    
    # Step 2: Apply winsorization to reduce influence of outliers
    winsorized_returns <- winsorize_returns(returns, clean_mask, method = "shrinkage")
    
    # Step 3: Apply robust PCA
    rob_pca <- NULL
    
    if(requireNamespace("pcaPP", quietly = TRUE)) {
      # Use pcaPP for robust PCA if available
      rob_pca <- tryCatch({
        pca_result <- pcaPP::PCAgrid(winsorized_returns, k = k, method = "qn", center = TRUE)
        list(scores = pca_result$scores, loadings = pca_result$loadings)
      }, error = function(e) {
        warning("PCAgrid failed:", e$message)
        # Fall back to standard PCA on winsorized data
        NULL
      })
    }
    
    # Fall back to standard PCA on winsorized data if pcaPP fails or is not available
    if(is.null(rob_pca)) {
      pca_result <- prcomp(winsorized_returns, center = TRUE, scale. = FALSE)
      rob_pca <- list(
        scores = pca_result$x[, 1:k, drop = FALSE],
        loadings = pca_result$rotation[, 1:k, drop = FALSE]
      )
    }
    
    # Extract factors and loadings
    factors <- rob_pca$scores
    loadings <- rob_pca$loadings
    
    # Step 4: Apply robust covariance estimation on factors
    if(n > 2*k) {  # Check if we have enough observations for MCD
      mcd_result <- tryCatch({
        rrcov::CovMcd(factors, alpha = alpha)
      }, error = function(e) {
        warning("MCD on factors failed:", e$message)
        # Fall back to standard covariance
        list(center = colMeans(factors), cov = cov(factors))
      })
      
      # Extract factor covariance matrix
      if("cov" %in% names(mcd_result)) {
        factor_cov <- mcd_result$cov
      } else {
        factor_cov <- mcd_result@cov
      }
    } else {
      # If not enough data for MCD, use robust pairwise correlations
      robust_var <- apply(factors, 2, function(x) {
        mad_x <- median(abs(x - median(x)), na.rm = TRUE)
        if(mad_x <= 0) mad_x <- sd(x, na.rm = TRUE) / 1.4826
        if(mad_x <= 0) mad_x <- 1e-8
        return((mad_x * 1.4826)^2)
      })
      
      # Use Spearman correlations for robustness
      cors <- cor(factors, method = "spearman")
      factor_cov <- diag(sqrt(robust_var)) %*% cors %*% diag(sqrt(robust_var))
    }
    
    # Step 5: Calculate residuals
    reconstructed <- factors %*% t(loadings)
    centered_returns <- scale(returns, center = TRUE, scale = FALSE)
    residuals <- centered_returns - reconstructed
    
    # Step 6: Detect outliers in residuals with fixed threshold
    residual_mask <- matrix(1, nrow = n, ncol = p)
    
    for(j in 1:p) {
      res_j <- residuals[, j]
      med_j <- median(res_j, na.rm = TRUE)
      mad_j <- mad(res_j, na.rm = TRUE)
      if(mad_j <= 0) mad_j <- sd(res_j, na.rm = TRUE) / 1.4826
      if(mad_j <= 0) mad_j <- 1e-8
      
      # Use fixed threshold for residuals
      outliers_j <- abs(res_j - med_j) > threshold * mad_j
      residual_mask[outliers_j, j] <- 0
    }
    
    # Step 7: Apply winsorization to residuals 
    winsorized_residuals <- winsorize_returns(residuals, residual_mask, method = "shrinkage")
    
    # Step 8: Compute robust residual variances
    D <- diag(p)
    for(j in 1:p) {
      res_j <- winsorized_residuals[, j]
      
      # Robust variance estimate
      mad_j <- median(abs(res_j - median(res_j)), na.rm = TRUE)
      if(mad_j <= 0) mad_j <- sd(res_j, na.rm = TRUE) / 1.4826
      if(mad_j <= 0) mad_j <- 1e-8
      
      D[j, j] <- (mad_j * 1.4826)^2
    }
    
    # Step 9: Reconstruct covariance matrix
    cov_matrix <- loadings %*% factor_cov %*% t(loadings) + D
    
    # Step 10: Ensure positive definiteness
    eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
    values <- pmax(eigen_decomp$values, 1e-8)
    cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
    
    return(cov_matrix)
  }, error = function(e) {
    warning("Hybrid MCD+RobPCA failed:", e$message)
    # Fall back to simpler hybrid approach
    return(hybrid_mcd_pca(returns, k, alpha))
  })
}

# Hybrid MCD + Robust PCA + GLasso covariance estimation
hybrid_mcd_robpca_glasso <- function(returns, k = NULL, alpha = NULL, lambda = NULL, threshold = 3.0) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  # Set default lambda if not provided
  if(is.null(lambda)) {
    lambda <- 0.1 * sqrt(log(p) / n)
    cat("Using default lambda =", lambda, "\n")
  }
  
  cat("Computing Hybrid MCD+RobPCA+GLasso\n")
  
  tryCatch({
    # Step 1: Get initial robust covariance estimate using hybrid approach
    initial_cov <- hybrid_mcd_robpca(returns, k, alpha, threshold)
    
    # Step 2: Apply graphical lasso regularization
    if(requireNamespace("glasso", quietly = TRUE)) {
      # Add small ridge to diagonal for stability before glasso
      diag(initial_cov) <- diag(initial_cov) * 1.01
      
      # Run glasso
      glasso_result <- glasso::glasso(initial_cov, rho = lambda, penalize.diagonal = FALSE)
      cov_matrix <- glasso_result$w
      
      # Check result validity
      if(any(is.na(cov_matrix)) || any(is.infinite(cov_matrix))) {
        warning("GLasso produced invalid results. Using initial estimate instead.")
        cov_matrix <- initial_cov
      } else {
        # Check condition number
        eig_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
        min_eig <- min(eig_vals)
        
        if(min_eig <= 0) {
          warning("GLasso result is not positive definite. Applying correction.")
          
          # Ensure positive definiteness
          eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
          values <- pmax(eigen_decomp$values, 1e-8)
          cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
        }
      }
      
      return(cov_matrix)
    } else {
      warning("glasso package not available, using initial estimate")
      return(initial_cov)
    }
  }, error = function(e) {
    warning("Hybrid MCD+RobPCA+GLasso failed:", e$message)
    # Fall back to MCD+RobPCA without GLasso
    return(hybrid_mcd_robpca(returns, k, alpha, threshold))
  })
}

# Hybrid Factor-Based Robust Estimator (HFBRE)
hybrid_factor_robust_cov <- function(returns, k = NULL, alpha = NULL, threshold = 3.0) {
  n <- nrow(returns)
  p <- ncol(returns)
  
  # Determine number of factors adaptively if not provided
  if(is.null(k)) {
    tryCatch({
      # Pre-process data to handle outliers
      outlier_result <- detect_outliers(returns, method = "MAD", threshold = threshold)
      winsorized_returns <- winsorize_returns(returns, outlier_result$outlier_mask, method = "shrinkage")
      
      # Use robust parallel analysis
      cov_mat <- cov(winsorized_returns, use = "pairwise.complete.obs")
      orig_eig <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
      
      # Generate random data with same dimensions
      n_simulations <- 10
      random_eigs <- matrix(0, n_simulations, p)
      
      for(i in 1:n_simulations) {
        random_data <- matrix(rnorm(n*p), nrow=n, ncol=p)
        random_cov <- cov(random_data)
        random_eigs[i,] <- eigen(random_cov, symmetric = TRUE, only.values = TRUE)$values
      }
      
      # Average the random eigenvalues
      mean_random_eigs <- colMeans(random_eigs)
      
      # Find where original eigenvalues > random eigenvalues
      k <- sum(orig_eig > mean_random_eigs)
      
      # Ensure k is reasonable
      k <- min(k, ceiling(p/4), floor(n/5))
      k <- max(k, 3)  # At least 3 factors
    }, error = function(e) {
      # Fallback option if factor selection fails
      warning("Factor number selection failed:", e$message)
      k <- min(ceiling(p/10), floor(n/5), 15)
      k <- max(k, 3)
    })
    
    cat("Selected optimal number of factors k =", k, "\n")
  }
  
  # Determine appropriate alpha based on factor dimensionality
  if(is.null(alpha)) {
    alpha <- min(0.875, max(0.5, 0.75 + 0.125 * (k/n)))
  }
  
  cat("Computing Hybrid Factor-Based Robust Estimator with k =", k, "and alpha =", alpha, "\n")
  
  # Step 1: Robust factor identification
  tryCatch({
    # Use robust PCA to identify factors
    rob_pca <- NULL
    
    if(requireNamespace("pcaPP", quietly = TRUE)) {
      # Try to use PCAgrid
      rob_pca <- tryCatch({
        pca_result <- pcaPP::PCAgrid(returns, k = k, method = "qn", center = TRUE)
        list(scores = pca_result$scores, loadings = pca_result$loadings)
      }, error = function(e) {
        warning("PCAgrid failed:", e$message)
        NULL
      })
    }
    
    # Fall back to standard PCA on winsorized data if PCAgrid fails
    if(is.null(rob_pca)) {
      # Detect and handle outliers
      outlier_result <- detect_outliers(returns, method = "MAD", threshold = threshold)
      winsorized_returns <- winsorize_returns(returns, outlier_result$outlier_mask, method = "shrinkage")
      
      # Apply PCA on winsorized data
      pca_result <- prcomp(winsorized_returns, center = TRUE, scale. = FALSE)
      rob_pca <- list(
        scores = pca_result$x[, 1:k, drop = FALSE],
        loadings = pca_result$rotation[, 1:k, drop = FALSE]
      )
    }
    
    # Extract factors and loadings
    factors <- rob_pca$scores
    loadings <- rob_pca$loadings
    
    # Verify that factors are valid
    if(is.null(factors) || !is.matrix(factors) || nrow(factors) < 2 || ncol(factors) < 1) {
      stop("Invalid factors matrix")
    }
    
    # Step 2: Robust estimation of factor covariance
    factor_cov <- tryCatch({
      if(n > 2*k) {  # Check if we have enough observations for MCD
        mcd_result <- rrcov::CovMcd(factors, alpha = alpha)
        mcd_result@cov
      } else {
        # If not enough observations, use robust pairwise approach
        diag_var <- apply(factors, 2, function(x) {
          mad_val <- mad(x, na.rm = TRUE)
          if(mad_val <= 0) mad_val <- sd(x, na.rm = TRUE) / 1.4826
          if(mad_val <= 0) mad_val <- 1e-8
          return((mad_val * 1.4826)^2)
        })
        
        # Use Spearman correlations for robustness
        cors <- cor(factors, method = "spearman")
        diag(sqrt(diag_var)) %*% cors %*% diag(sqrt(diag_var))
      }
    }, error = function(e) {
      warning("Robust factor covariance failed:", e$message)
      # Use standard covariance as fallback
      cov(factors)
    })
    
    # Step 3: Calculate residuals
    centered_returns <- scale(returns, center = TRUE, scale = FALSE)
    reconstructed <- factors %*% t(loadings)
    residuals <- centered_returns - reconstructed
    
    # Step 4: Detect outliers in residuals with fixed threshold
    residual_mask <- matrix(1, nrow = n, ncol = p)
    
    for(j in 1:p) {
      res_j <- residuals[, j]
      med_j <- median(res_j, na.rm = TRUE)
      mad_j <- mad(res_j, na.rm = TRUE)
      if(mad_j <= 0) mad_j <- sd(res_j, na.rm = TRUE) / 1.4826
      if(mad_j <= 0) mad_j <- 1e-8
      
      # Use fixed threshold
      outliers_j <- abs(res_j - med_j) > threshold * mad_j
      residual_mask[outliers_j, j] <- 0
    }
    
    # Step 5: Apply winsorization to residuals
    winsorized_residuals <- winsorize_returns(residuals, residual_mask, method = "shrinkage")
    
    # Step 6: Compute robust residual variances
    residual_var <- diag(p)
    for(j in 1:p) {
      res_j <- winsorized_residuals[, j]
      
      # Robust variance estimate
      mad_j <- median(abs(res_j - median(res_j)), na.rm = TRUE)
      if(mad_j <= 0) mad_j <- sd(res_j, na.rm = TRUE) / 1.4826
      if(mad_j <= 0) mad_j <- 1e-8
      
      residual_var[j, j] <- (mad_j * 1.4826)^2
    }
    
    # Step 7: Reconstruct covariance matrix
    cov_matrix <- loadings %*% factor_cov %*% t(loadings) + residual_var
    
    # Step 8: Apply regularization to ensure positive definiteness
    eigen_decomp <- eigen(cov_matrix, symmetric = TRUE)
    values <- pmax(eigen_decomp$values, 1e-8)
    cov_matrix <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
    
    return(cov_matrix)
  }, error = function(e) {
    warning("Hybrid factor robust covariance estimation failed:", e$message, ". Using sample covariance.")
    
    # Fall back to MCD if sufficient data
    if(n > 2*p) {
      return(mcd_cov(returns, alpha))
    } else {
      # Diagonal robust covariance as last resort
      diag_var <- apply(returns, 2, function(x) {
        mad_val <- mad(x, na.rm = TRUE)
        if(mad_val <= 0) mad_val <- sd(x, na.rm = TRUE) / 1.4826
        if(mad_val <= 0) mad_val <- 1e-8
        return((mad_val * 1.4826)^2)
      })
      return(diag(diag_var))
    }
  })
}

# Tyler M-estimator implementation 
tyler_m_cov <- function(returns, tol = 1e-6, max_iter = 100) {
  p <- ncol(returns)
  n <- nrow(returns)
  
  tryCatch({
    V <- diag(p)  # Initial estimate
    
    for (iter in 1:max_iter) {
      V_old <- V
      
      # Calculate distances with stabilization
      d <- sapply(1:n, function(i) {
        x <- returns[i,]
        # Ensure numeric vector
        if(is.null(dim(x))) {
          x <- as.numeric(x)
        } else {
          x <- as.numeric(x[1,])
        }
        
        # Compute Mahalanobis distance with safeguards
        tryCatch({
          V_inv <- solve(V)
          distance <- sqrt(max(1e-10, t(x) %*% V_inv %*% x))
          return(distance)
        }, error = function(e) {
          return(sqrt(sum(x^2)))  # Fallback to Euclidean distance
        })
      })
      
      # Update V with stabilization
      V_new <- matrix(0, p, p)
      for(i in 1:n) {
        x <- as.numeric(returns[i,])
        weight <- p / (max(d[i], 1e-6) * n)
        V_new <- V_new + weight * (x %*% t(x))
      }
      
      # Normalize to have determinant 1
      det_V <- det(V_new)
      if(det_V > 0) {
        V_new <- V_new / det_V^(1/p)
      } else {
        # If determinant is zero or negative, use near-PD correction
        V_new <- as.matrix(Matrix::nearPD(V_new)$mat)
        det_V <- det(V_new)
        V_new <- V_new / det_V^(1/p)
      }
      
      # Check convergence
      rel_diff <- norm(V_new - V_old, "F") / max(1e-10, norm(V_old, "F"))
      V <- V_new
      
      if(rel_diff < tol) {
        cat("Tyler M-estimator converged after", iter, "iterations\n")
        break
      }
    }
    
    # Ensure result is positive definite
    return(as.matrix(Matrix::nearPD(V)$mat))
  }, error = function(e) {
    warning("Tyler M-estimator failed: ", e$message, ". Using sample covariance.")
    return(cov(returns))
  })
}

# Adaptive Ledoit-Wolf shrinkage estimator
adaptive_lw_cov <- function(returns) {
  p <- ncol(returns)
  n <- nrow(returns)
  
  tryCatch({
    # Start with robust initial estimate
    robust_cov <- tryCatch({
      mcd_cov(returns)
    }, error = function(e) {
      # Fall back to robust diagonal covariance
      diag_var <- apply(returns, 2, function(x) {
        mad_val <- mad(x, na.rm = TRUE)
        if(mad_val <= 0) mad_val <- sd(x, na.rm = TRUE) / 1.4826
        if(mad_val <= 0) mad_val <- 1e-8
        return((mad_val * 1.4826)^2)
      })
      diag(diag_var)
    })
    
    # Calculate optimal shrinkage target and intensity
    if(requireNamespace("corpcor", quietly = TRUE)) {
      # Use corpcor package for linear shrinkage
      shrink_result <- corpcor::cov.shrink(returns, verbose = FALSE)
      return(shrink_result)
    } else {
      # Manual implementation of adaptive shrinkage
      # Compute sample covariance
      S <- cov(returns)
      
      # Use robust diagonal as target
      diag_var <- diag(robust_cov)
      F <- diag(diag_var)
      
      # Calculate adaptive shrinkage intensity
      gamma <- min(1, max(0, p/n))
      # Blend between sample and target
      shrinkage_cov <- (1-gamma) * S + gamma * F
      
      # Ensure positive definiteness
      eigen_vals <- eigen(shrinkage_cov, symmetric = TRUE, only.values = TRUE)$values
      if(min(eigen_vals) <= 0) {
        eigen_decomp <- eigen(shrinkage_cov, symmetric = TRUE)
        values <- pmax(eigen_decomp$values, 1e-8)
        shrinkage_cov <- eigen_decomp$vectors %*% diag(values) %*% t(eigen_decomp$vectors)
      }
      
      return(shrinkage_cov)
    }
  }, error = function(e) {
    warning("Adaptive LW shrinkage failed: ", e$message, ". Using sample covariance.")
    return(cov(returns))
  })
}
