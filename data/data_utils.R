# ========================================================================
# Data Acquisition and Synthetic Data Generation
# ========================================================================

# Download S&P 500 constituents and price data
get_sp500_data <- function(start_date, end_date, n_stocks = 100) {
  cat("Getting S&P 500 constituents...\n")
  
  # Get S&P 500 constituents
  sp500_stocks <- tryCatch({
    tq_index("SP500") %>%
      arrange(desc(weight)) %>%  # Sort by market cap weight
      head(n_stocks + 50)  # Get a few extra to account for missing data
  }, error = function(e) {
    # Fallback to a list of major S&P 500 stocks
    cat("Error retrieving S&P 500 components:", e$message, "\n")
    cat("Using a fallback list of major S&P 500 stocks\n")
    
    top_stocks <- c("AAPL", "MSFT", "AMZN", "NVDA", "GOOGL", "META", "GOOG", "TSLA", 
                    "BRK.B", "UNH", "JPM", "V", "PG", "MA", "HD", "CVX", "MRK", "LLY", 
                    "AVGO", "PEP", "KO", "ABBV", "COST", "BAC", "PFE", "TMO", "MCD", 
                    "CSCO", "CRM", "ABT", "ADBE", "WMT", "TXN", "CMCSA", "ACN", "NKE",
                    "DHR", "LIN", "PM", "VZ", "NEE", "WFC", "BMY", "RTX", "INTC")
    
    data.frame(
      symbol = top_stocks,
      company = top_stocks,
      weight = rev(seq_along(top_stocks))
    )
  })
  
  # Get symbols
  symbols <- sp500_stocks$symbol
  
  cat("Downloading price data for", length(symbols), "stocks...\n")
  
  # Download stock prices
  stock_prices <- tq_get(symbols, 
                         from = start_date, 
                         to = end_date, 
                         get = "stock.prices")
  
  # Also get S&P 500 index data for comparison
  sp500_index <- tq_get("^GSPC", 
                        from = start_date, 
                        to = end_date, 
                        get = "stock.prices")
  
  # Calculate daily returns
  stock_returns <- stock_prices %>%
    group_by(symbol) %>%
    tq_transmute(select = adjusted, 
                 mutate_fun = periodReturn, 
                 period = "daily", 
                 col_rename = "return")
  
  # Calculate S&P 500 index returns
  sp500_returns <- sp500_index %>%
    tq_transmute(select = adjusted, 
                 mutate_fun = periodReturn, 
                 period = "daily", 
                 col_rename = "return")
  
  # Reshape to wide format (date × stocks)
  returns_matrix <- stock_returns %>%
    pivot_wider(id_cols = date, names_from = symbol, values_from = return) %>%
    as.data.frame()
  
  # Extract dates
  dates <- returns_matrix$date
  returns_matrix <- as.matrix(returns_matrix[, -1])  # Remove date column
  
  # Handle missing values and select top n_stocks
  returns_matrix <- preprocess_returns(returns_matrix, n_stocks)
  
  # Create xts object for returns
  returns_xts <- xts(returns_matrix, order.by = dates)
  
  # Align market returns with stock returns dates
  market_returns <- sp500_returns$return
  names(market_returns) <- sp500_returns$date
  market_returns <- market_returns[as.character(dates)]
  market_returns_xts <- xts(market_returns, order.by = dates)
  
  # Return the processed data
  list(
    returns = returns_matrix,
    returns_xts = returns_xts,
    dates = dates,
    symbols = colnames(returns_matrix),
    market_returns = market_returns,
    market_returns_xts = market_returns_xts,
    stock_prices = stock_prices,
    sp500_index = sp500_index
  )
}

# Preprocess returns: handle missing values and select stocks
preprocess_returns <- function(returns_matrix, n_stocks = 100) {
  # Check dimensions
  n <- nrow(returns_matrix)
  p <- ncol(returns_matrix)
  
  # Calculate NA statistics
  na_counts <- colSums(is.na(returns_matrix))
  na_prop <- na_counts / n
  
  # Order columns by proportion of NAs
  ordered_cols <- order(na_prop)
  
  # Select top n_stocks with least missing data
  if(length(ordered_cols) >= n_stocks) {
    returns_matrix <- returns_matrix[, ordered_cols[1:n_stocks], drop = FALSE]
    cat("Selected", ncol(returns_matrix), "stocks with least missing data\n")
  } else {
    warning("Could not find", n_stocks, "stocks with sufficient data, using", length(ordered_cols), "stocks")
    returns_matrix <- returns_matrix[, ordered_cols, drop = FALSE]
  }
  
  # For remaining NA values, impute with column median
  p <- ncol(returns_matrix)
  for(j in 1:p) {
    na_idx <- which(is.na(returns_matrix[, j]))
    if(length(na_idx) > 0) {
      returns_matrix[na_idx, j] <- median(returns_matrix[, j], na.rm = TRUE)
    }
  }
  
  # Replace Inf values with NA then with column max/min
  for(j in 1:p) {
    inf_idx <- which(is.infinite(returns_matrix[, j]))
    if(length(inf_idx) > 0) {
      # Replace +Inf with max and -Inf with min
      pos_inf <- returns_matrix[inf_idx, j] == Inf
      neg_inf <- returns_matrix[inf_idx, j] == -Inf
      
      if(any(pos_inf)) {
        returns_matrix[inf_idx[pos_inf], j] <- max(returns_matrix[, j][is.finite(returns_matrix[, j])])
      }
      
      if(any(neg_inf)) {
        returns_matrix[inf_idx[neg_inf], j] <- min(returns_matrix[, j][is.finite(returns_matrix[, j])])
      }
    }
  }
  
  return(returns_matrix)
}

# Generate multivariate normal data with contamination
generate_data <- function(n, p, delta = 0.05, balanced = FALSE, heavy_tailed = FALSE) {
  # Create mean vector with sector structure
  mu <- rep(0, p)
  for (g in 1:ceiling(p/5)) {
    start_idx <- (g-1)*5 + 1
    end_idx <- min(g*5, p)
    mu[start_idx:end_idx] <- 0.001 + (g-1)*0.0002  # Realistic small values
  }
  
  # Create covariance matrix with factor structure
  sigma <- matrix(0.0005, p, p)  # Base correlation
  diag(sigma) <- 0.001 + mu/10   # Asset-specific variance
  
  # Generate clean data
  if (heavy_tailed) {
    # t-distribution with 5 degrees of freedom
    clean_data <- rmvt(n, sigma = sigma, df = 5)
    # Scale and shift to match the desired mean and variance
    clean_data <- scale(clean_data) * sqrt(diag(sigma))
    clean_data <- sweep(clean_data, 2, mu, "+")
  } else {
    clean_data <- rmvnorm(n, mu, sigma)
  }
  
  # Create contamination mask (1 = clean, 0 = contaminated)
  contamination_mask <- matrix(1, nrow = n, ncol = p)
  
  # Apply contamination
  if (delta > 0) {
    contaminated_data <- clean_data
    # Contamination parameters
    mu_outlier <- -1 * mu - 0.01  # Opposite direction
    
    if (balanced) {
      # Balanced contamination
      for (j in 1:p) {
        outlier_count <- round(delta * n)
        outlier_indices <- sample(n, outlier_count)
        contaminated_data[outlier_indices, j] <- rnorm(
          outlier_count, 
          mu_outlier[j], 
          sqrt(diag(sigma)[j]) * 3
        )
        contamination_mask[outlier_indices, j] <- 0
      }
    } else {
      # Unbalanced contamination - more common in financial data
      sp <- (p*(p+1))/2  # Factor for increasing contamination
      for (j in 1:p) {
        outlier_count <- round((delta * n * p) * j / sp)
        outlier_indices <- sample(n, outlier_count)
        contaminated_data[outlier_indices, j] <- rnorm(
          outlier_count, 
          mu_outlier[j], 
          sqrt(diag(sigma)[j]) * 3
        )
        contamination_mask[outlier_indices, j] <- 0
      }
    }
    return(list(
      X = contaminated_data,
      X_clean = clean_data,
      mu = mu,
      Sigma = sigma,
      contamination_mask = contamination_mask,
      contamination_level = delta
    ))
  } else {
    # No contamination
    return(list(
      X = clean_data,
      X_clean = clean_data,
      mu = mu,
      Sigma = sigma,
      contamination_mask = contamination_mask,
      contamination_level = 0
    ))
  }
}
