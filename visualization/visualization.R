# ========================================================================
# Visualization Functions for Robust Portfolio Optimization
# ========================================================================

# Create performance visualizations
create_performance_visualizations <- function(results) {
  plots <- list()
  
  # Extract performance summary
  perf <- results$performance_summary
  
  # Order methods
  perf$Method <- factor(perf$Method, levels = perf$Method[order(perf$OutOfSample_Sharpe, decreasing = TRUE)])
  
  # Identify method types for better visualization
  perf$OptimizationType <- ifelse(grepl("MaxSharpe", perf$Method), "Max Sharpe", 
                                  ifelse(grepl("MinVar", perf$Method), "Min Variance", "Equal Weight"))
  
  perf$CovMethod <- sapply(strsplit(as.character(perf$Method), "_"), function(x) x[1])
  
  # Sharpe ratio comparison
  plots$sharpe_ratio <- ggplot(perf, aes(x = Method, y = OutOfSample_Sharpe, fill = OptimizationType)) +
    geom_bar(stat = "identity") +
    labs(title = "Out-of-Sample Sharpe Ratio by Method",
         x = "", y = "Sharpe Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Return vs Volatility scatter plot
  plots$risk_return <- ggplot(perf, aes(x = OutOfSample_Volatility, y = OutOfSample_Return, 
                                        color = CovMethod, shape = OptimizationType)) +
    geom_point(size = 3) +
    geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
    labs(title = "Risk-Return Profile",
         x = "Annualized Volatility", y = "Annualized Return") +
    theme_minimal()
  
  # In-sample vs Out-of-sample comparison
  perf_long <- reshape2::melt(perf, id.vars = c("Method", "OptimizationType", "CovMethod"), 
                              measure.vars = c("InSample_Sharpe", "OutOfSample_Sharpe"),
                              variable.name = "Period", value.name = "Sharpe")
  
  plots$in_vs_out <- ggplot(perf_long, aes(x = Method, y = Sharpe, fill = Period)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "In-Sample vs Out-of-Sample Sharpe Ratio",
         x = "", y = "Sharpe Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Method comparison by covariance type
  cov_comparison <- ggplot(perf, aes(x = CovMethod, y = OutOfSample_Sharpe, fill = OptimizationType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Sharpe Ratio by Covariance Method and Optimization Type",
         x = "Covariance Method", y = "Out-of-Sample Sharpe Ratio") +
    theme_minimal()
  
  plots$cov_comparison <- cov_comparison
  
  # Performance heatmap
  heatmap_data <- perf %>%
    select(Method, CovMethod, OptimizationType, OutOfSample_Sharpe) %>%
    arrange(CovMethod, OptimizationType)
  
  performance_heatmap <- ggplot(heatmap_data, 
                                aes(x = OptimizationType, y = CovMethod, fill = OutOfSample_Sharpe)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = mean(heatmap_data$OutOfSample_Sharpe)) +
    labs(title = "Sharpe Ratio Heatmap",
         x = "Optimization Type", y = "Covariance Method", fill = "Sharpe Ratio") +
    theme_minimal()
  
  plots$performance_heatmap <- performance_heatmap
  
  # Maximum drawdown comparison
  plots$max_drawdown <- ggplot(perf, aes(x = Method, y = MaxDrawdown, fill = OptimizationType)) +
    geom_bar(stat = "identity") +
    labs(title = "Maximum Drawdown by Method",
         x = "", y = "Maximum Drawdown") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Condition number comparison (log scale)
  if("Condition_Number" %in% colnames(perf) && !all(is.na(perf$Condition_Number))) {
    plots$condition_number <- ggplot(perf, aes(x = Method, y = log10(Condition_Number), fill = CovMethod)) +
      geom_bar(stat = "identity") +
      labs(title = "Log10 Condition Number by Method",
           x = "", y = "Log10 Condition Number") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # Add regime-specific plots if available
  if(!is.null(results$regime_plots)) {
    for(name in names(results$regime_plots)) {
      plots[[paste0("regime_", name)]] <- results$regime_plots[[name]]
    }
  }
  return(plots)
}

# Create weight allocation visualizations
create_weight_visualizations <- function(results) {
  plots <- list()
  
  # Extract weights from all portfolios
  all_weights <- lapply(results$portfolio_results, function(x) x$weights)
  
  # Create data frame for plotting
  weights_df <- data.frame(
    Asset = rep(results$data$symbols, length(all_weights)),
    Weight = unlist(all_weights),
    Method = rep(names(all_weights), each = length(results$data$symbols))
  )
  
  # Order methods by out-of-sample performance
  weights_df$Method <- factor(weights_df$Method, 
                              levels = results$performance_summary$Method[order(results$performance_summary$OutOfSample_Sharpe, decreasing = TRUE)])
  
  # Identify method types for better visualization
  weights_df$OptimizationType <- ifelse(grepl("MaxSharpe", weights_df$Method), "Max Sharpe", 
                                        ifelse(grepl("MinVar", weights_df$Method), "Min Variance", "Equal Weight"))
  
  weights_df$CovMethod <- sapply(strsplit(as.character(weights_df$Method), "_"), function(x) x[1])
  
  # Weight allocation barplot (for top methods)
  top_methods <- head(levels(weights_df$Method), min(4, length(levels(weights_df$Method))))
  top_weights_df <- weights_df[weights_df$Method %in% top_methods, ]
  
  # For each top method, show top 15 weights
  top_weights_plots <- list()
  
  for(method in top_methods) {
    method_weights <- all_weights[[method]]
    names(method_weights) <- results$data$symbols
    
    # Sort weights and get top 15
    sorted_idx <- order(method_weights, decreasing = TRUE)
    top_15_assets <- results$data$symbols[sorted_idx[1:min(15, length(sorted_idx))]]
    top_15_weights <- method_weights[sorted_idx[1:min(15, length(sorted_idx))]]
    
    top_15_df <- data.frame(
      Asset = factor(top_15_assets, levels = top_15_assets),
      Weight = top_15_weights
    )
    
    # Create plot
    p <- ggplot(top_15_df, aes(x = Asset, y = Weight)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = paste("Top 15 Weights for Method:", method),
           x = "Asset", y = "Weight") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim(0, max(top_15_weights) * 1.1)  # Add some headroom
    
    top_weights_plots[[method]] <- p
  }
  
  plots$top_weights <- top_weights_plots
  
  # Calculate weight concentration (top 10 weights)
  concentration <- lapply(all_weights, function(w) {
    sorted_w <- sort(w, decreasing = TRUE)
    sum(sorted_w[1:min(10, length(sorted_w))])
  })
  
  conc_df <- data.frame(
    Method = names(concentration),
    Top10_Weight = unlist(concentration)
  )
  
  conc_df$Method <- factor(conc_df$Method, 
                           levels = results$performance_summary$Method[order(results$performance_summary$OutOfSample_Sharpe, decreasing = TRUE)])
  
  # Add method type information
  conc_df$OptimizationType <- ifelse(grepl("MaxSharpe", conc_df$Method), "Max Sharpe", 
                                     ifelse(grepl("MinVar", conc_df$Method), "Min Variance", "Equal Weight"))
  
  conc_df$CovMethod <- sapply(strsplit(as.character(conc_df$Method), "_"), function(x) x[1])
  
  # Weight concentration plot
  plots$weight_concentration <- ggplot(conc_df, aes(x = Method, y = Top10_Weight, fill = OptimizationType)) +
    geom_bar(stat = "identity") +
    labs(title = "Weight Concentration (Top 10 Assets)",
         x = "", y = "Cumulative Weight") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Concentration by covariance method
  conc_by_cov <- ggplot(conc_df, aes(x = CovMethod, y = Top10_Weight, fill = OptimizationType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Top 10 Weight Concentration by Covariance Method",
         x = "Covariance Method", y = "Top 10 Weight Concentration") +
    theme_minimal()
  
  plots$conc_by_cov <- conc_by_cov
  
  # Calculate standard deviation of weights (measure of diversification)
  diversification <- lapply(all_weights, function(w) {
    sd(w)
  })
  
  div_df <- data.frame(
    Method = names(diversification),
    WeightSD = unlist(diversification)
  )
  
  div_df$Method <- factor(div_df$Method, 
                          levels = results$performance_summary$Method[order(results$performance_summary$OutOfSample_Sharpe, decreasing = TRUE)])
  
  # Add method type information
  div_df$OptimizationType <- ifelse(grepl("MaxSharpe", div_df$Method), "Max Sharpe", 
                                    ifelse(grepl("MinVar", div_df$Method), "Min Variance", "Equal Weight"))
  
  div_df$CovMethod <- sapply(strsplit(as.character(div_df$Method), "_"), function(x) x[1])
  
  # Diversification plot
  plots$diversification <- ggplot(div_df, aes(x = Method, y = WeightSD, fill = OptimizationType)) +
    geom_bar(stat = "identity") +
    labs(title = "Weight Dispersion (Standard Deviation)",
         x = "", y = "Weight SD") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plots)
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

# Helper function to capitalize first letter (for visualizations)
capitalize <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
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

# Create visualizations for rolling window analysis results
create_rolling_visualizations <- function(rolling_results, metric = "sharpe") {
  # Extract rolling data
  sharpe_analysis <- rolling_results$sharpe_analysis
  
  # Create line plot for rolling metric
  rolling_plot <- ggplot(sharpe_analysis$plot_data, 
                         aes(x = Window, y = get(metric), color = Method, group = Method)) +
    geom_line() +
    geom_point(size = 1) +
    labs(title = paste("Rolling", capitalize(metric)),
         x = "Window", y = capitalize(metric)) +
    theme_minimal()
  
  # Create rank stability plot
  rank_plot <- ggplot(sharpe_analysis$rank_plot_data, 
                      aes(x = Window, y = Rank, color = Method, group = Method)) +
    geom_line() +
    geom_point(size = 1) +
    scale_y_reverse() +
    labs(title = paste("Rolling Rank Based on", capitalize(metric)),
         x = "Window", y = "Rank (1 = best)") +
    theme_minimal()
  
  # Create heatmap
  heatmap_plot <- ggplot(sharpe_analysis$plot_data, 
                         aes(x = Window, y = Method, fill = get(metric))) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = median(sharpe_analysis$plot_data[[metric]], na.rm = TRUE)) +
    labs(title = paste(capitalize(metric), "Heatmap Across Rolling Windows"),
         x = "Window", y = "Method") +
    theme_minimal()
  
  # Create summary barplot
  summary_plot <- ggplot(sharpe_analysis$mean_performance, 
                         aes(x = reorder(Method, Mean), y = Mean, fill = Method)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Average", capitalize(metric), "Across All Windows"),
         x = "Method", y = paste("Average", capitalize(metric))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Return visualizations
  list(
    rolling_plot = rolling_plot,
    rank_plot = rank_plot,
    heatmap_plot = heatmap_plot,
    summary_plot = summary_plot
  )
}

# Create interactive dashboard
create_portfolio_dashboard <- function(results) {
  library(shiny)
  library(shinydashboard)
  
  # Extract performance summary
  performance_summary <- results$performance_summary
  if(is.null(performance_summary)) {
    stop("No performance summary available in results")
  }
  
  # Initialize method types
  performance_summary$OptimizationType <- ifelse(grepl("MaxSharpe", performance_summary$Method), "Max Sharpe", 
                                                 ifelse(grepl("MinVar", performance_summary$Method), "Min Variance", "Equal Weight"))
  
  performance_summary$CovMethod <- sapply(strsplit(performance_summary$Method, "_"), function(x) {
    if(length(x) > 0) x[1] else "Unknown"
  })
  
  # Extract unique method types for filters
  cov_methods <- sort(unique(performance_summary$CovMethod))
  opt_methods <- sort(unique(performance_summary$OptimizationType))
  
  # UI definition
  ui <- dashboardPage(
    dashboardHeader(title = "Portfolio Optimization Analysis"),
    
    dashboardSidebar(
      sidebarMenu(
        menuItem("Performance Summary", tabName = "performance", icon = icon("dashboard")),
        menuItem("Risk Analysis", tabName = "risk", icon = icon("chart-line")),
        menuItem("Weight Analysis", tabName = "weights", icon = icon("balance-scale")),
        menuItem("Regime Analysis", tabName = "regime", icon = icon("random"))
      ),
      
      # Filter controls
      selectInput("covMethods", "Covariance Methods:",
                  choices = c("All", cov_methods),
                  selected = "All"),
      
      selectInput("optMethods", "Optimization Methods:",
                  choices = c("All", opt_methods),
                  selected = "All"),
      
      # Download buttons
      downloadButton("downloadData", "Download Results")
    ),
    
    dashboardBody(
      tabItems(
        # Performance tab
        tabItem(tabName = "performance",
                fluidRow(
                  box(plotOutput("sharpePlot"), width = 6),
                  box(plotOutput("riskReturnPlot"), width = 6)
                ),
                fluidRow(
                  box(plotOutput("inVsOutPlot"), width = 6),
                  box(plotOutput("covComparisonPlot"), width = 6)
                ),
                fluidRow(
                  box(DT::dataTableOutput("perfTable"), width = 12)
                )
        ),
        
        # Risk tab
        tabItem(tabName = "risk",
                fluidRow(
                  box(plotOutput("maxDrawdownPlot"), width = 6),
                  box(plotOutput("conditionNumberPlot"), width = 6)
                ),
                fluidRow(
                  box(plotOutput("sortinoPlot"), width = 6),
                  box(plotOutput("betaPlot"), width = 6)
                )
        ),
        
        # Weights tab
        tabItem(tabName = "weights",
                fluidRow(
                  box(plotOutput("weightConcentrationPlot"), width = 6),
                  box(plotOutput("diversificationPlot"), width = 6)
                ),
                fluidRow(
                  box(selectInput("selectedMethod", "Select Method:",
                                  choices = performance_summary$Method), width = 3),
                  box(plotOutput("weightDistributionPlot"), width = 9)
                )
        ),
        
        # Regime tab
        tabItem(tabName = "regime",
                fluidRow(
                  box(
                    selectInput("selectedRegime", "Select Regime:",
                                choices = c("overall", setdiff(names(results$regime_performance), "overall"))),
                    width = 3
                  ),
                  box(plotOutput("regimeSharpeBar"), width = 9)
                ),
                fluidRow(
                  box(plotOutput("regimeHeatmap"), width = 12)
                )
        )
      )
    )
  )
  
  # Server logic
  server <- function(input, output) {
    # Filtered data reactive
    filtered_data <- reactive({
      data <- performance_summary
      
      if(input$covMethods != "All") {
        data <- data[data$CovMethod == input$covMethods, ]
      }
      
      if(input$optMethods != "All") {
        data <- data[data$OptimizationType == input$optMethods, ]
      }
      
      data
    })
    
    # Sharpe ratio plot
    output$sharpePlot <- renderPlot({
      data <- filtered_data()
      data$Method <- factor(data$Method, levels = data$Method[order(-data$OutOfSample_Sharpe)])
      
      ggplot(data, aes(x = Method, y = OutOfSample_Sharpe, fill = OptimizationType)) +
        geom_bar(stat = "identity") +
        labs(title = "Out-of-Sample Sharpe Ratio by Method",
             x = "", y = "Sharpe Ratio") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    # Risk-Return scatter plot
    output$riskReturnPlot <- renderPlot({
      data <- filtered_data()
      
      ggplot(data, aes(x = OutOfSample_Volatility, y = OutOfSample_Return, 
                       color = CovMethod, shape = OptimizationType)) +
        geom_point(size = 3) +
        geom_text(aes(label = Method), vjust = -1, size = 3, check_overlap = TRUE) +
        labs(title = "Risk-Return Profile",
             x = "Annualized Volatility", y = "Annualized Return") +
        theme_minimal()
    })
    
    # In-sample vs Out-of-sample comparison
    output$inVsOutPlot <- renderPlot({
      data <- filtered_data()
      data_long <- reshape2::melt(data, id.vars = c("Method", "OptimizationType", "CovMethod"), 
                                  measure.vars = c("InSample_Sharpe", "OutOfSample_Sharpe"),
                                  variable.name = "Period", value.name = "Sharpe")
      
      # Order by out-of-sample Sharpe
      data_long$Method <- factor(data_long$Method, 
                                 levels = data$Method[order(-data$OutOfSample_Sharpe)])
      
      ggplot(data_long, aes(x = Method, y = Sharpe, fill = Period)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "In-Sample vs Out-of-Sample Sharpe Ratio",
             x = "", y = "Sharpe Ratio") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    # Covariance method comparison
    output$covComparisonPlot <- renderPlot({
      data <- filtered_data()
      
      ggplot(data, aes(x = CovMethod, y = OutOfSample_Sharpe, fill = OptimizationType)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Sharpe Ratio by Covariance Method and Optimization Type",
             x = "Covariance Method", y = "Out-of-Sample Sharpe Ratio") +
        theme_minimal()
    })
    
    # Max drawdown plot
    output$maxDrawdownPlot <- renderPlot({
      data <- filtered_data()
      data$Method <- factor(data$Method, levels = data$Method[order(data$MaxDrawdown)])
      
      ggplot(data, aes(x = Method, y = MaxDrawdown, fill = OptimizationType)) +
        geom_bar(stat = "identity") +
        labs(title = "Maximum Drawdown by Method",
             x = "", y = "Maximum Drawdown") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    # Condition number plot
    output$conditionNumberPlot <- renderPlot({
      data <- filtered_data()
      
      if("Condition_Number" %in% colnames(data) && !all(is.na(data$Condition_Number))) {
        ggplot(data, aes(x = Method, y = log10(Condition_Number), fill = CovMethod)) +
          geom_bar(stat = "identity") +
          labs(title = "Log10 Condition Number by Method",
               x = "", y = "Log10 Condition Number") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        # Show message if no condition numbers available
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Condition number data not available") +
          theme_void()
      }
    })
    
    # Sortino ratio plot
    output$sortinoPlot <- renderPlot({
      data <- filtered_data()
      
      if("Sortino" %in% colnames(data) && !all(is.na(data$Sortino))) {
        data$Method <- factor(data$Method, levels = data$Method[order(-data$Sortino)])
        
        ggplot(data, aes(x = Method, y = Sortino, fill = OptimizationType)) +
          geom_bar(stat = "identity") +
          labs(title = "Sortino Ratio by Method",
               x = "", y = "Sortino Ratio") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        # Show message if no Sortino ratios available
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Sortino ratio data not available") +
          theme_void()
      }
    })
    
    # Beta plot
    output$betaPlot <- renderPlot({
      data <- filtered_data()
      
      if("Beta" %in% colnames(data) && !all(is.na(data$Beta))) {
        ggplot(data, aes(x = Method, y = Beta, fill = OptimizationType)) +
          geom_bar(stat = "identity") +
          geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
          labs(title = "Market Beta by Method",
               x = "", y = "Beta") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        # Show message if no Beta data available
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Beta data not available") +
          theme_void()
      }
    })
    
    # Weight distribution plot
    output$weightDistributionPlot <- renderPlot({
      if(is.null(input$selectedMethod) || !input$selectedMethod %in% names(results$portfolio_results)) {
        # Show message if no method selected
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Please select a method") +
          theme_void()
      } else {
        # Get weights for selected method
        weights <- results$portfolio_results[[input$selectedMethod]]$weights
        
        # Check if weights are available
        if(is.null(weights)) {
          ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = "Weight data not available for this method") +
            theme_void()
        } else {
          # Create weight data frame
          weight_df <- data.frame(
            Asset = results$data$symbols,
            Weight = weights
          )
          
          # Order by weight
          weight_df <- weight_df[order(-weight_df$Weight), ]
          
          # Keep top 20 weights
          if(nrow(weight_df) > 20) {
            top_weights <- weight_df[1:20, ]
            
            # Sum the rest
            other_weight <- sum(weight_df$Weight[21:nrow(weight_df)])
            
            # Add "Other" category
            weight_df <- rbind(top_weights, data.frame(Asset = "Other", Weight = other_weight))
          }
          
          # Set factor levels for ordering
          weight_df$Asset <- factor(weight_df$Asset, levels = weight_df$Asset)
          
          # Create bar plot
          ggplot(weight_df, aes(x = Asset, y = Weight)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            labs(title = paste("Weight Distribution -", input$selectedMethod),
                 x = "Asset", y = "Weight") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      }
    })
    
    # Weight concentration plot
    output$weightConcentrationPlot <- renderPlot({
      # Extract weights from filtered methods
      filtered_methods <- filtered_data()$Method
      
      all_weights <- lapply(results$portfolio_results[filtered_methods], function(x) x$weights)
      
      # Calculate concentration (top 10 weights)
      concentration <- lapply(all_weights, function(w) {
        if(is.null(w)) return(NA)
        sorted_w <- sort(w, decreasing = TRUE)
        sum(sorted_w[1:min(10, length(sorted_w))])
      })
      
      conc_df <- data.frame(
        Method = names(concentration),
        Top10_Weight = unlist(concentration)
      )
      
      # Handle NA values
      conc_df <- conc_df[!is.na(conc_df$Top10_Weight), ]
      
      if(nrow(conc_df) == 0) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Weight data not available for filtered methods") +
          theme_void()
      } else {
        # Order by concentration
        conc_df$Method <- factor(conc_df$Method, levels = conc_df$Method[order(-conc_df$Top10_Weight)])
        
        # Add method type information
        conc_df$OptimizationType <- ifelse(grepl("MaxSharpe", conc_df$Method), "Max Sharpe", 
                                           ifelse(grepl("MinVar", conc_df$Method), "Min Variance", "Equal Weight"))
        
        ggplot(conc_df, aes(x = Method, y = Top10_Weight, fill = OptimizationType)) +
          geom_bar(stat = "identity") +
          labs(title = "Weight Concentration (Top 10 Assets)",
               x = "", y = "Cumulative Weight") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    })
    
    # Diversification plot
    output$diversificationPlot <- renderPlot({
      # Extract weights from filtered methods
      filtered_methods <- filtered_data()$Method
      
      all_weights <- lapply(results$portfolio_results[filtered_methods], function(x) x$weights)
      
      # Calculate standard deviation of weights
      diversification <- lapply(all_weights, function(w) {
        if(is.null(w)) return(NA)
        sd(w)
      })
      
      div_df <- data.frame(
        Method = names(diversification),
        WeightSD = unlist(diversification)
      )
      
      # Handle NA values
      div_df <- div_df[!is.na(div_df$WeightSD), ]
      
      if(nrow(div_df) == 0) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Weight data not available for filtered methods") +
          theme_void()
      } else {
        # Order by weight SD
        div_df$Method <- factor(div_df$Method, levels = div_df$Method[order(-div_df$WeightSD)])
        
        # Add method type information
        div_df$OptimizationType <- ifelse(grepl("MaxSharpe", div_df$Method), "Max Sharpe", 
                                          ifelse(grepl("MinVar", div_df$Method), "Min Variance", "Equal Weight"))
        
        ggplot(div_df, aes(x = Method, y = WeightSD, fill = OptimizationType)) +
          geom_bar(stat = "identity") +
          labs(title = "Weight Dispersion (Standard Deviation)",
               x = "", y = "Weight SD") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    })
    
    # Regime analysis plots
    output$regimeSharpeBar <- renderPlot({
      if(is.null(input$selectedRegime) || 
         is.null(results$regime_performance) || 
         !input$selectedRegime %in% names(results$regime_performance)) {
        # Show message if no regime data available
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Regime data not available") +
          theme_void()
      } else {
        # Get regime data
        regime_df <- results$regime_performance[[input$selectedRegime]]
        
        # Filter methods
        if(input$covMethods != "All") {
          regime_df <- regime_df[grepl(paste0("^", input$covMethods), regime_df$Method), ]
        }
        
        if(input$optMethods != "All") {
          if(input$optMethods == "Max Sharpe") {
            regime_df <- regime_df[grepl("MaxSharpe", regime_df$Method), ]
          } else if(input$optMethods == "Min Variance") {
            regime_df <- regime_df[grepl("MinVar", regime_df$Method), ]
          } else if(input$optMethods == "Equal Weight") {
            regime_df <- regime_df[grepl("EqualWeights", regime_df$Method), ]
          }
        }
        
        if(nrow(regime_df) == 0) {
          ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = "No data available for selected filters") +
            theme_void()
        } else {
          # Set method factor levels for ordering
          regime_df$Method <- factor(regime_df$Method, 
                                     levels = regime_df$Method[order(-regime_df$Sharpe)])
          
          # Identify method types for coloring
          regime_df$OptimizationType <- ifelse(grepl("MaxSharpe", regime_df$Method), "Max Sharpe", 
                                               ifelse(grepl("MinVar", regime_df$Method), "Min Variance", "Equal Weight"))
          
          ggplot(regime_df, aes(x = Method, y = Sharpe, fill = OptimizationType)) +
            geom_bar(stat = "identity") +
            labs(title = paste("Sharpe Ratio by Method -", capitalize(input$selectedRegime), "Regime"),
                 x = "", y = "Sharpe Ratio") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      }
    })
    
    # Regime heatmap
    output$regimeHeatmap <- renderPlot({
      if(is.null(results$regime_performance)) {
        # Show message if no regime data available
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Regime data not available") +
          theme_void()
      } else {
        heatmap_data <- data.frame()
        
        for(regime in setdiff(names(results$regime_performance), "overall")) {
          regime_df <- results$regime_performance[[regime]]
          regime_df$Regime <- regime
          
          # Apply filters
          if(input$covMethods != "All") {
            regime_df <- regime_df[grepl(paste0("^", input$covMethods), regime_df$Method), ]
          }
          
          if(input$optMethods != "All") {
            if(input$optMethods == "Max Sharpe") {
              regime_df <- regime_df[grepl("MaxSharpe", regime_df$Method), ]
            } else if(input$optMethods == "Min Variance") {
              regime_df <- regime_df[grepl("MinVar", regime_df$Method), ]
            } else if(input$optMethods == "Equal Weight") {
              regime_df <- regime_df[grepl("EqualWeights", regime_df$Method), ]
            }
          }
          
          heatmap_data <- rbind(heatmap_data, regime_df[, c("Method", "Regime", "Sharpe")])
        }
        
        if(nrow(heatmap_data) == 0) {
          ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = "No data available for selected filters") +
            theme_void()
        } else {
          # Convert Method to factor to control ordering
          methods_order <- results$regime_performance[["overall"]]$Method
          methods_order <- methods_order[methods_order %in% unique(heatmap_data$Method)]
          heatmap_data$Method <- factor(heatmap_data$Method, levels = methods_order)
          
          ggplot(heatmap_data, aes(x = Regime, y = Method, fill = Sharpe)) +
            geom_tile() +
            scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                                 midpoint = median(heatmap_data$Sharpe, na.rm = TRUE)) +
            labs(title = "Sharpe Ratio Heatmap Across Market Regimes",
                 x = "Market Regime", y = "Method") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 8))
        }
      }
    })
    
    # Performance table
    output$perfTable <- DT::renderDataTable({
      filtered_data() %>%
        select(Method, OutOfSample_Sharpe, OutOfSample_Return, OutOfSample_Volatility, 
               MaxDrawdown, Sortino, Beta) %>%
        DT::datatable(options = list(pageLength = 10),
                      rownames = FALSE) %>%
        DT::formatRound(columns = c("OutOfSample_Sharpe", "OutOfSample_Return", 
                                    "OutOfSample_Volatility", "MaxDrawdown", 
                                    "Sortino", "Beta"), digits = 4)
    })
    
    # Download handler
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("portfolio_analysis_results_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(performance_summary, file, row.names = FALSE)
      }
    )
  }
  
  # Run Shiny app
  shinyApp(ui, server)
}
