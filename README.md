# Robust Portfolio Optimization

This repository contains the code and implementation for "High-Dimensional Robust Portfolio Optimization Under Sparse Contamination: A Factor-Analytic Approach." The project develops and analyzes innovative robust covariance estimation methods designed to maintain portfolio performance under market stress and contamination.

## Overview

Modern portfolio theory relies on accurate covariance estimation to construct optimal portfolios. However, traditional methods often fail during market stress due to their sensitivity to outliers and violation of distributional assumptions. This research addresses these challenges through specialized robust estimation techniques with a focus on computational feasibility for high-dimensional portfolios.

The repository implements several key methodological innovations:

- **Parallel Factor Space Estimator (PFSE)**: A novel approach combining parallel analysis for factor dimension selection with robust estimation in the reduced factor space
- **Sequential Screening Robust Estimator (SSRE)**: A multi-stage approach with preliminary outlier detection followed by robust factor analysis
- **Hybrid robust methods**: Various combinations of dimension reduction and robust estimation techniques

These are compared against traditional methods including:

- Minimum Covariance Determinant (MCD)
- Tyler's M-estimator
- Adaptive Ledoit-Wolf shrinkage
- Sample covariance matrix
- Equal-weight portfolios

## Repository Structure

```
robust-portfolio-optimization/
├── data/                  # Data acquisition and preprocessing
│   ├── sp500_data.R       # S&P 500 data acquisition functions
│   └── synthetic_data.R   # Synthetic data generation functions
├── src/                   # Source code for all estimation methods
│   ├── estimation/        # Covariance estimation methods
│   │   ├── pfse.R         # Parallel Factor Space Estimator
│   │   ├── ssre.R         # Sequential Screening Robust Estimator
│   │   ├── mcd.R          # Minimum Covariance Determinant
│   │   ├── tyler.R        # Tyler's M-estimator
│   │   └── adaptive_lw.R  # Adaptive Ledoit-Wolf
│   ├── portfolio/         # Portfolio optimization methods
│   │   ├── min_variance.R # Minimum variance optimization
│   │   ├── max_sharpe.R   # Maximum Sharpe ratio optimization
│   │   └── constraints.R  # Portfolio constraints
│   └── utils/             # Utility functions
│       ├── outlier.R      # Outlier detection functions
│       ├── factor.R       # Factor analysis utilities
│       └── parallel.R     # Parallel processing utilities
├── analysis/              # Analysis modules
│   ├── simulation/        # Simulation study 
│   │   ├── contamination.R # Contamination generation
│   │   └── sim_study.R    # Full simulation pipeline
│   ├── empirical/         # Empirical analysis
│   │   ├── market_regimes.R # Market regime identification
│   │   └── sp500_analysis.R # S&P 500 analysis
│   ├── stress_testing/    # Stress testing framework
│   │   ├── scenarios.R    # Stress scenario definitions
│   │   └── stress_test.R  # Stress test implementation
│   └── validation/        # Statistical validation
│       ├── bootstrap.R    # Bootstrap analysis
│       └── hypothesis.R   # Hypothesis testing
├── visualization/         # Visualization modules
│   ├── performance.R      # Performance visualization
│   ├── condition.R        # Condition number analysis
│   └── stability.R        # Weight stability visualization
├── examples/              # Example scripts
│   ├── run_simulation.R   # Run complete simulation study
│   ├── run_empirical.R    # Run S&P 500 analysis
│   └── run_stress_test.R  # Run stress testing
└── README.md              # This file
```

## Key Features

### 1. Novel Robust Estimation Methods

- **PFSE (Parallel Factor Space Estimator)**: Combines robust factor dimension selection via parallel analysis with targeted robustification in the factor space, achieving superior computational efficiency while maintaining statistical robustness
- **SSRE (Sequential Screening Robust Estimator)**: Implements a multi-stage approach with preliminary outlier detection, robust dimension reduction, and structured residual handling

### 2. Comprehensive Analysis Framework

- **Simulation Study**: Evaluates performance under controlled contamination with various distribution types
- **Empirical Analysis**: Tests methods on S&P 500 data across multiple market regimes (2015-2025)
- **Stress Testing**: Assesses resilience under systematically generated stress scenarios
- **Statistical Validation**: Provides rigorous statistical assessment of performance differences

### 3. Performance Metrics

- Risk-adjusted returns (Sharpe ratio, Sortino ratio)
- Maximum drawdown and recovery time
- Value-at-Risk (VaR) and Conditional Value-at-Risk (CVaR)
- Portfolio turnover and weight stability
- Computational efficiency

## Getting Started

### Prerequisites

The code requires R with the following packages:

```r
required_packages <- c(
  "mvtnorm", "robustbase", "rrcov", "pcaPP", "DEoptim", "glasso", "corpcor",
  "GA", "Matrix", "parallel", "doParallel", "foreach", "ggplot2", "reshape2",
  "dplyr", "plotly", "tidyquant", "quantmod", "xts", "PerformanceAnalytics", 
  "tidyr", "quadprog", "tseries", "boot", "moments", "nortest"
)

# Install required packages
install.packages(required_packages)
```

### Running the Code

#### Simulation Study

To run the simulation study with default parameters:

```r
source("examples/run_simulation.R")
```

This will execute the full simulation pipeline, including:
- Data generation with varying contamination levels
- Estimation using all methods
- Performance evaluation
- Results visualization

#### Empirical Analysis

To analyze S&P 500 data:

```r
source("examples/run_empirical.R")
```

This runs the analysis on S&P 500 constituents, including:
- Data acquisition and preprocessing
- Market regime identification
- Portfolio optimization with all methods
- Performance evaluation across regimes
- Statistical validation

#### Stress Testing

To run the stress testing framework:

```r
source("examples/run_stress_test.R")
```

This executes the multi-scenario stress testing framework to evaluate method resilience.

## Key Results

The repository implements the methods that achieved the following results in our research:

1. **Performance Improvement**: PFSE achieves approximately 16% higher Sharpe ratios than conventional approaches while requiring only 10% of the computational resources of traditional robust estimators

2. **Weight Stability**: Robust methods reduce portfolio turnover by 30-40% compared to conventional approaches, substantially lowering transaction costs

3. **Stress Resilience**: Under stress conditions, PFSE maintains over 95% of its performance while conventional methods deteriorate by 30-50%

4. **Computational Efficiency**: Hybrid methods achieve robust estimation with computation times reduced by 80-90% compared to traditional robust approaches

## Usage Examples

### Basic Usage

```r
# Load libraries and source files
source("src/estimation/pfse.R")
source("src/portfolio/min_variance.R")

# Generate or load data
returns <- read.csv("data/returns.csv")

# Apply PFSE estimation
pfse_result <- pfse(returns, k = 5, threshold = 3.0)

# Construct minimum variance portfolio
weights <- min_variance_portfolio(pfse_result$covariance, max_weight = 0.2)

# Evaluate performance
performance <- evaluate_portfolio(returns, weights)
print(performance)
```

### Advanced Usage: Custom Contamination

```r
# Generate clean data
n <- 500  # observations
p <- 100  # assets
clean_returns <- generate_returns(n, p, distribution = "t", df = 5)

# Apply contamination
contaminated_returns <- apply_contamination(
  clean_returns, 
  type = "tail", 
  proportion = 0.1, 
  intensity = 0.3
)

# Compare methods
methods <- c("Sample", "MCD", "Tyler", "PFSE", "SSRE")
results <- compare_methods(contaminated_returns, methods, 
                          portfolio_type = "MinVariance")

# Visualize results
plot_performance(results, metric = "sharpe_ratio")
```

## Contributing

We welcome contributions to this project. Please feel free to submit a pull request or open an issue to discuss potential improvements.

## License

This project is licensed under the MIT License - see the LICENSE file for details.



## Acknowledgments

- The implementation of robust covariance estimation methods builds upon work in the robustbase and rrcov packages
- Thanks to all contributors to the statistical and financial packages used in this project
