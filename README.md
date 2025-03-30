# Robust Portfolio Optimization

This repository contains the code and implementation for "High-Dimensional Robust Portfolio Optimization Under Sparse Contamination: A Factor-Analytic Approach." The project develops and analyzes innovative robust covariance estimation methods designed to maintain portfolio performance under market stress and contamination.

## Overview

Modern portfolio theory relies on accurate covariance estimation to construct optimal portfolios. However, traditional methods often fail during market stress due to their sensitivity to outliers and violation of distributional assumptions. This research addresses these challenges through specialized robust estimation techniques with a focus on computational feasibility for high-dimensional portfolios.

The repository implements several key methodological innovations:

- **Parallel Factor Space Estimator (PFSE)**: A novel approach combining robust factor dimension selection with targeted robustification in the factor space, achieving superior computational efficiency while maintaining statistical robustness
- **Sequential Screening Robust Estimator (SSRE)**: A multi-stage approach with preliminary outlier detection followed by robust factor analysis
- **SSRE with Graphical Lasso (SSRE_GLasso)**: Enhanced estimation with graphical lasso regularization for better conditioning

These are compared against traditional methods including:

- Minimum Covariance Determinant (MCD)
- Tyler's M-estimator
- Adaptive Ledoit-Wolf shrinkage
- Sample covariance matrix
- Equal-weight portfolios

## Repository Structure

```
robust-portfolio-optimization/
├── data/                 
│   └── data_utils.R       # Data acquisition and synthetic data generation
├── src/                   
│   ├── estimation.R       # All covariance estimation methods in one file
│   ├── portfolio.R        # All portfolio optimization methods in one file
│   └── utils.R            # Utility functions for outliers, factors, etc.
├── analysis/              
│   ├── simulation.R       # Simulation study in a single file
│   ├── empirical.R        # S&P 500 and empirical analysis
│   ├── stress_testing.R   # Stress testing framework
│   └── validation.R       # Statistical validation methods
├── visualization/         
│   └── visualization.R    # All visualization functions in one file
├── examples/              
│   ├── run_simulation.R   # Run simulation examples
│   ├── run_empirical.R    # Run empirical analysis examples
│   └── run_stress_test.R  # Run stress test examples
└── README.md             
```

## Key Features

### 1. Novel Robust Estimation Methods

- **PFSE (Parallel Factor Space Estimator)**: Combines robust factor dimension selection via parallel analysis with targeted robustification in the factor space, achieving superior computational efficiency while maintaining statistical robustness
- **SSRE (Sequential Screening Robust Estimator)**: Implements a multi-stage approach with preliminary outlier detection, robust dimension reduction, and structured residual handling
- **SSRE_GLasso**: Enhances estimation with graphical lasso regularization for better conditioning

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

1. **Performance Improvement**: Robust methods (PFSE, SSRE) achieve approximately 15-20% higher Sharpe ratios than conventional approaches while requiring only 10-15% of the computational resources of traditional robust estimators

2. **Weight Stability**: Robust methods reduce portfolio turnover by 30-40% compared to conventional approaches, substantially lowering transaction costs

3. **Stress Resilience**: Under stress conditions, hybrid methods maintain over 90% of their performance while conventional methods deteriorate by 30-50%

4. **Computational Efficiency**: PFSE and SSRE achieve robust estimation with computation times reduced by 80-90% compared to traditional robust approaches

## Usage Examples

### Basic Usage

```r
# Load modules
source("src/estimation.R")
source("src/portfolio.R")
source("src/utils.R")

# Generate or load data
returns <- read.csv("your_returns_data.csv")

# Apply PFSE estimation
cov_matrix <- pfse(returns, k = 5, threshold = 3.0)

# Construct minimum variance portfolio
weights <- min_var_portfolio(returns, cov_matrix)

# Evaluate performance
risk_measures <- calculate_risk_measures(calculate_portfolio_returns(weights, returns))
print(risk_measures)
```

### Advanced Usage: Stress Testing

```r
# Load modules
source("src/estimation.R")
source("src/portfolio.R")
source("analysis/stress_testing.R")

# Load or generate data
returns <- read.csv("your_returns_data.csv")

# Apply stress testing
stress_results <- run_stress_test_analysis(
  returns, 
  stress_method = "tail_contamination", 
  stress_intensity = 0.3
)

# Compare method performance
print(stress_results$performance_summary)
```

## Contributing

We welcome contributions to this project. Please feel free to submit a pull request or open an issue to discuss potential improvements.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The implementation of robust covariance estimation methods builds upon work in the robustbase and rrcov packages
- Thanks to all contributors to the statistical and financial packages used in this project
