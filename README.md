# Robust Portfolio Optimization

This repository contains the code used in the thesis "High-Dimensional Robust Portfolio Optimization Under Sparse Contamination: A Factor-Analytic Approach." The project implements and analyzes various robust covariance estimation methods for portfolio optimization under different market conditions, including stressed environments.

## Overview

Modern portfolio theory relies heavily on accurate covariance estimation to build optimal portfolios. However, traditional methods can be sensitive to outliers and violate distributional assumptions, particularly during market stress. This project implements and tests several robust covariance estimation techniques, with a focus on:

- Minimum Covariance Determinant (MCD)
- Hybrid MCD with Principal Component Analysis (PCA)
- Hybrid MCD with Robust PCA
- Hybrid MCD with Robust PCA and Graphical Lasso regularization
- Hybrid Factor-Based Robust Estimator (HFBRE)
- Tyler M-estimator
- Adaptive Ledoit-Wolf shrinkage

The code includes a comprehensive simulation study, analysis of real S&P 500 data, stress testing, and statistical validation of the results.

## Repository Structure

- `src/`: Source code implementing all methods and analysis techniques
  - `simulation/`: Code for simulation studies
  - `analysis/`: Code for real data analysis
  - `stress_testing/`: Code for stress testing market conditions
  - `statistics/`: Statistical testing and validation code
- `examples/`: Example scripts showing how to use the codebase
- `data/`: Data loading functions and sample data
- `docs/`: Additional documentation
- `results/`: Directory for storing analysis results (not tracked by git)

## Key Features

- **Simulation Study**: Generate synthetic financial data with controlled contamination and evaluate the performance of different covariance estimation methods
- **S&P 500 Analysis**: Apply robust portfolio optimization to real S&P 500 data
- **Market Regime Analysis**: Performance evaluation under different market regimes
- **Stress Testing**: Test portfolio robustness under simulated market stress conditions
- **Statistical Analysis**: Comprehensive statistical testing framework to validate results

## Getting Started

### Prerequisites

The code requires R with the following packages:
```r
required_packages <- c(
  "mvtnorm", "robustbase", "rrcov", "pcaPP", "DEoptim", "glasso", "corpcor",
  "GA", "Matrix", "parallel", "doParallel", "foreach", "ggplot2", "reshape2",
  "dplyr", "plotly", "shiny", "shinydashboard", "tidyquant", "quantmod", "xts",
  "PerformanceAnalytics", "tidyr", "quadprog", "tseries", "boot", "car",
  "moments", "nortest"
)

# Install required packages
install.packages(required_packages)
```

### Running the Code

#### Simulation Study

To run the simulation study:

```r
source("examples/run_simulation.R")
```

#### S&P 500 Analysis

To analyze S&P 500 data:

```r
source("examples/run_sp500_analysis.R")
```

#### Stress Testing

To run stress tests:

```r
source("examples/run_stress_test.R")
```

## Results

The code produces various outputs including:
- Performance metrics for different portfolio methods
- Visualizations comparing methods across different conditions
- Statistical test results
- Regime-specific performance summaries

Results are saved in the `results/` directory by default.

## Citation

If you use this code in your research, please cite:

```
@thesis{author2023robust,
  title={Robust Portfolio Optimization with Applications to Financial Data},
  author={Author Name},
  year={2023},
  school={University Name}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The implementation of robust covariance estimation methods builds upon the excellent work in the robustbase and rrcov packages
- Thanks to all contributors to the statistical and financial packages used in this project
