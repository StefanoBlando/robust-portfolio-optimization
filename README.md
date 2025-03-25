# High-Dimensional Robust Portfolio Optimization

This repository contains the code and resources for the Master's thesis "High-Dimensional Robust Portfolio Optimization Under Sparse Contamination: A Factor-Analytic Approach" by Stefano Blando, completed at Università degli studi di Roma "Tor Vergata".

## Overview

This research introduces novel robust estimation methods for portfolio optimization in high-dimensional settings with contaminated financial data. Key contributions include:

- The Parallel Factor Space Estimator (PFSE) method
- The Sequential Screening Robust Estimator (SSRE) method
- Multi-scenario stress testing framework for portfolio resilience
- Interactive dashboard for method exploration and visualization

## Repository Structure

- `data/`: Financial datasets and processing scripts
- `src/`: Source code organized by function
  - `data_processing/`: Data cleaning and preparation
  - `estimation_methods/`: Robust covariance estimation implementations
  - `portfolio_construction/`: Portfolio optimization algorithms
  - `stress_testing/`: Stress scenario generators and analysis
  - `visualization/`: Plotting functions and utilities
- `notebooks/`: Jupyter notebooks with examples and tutorials
- `results/`: Generated figures and tables
- `shiny-app/`: Interactive dashboard application
- `docs/`: Additional documentation

## Installation

### Prerequisites
- R 4.1.0 or higher
- Required packages listed in `requirements.txt`

### Setup
```bash
# Clone the repository
git clone https://github.com/yourusername/robust-portfolio-optimization.git
cd robust-portfolio-optimization

# Install R dependencies
Rscript -e "install.packages(c('quadprog', 'robustbase', 'rrcov', 'tidyverse', 'shiny'))"


