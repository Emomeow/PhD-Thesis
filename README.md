# PhD-Thesis
Joint analysis for multivariate biomarkers with a change point anchored at interval-censored event time, with applications to Huntington's Disease

## Overview

This repository contains R code for a joint modeling framework that analyzes the relationship between longitudinal biomarkers and interval-censored event times. The study focuses on modeling how biomarker trajectories change at an unobserved event time, with applications to Huntington's Disease progression.

## Study Aims

The primary objectives of this research are to:

1. **Develop a joint model** that simultaneously analyzes:
   - Multivariate longitudinal biomarker trajectories
   - Interval-censored event times (diagnosis times)
   - A change point in biomarker trajectories occurring at the event time

2. **Estimate likelihood of flexible baseline hazard functions** using:
   - Constant hazard (mode 1)
   - Piecewise constant hazard (mode 2)
   - B-spline representation (mode 3)

3. **Apply the methodology** to Huntington's Disease data to understand how cognitive and motor biomarkers change before and after diagnosis

4. **Study causal mediation effect** of covariate on survival probability through biomarkers and effect of covariate on biomarkers through event time

## Methodology

The joint model consists of:

- **Longitudinal sub-model**: Mixed-effects model for biomarker trajectories with a change point at the event time
  - Before event: `M*(t) = β₀ + β₁t + β₂X + a + ε`
  - After event: `M*(t) = β₀ + β₁t + β₂X + γ(t-T) + a + ε`

- **Survival sub-model**: Proportional hazards model for the interval-censored event time
  - Hazard: `λ(t|X,a,M*) = λ₀(t)exp(θₓX + θₐa + θₘM*(t))`

- **Estimation**: Maximum likelihood using Gauss-Hermite quadrature for random effects integration and Fisher-scoring optimization

Causal mediation analysis:

- **g-formula** is used to approximately calculate the causal mediation effect, i.e. natural direct effect (NDE) and natural indirect effect (NIE) of covariates through mediators (biomarkers/event time)
  
- **Simulations with repeated trials** are performed to demonstrate the correctness of our method

- **Application on HD data** reveals the mediation effect of CAP score on survival probability through two biomarkers sydigtot and stroopwo

## R Programs

### Core Analysis Files

#### 1. `Initial parameter.R`
**Purpose**: Initialize all model parameters and set up the analysis framework.

**Key Functions**:
- Defines sample size (`n.id`), number of visits (`J`), and biomarker dimensions (`K`)
- Initializes fixed effects parameters (`beta.0`, `beta.1`, `beta.2`)
- Initializes survival model parameters (`theta.x`, `theta.a`, `theta.m`)
- Sets up variance-covariance matrices (`Sigma.a`, `Sigma.e`)
- Configures Gauss-Hermite quadrature and Gauss-Legendre quadrature nodes and weights
- Defines different hazard modes (constant, piecewise, B-spline)

**Parameters Set**:
- `mode`: 1 (constant), 2 (piecewise), or 3 (B-spline) hazard
- `K`: Number of biomarker domains (default: 2)
- `p`: Number of covariates (default: 1)
- `J`: Number of visits (default: 11)
- `n.id`: Sample size (default: 400)

#### 2. `Data generation.R`
**Purpose**: Generate simulated longitudinal and survival data under the joint model framework.

**Key Functions**:
- `data.gen.vec()`: Generalized version supporting multiple biomarkers and covariates (updated 11/05/2025 for time-dependent covariates)
  - Main data generation function with SIMEX capability
  - Generates visit times with random noise
  - Simulates covariates, random effects, and biomarkers
  - Generates event times from proportional hazards model
  - Implements change point in biomarker trajectories after event
- `Lambda.0(t, par)`: Baseline cumulative hazard function (Weibull)
  
- `inv.Lambda.0(lam, par)`: Inverse cumulative hazard for event time generation

**Output**: Data frame with longitudinal biomarker measurements and interval-censored event times

#### 3. `Numeric Diff.R`
**Purpose**: Compute likelihoods, score functions and Hessian matrices using numerical differentiation.

**Key Functions**:
- `likelihood.vec2()`: Calculate likelihood with constant baseline hazard, using vectorized computation to save memory and time
- `likelihood.piecewise()`: Calculate likelihood with piecewise constant baseline hazard, using vectorized computation to save memory and time
- `spline_cumulative_hazard()`: Calculate B-spline's value as cumulative hazard given `x`, `knots`, `alpha`, `boundary_knots` and `degree`
- `spline_hazard()`: Calculate derivatives of B-spline's value as hazard value 
- `spline_hazard_matrix()`: Calculate derivatives of B-spline's value as hazard value when `x` input is matrix form
- `likelihood.spline2()`: Calculate likelihood with B-spline cumulative hazard given `parameters`, `data` and `knots`, using vectorized computation to save memory and time
- `diff.likelihood.vec()`: Compute score, and Hessian for vectorized parameters

**Method**: Uses finite differences with delta = 1e-6 for numerical stability

#### 4. `NR.R`
**Purpose**: Implement Fisher-scoring algorithms for parameter estimation, create Bootstrap samples of data to estimate standard deviation of all parameters.

**Key Functions**:
- `NR()`: Fisher-scoring algorithm for piecewise constant hazard model
  - Iterative optimization using score and Hessian
  - Line search with backtracking for step size selection
  - Ensures positive hazard and variance parameter constraints
  - Relative Convergence tolerance: 1e-3, max iterations: 200

- `NR_spline()`: Fisher-scoring algorithm specifically for B-spline hazard model
  - Similar structure to `NR()` but adapted for spline parameters
  - Uses `likelihood.spline2()` for B-spline evaluation
  - Enforces positive variance constraints

- `bootstrap_sample()`: Creates bootstrap samples by resampling subjects with replacement

**Convergence Criteria**:
- Maximum parameter relative change < 1e-3
- Step size * max(abs(step)) < 1e-6 triggers convergence
- Maximum 200 iterations

#### 5. `Causal Mediation.R`
**Purpose**: Conduct causal mediation analysis using piecewise constant hazard 

**Key Functions**:
-`With_cov()`: Calculate causal mediation effect of covariate on survival function through biomarkers given `x1`, `x2`, `t`, `parameters` and `data`
-`With_cov_true()`: Calculate true causal mediation effect of covariate on survival function through biomarkers using true random effects
-`NE_M2()`: Calculate causal mediation effect of covariate on biomarkers through event time given `x1`, `x2`, `t`, `parameters` and `data`
-`NE_M2_true()`: Calculate true causal mediation effect of covariate on biomarkers through event time using true random effects

#### 6. `One trial spline.R`
**Purpose**: Execute a single simulation trial with B-spline hazard estimation.

**Key Functions**:
- Complete workflow from data generation to parameter estimation
- `estimate_bspline_hazard()`: Maximum likelihood estimation of B-spline coefficients
  - Uses interval-censored survival data
  - Optimizes B-spline coefficients (alpha) via L-BFGS-B
  - Computes standard errors from Hessian matrix

**Workflow**:
1. Generate data using `data.gen.vec()`
2. Set up B-spline knots based on event time quantiles
3. Get initial estimates for alpha coefficients
4. Run Newton-Raphson optimization
5. Compute standard errors (via bootstrap or Hessian)
6. Save results: coefficients, hazard estimates, coverage probabilities

**Output Files** (for each trial i):
- `est_coef_i={i}_K={K}_p={p}_n={n}_seed={seed}.csv`
- `hazard_i={i}_K={K}_p={p}_n={n}_seed={seed}.csv`
- `est_std_i={i}_K={K}_p={p}_n={n}_seed={seed}.csv`
- `CP_i={i}_K={K}_p={p}_n={n}_seed={seed}.csv`

**HCC parallel computation**:
- Since repeated trials are needed to demonstrate the performance of our joint model, we use HCC server to run repeated trials simutaneously to save the time
- **Simulation.submit** specifies the repeated times, memory-per-computer and number of cores for CPU

### Data Analysis Files

#### 7. `Processing data.R`
**Purpose**: Aggregate and analyze results from multiple simulation trials or bootstrap samples.

**Key Functions**:
- Reads multiple CSV files for coefficients, standard errors, coverage probabilities
- Computes summary statistics across trials:
  - Mean estimates
  - Bias (mean estimate - true value)
  - Standard deviation of estimates
  - Mean of standard errors
  - Coverage probability of 95% confidence intervals
  - 2.5% and 97.5% percentiles

- Generates plots for piecewise and B-spline hazard estimates:
  - Cumulative hazard curves with confidence bands
  - Comparison with true hazard function

- Summarizes causal mediation analysis results (NDE, NIE)

**Output**: 
- Combined results matrix with rows: Mean_coef, Mean_Bias, Std_coef, Mean_std, CP, 2.5%pct, 97.5%pct
- Saved as CSV files in Results directory
- Deletes intermediate CSV files after processing

#### 8. `HDanalysis.R`
**Purpose**: Analyze real Huntington's Disease data from the PREDICT-HD study.

**Key Sections**:

**Data Cleaning**:
- Filters for cases (diagnosed patients)
- Removes subjects with missing education years
- Excludes subjects diagnosed at first observation
- Imputes missing biomarker values (SDMT: `sydigtot`, Stroop: `stroopwo`)
- Removes subjects with completely missing biomarkers
- Handles missing diagnosis indicators (`dx17`)

**Variable Construction**:
- `status`: Event indicator (1 = diagnosed, 0 = censored)
- `visittime`: Time from baseline (age - age1)
- `preetime`: Last visit before diagnosis (V)
- `etime`: First visit after diagnosis (U)
- `befdiag`: Number of visits before/at diagnosis
- `afterdiag`: Number of visits after diagnosis
- `base_cap`: Baseline CAG-Age Product score (normalized)

**Initial Parameter Estimation**:
1. Fits separate linear mixed models for each biomarker using `lme()`
2. Estimates variance components (Sigma.e, Sigma.a)
3. Fits Cox proportional hazards model to estimate theta parameters
4. Fits joint model using `JM` package for refined initial values

**Main Analysis**:
- Runs B-spline hazard model (mode 3)
- Uses Newton-Raphson optimization via `NR_spline()`
- Computes standard errors using bootstrap (B=50 samples) or Hessian
- Generates hazard plots with confidence intervals
- Performs causal mediation analysis

**Biomarkers Analyzed**:
- `sydigtot`: Symbol Digit Modalities Test (cognitive function)
- `stroopwo`: Stroop Word Reading Test (cognitive processing)

**Covariates**:
- `age1`: Baseline age
- `educ_yrs`: Years of education
- `gender`: 0=female, 1=male
- `base_cap`: Baseline CAG-Age Product

**Output**:
- `cleaned data.csv`: Processed data
- `coefficients_spline_knots=6.csv`: Estimated parameters
- Confidence intervals and p-values
- Cumulative hazard plots

## Data Files

- `cleaned data.csv`: Processed Huntington's Disease data from PREDICT-HD study
  - Contains longitudinal biomarker measurements
  - Interval-censored diagnosis times
  - Subject-level covariates

## Dependencies

Required R packages:
```r
library(MASS)          # For multivariate normal distribution
library(dplyr)         # Data manipulation
library(tidyr)         # Data tidying
library(statmod)       # Gauss-Hermite quadrature
library(Matrix)        # Matrix operations
library(matrixcalc)    # Matrix calculations
library(expm)          # Matrix exponential
library(splines)       # Spline functions
library(splines2)      # Additional spline functions
library(ggplot2)       # Plotting
library(lme4)          # Linear mixed models (for HD analysis)
library(JMbayes2)      # Joint modeling (for HD analysis)
library(survival)      # Survival analysis
library(tibble)        # Modern data frames
```

## Usage

### Running a Simulation Study

1. **Set parameters** in `Initial parameter.R`:
```r
mode = 3              # 1=constant, 2=piecewise, 3=B-spline
n.id = 400           # Sample size
K = 2                 # Number of biomarkers
p = 1                 # Number of covariates
J = 11                # Number of visits
```

2. **Run a single trial**:
```r
source('One trial spline.R')
```

3. **Process multiple trials**:
```r
source('Processing data.R')
```

### Analyzing Huntington's Disease Data

```r
source('HDanalysis.R')
```

This will:
- Clean and process the PREDICT-HD data
- Fit the joint model with B-spline hazard
- Compute parameter estimates and standard errors
- Generate hazard plots and inference

## Model Modes

The code supports three baseline hazard specifications:

1. **Mode 1 (Constant)**: `λ₀(t) = λ₀` (constant hazard rate)
2. **Mode 2 (Piecewise)**: Hazard is constant within predefined intervals
3. **Mode 3 (B-spline)**: Flexible hazard using B-spline basis functions with specified degree and knots

## Output

Simulation results include:
- Parameter estimates
- Standard errors
- Coverage probabilities
- Bias
- Hazard function estimates

Real data analysis produces:
- Parameter estimates with confidence intervals
- P-values for hypothesis tests
- Cumulative hazard plots
- Model diagnostics

## References

This code implements methods for joint modeling of longitudinal and interval-censored survival data with a change point, extending classical joint modeling frameworks to handle:
- Interval-censored event times
- Change points in biomarker trajectories at unobserved event times
- Measurement error in biomarkers
- Flexible baseline hazard specifications

## Contact

For questions about this code, please contact the repository owner.

## License

[Specify license if applicable]
