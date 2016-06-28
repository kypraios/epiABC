# epiABC

This repository collects code to perform ABC inference of infectious disease epidemics, used in the paper "A tutorial introduction to Bayesian inference for stochastic epidemic models using Approximate Bayesian Computation" by Kypraios, Neal and Prangle.

The code is presented as an R package for solely for convenience. It is not our intention to develop or support this as a fully featured package.

## Organisation

The package contains functions to run various ABC analyses.
These functions fall into several groups, indicated by a prefix.
More details of some groups in given in sections below.

Several scripts are available to run analyses as in the paper.
These are contained within the package as demos.

## "coupled" prefix: homogeneous SIR epidemics

### Vanilla/SMC ABC

#### Demo

* `coupled-Toni_SMC` - Script for running SMC ABC and estimating the infection rate. Setting `epsil` (thresholds) to
a single number gives the Vanilla ABC algorithm.

#### Function files

* `coupled-ABC_reject.R` - Implements Vanilla ABC (sampling parameters from the prior).
* `coupled-ABC_importance.R` - Implements ABC using importance samples from previous runs.
* `coupled-sim_homo_SIR.R` - Simulates homogeneous epidemics.

### Coupled ABC

#### Demo

* `coupled-homo_SIR_run_F` - Script for running coupled ABC and estimating the infection rate.

#### Function files

* `coupled-threshold.R` - Sets thresholds.
* `coupled-homo_gamma_epi_F.R` - Simulates epidemics with Gamma infectious periods. (Constant infectious period k=0).
* `coupled_moment_calc.R` - Function to compute moments for parameters.

## "house" prefix: household SIR epidemics

Note that priors for lambda_G and lambda_L have been semi-hard coded to be Exp (1).

Seattle and Tecumseh datasets are included in the package.

### Partially coupled ABC

#### Demo

* `house-Main_coupled.R` - Main file for running partially coupled ABC for household epidemics
and computing estimates of the posterior mean and variance of the parameters.

#### Function files

* `house-thres.R` - Sets thresholds within a household.
* `house-sim.R` - Simulates household epidemics
* `house-coupled.R` - Implements partially coupled SIR epidemic for household epidemics.
* `house-process.R` - Extracts key information from dataset chosen.
* `house-couple_run.R` - Runs partially coupled ABC to obtain sample of desired size.
* `house-moment_exp.R` - Function to assist calculating moments

### Vanilla ABC

#### Demo

* `house-Main_vanilla.R` - Main file for running vanilla ABC for household epidemics
and computing estimates of the posterior mean and variance of the parameters.

#### Function files

* `house-thres.R` - Sets thresholds within a household.
* `house-sim.R` - Simulates household epidemics
* `house-vanilla.R` - Implements ABC SIR epidemic for household epidemics.
* `house-process.R` - Extracts key information from dataset chosen.
* `house-vanilla_run.R` - Runs vanilla ABC to obtain sample of desired size.

### SMC-ABC

#### Demo

* `house-Main_Toni.R` - Main file for running SMC-ABC for household epidemics
and computing estimates of the posterior mean and variance of the parameters.

#### Function files

* `house-thres.R` - Sets thresholds within a household.
* `house-sim.R` - Simulates household epidemics
* `house-vanilla.R` - Implements ABC SIR epidemic for household epidemics.
* `house-process.R` - Extracts key information from dataset chosen.
* `house-vanilla_run.R` - Runs vanilla ABC to obtain sample of desired size. (Initial run).
* `house-import_run.R` - Runs importance sampling ABC for subsequent runs to obtain samples of the the desired size.

### Partially coupled ABC

#### Demo

* `house-Main_coupled_SMC.R` - Main file for running partially coupled ABC for household epidemics with SMC on \lambda_L parameter.
Computes estimates of the posterior mean and variance of the parameters.

#### Function files

* `house-thres.R` - Sets thresholds within a household.
* `house-sim.R` - Simulates household epidemics
* `house-coupled.R` - Implements partially coupled SIR epidemic for household epidemics.
* `house-process.R` - Extracts key information from dataset chosen.
* `house-couple_run.R` - Runs partially coupled ABC to obtain sample of desired size.
(Initial run).
* `house-coupSMC_run.R` - Runs importance sampling partially coupled ABC for subsequent runs to obtain samples of the the desired size.
* `house-moment_exp.R` - Function to assist calculating moments
