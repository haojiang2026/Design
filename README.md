# Reproducibility Code for the Simulation Study

This repository contains the R code required to reproduce the simulation results presented in the manuscript.

## Repository Contents

| File | Description |
|---|---|
| `Scenarios.R` | Simulation scenarios and sensitivity analysis settings. |
| `Functions_for_Design.R` | Functions for the proposed design. |
| `Functions_for_Benchmark.R` | Functions for the benchmark design. |
| `Functions_for_MCPMod.R` | Functions for the MCP-Mod design. |
| `Functions_for_spike_and_slab_prior.R` | Functions for the proposed design integrating BOIN design with a spike-and-slab prior. |
| `Functions_for_uniform_prior.R` | Functions for the proposed design integrating BOIN design with a uniform prior. |

---

## Arguments

The `get.oc()` function in these files requires some or all of the following arguments, depending on the design:

| Argument | Description |
|---|---|
| `pI.true` | True immune response probability for each dose level. |
| `pT0.true` | True toxicity probability for each dose level without immune response. |
| `pT1.true` | True toxicity probability for each dose level with immune response. |
| `pE0.true` | True efficacy probability for each dose level without immune response. |
| `pE1.true` | True efficacy probability for each dose level with immune response. |
| `rho0` | Correlation between efficacy and toxicity without immune response. |
| `rho1` | Correlation between efficacy and toxicity with immune response. |
| `phi.pT` | Highest acceptable toxicity probability. |
| `phi.pE` | Lowest acceptable efficacy probability. |
| `cf.pT` | Posterior probability cutoff for determining the admissible set of toxicity. |
| `cf.pE` | Posterior probability cutoff for determining the admissible set of efficacy. |
| `ndose` | Number of dose levels in the trial. |
| `ntrial` | Number of simulated trials. |
| `ndraw` | Number of posterior samples drawn for the parameters of interest. |
| `nstage` | Number of stages in the trial. |
| `nsample` | Prespecified sample size for each stage. |
| `seed_number` | Random seed used to reproduce the simulation results. |
| `utable` | Utility table specifying the utility value for each possible outcome. |
| `target.pT` | Target toxicity probability for BOIN design (Required only when integrating BOIN design).|
| `ncohort` | Number of cohorts in the BOIN design (Required only when integrating BOIN design).|
| `cohortsize` | Number of patients in each cohort (Required only when integrating BOIN design).|
| `sigma.spike` | Standard deviation of the spike component (Required only for the spike-and-slab prior).|
| `sigma.slab` | Standard deviation of the slab component (Required only for the spike-and-slab prior).|
| `nchain` | Number of MCMC chains (Required only for MCMC). |
| `burn` | Number of MCMC iterations discarded during the burn-in period (Required only for MCMC). |
| `thin` | Number of MCMC iterations per retained sample (Required only for MCMC). |

---

## Example

The following example illustrates how to reproduce the simulation results for Scenario 1 using the proposed design implemented in `Functions_for_Design.R`.

```r
library(dplyr)
library(mvtnorm)
library(Iso)
library(BOIN)
library(DoseFinding)
library(rjags)

#scenario 1
pI.true =  c(0.70,0.70,0.70,0.70)
pT0.true = c(0.20,0.40,0.45,0.50)
pT1.true = c(0.22,0.45,0.50,0.55)
pE0.true = c(0.40,0.55,0.58,0.60)
pE1.true = c(0.50,0.60,0.65,0.68)

#set arguments
rho0=0
rho1=0
phi.pT=0.30
phi.pE=0.25
cf.pT=0.80
cf.pE=0.80
ndose=4
ntrial=10000
ndraw=1000
nstage=5
nsample=c(40,5,5,5,5)
seed_number=1234
utable=data.frame(yI=integer(8),yT=integer(8),yE=integer(8),omega=integer(8))
utable$yI=c(0,0,0,0,1,1,1,1)
utable$yT=c(0,0,1,1,0,0,1,1)
utable$yE=c(0,1,0,1,0,1,0,1)
utable$omega=c(0,80,0,35,5,100,0,45)

#source functions
source('Functions_for_Design.R')

#run get.oc()
get.oc(pI.true,pT0.true,pT1.true,pE0.true,pE1.true,rho0,rho1,
       phi.pT,phi.pE,cf.pT,cf.pE,ndose,ntrial,ndraw,nstage,nsample,
       utable,seed_number)






