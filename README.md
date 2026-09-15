# A Bayesian adaptive randomized phase I/II design for immunotherapy trials

This repository contains the R code for reproducing the simulation results presented in the manuscript and for implementing the proposed design in practice.

## Repository Contents

| File | Description |
|---|---|
| `Scenarios.R` | Simulation scenarios and sensitivity analysis settings. |
| `Functions_for_Design.R` | Functions for the proposed design. |
| `Functions_for_Benchmark.R` | Functions for the benchmark design. |
| `Functions_for_MCPMod.R` | Functions for the MCP-Mod design. |
| `Functions_for_spike_and_slab_prior.R` | Functions for proposed design integrating BOIN design with spike-and-slab prior. |
| `Functions_for_uniform_prior.R` | Functions for proposed design integrating BOIN design with uniform prior. |
| `Functions_for_user.R` | Functions for implementing the proposed design in practice. |

---

## Arguments

The functions in these files requires some or all of the following arguments, depending on the design:

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
| `nchain` | Number of MCMC chains (Required only for the spike-and-slab prior). |
| `burn` | Number of MCMC iterations discarded during burn-in period (Required only for the spike-and-slab prior). |
| `thin` | Number of MCMC iterations per retained sample (Required only for the spike-and-slab prior). |

---

## Simulation example

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
#target.pT=0.30
#ncohort=12
#cohortsize=3
#sigma.spike=0.01
#sigma.slab=1.1
#nchain=1
#burn=5000
#thin=1


#source functions
source('Functions_for_Design.R')


#run get.oc()
get.oc(pI.true,pT0.true,pT1.true,pE0.true,pE1.true,rho0,rho1,
       phi.pT,phi.pE,cf.pT,cf.pE,ndose,ntrial,ndraw,nstage,nsample,
       utable,seed_number)


#result
"
$pI.true (immune response)
[1] 0.7 0.7 0.7 0.7

$pT.true (marginal toxicity)
[1] 0.214 0.435 0.485 0.535

$pE.true (marginal efficacy)
[1] 0.470 0.585 0.629 0.656

$U.true (utility)
[1] 40.6500 42.6050 43.9965 44.0550

$mud.per (MUD selection percentage)
[1] 0.7709 0.1100 0.0082 0.0003 0.1106

$pat.per (patients allocation percentage)
[1] 0.4290586 0.2212824 0.1774359 0.1722231

$pat.mean (patients allocation)
[1] 24.9725 12.8793 10.3273 10.0239

$pat (mean sample size)
[1] 58.203
"
```

## Implementation Example

The following example illustrates how to implement the proposed design in practice with the functions provided in `Functions_for_user.R`.

```r
library(dplyr)
library(mvtnorm)
library(Iso)
library(BOIN)
library(DoseFinding)
library(rjags)

#set arguments
phi.pT=0.30
phi.pE=0.25
cf.pT=0.80
cf.pE=0.80
ndose=4
nstage=5
ndraw=1000
nsample=c(40,5,5,5,5)
utable=data.frame(yI=integer(8),yT=integer(8),yE=integer(8),omega=integer(8))
utable$yI=c(0,0,0,0,1,1,1,1)
utable$yT=c(0,0,1,1,0,0,1,1)
utable$yE=c(0,1,0,1,0,1,0,1)
utable$omega=c(0,80,0,35,5,100,0,45)


#observed data in stage 1
set.seed(1234)
ob_data=data.frame(ID=1:40,stage=rep(1,40),
  dose=c(rep(1,10),rep(2,10),rep(3,10),rep(4,10)),
  yI=c(rbinom(5,1,0.30),rbinom(5,1,0.60),rbinom(5,1,0.60),rbinom(5,1,0.60)),
  yT=c(rbinom(5,1,0.09),rbinom(5,1,0.14),rbinom(5,1,0.20),rbinom(5,1,0.49)),
  yE=c(rbinom(5,1,0.33),rbinom(5,1,0.46),rbinom(5,1,0.59),rbinom(5,1,0.59)))

#source functions
source('C:/Users/14198/Desktop/TrialDesign/package/code/Functions_for_user.R')

#get immunity probabilities
immunity=monitor.I(ob_data,ndose,ndraw)

#get toxicity probabilities
toxicity=monitor.T(ob_data,ndose,ndraw,immunity)

#get efficacy probabilities
efficacy=monitor.E(ob_data,ndose,ndraw,immunity)

#get utility scores
utility=monitor.U(immunity,toxicity,efficacy,utable)

#get admissible dose set
adm.set=get.adm.set(ndose,phi.pT,phi.pE,cf.pT,cf.pE,toxicity,efficacy)

#get patients allocation in next stage
stage=2 #set the next stage; nsample[2]=5 patients
allocation=get.allocation(adm.set,stage,nsample,ndraw,utility)


#result

"
adm.set
[1] 1 1 1 0

allocation$pmud (allocation probabilities)
[1] 0.018 0.589 0.393 0.000

allocation$n (number of patients allocated to each dose level in the next stage)
[1] 0 2 3 0
"


```




