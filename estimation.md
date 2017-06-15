# Estimation and Diagnostic
Tommaso Rigon  

## Description

In this document we describe how to estimate the model we described in Section 3 of our paper. We implemented our algorithm in R, making use of the `Rcpp` and `RcppArmadillo` R packages for the most computationally intensive steps. The [R code](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/core_functions.R) is made available as well as the [C++](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/core_functions.cpp) code. 

**Important**. If you encounter error compiling the C++ code on a Mac Os X, this could be due to recent updates of the software `R` (version 3.4.0). Please refer to the [official R documentation](https://cloud.r-project.org/bin/macosx/tools/) for more details on this issue.

The Gibbs sampler for our model is contained in the function called `fit_logit`, which we will use extensively later on. As a first step, we load in memory all the required libraries and we compile the C++ code. Moreover, we load in memory also the previously obtained [`dataset`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.md).

The functions `logit_ranef`, `logit_ranef_spline` and `logit_ranefDP` allow to estimate the submodels described in Section 4.1 of the paper. The function `logit_ranefDP_spline` estimates the full model. All these functions  wrapped by the global `fit_logit` function.


```r
library(dplyr)
library(BayesLogit)
library(splines)
library(Rcpp)
library(RcppArmadillo)

rm(list=ls())

# Load the clean dataset
load("dataset.RData")

# Compile the C++ functions
sourceCpp("core_functions.cpp")

# Load the R functions
source("core_functions.R")
```

## Prior specification and estimation

The prior hyperparameters for all our models, discussed in Section 3 and Section 4.1 for the submodels, is specified through a list. Fixing `P_Fix_const = 1e-2` is equivalent to set ${\bf B} = \text{diag}(100,\dots,100)$.


```r
prior = list(P_Fix_const=1e-2, 
             H=32,
             a_lambda=1e-3, b_lambda=1e-3,
             a_tau=1e-4, b_tau=1e-4, 
             a_alpha=.1, b_alpha=.1)
```

Then, we set the number of MCMC iterations equal to `R= 20000` with a burn-in equal to `burn_in=2000`. We also define two `formula`: one for the specification in which `age` enters in the predictor linearly and one in which `age` is modeled using Bayesian penalized splines.


```r
R       <- 20000
burn_in <- 2000

# Age enters in the predictor linearly
f   <- as.formula(target ~ age + child + area + religion + education)
# Age is absent, since it is modeled flexibly through Bayesian P-splines.
f_s <- as.formula(target ~ child + area + religion + education)
```

The estimation process **requires a non-negligible amount of time** to be completed. On standard laptop, this will need about 4-5 hours. We made available the results of the MCMC chain in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/workspaces) folder, which can be loaded in memory without running the following steps.

#### 1. Usage choice

The model and the corresponding submodels for the usage choice are estimated below. 


```r
set.seed(123) # We set a seed so that our results are fully reproducible.

# We define the target variable. 
dataset$target     <- factor(dataset$method!="1. No contraceptive method")

# Estimate the submodels
fit1_ranef         <- fit_logit(f,dataset$state,dataset$age,dataset,method="ranef",prior,R,burn_in)
fit1_ranef_s       <- fit_logit(f_s,dataset$state,dataset$age,dataset,method="ranef_s",prior,R,burn_in)
fit1_dp_ranef      <- fit_logit(f,dataset$state,dataset$age,dataset,method="dp_ranef",prior,R,burn_in)

# Estimate the full model
fit1_dp_ranef_s    <- fit_logit(f_s,dataset$state,dataset$age,dataset,method="dp_ranef_s",prior,R,burn_in)
```

#### 2. Reversible choice

We condition to a smaller dataset of observations and therefore we select only those women who are making use of contraception.


```r
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
dataset2            <- dataset[dataset$method != "1. No contraceptive method",]

# We define the new target variable. 
dataset2$target     <- factor(dataset2$method != "2. Sterilization") # table(dataset2$target, dataset2$method)

# Estimate the submodels
fit2_ranef         <- fit_logit(f,dataset2$state,dataset2$age,dataset2,method="ranef",prior,R,burn_in)
fit2_ranef_s       <- fit_logit(f_s,dataset2$state,dataset2$age,dataset2,method="ranef_s",prior,R,burn_in)
fit2_dp_ranef      <- fit_logit(f,dataset2$state,dataset2$age,dataset2,method="dp_ranef",prior,R,burn_in)

# Estimate the full model
fit2_dp_ranef_s    <- fit_logit(f_s,dataset2$state,dataset2$age,dataset2,method="dp_ranef_s",prior,R,burn_in)
```


#### 3. Method choice

Finally for the method choice model we perform similar steps as before.


```r
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
dataset3            <- dataset2[dataset2$method != "2. Sterilization",]

# We define the new target variable. 
dataset3$target     <- factor(dataset3$method != "3. Natural methods") # table(dataset3$target,dataset3$method)

# Estimate the submodels
fit3_ranef         <- fit_logit(f,dataset3$state,dataset3$age,dataset3,method="ranef",prior,R,burn_in)
fit3_ranef_s       <- fit_logit(f_s,dataset3$state,dataset3$age,dataset3,method="ranef_s",prior,R,burn_in)
fit3_dp_ranef      <- fit_logit(f,dataset3$state,dataset3$age,dataset3,method="dp_ranef",prior,R,burn_in)

# Estimate the full model
fit3_dp_ranef_s    <- fit_logit(f_s,dataset3$state,dataset3$age,dataset3,method="dp_ranef_s",prior,R,burn_in)
```

We end the estimation step by thinning the chains, cleaning the workspace and saving the results.


```r
thinning <- 5*(1:4000)

# DP
fit1_dp_ranef$beta_RF     <- fit1_dp_ranef$beta_RF[thinning,]
fit2_dp_ranef$beta_RF     <- fit2_dp_ranef$beta_RF[thinning,]
fit3_dp_ranef$beta_RF     <- fit3_dp_ranef$beta_RF[thinning,]
fit1_dp_ranef$beta_Fix    <- fit1_dp_ranef$beta_Fix[thinning,]
fit2_dp_ranef$beta_Fix    <- fit2_dp_ranef$beta_Fix[thinning,]
fit3_dp_ranef$beta_Fix    <- fit3_dp_ranef$beta_Fix[thinning,]
fit1_dp_ranef$beta_spline <- fit1_dp_ranef$beta_spline[thinning,]
fit2_dp_ranef$beta_spline <- fit2_dp_ranef$beta_spline[thinning,]
fit3_dp_ranef$beta_spline <- fit3_dp_ranef$beta_spline[thinning,]

# Gaussian 
fit1_ranef$beta_RF <- fit1_ranef$beta_RF[thinning,]
fit2_ranef$beta_RF <- fit2_ranef$beta_RF[thinning,]
fit3_ranef$beta_RF <- fit3_ranef$beta_RF[thinning,]
fit1_ranef$beta_Fix <- fit1_ranef$beta_Fix[thinning,]
fit2_ranef$beta_Fix <- fit2_ranef$beta_Fix[thinning,]
fit3_ranef$beta_Fix <- fit3_ranef$beta_Fix[thinning,]
fit1_ranef$beta_spline <- fit1_ranef$beta_spline[thinning,]
fit2_ranef$beta_spline <- fit2_ranef$beta_spline[thinning,]
fit3_ranef$beta_spline <- fit3_ranef$beta_spline[thinning,]

# Gaussian + Splines
fit1_ranef_s$beta_RF <- fit1_ranef_s$beta_RF[thinning,]
fit2_ranef_s$beta_RF <- fit2_ranef_s$beta_RF[thinning,]
fit3_ranef_s$beta_RF <- fit3_ranef_s$beta_RF[thinning,]
fit1_ranef_s$beta_Fix <- fit1_ranef_s$beta_Fix[thinning,]
fit2_ranef_s$beta_Fix <- fit2_ranef_s$beta_Fix[thinning,]
fit3_ranef_s$beta_Fix <- fit3_ranef_s$beta_Fix[thinning,]
fit1_ranef_s$beta_spline <- fit1_ranef_s$beta_spline[thinning,]
fit2_ranef_s$beta_spline <- fit2_ranef_s$beta_spline[thinning,]
fit3_ranef_s$beta_spline <- fit3_ranef_s$beta_spline[thinning,]

# DP + splines
fit1_dp_ranef_s$beta_RF <- fit1_dp_ranef_s$beta_RF[thinning,]
fit2_dp_ranef_s$beta_RF <- fit2_dp_ranef_s$beta_RF[thinning,]
fit3_dp_ranef_s$beta_RF <- fit3_dp_ranef_s$beta_RF[thinning,]
fit1_dp_ranef_s$beta_Fix <- fit1_dp_ranef_s$beta_Fix[thinning,]
fit2_dp_ranef_s$beta_Fix <- fit2_dp_ranef_s$beta_Fix[thinning,]
fit3_dp_ranef_s$beta_Fix <- fit3_dp_ranef_s$beta_Fix[thinning,]
fit1_dp_ranef_s$beta_spline <- fit1_dp_ranef_s$beta_spline[thinning,]
fit2_dp_ranef_s$beta_spline <- fit2_dp_ranef_s$beta_spline[thinning,]
fit3_dp_ranef_s$beta_spline <- fit3_dp_ranef_s$beta_spline[thinning,]

# Save relevants files

# Baseline
save(fit1_ranef,
     fit2_ranef,
     fit3_ranef,
     file="workspaces/ranef.RData")

# Splines
save(fit1_ranef_s,file="workspaces/ranef_s_part1.RData")
save(fit2_ranef_s,file="workspaces/ranef_s_part2.RData")
save(fit3_ranef_s,file="workspaces/ranef_s_part3.RData")

# DP
save(fit1_dp_ranef,
     fit2_dp_ranef,
     fit3_dp_ranef,
     file="workspaces/dp_ranef.RData")

# Full model
save(fit1_dp_ranef_s,
     fit2_dp_ranef_s,
     fit3_dp_ranef_s,
     file="workspaces/dp_ranef_s.RData")
```


## Convergence diagnostic

We load again everything in memory, on a clean workspace. Notice that the results are available also in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/workspaces) folder, if the computations are excessive. We uploaded the workspaces separately due to limitation of GitHub in file size (max 25Mb).


```r
rm(list=ls())

# Load the clean dataset
load("dataset.RData")

# Load the results of the MCMC chain
load("workspaces/ranef.RData")
load("workspaces/ranef_s_part1.RData")
load("workspaces/ranef_s_part2.RData")
load("workspaces/ranef_s_part3.RData")
load("workspaces/dp_ranef.RData")
load("workspaces/dp_ranef_s.RData")
```


```r
library(coda)

beta_RF1  <- as.mcmc(as.matrix(fit1_dp_ranef_s$beta_RF))
beta_RF2  <- as.mcmc(as.matrix(fit2_dp_ranef_s$beta_RF))
beta_RF3  <- as.mcmc(as.matrix(fit3_dp_ranef_s$beta_RF))

beta_spline1 <- as.mcmc(as.matrix(fit1_dp_ranef_s$beta_spline))
beta_spline2 <- as.mcmc(as.matrix(fit2_dp_ranef_s$beta_spline))
beta_spline3 <- as.mcmc(as.matrix(fit3_dp_ranef_s$beta_spline))

beta_Fix1  <- as.mcmc(as.matrix(fit1_dp_ranef_s$beta_Fix))
beta_Fix2  <- as.mcmc(as.matrix(fit2_dp_ranef_s$beta_Fix))
beta_Fix3  <- as.mcmc(as.matrix(fit3_dp_ranef_s$beta_Fix))
```

#### Effective sample sizes

The effective sample size is monitored using the function `effectiveSize` of the [`coda`](https://cran.r-project.org/web/packages/coda/index.html) R package.


```r
# Effective sample size of random effects
tab1 <- round(rbind(quantile(effectiveSize(beta_RF1),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_RF2),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_RF3),c(0.25,0.5,0.75))))
rownames(tab1) <- c("Usage choice","Reversibility choice","Method choice")
knitr::kable(tab1,format='markdown')
```



|                     |  25%|  50%|  75%|
|:--------------------|----:|----:|----:|
|Usage choice         | 1584| 2305| 3405|
|Reversibility choice | 1738| 2453| 3621|
|Method choice        | 2724| 3121| 3781|

```r
# Effective sample size of splines coefficients
tab2 <- round(rbind(
quantile(effectiveSize(beta_spline1),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_spline2),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_spline3),c(0.25,0.5,0.75))))
rownames(tab2) <- c("Usage choice","Reversibility choice","Method choice")
knitr::kable(tab2,format='markdown')
```



|                     |  25%|  50%|  75%|
|:--------------------|----:|----:|----:|
|Usage choice         | 1138| 1261| 1431|
|Reversibility choice | 1111| 1252| 1846|
|Method choice        | 1733| 1838| 2135|

```r
# Effective sample size of fixed effects
tab3 <- round(rbind(
quantile(effectiveSize(beta_Fix1),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_Fix2),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_Fix3),c(0.25,0.5,0.75))))
rownames(tab3) <- c("Usage choice","Reversibility choice","Method choice")
knitr::kable(tab3,format='markdown')
```



|                     |  25%|  50%|  75%|
|:--------------------|----:|----:|----:|
|Usage choice         | 3342| 3413| 3617|
|Reversibility choice | 2898| 3263| 3687|
|Method choice        | 2983| 3217| 3414|

