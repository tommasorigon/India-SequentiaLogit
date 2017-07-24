# Estimation and convergence diagnostic

## Description

In this document we perform posterior computation for the model described in Section 2, via the Gibbs sampler outlined in Section 3 of the paper. We implemented the algorithms in R: the source [R code](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/core_functions.R) is made available in the github repository. 

The Gibbs sampler for our model can be called using the function `fit_logit()`, which we will use extensively during this document. As a first step, we load in the memory all the required libraries,  and the previously obtained [`dataset`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.md), which **is not made available** in this repository. See the [`data-cleaning.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.md) documentation for guidelines on how to obtain the final dataset under analysis.

The function `fit_logit()` also allows the estimation of the three sub-models described in Section 4.1.1 of the paper. Notice that the implementation of the simpler Gaussian random effects model, can be easily performed by fixing the number of mixture components `H = 1` in the code for our more general formulation.


```r
library(dplyr)
library(BayesLogit)
library(splines)

rm(list=ls())

# Load the cleaned dataset
load("dataset.RData")

# Load the R functions
source("core_functions.R")
```

## Prior specification and estimation

Posterior inference relies on `R=4000` MCMC samples, with a burn-in equal to `burn_in=2000`. These `R=4000` MCMC samples are obtained by running the Gibbs sampler for `20000` iterations (after burn-in), and thinning the chain every `5` iterations. We also define two `formula`. One for the case in which `age` enters the predictor linearly, and one in which `age` is modeled using Bayesian penalized splines.


```r
R       <- 4000
burn_in <- 2000

# Age enters the predictor linearly
f   <- as.formula(target ~ age + child + area + religion + education)
# Age is absent, since it is modeled flexibly through Bayesian P-splines.
f_s <- as.formula(target ~ child + area + religion + education)
```

As discussed in Section 4 of the paper, we fix the hyperparameters for the Gaussian kernels in the mixture of Gaussians prior for the state-specific effects, via a data driven approach. In the following code, provides details associated with this data driven prior calibration procedure. 

In particular, we estimate a classical generalized linear model for the `usage choice`, the `reversibility choice`, and the `method choice`, respectively. Then, we clustered together the `state`-specific parameters, treated here as fixed effect, via hieriarchical clustering with complete linkage. The number of clusters is selected via graphical inspection of the dendrograms. Then, we compute the average variance within cluster and the average of the squared cluster means, for all the models. These quantities will be helpful in specifying hyperparameters. See Section 4 of the paper for a detailed discussion.


```r
tab <- matrix(0,2,3)

# Usage choice - 3 clusters are selected
dataset$target     <- factor(dataset$method!="1. No contraceptive method")
coef1_glm          <- coef(glm(update(f,.~ state +. ),data=dataset,family="binomial"))[2:33]
cl       <- hclust(dist(coef1_glm),method = "complete")
cl_means <- tapply(coef1_glm,cutree(cl,3),mean) 
tab[1,1] <- 1/mean(cl_means^2)
tab[2,1] <- 1/mean(tapply(coef1_glm,cutree(cl,3),var))

# Reversibility choice - 4 clusters are selected
dataset2            <- dataset[dataset$method != "1. No contraceptive method",]
dataset2$target     <- factor(dataset2$method != "2. Sterilization") 
coef2_glm           <- coef(glm(update(f,.~ state +. ),data=dataset2,family="binomial"))[2:33]
cl       <- hclust(dist(coef2_glm),method = "complete")
cl_means <- tapply(coef2_glm,cutree(cl,4),mean)
tab[1,2] <- 1/mean(cl_means^2)
tab[2,2] <- 1/mean(tapply(coef2_glm,cutree(cl,3),var),na.rm=TRUE)

# Method choice - 3 cluster are selected
dataset3            <- dataset2[dataset2$method != "2. Sterilization",]
dataset3$target     <- factor(dataset3$method != "3. Natural methods") 
coef3_glm <- coef(glm(update(f,.~ state +. ),data=dataset3,family="binomial"))[2:33]
cl       <- hclust(dist(coef3_glm),method = "complete")
cl_means <- tapply(coef3_glm,cutree(cl,3),mean)

tab[1,3] <- 1/mean(cl_means^2)
tab[2,3] <- 1/mean(tapply(coef3_glm,cutree(cl,3),var))

colnames(tab) <- c("Usage choice","Reversibility choice","Method choice")
rownames(tab) <- c("Precision of the cluster means","Precision within the cluster")
knitr::kable(tab,format='markdown')
```



|                               | Usage choice| Reversibility choice| Method choice|
|:------------------------------|------------:|--------------------:|-------------:|
|Precision of the cluster means |    0.2333634|            0.0215313|     0.0111064|
|Precision within the cluster   |    5.0992340|            1.1460989|     1.6656611|

The prior hyperparameters for our model and the sub-models, discussed in Section 4, are specified through a list. Note that fixing `P_Fix_const = 1e-2` is equivalent to set **B** = (100,...,100). We derived different prior distribution for some parameters in the `usage choice` model, the `reversibility choice`, and the `method choice` model, basing this decision on the table reported above.


```r
prior1 <- list(P_Fix_const=1e-2, 
             H=32,
             a_lambda=1.5, b_lambda=5e-4,
             a_tau=0.5, b_tau=0.1,
             tau_mu = 0.2)

prior2 <- list(P_Fix_const=1e-2, 
             H=32,
             a_lambda=1.5, b_lambda=5e-4,
             a_tau=0.1, b_tau=0.1,
             tau_mu = 0.02)

prior3 <- list(P_Fix_const=1e-2, 
             H=32,
             a_lambda=1.5, b_lambda=5e-4,
             a_tau=0.15, b_tau=0.1,
             tau_mu = 0.01)
```


The estimation process **requires a non-negligible amount of time** to be completed. On a standard laptop, this will need about 4-5 hours. We made available the results of the MCMC chain in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/workspaces) folder, which can be loaded in the memory directly, without running the next steps. If the reader is not interested in obtaining by itself these workspaces, he / she  can skip to the convergence diagnostic step.

#### 1. Usage choice

The model and the corresponding sub-models for the `usage choice` are estimated below. Notice that the binary outcome is referred as `target`. For the `usage choice` model this target is `1` if the woman is using a contraceptive method, and `0` otherwise.


```r
set.seed(123) # We set a seed so that our results are fully reproducible.

# We define the target variable. 
dataset$target     <- factor(dataset$method!="1. No contraceptive method")

# Estimate the submodels
fit1_ranef         <- fit_logit(f,dataset$state,dataset$age,dataset,method="ranef",prior1,R,burn_in)
fit1_ranef_s       <- fit_logit(f_s,dataset$state,dataset$age,dataset,method="ranef_s",prior1,R,burn_in)
fit1_dp_ranef      <- fit_logit(f,dataset$state,dataset$age,dataset,method="dp_ranef",prior1,R,burn_in)

# Estimate the full model
fit1_dp_ranef_s    <- fit_logit(f_s,dataset$state,dataset$age,dataset,method="dp_ranef_s",prior1,R=R,burn_in=burn_in)
```

#### 2. Reversibility choice

For the `reversibility choice` model we restrict the analysis to the smaller dataset (called `dataset2`) comprising only women who are making use of contraception. See discussion in Section 3 of the paper. In this case the target is `1` if the woman is using a temporary method, and `0` if the woman is using sterilization.


```r
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
dataset2            <- dataset[dataset$method != "1. No contraceptive method",]

# We define the new target variable. 
dataset2$target     <- factor(dataset2$method != "2. Sterilization") 

# Estimate the submodels
fit2_ranef         <- fit_logit(f,dataset2$state,dataset2$age,dataset2,method="ranef",prior2,R,burn_in)
fit2_ranef_s       <- fit_logit(f_s,dataset2$state,dataset2$age,dataset2,method="ranef_s",prior2,R,burn_in)
fit2_dp_ranef      <- fit_logit(f,dataset2$state,dataset2$age,dataset2,method="dp_ranef",prior2,R,burn_in)

# Estimate the full model
fit2_dp_ranef_s    <- fit_logit(f_s,dataset2$state,dataset2$age,dataset2,method="dp_ranef_s",prior2,R,burn_in)
```


#### 3. Method choice

Finally for the `method choice` model we perform similar steps as before. In particular, we further restrict the focus on those women using reversible contraceptive (comprising `dataset3`). In this case the target is `1` if the woman is using a modern temporary method, and `0` if the woman is a natural temporary method.



```r
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
dataset3            <- dataset2[dataset2$method != "2. Sterilization",]

# We define the new target variable. 
dataset3$target     <- factor(dataset3$method != "3. Natural methods") # table(dataset3$target,dataset3$method)

# Estimate the submodels
fit3_ranef         <- fit_logit(f,dataset3$state,dataset3$age,dataset3,method="ranef",prior3,R,burn_in)
fit3_ranef_s       <- fit_logit(f_s,dataset3$state,dataset3$age,dataset3,method="ranef_s",prior3,R,burn_in)
fit3_dp_ranef      <- fit_logit(f,dataset3$state,dataset3$age,dataset3,method="dp_ranef",prior3,R,burn_in)

# Estimate the full model
fit3_dp_ranef_s    <- fit_logit(f_s,dataset3$state,dataset3$age,dataset3,method="dp_ranef_s",prior3,R,burn_in)
```

We end the estimation step by thinning the beta coefficients, cleaning the workspace and finally saving the results in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/workspaces) folder.


```r
# Save relevants files

# Baseline
save(fit1_ranef,
     fit2_ranef,
     fit3_ranef,
     file="workspaces/ranef.RData")

# Splines
save(fit1_ranef_s,
     fit2_ranef_s,
     fit3_ranef_s,file="workspaces/ranef_s.RData")

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

## Convergence diagnostic for the Bayesian semiparametric model

We load again everything in memory, on a clean envirnoment. As previously mentioned, the results are available in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/workspaces) folder, if the computations are excessive. We uploaded the workspaces separately due to limitation of GitHub in file size.


```r
rm(list=ls())

# Load the clean dataset
load("dataset.RData")

# Load the results of the MCMC chain
load("workspaces/ranef.RData")
load("workspaces/ranef_s.RData")
load("workspaces/dp_ranef.RData")
load("workspaces/dp_ranef_s.RData")
```

The chain are converted according the the standards of the [`coda`](https://cran.r-project.org/web/packages/coda/index.html) R package.


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

The effective sample size is monitored using the function `effectiveSize` of the [`coda`](https://cran.r-project.org/web/packages/coda/index.html) R package. We provide a summary of these quantities, reporting the quartiles and the median for each group of coefficients.


```r
# Effective sample size of RANDOM EFFECTS
tab1 <- round(rbind(quantile(effectiveSize(beta_RF1),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_RF2),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_RF3),c(0.25,0.5,0.75))))
rownames(tab1) <- c("Usage choice","Reversibility choice","Method choice")
knitr::kable(tab1,format='markdown')
```



|                     |  25%|  50%|  75%|
|:--------------------|----:|----:|----:|
|Usage choice         | 1887| 2363| 3074|
|Reversibility choice | 1629| 1992| 3075|
|Method choice        | 1244| 2822| 3401|

```r
# Effective sample size of SPLINES COEFFICIENTS
tab2 <- round(rbind(
quantile(effectiveSize(beta_spline1),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_spline2),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_spline3),c(0.25,0.5,0.75))))
rownames(tab2) <- c("Usage choice","Reversibility choice","Method choice")
knitr::kable(tab2,format='markdown')
```



|                     |  25%|  50%|  75%|
|:--------------------|----:|----:|----:|
|Usage choice         | 1122| 1242| 1532|
|Reversibility choice | 1088| 1152| 1710|
|Method choice        | 1460| 1549| 1938|

```r
# Effective sample size of FIXED EFFECTS
tab3 <- round(rbind(
quantile(effectiveSize(beta_Fix1),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_Fix2),c(0.25,0.5,0.75)),
quantile(effectiveSize(beta_Fix3),c(0.25,0.5,0.75))))
rownames(tab3) <- c("Usage choice","Reversibility choice","Method choice")
knitr::kable(tab3,format='markdown')
```



|                     |  25%|  50%|  75%|
|:--------------------|----:|----:|----:|
|Usage choice         | 3327| 3435| 3639|
|Reversibility choice | 3158| 3564| 3630|
|Method choice        | 3007| 3532| 3728|

