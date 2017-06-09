# Estimation and Results
Tommaso Rigon  

## Description

In this tutorial we describe how to estimate the model we described in Section~3 of our paper. We implemented our algorithm in R, making use of the `Rcpp` and `RcppArmadillo` for the most computationally intensive steps. The [R code]() is made available as well as the [C++]() code.

The Gibbs sampler for our model is contained in the function called `logit_ranefDP_spline`, which we will use extensively later on. As a first step, we load in memory all the required libraries and we compile the C++ code. Moreover, we load in memory also the [previously cleaned](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/data-cleaning.md) dataset.


```r
print(getwd())
```

```
## [1] "C:/Users/3005213/GoogleDrive/Università/Lavori/INDIA/GitHub"
```


```r
library(BayesLogit)
library(splines)
library(Rcpp)
library(RcppArmadillo)

rm(list=ls())
load("dataset.RData")
#sourceCpp("core_functions.cpp")
source("core_functions.R")
```

Some unnecessary functions are loaded in memory as well: they will be useful during the simulation.

## Prior specification and estimation

## Convergence diagnostic

## Results

#### Information criteria (DIC and WAIC)

#### Fixed effects

#### Age effect

#### Random effects: clustering and graphs



