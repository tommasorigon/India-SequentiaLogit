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

# Prior distribution
prior = list(P_Fix_const=1e-2, 
             H=32,
             a_lambda=1e-3, b_lambda=1e-3,
             a_tau=1e-4, b_tau=1e-4, 
             a_alpha=.1, b_alpha=.1)

# Iterations and burn_in
R       <- 20000
burn_in <- 2000

# Age enters in the predictor linearly
f   <- as.formula(target ~ age + child + area + religion + education)
# Age is absent, since it is modeled flexibly through Bayesian P-splines.
f_s <- as.formula(target ~ child + area + religion + education)


# Dataset splitting: training and test set.
set.seed(309)
idTrain           <- sample(nrow(dataset),floor(0.75*nrow(dataset)))
data_training     <- dataset[idTrain,]
data_validation   <- dataset[-idTrain,]

# Usage choice ----------------------------------------------------
set.seed(123) # We set a seed so that our results are fully reproducible.

# We define the target variable. 
data_training$target     <- factor(data_training$method!="1. No contraceptive method")

# Estimate the submodels
fit1_ranef         <- fit_logit(f,data_training$state,data_training$age,data_training,method="ranef",prior,R,burn_in)
fit1_ranef_s       <- fit_logit(f_s,data_training$state,data_training$age,data_training,method="ranef_s",prior,R,burn_in)
fit1_dp_ranef      <- fit_logit(f,data_training$state,data_training$age,data_training,method="dp_ranef",prior,R,burn_in)

# Estimate the full model
fit1_dp_ranef_s    <- fit_logit(f_s,data_training$state,data_training$age,data_training,method="dp_ranef_s",prior,R,burn_in)

# Reversible choice ----------------------------------------------------
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
data_training2            <- data_training[data_training$method != "1. No contraceptive method",]

# We define the new target variable. 
data_training2$target     <- factor(data_training2$method != "2. Sterilization") # table(data_training2$target, data_training2$method)

# Estimate the submodels
fit2_ranef         <- fit_logit(f,data_training2$state,data_training2$age,data_training2,method="ranef",prior,R,burn_in)
fit2_ranef_s       <- fit_logit(f_s,data_training2$state,data_training2$age,data_training2,method="ranef_s",prior,R,burn_in)
fit2_dp_ranef      <- fit_logit(f,data_training2$state,data_training2$age,data_training2,method="dp_ranef",prior,R,burn_in)

# Estimate the full model
fit2_dp_ranef_s    <- fit_logit(f_s,data_training2$state,data_training2$age,data_training2,method="dp_ranef_s",prior,R,burn_in)

# Method choice ----------------------------------------------------
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
data_training3            <- data_training2[data_training2$method != "2. Sterilization",]

# We define the new target variable. 
data_training3$target     <- factor(data_training3$method != "3. Natural methods") # table(data_training3$target,data_training3$method)

# Estimate the submodels
fit3_ranef         <- fit_logit(f,data_training3$state,data_training3$age,data_training3,method="ranef",prior,R,burn_in)
fit3_ranef_s       <- fit_logit(f_s,data_training3$state,data_training3$age,data_training3,method="ranef_s",prior,R,burn_in)
fit3_dp_ranef      <- fit_logit(f,data_training3$state,data_training3$age,data_training3,method="dp_ranef",prior,R,burn_in)

# Estimate the full model
fit3_dp_ranef_s    <- fit_logit(f_s,data_training3$state,data_training3$age,data_training3,method="dp_ranef_s",prior,R,burn_in)

# Thinning
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
     file="workspaces/pred_ranef.RData")

# Splines
save(fit1_ranef_s,file="workspaces/pred_ranef_s_part1.RData")
save(fit2_ranef_s,file="workspaces/pred_ranef_s_part2.RData")
save(fit3_ranef_s,file="workspaces/pred_ranef_s_part3.RData")

# DP
save(fit1_dp_ranef,
     fit2_dp_ranef,
     fit3_dp_ranef,
     file="workspaces/pred_dp_ranef.RData")

# Full model
save(fit1_dp_ranef_s,
     fit2_dp_ranef_s,
     fit3_dp_ranef_s,
     file="workspaces/pred_dp_ranef_s.RData")

# Train and validation sets
save(data_training,
     data_validation,
     file="workspaces/train_validation.RData")