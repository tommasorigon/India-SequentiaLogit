# Predictive performance

## Description

In this document we describe how to reproduce the predictive performance discussed in Section 4 of our paper.

The first part of this document follows closely the estimation of the full model explained in the [`estimation.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/estimation.md) document. However, we will estimate all the models and submodels only on a fraction of the data (about the 75%).

The estimation process **requires a non-negligible amount of time** to be completed. On standard laptop, this will need some hours. We made available the results of the MCMC chain in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/workspaces) folder, which can be loaded in memory without running the following step. 

## Model estimation

We do not enters in the detail of the following code chunk, since it is almost an exact copy of the one used in the [`estimation.md`]((https://github.com/tommasorigon/India-SequentiaLogit/blob/master/estimation.md)) file. However, please note that the models are estimated using a `data_training` dataset. The validation set is stored, instead, the `data_validation` dataset.


```r
library(dplyr)
library(BayesLogit)
library(splines)

rm(list=ls())

# Load the clean dataset
load("dataset.RData")

# Load the R functions
source("core_functions.R")

# Prior distribution
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

# Iterations and burn_in
R       <- 40
burn_in <- 20

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
fit1_ranef         <- fit_logit(f,data_training$state,data_training$age,data_training,method="ranef",prior1,R,burn_in)
fit1_ranef_s       <- fit_logit(f_s,data_training$state,data_training$age,data_training,method="ranef_s",prior1,R,burn_in)
fit1_dp_ranef      <- fit_logit(f,data_training$state,data_training$age,data_training,method="dp_ranef",prior1,R,burn_in)

# Estimate the full model
fit1_dp_ranef_s    <- fit_logit(f_s,data_training$state,data_training$age,data_training,method="dp_ranef_s",prior1,R,burn_in)

# Reversible choice ----------------------------------------------------
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
data_training2            <- data_training[data_training$method != "1. No contraceptive method",]

# We define the new target variable. 
data_training2$target     <- factor(data_training2$method != "2. Sterilization") # table(data_training2$target, data_training2$method)

# Estimate the submodels
fit2_ranef         <- fit_logit(f,data_training2$state,data_training2$age,data_training2,method="ranef",prior2,R,burn_in)
fit2_ranef_s       <- fit_logit(f_s,data_training2$state,data_training2$age,data_training2,method="ranef_s",prior2,R,burn_in)
fit2_dp_ranef      <- fit_logit(f,data_training2$state,data_training2$age,data_training2,method="dp_ranef",prior2,R,burn_in)

# Estimate the full model
fit2_dp_ranef_s    <- fit_logit(f_s,data_training2$state,data_training2$age,data_training2,method="dp_ranef_s",prior2,R,burn_in)

# Method choice ----------------------------------------------------
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
data_training3            <- data_training2[data_training2$method != "2. Sterilization",]

# We define the new target variable. 
data_training3$target     <- factor(data_training3$method != "3. Natural methods") # table(data_training3$target,data_training3$method)

# Estimate the submodels
fit3_ranef         <- fit_logit(f,data_training3$state,data_training3$age,data_training3,method="ranef",prior3,R,burn_in)
fit3_ranef_s       <- fit_logit(f_s,data_training3$state,data_training3$age,data_training3,method="ranef_s",prior3,R,burn_in)
fit3_dp_ranef      <- fit_logit(f,data_training3$state,data_training3$age,data_training3,method="dp_ranef",prior3,R,burn_in)

# Estimate the full model
fit3_dp_ranef_s    <- fit_logit(f_s,data_training3$state,data_training3$age,data_training3,method="dp_ranef_s",prior3,R,burn_in)

# Save relevants files

# Baseline
save(fit1_ranef,
     fit2_ranef,
     fit3_ranef,
     file="workspaces/pred_ranef.RData")

# Splines
save(fit1_ranef_s,
     fit2_ranef_s,
     fit3_ranef_s,
     file="workspaces/pred_ranef_s.RData")

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
```

## Predictive performance

As explained in our paper, the predictive performance are obtained in two steps. We load in memory the required workspaces, as well as the required libraries. Moreover, we define the functions for computing the misclassification rate and other relevant measures.


```r
library(dplyr)     # Required to manipulate the dataset
library(splines)   # Required for computing the B-spline basis
library(reshape2)  # Manipulating data
library(ggplot2)   # Graphical library
library(dummies)   # Create dummy variables
library(nnet)      # For computing the multinomial model
library(ranger)    # Fast implementation of Random Forest
library(plotROC)   # Package for computing the ROC curve

# Clean the workspace
rm(list=ls())

# Load the clean dataset
load("dataset.RData")

# Load the results of the MCMC chain
load("workspaces/pred_ranef.RData")
load("workspaces/pred_ranef_s.RData")
load("workspaces/pred_dp_ranef.RData")
load("workspaces/pred_dp_ranef_s.RData")
load("workspaces/train_validation.RData")

# Binary misclassification rate, with a cut-off 
missclass<- function(pred,target,cutoff=0.5){
  t <- table(pred > cutoff,target)
  1 - sum(diag(t))/sum(t)
}

# Misclassification rate
missclassM<- function(pred,target){
  t <- table(apply(pred,1,which.max),target)
  1 - sum(diag(t))/sum(t)
}

# False positive rate
fpr <- function(pred,target,cutoff=0.5){
  t <- table(pred > cutoff,target)
  t[2,1]/sum(t[,1])
}

# False negative rate
fnr <- function(pred,target,cutoff=0.5){
  t <- table(pred > cutoff,target)
  t[1,2]/sum(t[,2])
}
```

Finally, we compute some other quantity of interest (e.g. the design matrix on the validation set). The `data_validation2` object defined below is a `data.frame` of women that are using a contraceptive method and that are a subset of the validation set contained in the `data_validation` dataframe.


```r
inner_knots <- 40; degree <- 3
xl          <- min(data_training$age); xr <- max(data_training$age); dx <- (xr - xl) / (inner_knots-1)
knots       <- seq(xl - degree * dx, xr + degree * dx, by = dx)

##############
# First step
##############

# DP + splines and Gaussian + splines.
B_val        <- spline.des(knots, data_validation$age, degree + 1, 0 * data_validation$age, outer.ok=TRUE)$design
X_val_RF     <- dummy(data_validation$state)[,-1]
X_val_Fix    <- model.matrix(method ~ child + area + religion + education, data = data_validation)[,-1]

# DP + Gaussian.
X_val_RF2   <- dummy(data_validation$state)
X_val_Fix2  <- model.matrix(method ~ age + child + area + religion + education, data = data_validation)[,-1]

###############
# Second step 
###############
data_validation2            <- data_validation[data_validation$method != "1. No contraceptive method",]
data_validation2$method     <- factor(data_validation2$method)

# DP + splines and Gaussian + splines. 
B_val2        <- spline.des(knots, data_validation2$age, degree + 1, 0 * data_validation2$age, outer.ok=TRUE)$design
X_val2_RF   <- dummy(data_validation2$state,drop = FALSE)[,-1]
X_val2_Fix  <- model.matrix(method ~ child + area + religion + education, data = data_validation2)[,-1]

# DP + Gaussian. 
X_val2_RF2   <- dummy(data_validation2$state,drop=FALSE)
X_val2_Fix2  <- model.matrix(method ~ age + child + area + religion + education, data = data_validation2)[,-1]
```

## Predictive performance: first step

In the first step we evaluate the performance only of the `usage choice` model. The posterior samples of the probability of using a contraceptive are obtained. For each unit, we considered the sample mean of the posterior draws.


```r
# DP + splines
rho1_dp_s <- rowMeans(1/(1+exp(-(X_val_RF%*%t(fit1_dp_ranef_s$beta_RF) + B_val%*%t(fit1_dp_ranef_s$beta_spline) + X_val_Fix%*%t(fit1_dp_ranef_s$beta_Fix)))))
# Gaussian + splines
rho1_rf_s <- rowMeans(1/(1+exp(-(X_val_RF%*%t(fit1_ranef_s$beta_RF) + B_val%*%t(fit1_ranef_s$beta_spline) + X_val_Fix%*%t(fit1_ranef_s$beta_Fix)))))
# Gaussian
rho1_rf <- rowMeans(1/(1+exp(-(X_val_RF2%*%t(fit1_ranef$beta_RF) + X_val_Fix2%*%t(fit1_ranef$beta_Fix)))))
# DP
rho1_dp <- rowMeans(1/(1+exp(-(X_val_RF2%*%t(fit1_dp_ranef$beta_RF) + X_val_Fix2%*%t(fit1_dp_ranef$beta_Fix)))))
```

#### Random forest

We used here the [`ranger`](https://cran.r-project.org/web/packages/ranger/index.html) R package, which is a fast implementation of random forest.


```r
set.seed(123)
fit_ranf1 <- ranger(target ~ state + age + child + area + religion + education, 
                    data=data_training,
                    probability=TRUE, 
                    num.trees = 1000)
rho1_ranger <- predict(fit_ranf1,data=data_validation,type="response")$prediction[,2]
```

#### ROC curve and performances

For each model, we computed the AUC, the misclassification rate, the false positive and negative rates. The ROC curves for each model are computed and the performance indexes are reported in the table below. For the misclassification, we used a `cutoff=0.75` in order to obtain a balance between false positive and false negative.


```r
target_val  <- as.numeric(data_validation$method != "1. No contraceptive method")
n_val       <- length(target_val)

# ROC curve
rocdata <- data.frame(out = rep(target_val,5), 
                      prediction = c(rho1_rf,rho1_dp, rho1_rf_s,rho1_dp_s, rho1_ranger), 
                      Model = c(
                        rep("baseline",n_val),
                        rep("DP",n_val),
                        rep("Splines",n_val),
                        rep("DP + splines",n_val), 
                        rep("Random Forest",n_val)),
                        k = "Usage choice")

p7 <- ggplot(rocdata, aes(m = prediction, d = out, col = Model)) + geom_roc(n.cuts=0) + xlab("False Positive Rate") + ylab("True Positive Rate") + theme_bw()+ facet_grid(~k)

# AUC
part1_AUC <- calc_auc(p7)$AUC

cutoff    <- 0.75 
part1_misclass <- c(missclass(rho1_rf,target_val,cutoff),
                    missclass(rho1_dp,target_val,cutoff),
                    missclass(rho1_rf_s,target_val,cutoff),
                    missclass(rho1_dp_s,target_val,cutoff),
                    missclass(rho1_ranger,target_val,cutoff))
part1_fpr <- c(fpr(rho1_rf,target_val,cutoff),
                    fpr(rho1_dp,target_val,cutoff),
                    fpr(rho1_rf_s,target_val,cutoff),
                    fpr(rho1_dp_s,target_val,cutoff),
                    fpr(rho1_ranger,target_val,cutoff))
part1_fnr <- c(fnr(rho1_rf,target_val,cutoff),
                    fnr(rho1_dp,target_val,cutoff),
                    fnr(rho1_rf_s,target_val,cutoff),
                    fnr(rho1_dp_s,target_val,cutoff),
                    fnr(rho1_ranger,target_val,cutoff))
# Output
tab_part1 <- cbind(AUC=part1_AUC,
                   Misclassification=part1_misclass,
                   FPR=part1_fpr,
                   FNR=part1_misclass)
rownames(tab_part1)<- c("baseline","DP","Splines", "DP + Splines","Random Forest")
knitr::kable(round(t(tab_part1),digits = 3),format="markdown")
```



|                  | baseline|    DP| Splines| DP + Splines| Random Forest|
|:-----------------|--------:|-----:|-------:|------------:|-------------:|
|AUC               |    0.790| 0.790|   0.799|        0.797|         0.799|
|Misclassification |    0.261| 0.261|   0.254|        0.254|         0.253|
|FPR               |    0.335| 0.335|   0.330|        0.333|         0.326|
|FNR               |    0.261| 0.261|   0.254|        0.254|         0.253|

```r
ggsave("img/Roc_curve.pdf",p7,device="pdf",width=11,height=9)
ggsave("img/Roc_curve.jpg",p7,device="jpg",width=11,height=9)
```

<!-- ![](https://raw.githubusercontent.com/tommasorigon/India-SequentiaLogit/master/img/Roc_curve.jpg) -->

## Predictive performance: second step

For the second step, we evaluated the remaining models jointly, so that we can make a fair comparison with the model of [De Oliveira et al. (2014)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0086654). Women not using a contraceptive methods were excluded also from the validation set.


```r
rho2 <- 1/(1+exp(-(X_val2_RF%*%t(fit2_dp_ranef_s$beta_RF) + B_val2%*%t(fit2_dp_ranef_s$beta_spline) + X_val2_Fix%*%t(fit2_dp_ranef_s$beta_Fix))))
rho3 <- 1/(1+exp(-(X_val2_RF%*%t(fit3_dp_ranef_s$beta_RF) + B_val2%*%t(fit3_dp_ranef_s$beta_spline) + X_val2_Fix%*%t(fit3_dp_ranef_s$beta_Fix))))
rho4 <- 1 - rho3

# DP + splines
prob_dp_s <- data.frame(Sterilization    =  rowMeans(1 - rho2),
                   TraditionalMethods    =  rowMeans(rho2*rho4), 
                   ModernMethods         =  rowMeans(rho2*rho3))

# Gaussian + splines
rho2 <- 1/(1+exp(-(X_val2_RF%*%t(fit2_ranef_s$beta_RF) + B_val2%*%t(fit2_ranef_s$beta_spline) + X_val2_Fix%*%t(fit2_ranef_s$beta_Fix))))
rho3 <- 1/(1+exp(-(X_val2_RF%*%t(fit3_ranef_s$beta_RF) + B_val2%*%t(fit3_ranef_s$beta_spline) + X_val2_Fix%*%t(fit3_ranef_s$beta_Fix))))
rho4 <- 1 - rho3
prob_rf_s <- data.frame(Sterilization         =  rowMeans(1 - rho2),
                        TraditionalMethods    =  rowMeans(rho2*rho4), 
                        ModernMethods         =  rowMeans(rho2*rho3))

# Gaussian
rho2 <- 1/(1+exp(-(X_val2_RF2%*%t(fit2_ranef$beta_RF) + X_val2_Fix2%*%t(fit2_ranef$beta_Fix))))
rho3 <- 1/(1+exp(-(X_val2_RF2%*%t(fit3_ranef$beta_RF) + X_val2_Fix2%*%t(fit3_ranef$beta_Fix))))
rho4 <- 1-rho3
prob_rf <- data.frame(Sterilization         =  rowMeans(1 - rho2),
                      TraditionalMethods    =  rowMeans(rho2*rho4), 
                      ModernMethods         =  rowMeans(rho2*rho3))

# DP
rho2 <- 1/(1+exp(-(X_val2_RF2%*%t(fit2_dp_ranef$beta_RF) + X_val2_Fix2%*%t(fit2_dp_ranef$beta_Fix))))
rho3 <- 1/(1+exp(-(X_val2_RF2%*%t(fit3_dp_ranef$beta_RF) + X_val2_Fix2%*%t(fit3_dp_ranef$beta_Fix))))
rho4 <- 1-rho3
prob_dp <- data.frame(Sterilization         =  rowMeans(1 - rho2),
                      TraditionalMethods    =  rowMeans(rho2*rho4), 
                      ModernMethods         =  rowMeans(rho2*rho3))
```

#### De Oliveira et al. (2014) model

A multinomial model is estimated in order to reproduce the model of [De Oliveira et al. (2014)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0086654). The variable `state` enters in the multinomial specification a a fixed effect. As explained in the paper, the variable `age` enters in the model via a piecewise constant specification.


```r
# The datatraing2 dataset has to be redefined.
data_training2           <- data_training[data_training$method != "1. No contraceptive method",]
data_training2$method    <- factor(data_training2$method)

# Threshold for age
data_training2$age_cut   <- cut(data_training2$age,c(14,25,34,49))
data_validation2$age_cut <- cut(data_validation2$age,c(14,25,34,49))

fit_multinom <- multinom(method ~  state + age_cut + area + religion + education + child,
                         data=data_training2,maxit=10000)
```

```
## # weights:  135 (88 variable)
## initial  value 18533.589310 
## iter  10 value 10542.411193
## iter  20 value 10126.713910
## iter  30 value 9926.304238
## iter  40 value 9770.871433
## iter  50 value 9736.944169
## iter  60 value 9728.673558
## iter  70 value 9725.030136
## iter  80 value 9724.068962
## iter  90 value 9723.720561
## final  value 9723.717673 
## converged
```

```r
prob_multinom <- predict(fit_multinom,newdata=data_validation2,type="probs")
```

#### Random forest

As before, the random forest are obtained using the `ranger` R package.


```r
set.seed(123)
fit_ranf2 <- ranger(method ~ state + age + child + area + religion + education, 
                    data=data_training2,
                    probability=TRUE, 
                    num.trees = 1000)
prob_ranger <- predict(fit_ranf2,data=data_validation2,type="response")$prediction
```

#### Misclassification rates

For this step, we computed only the misclassification rate for all the models.


```r
part2_misclass <- c(missclassM(prob_rf,data_validation2$method),
                    missclassM(prob_rf_s,data_validation2$method),
                    missclassM(prob_dp,data_validation2$method),
                    missclassM(prob_dp_s,data_validation2$method),
                    missclassM(prob_ranger,data_validation2$method),
                    missclassM(prob_multinom,data_validation2$method)
                    )

tab_part2 <- cbind(Misclassification=part2_misclass)
rownames(tab_part2)<- c("baseline","splines","DP", "DP + splines","Random Forest","Multinomial")

knitr::kable(round(t(tab_part2),digits = 3),format="markdown")
```



|                  | baseline| splines|    DP| DP + splines| Random Forest| Multinomial|
|:-----------------|--------:|-------:|-----:|------------:|-------------:|-----------:|
|Misclassification |    0.234|   0.229| 0.234|        0.229|         0.237|       0.233|

