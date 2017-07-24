# Out-of-Sample Predictive performance

## Description

In this document we reproduce the results associated with the out-of-sample predictive assessments discussed in Section 4.1.2 of the paper.

The first part of this document follows closely the estimation procedure of the full model explained in the [`estimation.md`](https://github.com/tommasorigon/India-SequentiaLogit/blob/master/estimation.md) document. However, estimation is based only on a fraction of training data (about the 75%). The remaining 25% is used for evaluating out-of-sample predictive performance.

The estimation process **requires a non-negligible amount of time** to be completed. On a standard laptop, this will need a couple of hours. We made available the results of the posterior computation in the [`workspaces`](https://github.com/tommasorigon/India-SequentiaLogit/tree/master/workspaces) folder, which can be loaded in the memory without running the following posterior computation step. 

## Posterior computation for the semiparametric Bayesian model 

We do not enter in the details of the following code, since it is almost the same of the one used in the [`estimation.md`]((https://github.com/tommasorigon/India-SequentiaLogit/blob/master/estimation.md)) document. However, please note that now the model is estimated using a `data_training` dataset. The validation set is stored, instead, in the `data_validation` dataset.


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
R       <- 4000
burn_in <- 2000

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

# Estimate the full model
fit1_dp_ranef_s    <- fit_logit(f_s,data_training$state,data_training$age,data_training,method="dp_ranef_s",prior1,R,burn_in)

# Reversible choice ----------------------------------------------------
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
data_training2            <- data_training[data_training$method != "1. No contraceptive method",]

# We define the new target variable. 
data_training2$target     <- factor(data_training2$method != "2. Sterilization") 
# Estimate the full model
fit2_dp_ranef_s    <- fit_logit(f_s,data_training2$state,data_training2$age,data_training2,method="dp_ranef_s",prior2,R,burn_in)

# Method choice ----------------------------------------------------
set.seed(123) # We set a seed so that our results are fully reproducible.

# Subsetting observations
data_training3            <- data_training2[data_training2$method != "2. Sterilization",]

# We define the new target variable. 
data_training3$target     <- factor(data_training3$method != "3. Natural methods") 
# Estimate the full model
fit3_dp_ranef_s    <- fit_logit(f_s,data_training3$state,data_training3$age,data_training3,method="dp_ranef_s",prior3,R,burn_in)

# Save relevants files
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

As explained in Section 4.1.2 of the paper, the predictive performance is assessed in two steps. We load in the memory the required workspaces, as well as the required libraries. Moreover, we define the functions for computing the misclassification rate and other relevant measures.


```r
library(dplyr)     # Required to manipulate the dataset
library(splines)   # Required for computing the B-spline basis
library(reshape2)  # Manipulating data
library(ggplot2)   # Graphical library
library(dummies)   # Create dummy variables
library(nnet)      # For computing the multinomial model
library(ranger)    # Fast implementation of random forest
library(xgboost)   # Fast implementation of gradient boosting
library(MASS)      # Will be used for the lda function
library(ROCR)      # Package for computing the ROC curve

# Clean the workspace
rm(list=ls())

# Load the clean dataset
load("dataset.RData")

# Load the results of the MCMC chain
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

# Area under the ROC curve (using the ROCR package)
AUC <- function(pred,target){
  auc <- performance(prediction(pred,target), measure = "auc")  
  auc@y.values[[1]]  
}
```

Finally, we compute the design matrix of the women in the validation set to perform out-of-sample prediction. The `data_validation2` object defined below is the `data.frame` comprising the women in the `data_validation` dataset that are using a contraceptive method.


```r
inner_knots <- 40; degree <- 3
xl          <- min(data_training$age); xr <- max(data_training$age); dx <- (xr - xl) / (inner_knots-1)
knots       <- seq(xl - degree * dx, xr + degree * dx, by = dx)

# --------------------
# First step
# --------------------

# Mixture + splines
B_val      <- spline.des(knots, data_validation$age, degree + 1, 0 * data_validation$age, outer.ok=TRUE)$design
X_val_RF   <- dummy(data_validation$state)[,-1]
X_val_Fix  <- model.matrix(method ~ child + area + religion + education, data = data_validation)[,-1]

# -------------------
# Second step 
# -------------------

data_validation2        <- data_validation[data_validation$method != "1. No contraceptive method",]
data_validation2$method <- factor(data_validation2$method)

# DP + splines 
B_val2      <- spline.des(knots, data_validation2$age, degree + 1, 0 * data_validation2$age, outer.ok=TRUE)$design
X_val2_RF   <- dummy(data_validation2$state,drop = FALSE)[,-1]
X_val2_Fix  <- model.matrix(method ~ child + area + religion + education, data = data_validation2)[,-1]
```

## Predictive performance: first step

In the first step we evaluate the performance only for the `usage choice` model. The posterior samples of the probability of using a contraceptive are obtained. For each unit, we considered the sample mean of the posterior draws.


```r
# Mixture + splines
rho1_dp_s <- rowMeans(1/(1+exp(-(X_val_RF%*%t(fit1_dp_ranef_s$beta_RF) + B_val%*%t(fit1_dp_ranef_s$beta_spline) + X_val_Fix%*%t(fit1_dp_ranef_s$beta_Fix)))))
```

#### Linear discriminant analysis (LDA)

In the code below, it is conducted a linear discriminant analysis for two categories, using the `lda` function of the `MASS` R package.


```r
set.seed(123)
fit_lda1 <- lda(target ~ state + age + child + area + religion + education, 
                    data=data_training)
rho1_lda <- predict(fit_lda1,newdata = data_validation)$posterior[,2]
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

#### Gradient Boosting

We used the [`xgboost`](http://xgboost.readthedocs.io/en/latest/) R package, which is a fast implementation for gradient boosting.


```r
set.seed(123)
# Matrix used for the estimation
xgb_training <- model.matrix(target ~ state + age + child + area + religion + education, data=data_training)[,-1] 
# Convert it into a sparse matrix - Indicate the target variable
xgb_training <- xgb.DMatrix(xgb_training, label = as.numeric(data_training$target)-1)

# List of tuning parameters
nround   <- 75
param    <- list(max_depth=3, eta=0.5, objective='binary:logistic')
fit1_xgb <- xgboost(data=xgb_training, nrounds = nround, params=param,verbose=0)

# Prediction on the validation dataset
xgb_validation <- model.matrix(method ~ state + age + child + area + religion + education, data=data_validation)[,-1] # Remove intercept
target_val     <- as.numeric(data_validation$method != "1. No contraceptive method")
xgb_validation <- xgb.DMatrix(xgb_validation, label = target_val)

# Output probabilities
rho1_xgb <- predict(fit1_xgb,newdata=xgb_validation)
```


#### ROC curve and performances

For each model, we computed the AUC and the misclassification rate. For the misclassification, we used a `cutoff=0.5`.


```r
target_val  <- as.numeric(data_validation$method != "1. No contraceptive method")
n_val       <- length(target_val)

cutoff    <- 0.5
part1_AUC <- c(AUC(rho1_dp_s,target_val),
               AUC(rho1_lda,target_val),
               AUC(rho1_ranger,target_val),
               AUC(rho1_xgb,target_val))
part1_misclass <- c(missclass(rho1_dp_s,target_val,cutoff),
                    missclass(rho1_lda,target_val,cutoff),
                    missclass(rho1_ranger,target_val,cutoff),
                    missclass(rho1_xgb,target_val,cutoff))
# Output
tab_part1 <- cbind(AUC=part1_AUC,
                   Misclassification=part1_misclass)
rownames(tab_part1)<- c("Mixture + Splines","LDA","Random Forest","Gradient Boosting")
knitr::kable(round(t(tab_part1),digits = 3),format="markdown")
```



|                  | Mixture + Splines|   LDA| Random Forest| Gradient Boosting|
|:-----------------|------------:|-----:|-------------:|-----------------:|
|AUC               |        0.799| 0.790|         0.797|             0.803|
|Misclassification |        0.197| 0.196|         0.196|             0.194|


## Predictive performance: second step

For the second step, we evaluated the remaining models jointly, so that we can make a fair comparison with the model of [De Oliveira et al. (2014)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0086654). Women not using a contraceptive methods were excluded also from the validation set. See discussion in Section 4.1.2.


```r
rho2 <- 1/(1+exp(-(X_val2_RF%*%t(fit2_dp_ranef_s$beta_RF) + B_val2%*%t(fit2_dp_ranef_s$beta_spline) + X_val2_Fix%*%t(fit2_dp_ranef_s$beta_Fix))))
rho3 <- 1/(1+exp(-(X_val2_RF%*%t(fit3_dp_ranef_s$beta_RF) + B_val2%*%t(fit3_dp_ranef_s$beta_spline) + X_val2_Fix%*%t(fit3_dp_ranef_s$beta_Fix))))
rho4 <- 1 - rho3

# Mixture + splines
prob_dp_s <- data.frame(Sterilization    =  rowMeans(1 - rho2),
                   TraditionalMethods    =  rowMeans(rho2*rho4), 
                   ModernMethods         =  rowMeans(rho2*rho3))
```

#### De Oliveira et al. (2014) model

A multinomial model is estimated in order to reproduce the model of [De Oliveira et al. (2014)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0086654). The variable `state` enters in the multinomial specification as a fixed effect. As explained in the paper, the variable `age` enters in the model via a piecewise constant specification.


```r
# The datatraing2 dataset has to be redefined.
data_training2           <- data_training[data_training$method != "1. No contraceptive method",]
data_training2$method    <- factor(data_training2$method)

# Threshold for age
data_training2$age_cut   <- cut(data_training2$age,c(14,25,34,49))
data_validation2$age_cut <- cut(data_validation2$age,c(14,25,34,49))

fit_multinom <- multinom(method ~  state + age_cut + area + religion + education + child,
                         data=data_training2,maxit=10000,trace=FALSE)

prob_multinom <- predict(fit_multinom,newdata=data_validation2,type="probs")
```

#### Linear discriminant analysis (LDA)

In the code below, it is conducted a linear discriminant analysis.


```r
set.seed(123)
fit_lda2 <- lda(method ~  state + age_cut + area + religion + education + child,
                         data=data_training2)
prob_lda <- predict(fit_lda2,newdata = data_validation2)$posterior
```

#### Random forest

As before, the random forest is  obtained using the `ranger` R package.


```r
set.seed(123)
fit_ranf2 <- ranger(method ~ state + age + child + area + religion + education, 
                    data=data_training2,
                    mtry = 3,
                    probability=TRUE, 
                    num.trees = 1000)
prob_ranger <- predict(fit_ranf2,data=data_validation2,type="response")$prediction
```

#### Gradient Boosting

As before, the gradient boosting is obtained using the `xgboost` R package.


```r
set.seed(123)
# Matrix used for the estimation
xgb_training <- model.matrix(method ~ state + age + child + area + religion + education, data=data_training2)[,-1] # Remove intercept
# Convert it into a sparse matrix - Indicate the target variable
xgb_training <- xgb.DMatrix(xgb_training, label = as.numeric(data_training2$method)-1)

# List of tuning parameters
nround   <- 50
param    <- list(max_depth=3, eta=1,num_class=3,objective='multi:softprob')
fit2_xgb <- xgboost(data=xgb_training, nrounds = nround, params=param,verbose=0)

# Prediction on the validation dataset
xgb_validation <- model.matrix(method ~ state + age + child + area + religion + education, data=data_validation2)[,-1] # Remove intercept
xgb_validation <- xgb.DMatrix(xgb_validation, label = data_validation2$method)
# Output probabilities
prob_xgb <- matrix(predict(fit2_xgb,newdata=xgb_validation),nrow(data_validation2),3,byrow = TRUE)
```


#### Misclassification rates

For this step, we computed only the misclassification rate for all the models.


```r
part2_misclass <- c(missclassM(prob_dp_s,data_validation2$method),
                    missclassM(prob_lda,data_validation2$method),
                    missclassM(prob_ranger,data_validation2$method),
                    missclassM(prob_xgb,data_validation2$method),
                    missclassM(prob_multinom,data_validation2$method)
                    )

tab_part2 <- cbind(Misclassification=part2_misclass)
rownames(tab_part2)<- c("Mixture + splines","LDA","Random Forest","Gradient Boosting","Multinomial (De Oliveira et al.)")

knitr::kable(round(t(tab_part2),digits = 3),format="markdown")
```



|                  | Mixture + splines|   LDA| Random Forest| Gradient Boosting| Multinomial (De Oliveira et al.)|
|:-----------------|------------:|-----:|-------------:|-----------------:|--------------------------------:|
|Misclassification |        0.229| 0.249|         0.234|             0.231|                            0.233|

