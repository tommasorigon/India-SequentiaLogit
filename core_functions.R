# -----------------------------------------------------------------------
# In this file it is available an implemention of the Gibbs sampler
# described in the paper, Section 3.
# ----------------------------------------------------------------------

# A function for sampling from a Dirichlet distribution
rdir <- function(n,alpha) {
  stopifnot(all(alpha> 0), n > 0)
  p   <- length(alpha)
  sim <- matrix(rgamma(n*p,alpha),ncol=p, byrow=TRUE)
  sum <- sim %*% rep(1,p)
  sim / as.vector(sum)
}

# Models without the spline component
logit_ranefDP <- function(formula, strata, data, prior, R, burn_in, thinning, verbose) {
  
  # Print the current iteration every "verbose_step" iteration. Ignore if verbose=FALSE
  verbose_step = ceiling(R/50)
  
  # Fixed quantities
  frame  <- model.frame(formula, data = data)     
  X_Fix  <- model.matrix(frame, data = data)[,-1] # Design matrix fixed effects
  y      <- as.numeric(model.response(frame))-1   # Response binary variabile
  strata <- as.factor(strata)                     # Categorical variable denoting the intercept
  
  n     <- length(y)    # Number of observations
  p_Fix <- NCOL(X_Fix)  # Number of parameters of the fixed coefficients
  p_RF  <- length(levels(strata)) # Number of parameters of the random effects. 
  
  # Hyperparameters settings
  P_Fix   <- diag(prior$P_Fix_const,p_Fix)  # Inverse of B, the prior covariance matrix of the fixed coefficients --- here called P_Fix.
  a_tau   <- prior$a_tau                    # Prior hyperparameters for the within cluster precision. See Equation (7) of the paper.
  b_tau   <- prior$b_tau                    # As above
  tau_mu  <- prior$tau_mu                   # Precision of the cluster mean. See equation (7) of the paper.
  H       <- prior$H                        # Number of mixture components. See equation (7) of the paper.
  
  # Initialization of the output
  beta_Fix_out <- matrix(0, R, p_Fix) # Fixed coefficients (Eq. 10)
  beta_RF_out  <- matrix(0, R, p_RF)  # Random intercepts (Eq. 6)
  mu_RF_bar_out<- matrix(0,R,H)       # Cluster specific means (Eq. 7)
  tau_out      <- numeric(R)          # Precision of each Gaussian mixture component (Eq. 7)
  G_out        <- matrix(0,R,p_RF)    # Cluster indicator 

  # Quantities useful for evaluating the DIC, WAIC indexes and for prediction
  loglik_out   <- numeric(R) # Log-likelihood of each iteration of the MCMC chain        - Necessary for DIC
  eta_hat      <- numeric(n) # Posterior mean of the each term of the linear predictor   - Necessary for DIC
  exp_lppd     <- numeric(n) # Posterior mean of each contribution to the likelihood     - Necessary for WAIC
  log_pdf_hat  <- numeric(n) # Posterior mean of each contribution to the log-likelihood - Necessary for WAIC 
  prob_hat     <- numeric(n) # Posterior mean of the probability                 - Necessary for prediction
  
  # Initialization
  beta_RF   <- numeric(p_RF)
  beta_Fix  <- numeric(p_Fix)
  eta_RF    <- numeric(n)                      # Linear predictor of the random effects
  eta_Fix   <- numeric(n)                      # Linear predictor of the fixed coefficients
  eta       <- numeric(n)                      # Linear predictor (sum of the above 3 terms)
  omega     <- numeric(n)                      # Omega Polya-gamma random variables
  G         <- factor(rep(1,p_RF),levels=1:H)  # Cluster indicator - All the values are allocated in the first cluster
  mu_RF_bar <- numeric(H)                     # Cluster specific means
  tau       <- a_tau/b_tau                     # Precision of each Gaussian mixture component

  # Compared to the paper, notice that some steps are in a slightly different order. Clearly, some other are omitted because there is not the spline component
    
  # Starting the Gibbs sampling
  for (r in 1:(R*thinning + burn_in)) {
    
    # Step 1] Updating the Polya gamma random variable omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)
    
    # Step 4] Updating the beta_Fix fixed coefficients
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*eta_RF))
    A1        <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix  <- mu_Fix + c(matrix(rnorm(1 * p_Fix), 1, p_Fix) %*% A1)
    eta_Fix   <- X_Fix %*% beta_Fix
    
    # Step 6] Updating the beta_RF random effects 
    reg       <- tapply(y - 0.5 - omega*eta_Fix, strata, sum) + mu_RF_bar[G]*tau
    Sigma_RF  <- 1/(tapply(omega, strata, sum) + tau)
    mu_RF     <- Sigma_RF * reg
    beta_RF   <- rnorm(p_RF, mu_RF, sqrt(Sigma_RF))
    eta_RF    <- beta_RF[as.numeric(strata)]
    
    # Update the global linear predictor
    eta <- as.numeric(eta_RF + eta_Fix)
    
    # Step 7] Updating the mixture precision
    tau <- rgamma(1, a_tau + p_RF/2, b_tau + sum((beta_RF - mu_RF_bar[G])^2)/2)
    
    # Step 8] Updating the mixture mean
    n_G              <- as.numeric(table(G)) # Number of observations for each cluster
    reg              <- tapply(beta_RF*tau, G, sum)
    reg[is.na(reg)]  <- 0 # If in some cluster there are no observation, set the mean to 0 (the hyperprior)
    Sigma_mu_RF_bar  <- 1/(tau_mu + tau*n_G) 
    mu_mu_RF_bar     <- Sigma_mu_RF_bar * reg
    mu_RF_bar        <- rnorm(H, mu_mu_RF_bar, sqrt(Sigma_mu_RF_bar))
    
    # Step 9] Updating the mixture weights
    nu  <- c(rdir(1, 1/H + n_G)) # If H=1 it returns nu = 1.
    
    # Step 5] Updating the cluster indicator
    if(H > 1){
      lprobs <- t(sapply(beta_RF, function(x) log(nu) + dnorm(x,mu_RF_bar,sqrt(1/tau),log=TRUE)))
      probs  <- exp(t(apply(lprobs,1, function(x) x - max(x)))) # Probabilities are rescaled for numerical stability
      # probs  <- probs/rowSums(probs). This step is not necessary: sample automatically perform this operation
      G      <- factor(apply(probs,1,function(x) sample(H,1,prob = x)),levels=1:H)
    }
    
    # Output
    if (r > burn_in & ((r - burn_in) %% thinning == 0)) {
      rr                  <- floor((r - burn_in)/thinning)
      beta_RF_out[rr, ]   <- beta_RF
      beta_Fix_out[rr, ]  <- beta_Fix
      mu_RF_bar_out[rr,]      <- mu_RF_bar
      tau_out[rr]         <- tau
      
      log_pdf            <-  y*eta - log(1 + exp(eta))
      log_pdf_hat        <- ((rr-1)*log_pdf_hat + log_pdf)/rr
      loglik_out[rr]     <- sum(log_pdf)
      G_out[rr,]         <- as.numeric(G)
      
      prob_hat          <-  ((rr-1)*prob_hat + plogis(eta))/rr
      eta_hat           <-  ((rr-1)*eta_hat + eta)/rr
      exp_lppd          <-  ((rr-1)*exp_lppd + exp(log_pdf))/rr
    }
    
    if (verbose) {
      if (r%%(verbose_step*thinning) == 0) 
        cat(paste("Sampling iteration: ", r, " out of ",R*thinning + burn_in, "\n",
                  sep = ""))
    }
  }
  
  # Likelihood evaluated in the Bayes estimator - Useful for computing the DIC
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  
  list(beta_Fix = beta_Fix_out, 
       beta_RF = beta_RF_out, 
       mu_RF_bar=mu_RF_bar_out,
       tau = tau_out, 
       S =  G_out,
       prob_hat = prob_hat,
       eta_hat = eta_hat,
       log_pdf = log_pdf_hat,
       lppd    = log(exp_lppd), 
       loglik  = loglik_out,
       loglik_hat=loglik_hat,
       parameters=(p_Fix+p_RF)
  )
}

# Models with the spline component
logit_ranefDP_spline <- function(formula, strata, x_spline, data, inner_knots, degree, dif, prior, R , burn_in, thinning, verbose) {
  
  # Print the current iteration every "verbose_step" iteration. Ignore if verbose=FALSE
  verbose_step = ceiling(R/50)
  
  # Fixed quantities
  frame  <- model.frame(formula, data = data)
  X_Fix  <- model.matrix(frame, data = data)[,-1] # Design matrix fixed effects
  y      <- as.numeric(model.response(frame))-1   # Response binary variabile
  strata <- as.factor(strata)                     # Categorical variable denoting the intercept
  
  n     <- length(y)    # Number of observations
  p_Fix <- NCOL(X_Fix)  # Number of parameters of the fixed coefficients
  p_RF  <- length(levels(strata)) - 1 # Number of parameters of the random effects. 
  
  # Hyperparameters settings
  P_Fix   <- diag(prior$P_Fix_const,p_Fix) # Inverse of B, the prior covariance matrix of the fixed coefficients --- here called P_Fix.
  a_tau   <- prior$a_tau                   # Prior hyperparameters for the within cluster precision. See Equation (7) of the paper.
  b_tau   <- prior$b_tau                   # As above
  tau_mu  <- prior$tau_mu                  # Precision of the cluster mean. See equation (7) of the paper.
  H       <- prior$H                       # Number of mixture components. See equation (7) of the paper.
  a_lambda      <- prior$a_lambda          # Prior hyperparameters for the smoothing parameter. See equation (9) of the paper.
  b_lambda      <- prior$b_lambda          # As above
  
  # Placements of the knots for the spline basis
  xl    <- min(x_spline); xr <- max(x_spline); dx <- (xr - xl) / (inner_knots-1)
  knots <- seq(xl - degree * dx, xr + degree * dx, by = dx)
  
  # Fixed quantities related to P-splines
  X_Spline <- spline.des(knots, x_spline, degree + 1, 0 * x_spline, outer.ok=TRUE)$design # B-spline design matrix 
  rankD    <- NCOL(X_Spline) - dif # Rank of the smoothing matrix D
  p_Spline <- NCOL(X_Spline)       # Number of parameters for the spline
  
  # Creation of the penalty matrix D as for Equation (9) of the paper
  if(dif==0) {D_root <- diag(p_Spline)} else{D_root <- diff(diag(p_Spline),dif=dif)}
  D <- crossprod(D_root) 
  
  # Initialization of the output
  beta_Fix_out <- matrix(0, R, p_Fix)        # Fixed coefficients (Eq. 10)
  beta_RF_out  <- matrix(0, R, p_RF)         # Random intercepts (Eq. 6)
  beta_spline_out  <- matrix(0, R, p_Spline) # Spline coefficients (Eq. 9)
  mu_RF_bar_out    <- matrix(0,R,H)          # Cluster specific means (Eq. 7)
  tau_out          <- numeric(R)             # Precision of each Gaussian mixture component (Eq. 7)
  G_out            <- matrix(0,R,p_RF)       # Cluster indicator 
  lambda_out       <- numeric(R)             # Smoothing parameter 
  
  # Quantities useful for evaluating the DIC, WAIC indexes and for prediction
  loglik_out   <- numeric(R) # Log-likelihood of each iteration of the MCMC chain        - Necessary for DIC
  eta_hat      <- numeric(n) # Posterior mean of the each term of the linear predictor   - Necessary for DIC
  exp_lppd     <- numeric(n) # Posterior mean of each contribution to the likelihood     - Necessary for WAIC
  log_pdf_hat  <- numeric(n) # Posterior mean of each contribution to the log-likelihood - Necessary for WAIC 
  prob_hat     <- numeric(n) # Posterior mean of the probability                 - Necessary for prediction
  
  # Initialization
  beta_RF   <- numeric(p_RF) 
  beta_Fix  <- numeric(p_Fix)
  beta_spline   <- numeric(p_Spline)
  eta_spline    <- numeric(n) # Linear predictor of the splines 
  eta_RF    <- numeric(n)     # Linear predictor of the random effects
  eta_Fix   <- numeric(n)     # Linear predictor of the fixed coefficients
  eta       <- numeric(n)     # Linear predictor (sum of the above 3 terms)
  omega     <- numeric(n)     # Omega Polya-gamma random variables
  G         <- factor(rep(1,p_RF),levels=1:H) # Cluster indicator - All the values are allocated in the first cluster
  mu_RF_bar <- numeric(H)        # Cluster specific means
  lambda    <- a_lambda/b_lambda # Smoothing parameters
  tau       <- a_tau/b_tau       # Precision of each Gaussian mixture component

  # Compared to the paper, notice that some steps are in a slightly different order.
  
  # Starting the Gibbs sampling
  for (r in 1:(R*thinning + burn_in)) {
    
    # Step 1] Updating the Polya gamma random variable omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)
    
    # Step 2] Updating the beta_splines splines coefficient
    eig          <- eigen(crossprod(X_Spline * sqrt(omega)) + lambda*D, symmetric = TRUE)
    Sigma_spline <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_spline    <- Sigma_spline %*% (crossprod(X_Spline, y - 0.5- omega*(eta_RF + eta_Fix)))
    A1           <- t(eig$vectors)/sqrt(eig$values)
    beta_spline  <- mu_spline + c(matrix(rnorm(p_Spline), 1, p_Spline) %*% A1)
    eta_spline   <- as.numeric(X_Spline%*%beta_spline) 
    
    # Step 3] Updating the smoothing parameter lambda
    mahalanob       <- as.numeric(crossprod(D_root%*%beta_spline)) 
    a_lambda_tilde  <- a_lambda + rankD/2
    b_lambda_tilde  <- b_lambda + mahalanob/2
    lambda          <- rgamma(1, a_lambda_tilde, b_lambda_tilde)
    
    # Step 4] Updating the beta_Fix fixed coefficients
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*(eta_RF + eta_spline)))
    A1        <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix  <- mu_Fix + c(matrix(rnorm(p_Fix), 1, p_Fix) %*% A1)
    eta_Fix   <- X_Fix %*% beta_Fix
    
    # Step 6] Updating the beta_RF random effects 
    reg       <- tapply(y - 0.5 - omega*(eta_Fix + eta_spline), strata, sum)[-1] + mu_RF_bar[G]*tau
    Sigma_RF  <- 1/(tapply(omega, strata, sum) + tau)[-1]
    mu_RF     <- Sigma_RF * reg
    beta_RF   <- rnorm(p_RF, mu_RF, sqrt(Sigma_RF))
    eta_RF    <- c(0,beta_RF)[as.numeric(strata)]
    
    # Update the global linear predictor
    eta <- as.numeric(eta_RF + eta_Fix + eta_spline)
    
    # Step 7] Updating the mixture precision
    tau <- rgamma(1, a_tau + p_RF/2, b_tau + sum((beta_RF - mu_RF_bar[G])^2)/2)
    
    # Step 8] Updating the mixture mean
    n_G              <- as.numeric(table(G))        # Number of observations for each cluster
    reg              <- tapply(beta_RF*tau, G, sum)
    reg[is.na(reg)]  <- 0   # If in some cluster there are no observation, set the mean to 0 (the hyperprior)
    Sigma_mu_RF_bar  <- 1/(tau_mu + tau*n_G) 
    mu_mu_RF_bar     <- Sigma_mu_RF_bar * reg
    mu_RF_bar        <- rnorm(H, mu_mu_RF_bar, sqrt(Sigma_mu_RF_bar))
    
    # Step 9] Updating the mixture weights
    nu  <- c(rdir(1, 1/H + n_G)) # If H = 1 it works and returns nu = 1.
    
    # Step 5] Updating the cluster indicator
    if(H > 1){
      lprobs <- t(sapply(beta_RF, function(x) log(nu) + dnorm(x,mu_RF_bar,sqrt(1/tau),log=TRUE)))
      probs  <- exp(t(apply(lprobs,1, function(x) x - max(x)))) # Probabilities are rescaled for numerical stability
      # probs  <- probs/rowSums(probs). This step is not necessary: sample automatically perform this operation
      G      <- factor(apply(probs,1,function(x) sample(H,1,prob = x)),levels=1:H)
    }
    
    # Output
    if (r > burn_in & ((r - burn_in) %% thinning == 0)) {
      rr                  <- floor((r - burn_in)/thinning) 
      beta_RF_out[rr, ]   <- beta_RF
      beta_Fix_out[rr, ]  <- beta_Fix
      beta_spline_out[rr,]<- beta_spline
      lambda_out[rr]      <- lambda
      mu_RF_bar_out[rr,]  <- mu_RF_bar
      tau_out[rr]         <- tau

      log_pdf            <-  y*eta - log(1 + exp(eta)) 
      log_pdf_hat        <- ((rr-1)*log_pdf_hat + log_pdf)/rr 
      loglik_out[rr]     <- sum(log_pdf)
      G_out[rr,]         <- as.numeric(G)
      
      prob_hat          <-  ((rr-1)*prob_hat + plogis(eta))/rr
      eta_hat           <-  ((rr-1)*eta_hat + eta)/rr
      exp_lppd          <-  ((rr-1)*exp_lppd + exp(log_pdf))/rr
    }
    
    if (verbose) {
      if (r%%(verbose_step*thinning) == 0) 
        cat(paste("Sampling iteration: ", r, " out of ",R*thinning + burn_in, "\n",
                  sep = ""))
    }
  }
  
  # Likelihood evaluated in the Bayes estimator - Useful for computing the DIC
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  
  list(beta_Fix = beta_Fix_out, 
       beta_RF = beta_RF_out, 
       beta_spline=beta_spline_out,
       mu_RF_bar=mu_RF_bar_out,
       tau = tau_out, 
       lambda=lambda_out,
       S =  G_out,
       prob_hat = prob_hat,
       eta_hat = eta_hat,
       log_pdf = log_pdf_hat,
       lppd    = log(exp_lppd), 
       loglik  = loglik_out,
       loglik_hat=loglik_hat,  
       parameters=(p_Fix+p_RF+p_Spline)
  )
}

# Global function that allows to call two functions above with a unique interface
fit_logit <- function(formula, strata, x_spline, data, method, prior,  R = 20000, burn_in = 2000, thinning=5, inner_knots = min(round(length(strata)/4),40), degree=3, dif = 2, verbose = TRUE){
  
  if(any(table(strata)==0)){stop("In the provided strata vector some levels have no observations")}
  if(method=="ranef") {
    prior$H <- 1
    return(logit_ranefDP(formula=formula, strata=strata, data=data, prior=prior, 
                       R=R, burn_in=burn_in,  thinning=thinning, verbose=verbose))
  } else if(method=="ranef_s"){ 
    prior$H <- 1
    return(logit_ranefDP_spline(formula=formula, x_spline=x_spline, strata=strata, data=data, prior=prior, 
                              inner_knots=inner_knots, degree=degree, dif=dif, 
                              R=R, burn_in=burn_in, thinning=thinning, verbose=verbose))
  } else if(method=="dp_ranef"){
    return(logit_ranefDP(formula=formula, strata=strata, data=data, 
                         prior=prior, R=R, burn_in=burn_in, thinning=thinning, verbose=verbose))
  } else if(method=="dp_ranef_s"){
      return(logit_ranefDP_spline(formula=formula, strata=strata, x_spline=x_spline, data=data, 
                                  inner_knots=inner_knots, degree=degree, dif=dif, 
                                  prior=prior, R =R, burn_in=burn_in, thinning=thinning, verbose=verbose))
    }
  else print("Provide a valid method.")
}

# Function for computing the Information Criteria DIC and WAIC
IC <- function(model){
  p      <- model$parameters
  # As in Gleman et. al (2014), equation (8)
  p_DIC  <- 2*(model$loglik_hat - mean(model$loglik))
  DIC    <- -2*model$loglik_hat + 2*p_DIC

  p_WAIC <- 2*sum( model$lppd - model$log_pdf)
  WAIC   <- sum(model$lppd) - p_WAIC

  return(cbind(DIC,WAIC,p,p_DIC,p_WAIC))
}
