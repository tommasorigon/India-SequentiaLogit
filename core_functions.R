
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
  
  verbose_step = ceiling(R/50)
  
  # Fixed quantities
  frame  <- model.frame(formula, data = data)
  X_Fix  <- model.matrix(frame, data = data)[,-1]
  y      <- as.numeric(model.response(frame))-1
  strata <- as.factor(strata)
  
  n     <- length(y)
  p_Fix <- NCOL(X_Fix)
  
  # Here we use all the levels, since no identifiability issues arises
  p_RF  <- length(levels(strata))
  
  # Hyperpriors
  P_Fix   <- diag(prior$P_Fix_const,p_Fix)
  a_tau   <- prior$a_tau
  b_tau   <- prior$b_tau
  tau_mu  <- prior$tau_mu
  H       <- prior$H
  
  # Output
  beta_Fix_out <- matrix(0, R, p_Fix)
  beta_RF_out  <- matrix(0, R, p_RF)
  theta_out    <- matrix(0,R,H)
  tau_out      <- numeric(R)
  S_out        <- matrix(0,R,p_RF)

  # Output for DIC and WAIC
  loglik_out   <- numeric(R)
  eta_hat      <- numeric(n)
  prob_hat     <- numeric(n)
  exp_lppd     <- numeric(n)
  log_pdf_hat  <- numeric(n)
  
  # Initialization
  beta_RF   <- numeric(p_RF)
  beta_Fix  <- numeric(p_Fix)
  eta_RF    <- numeric(n)
  eta_Fix   <- numeric(n)
  eta       <- numeric(n)
  omega     <- numeric(n) 
  S         <- factor(rep(1,p_RF),levels=1:H) # All the intercept are allocated in the first cluster
  theta_RF  <- numeric(H)
  tau       <- a_tau/b_tau
  
  # Gibbs sampling
  for (r in 1:(R*thinning + burn_in)) {
    
    # Step 1: Updating omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)
    
    # Step 3: Updating the fixed effects (beta_Fix coefficients)
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*eta_RF))
    A1        <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix  <- mu_Fix + c(matrix(rnorm(1 * p_Fix), 1, p_Fix) %*% A1)
    eta_Fix   <- X_Fix %*% beta_Fix
    
    # Step 5: Updating the random effects (beta_RF coefficients)
    reg       <- tapply(y - 0.5 - omega*eta_Fix, strata, sum) + theta_RF[S]*tau
    Sigma_RF  <- 1/(tapply(omega, strata, sum) + tau)
    mu_RF     <- Sigma_RF * reg
    beta_RF   <- rnorm(p_RF, mu_RF, sqrt(Sigma_RF))
    eta_RF    <- beta_RF[as.numeric(strata)]
    
    # Linear predictor
    eta <- as.numeric(eta_RF + eta_Fix)
    
    # Step 6: Updating the mixture variance
    tau <- rgamma(1, a_tau + p_RF/2, b_tau + sum((beta_RF - theta_RF[S])^2)/2)
    
    # Step 7: Updating the mixture mean
    n_S             <- as.numeric(table(S)) # Observations for each cluster
    reg             <- tapply(beta_RF*tau, S, sum)
    reg[is.na(reg)] <- 0 # If in some cluster there are no observation set the mean to 0.
    Sigma_theta_RF  <- 1/(tau_mu + tau*n_S) 
    mu_theta_RF     <- Sigma_theta_RF * reg
    theta_RF        <- rnorm(H, mu_theta_RF, sqrt(Sigma_theta_RF))
    
    # Step 4: Updating the mixture weights
    nu  <- c(rdir(1, 1/H + n_S)) # If H=1 it returns nu = 1.
    
    # Step 8: Updating the cluster indicator
    if(H > 1){
      lprobs <- t(sapply(beta_RF, function(x) log(nu) + dnorm(x,theta_RF,sqrt(1/tau),log=TRUE)))
      probs  <- exp(t(apply(lprobs,1, function(x) x - max(x)))) # Numerical stability!
      # probs  <- probs/rowSums(probs). Not necessary if the sample command is used.
      S      <- factor(apply(probs,1,function(x) sample(H,1,prob = x)),levels=1:H)
    }
    
    # Output
    if (r > burn_in & ((r - burn_in) %% thinning == 0)) {
      rr                  <- floor((r - burn_in)/thinning)
      beta_RF_out[rr, ]   <- beta_RF
      beta_Fix_out[rr, ]  <- beta_Fix
      theta_out[rr,]      <- theta_RF
      tau_out[rr]         <- tau
      
      log_pdf            <-  y*eta - log(1 + exp(eta))
      log_pdf_hat        <- ((rr-1)*log_pdf_hat + log_pdf)/rr
      loglik_out[rr]     <-  sum(log_pdf)
      S_out[rr,]         <- as.numeric(S)
      
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
  
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  
  list(beta_Fix = beta_Fix_out, 
       beta_RF = beta_RF_out, 
       theta=theta_out,
       tau = tau_out, 
       S =  S_out,
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
logit_ranefDP_spline <- function(formula, strata, x_spline, data,inner_knots, degree, dif,prior, R , burn_in, thinning, verbose) {
  
  verbose_step = ceiling(R/50)
  
  # Fixed quantities
  frame  <- model.frame(formula, data = data)
  X_Fix <- model.matrix(frame, data = data)[,-1]
  y      <- as.numeric(model.response(frame))-1
  strata <- as.factor(strata); 
  
  n     <- length(y)
  p_Fix <- NCOL(X_Fix)
  
  # Remove one element, since the first is used as baseline. 
  p_RF  <- length(levels(strata)) - 1; 
  
  # Hyperpriors
  P_Fix   <- diag(prior$P_Fix_const,p_Fix)
  a_tau   <- prior$a_tau
  b_tau   <- prior$b_tau
  tau_mu  <- prior$tau_mu
  H       <- prior$H
  a_lambda      <- prior$a_lambda
  b_lambda      <- prior$b_lambda
  
  # Knots placement
  xl    <- min(x_spline); xr <- max(x_spline); dx <- (xr - xl) / (inner_knots-1)
  knots <- seq(xl - degree * dx, xr + degree * dx, by = dx)
  
  # Fixed quantities
  B        <- spline.des(knots, x_spline, degree + 1, 0 * x_spline, outer.ok=TRUE)$design
  rankD    <- NCOL(B) - dif
  p_spline <- NCOL(B)
  
  if(dif==0) {D <- diag(p_spline)} else{D <- diff(diag(p_spline),dif=dif)}
  DtD <- crossprod(D) 
  
  # Output
  beta_Fix_out <- matrix(0, R, p_Fix)
  beta_RF_out  <- matrix(0, R, p_RF)
  theta_out    <- matrix(0,R,H)
  tau_out      <- numeric(R)
  S_out        <- matrix(0,R,p_RF)
  beta_spline_out <- matrix(0, R, p_spline)
  lambda_out      <- numeric(R)
  
  # Output for DIC and WAIC
  loglik_out   <- numeric(R)
  eta_hat      <- numeric(n)
  prob_hat     <- numeric(n)
  exp_lppd     <- numeric(n)
  log_pdf_hat  <- numeric(n)
  
  # Initialization
  beta_RF   <- numeric(p_RF)
  beta_Fix  <- numeric(p_Fix)
  beta_spline   <- numeric(p_spline)
  eta_spline    <- numeric(n)
  eta_RF    <- numeric(n)
  eta_Fix   <- numeric(n)
  eta       <- numeric(n)
  omega     <- numeric(n) 
  S         <- factor(rep(1,p_RF),levels=1:H) # All the intercept are allocated in the first cluster
  theta_RF  <- numeric(H)
  lambda    <- a_lambda/b_lambda
  tau       <- a_tau/b_tau

  # Gibbs sampling
  for (r in 1:(R*thinning + burn_in)) {
    
    # Step 1: Updating omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)
    
    # Step 2: Updating gamma (the beta_splines coefficient)
    eig          <- eigen(crossprod(B * sqrt(omega)) + lambda*DtD, symmetric = TRUE)
    Sigma_spline <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_spline    <- Sigma_spline %*% (crossprod(B, y - 0.5- omega*(eta_RF + eta_Fix)))
    A1           <- t(eig$vectors)/sqrt(eig$values)
    beta_spline  <- mu_spline + c(matrix(rnorm(1 * p_spline), 1, p_spline) %*% A1)
    eta_spline   <- as.numeric(B%*%beta_spline) 
    
    # Step 9: Updating the smoothing parameter lambda
    mahalanob       <- as.numeric(crossprod(D%*%beta_spline)) 
    a_lambda_tilde  <- a_lambda + rankD/2
    b_lambda_tilde  <- b_lambda + mahalanob/2
    lambda          <- rgamma(1, a_lambda_tilde, b_lambda_tilde)
    
    # Step 3: Updating the fixed effects (beta_Fix coefficients)
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*(eta_RF + eta_spline)))
    A1        <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix  <- mu_Fix + c(matrix(rnorm(1 * p_Fix), 1, p_Fix) %*% A1)
    eta_Fix   <- X_Fix %*% beta_Fix
    
    # Step 5: Updating the random effects (beta_RF coefficients)
    reg       <- tapply(y - 0.5 - omega*(eta_Fix + eta_spline), strata, sum)[-1] + theta_RF[S]*tau
    Sigma_RF  <- 1/(tapply(omega, strata, sum) + tau)[-1]
    mu_RF     <- Sigma_RF * reg
    beta_RF   <- rnorm(p_RF, mu_RF, sqrt(Sigma_RF))
    eta_RF    <- c(0,beta_RF)[as.numeric(strata)]
    
    # Linear predictor
    eta <- as.numeric(eta_RF + eta_Fix + eta_spline)
    
    # Step 6: Updating the mixture variance
    tau <- rgamma(1, a_tau + p_RF/2, b_tau + sum((beta_RF - theta_RF[S])^2)/2)
    
    # Step 7: Updating the mixture mean
    n_S             <- as.numeric(table(S)) # Observations for each cluster
    reg             <- tapply(beta_RF*tau, S, sum)
    reg[is.na(reg)] <- 0 # If in some cluster there are no observation set the mean to 0.
    Sigma_theta_RF  <- 1/(tau_mu + tau*n_S) 
    mu_theta_RF     <- Sigma_theta_RF * reg
    theta_RF        <- rnorm(H, mu_theta_RF, sqrt(Sigma_theta_RF))
    
    # Step 4: Updating the mixture weights
    nu  <- c(rdir(1, 1/H + n_S)) # If H=1 it returns nu = 1.
    
    # Step 8: Updating the cluster indicator
    if(H > 1){
      lprobs <- t(sapply(beta_RF, function(x) log(nu) + dnorm(x,theta_RF,sqrt(1/tau),log=TRUE)))
      probs  <- exp(t(apply(lprobs,1, function(x) x - max(x)))) # Numerical stability!
      # probs  <- probs/rowSums(probs). Not necessary if the sample command is used.
      S      <- factor(apply(probs,1,function(x) sample(H,1,prob = x)),levels=1:H)
    }
    
    # Output
    if (r > burn_in & ((r - burn_in) %% thinning == 0)) {
      rr                  <- floor((r - burn_in)/thinning)
      beta_RF_out[rr, ]   <- beta_RF
      beta_Fix_out[rr, ]  <- beta_Fix
      beta_spline_out[rr,]<- beta_spline
      lambda_out[rr]      <- lambda
      theta_out[rr,]       <- theta_RF
      tau_out[rr]         <- tau

      log_pdf            <-  y*eta - log(1 + exp(eta))
      log_pdf_hat        <- ((rr-1)*log_pdf_hat + log_pdf)/rr
      loglik_out[rr]     <-  sum(log_pdf)
      S_out[rr,]         <- as.numeric(S)
      
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
  
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  
  list(beta_Fix = beta_Fix_out, 
       beta_RF = beta_RF_out, 
       beta_spline=beta_spline_out,
       theta=theta_out,
       tau = tau_out, 
       lambda=lambda_out,
       S =  S_out,
       prob_hat = prob_hat,
       eta_hat = eta_hat,
       log_pdf = log_pdf_hat,
       lppd    = log(exp_lppd), 
       loglik  = loglik_out,
       loglik_hat=loglik_hat,
       parameters=(p_Fix+p_RF+p_spline)
  )
}

# Global wrapper
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
  p_DIC  <- 2*(model$loglik_hat - mean(model$loglik))
  DIC    <- -2*model$loglik_hat + 2*p_DIC

  p_WAIC <- 2*sum( model$lppd - model$log_pdf)
  WAIC   <- sum(model$lppd) - p_WAIC

  return(cbind(DIC,WAIC,p,p_DIC,p_WAIC))
}