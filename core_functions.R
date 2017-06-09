logit_ranef <- function(formula, strata, data, prior, R, burn_in,  verbose) {
  
  verbose_step = ceiling(R/50)
  
  # Fixed quantities
  frame <- model.frame(formula, data = data)
  X_Fix <- model.matrix(frame, data = data)[,-1] # REMOVE INTERCEPT.
  y     <- as.numeric(model.response(frame))-1
  strata <- as.factor(strata)
  
  n     <- length(y)
  p_Fix <- NCOL(X_Fix)
  p_RF  <- length(levels(strata))
  
  # Hyperpriors
  P_Fix <- diag(prior$P_Fix_const,p_Fix)
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Output
  beta_Fix_out <- matrix(0, R, p_Fix)
  beta_RF_out  <- matrix(0, R, p_RF)
  tau_out      <- numeric(R)
  
  # Output for DIC and WAIC
  loglik_out   <- numeric(R)
  eta_hat      <- numeric(n)
  prob_hat     <- numeric(n)
  exp_lppd     <- numeric(n)
  log_pdf_hat  <- numeric(n)
  
  # Initialization
  beta_RF   <- numeric(p_RF)
  beta_Fix  <- numeric(p_Fix)
  eta_RF  <- numeric(n)
  eta_Fix <- numeric(n)
  eta     <- numeric(n)
  omega <- numeric(n)
  
  # Gibbs sampling
  for (r in 1:(R + burn_in)) {
    
    # Updating random effects
    tau <- rgamma(1, a_tau + p_RF/2, b_tau + sum(beta_RF^2)/2)
    
    # Updating Omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)

    # Updating Random Effects
    reg       <- tapply(y - 0.5 - omega*eta_Fix, strata, sum)
    Sigma_RF  <- 1/(tapply(omega, strata, sum) + tau)
    mu_RF     <- Sigma_RF * reg
    beta_RF   <- rnorm(p_RF, mu_RF, sqrt(Sigma_RF))
    
    # Linear predictor random effects
    eta_RF <- beta_RF[as.numeric(strata)]
    
    # Updating Fixed Effects
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*eta_RF))
    
    A1       <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix <- mu_Fix + c(matrix(rnorm(1 * p_Fix), 1, p_Fix) %*% A1)
    
    # Linear predictor fixed effects
    eta_Fix <- X_Fix %*% beta_Fix
    
    # Linear predictor - updating
    eta <- as.numeric(eta_RF + eta_Fix)
    
    # Output
    if (r > burn_in) {
      rr                <- r - burn_in
      beta_RF_out[rr, ] <- beta_RF
      beta_Fix_out[rr, ]<- beta_Fix
      tau_out[rr]       <- tau
      
      log_pdf           <-  y*eta - log(1 + exp(eta))
      log_pdf_hat       <- ((rr-1)*log_pdf_hat + log_pdf)/rr
      loglik_out[rr]    <-  sum(log_pdf)
      eta_hat           <-  ((rr-1)*eta_hat + eta)/rr
      prob_hat          <-  ((rr-1)*prob_hat + plogis(eta))/rr
      exp_lppd          <-  ((rr-1)*exp_lppd + exp(log_pdf))/rr
    }
    
    if (verbose) {
      if (r%%verbose_step == 0) 
        cat(paste("Sampling iteration: ", r, ".\n", 
                  sep = ""))
    }
  }
  
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  list(beta_Fix = beta_Fix_out, beta_RF = beta_RF_out, tau = tau_out, 
       eta_hat = eta_hat,
       prob_hat = prob_hat,
       
       log_pdf = log_pdf_hat,
       lppd    = log(exp_lppd), 
       loglik  = loglik_out,
       loglik_hat=loglik_hat,
       parameters=(p_Fix+p_RF)
       )
}

logit_ranef_spline <- function(formula, x_spline, strata, data, prior, inner_knots, degree, dif, R, burn_in, verbose) {
  
  verbose_step = ceiling(R/50)
  # Fixed quantities
  
  frame <- model.frame(formula, data = data)
  X_Fix <- model.matrix(frame, data = data)[,-1] # REMOVE INTERCEPT.
  y     <- as.numeric(model.response(frame))-1
  strata <- as.factor(strata)
  
  n     <- length(y)
  p_Fix <- NCOL(X_Fix)
  p_RF  <- length(levels(strata)) - 1; levels(strata) <- 0:p_RF
  
  # Hyperpriors
  P_Fix <- diag(prior$P_Fix_const,p_Fix)
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  a_lambda      <- prior$a_lambda
  b_lambda      <- prior$b_lambda
  
  # Knots placement
  xl    <- min(x_spline); xr <- max(x_spline); dx <- (xr - xl) / (inner_knots-1)
  knots <- seq(xl - degree * dx, xr + degree * dx, by = dx)
  
  # Fixed quantities
  B        <- spline.des(knots, x_spline, degree + 1, 0 * x_spline, outer.ok=TRUE)$design
  rankB    <- NCOL(B) - dif
  p_spline <- NCOL(B)
  
  if(dif==0) {D <- diag(p_spline)} else{D <- diff(diag(p_spline),dif=dif)}
  DtD <- crossprod(D) 
  
  # Output
  beta_Fix_out <- matrix(0, R, p_Fix)
  beta_RF_out  <- matrix(0, R, p_RF)
  tau_out      <- numeric(R)
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
  eta_RF  <- numeric(n)
  eta_Fix <- numeric(n)
  eta     <- numeric(n)
  omega <- numeric(n)
  
  # Gibbs sampling
  for (r in 1:(R + burn_in)) {
    
    # Beta quantities
    mahalanob     <- as.numeric(crossprod(D%*%beta_spline)) 
    
    # Smoothing component
    a_lambda_tilde  <- a_lambda + rankB/2
    b_lambda_tilde  <- b_lambda + mahalanob/2
    lambda          <- rgamma(1, a_lambda_tilde, b_lambda_tilde)
    
    # Updating random effects
    tau <- rgamma(1, a_tau + p_RF/2, b_tau + sum(beta_RF^2)/2)
    
    # Updating Omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)
    
    # Updating Random Effects
    reg       <- tapply(y - 0.5 - omega*(eta_Fix + eta_spline), strata, sum)[-1]
    Sigma_RF  <- 1/(tapply(omega, strata, sum) + tau)[-1]
    mu_RF     <- Sigma_RF * reg
    beta_RF   <- rnorm(p_RF, mu_RF, sqrt(Sigma_RF))
    
    # Linear predictor random effects
    eta_RF <- c(0,beta_RF)[as.numeric(strata)]
    
    # Updating Fixed Effects
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*(eta_RF + eta_spline)))
    
    A1       <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix <- mu_Fix + c(matrix(rnorm(1 * p_Fix), 1, p_Fix) %*% A1)
    
    # Linear predictor fixed effects
    eta_Fix <- X_Fix %*% beta_Fix
    
    # Updating Splines Effects
    eig          <- eigen(crossprod(B * sqrt(omega)) + lambda*DtD, symmetric = TRUE)
    Sigma_spline <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_spline    <- Sigma_spline %*% (crossprod(B, y - 0.5- omega*(eta_RF + eta_Fix)))
    
    A1       <- t(eig$vectors)/sqrt(eig$values)
    beta_spline <- mu_spline + c(matrix(rnorm(1 * p_spline), 1, p_spline) %*% A1)
    
    # Linear predictor spline effect
    eta_spline    <- as.numeric(B%*%beta_spline) 
    
    # Linear predictor - updating
    eta <- as.numeric(eta_RF + eta_Fix + eta_spline)
    
    # Output
    if (r > burn_in) {
      rr                <- r - burn_in
      beta_RF_out[rr, ] <- beta_RF
      beta_Fix_out[rr, ]<- beta_Fix
      beta_spline_out[rr, ]<- beta_spline
      tau_out[rr]       <- tau
      lambda_out[rr]    <- lambda
      
      log_pdf           <-  y*eta - log(1 + exp(eta))
      log_pdf_hat       <- ((rr-1)*log_pdf_hat + log_pdf)/rr
      loglik_out[rr]    <-  sum(log_pdf)
      
      eta_hat           <-  ((rr-1)*eta_hat + eta)/rr
      prob_hat          <-  ((rr-1)*prob_hat + plogis(eta))/rr
      exp_lppd          <-  ((rr-1)*exp_lppd + exp(log_pdf))/rr
    }
    
    if (verbose) {
      if (r%%verbose_step == 0) 
        cat(paste("Sampling iteration: ", r, ".\n", 
                  sep = ""))
    }
  }
  
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  list(beta_Fix = beta_Fix_out, 
       beta_RF = beta_RF_out,
       beta_spline=beta_spline_out,
       tau = tau_out, 
       lambda=lambda_out,
       B = B,
       eta_hat = eta_hat,
       prob_hat = prob_hat,
       log_pdf = log_pdf_hat,
       lppd    = log(exp_lppd), 
       loglik  = loglik_out,
       loglik_hat=loglik_hat,
       parameters=(p_Fix+p_RF+p_spline)
  )
}

logit_ranefDP <- function(formula, strata, data, prior, R, burn_in, verbose) {
  
  verbose_step = ceiling(R/50)
  
  # Fixed quantities
  frame  <- model.frame(formula, data = data)
  X_Fix  <- model.matrix(frame, data = data)[,-1] # REMOVE INTERCEPT.
  y      <- as.numeric(model.response(frame))-1
  strata <- as.factor(strata); 
  
  n     <- length(y)
  p_Fix <- NCOL(X_Fix)
  p_RF  <- length(levels(strata)); levels(strata) <- 1:p_RF
  
  # Hyperpriors
  P_Fix <- diag(prior$P_Fix_const,p_Fix)
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  a_alpha <- prior$a_alpha
  b_alpha <- prior$b_alpha
  H       <- prior$H
  
  
  # Output
  beta_Fix_out <- matrix(0, R, p_Fix)
  beta_RF_out  <- matrix(0, R, p_RF)
  tau_out      <- numeric(R)
  alpha_out    <- numeric(R)
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
  S         <- sample(H,p_RF,replace=TRUE)
  theta_RF  <- numeric(H)
  alpha     <- 1
  
  # Gibbs sampling
  for (r in 1:(R + burn_in)) {
    
    # Updating random effects
    tau <- rgamma(1, a_tau + H/2, b_tau + sum(theta_RF^2)/2)
    
    # Updating Omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)
    
    # Updating the Dirichlet Process Random Effects
    
    # Updating the nu weights
    V     <- V_update(S, alpha, H); V[-H] <- pmin(V[-H],1-1e-16)
    nu    <- nu_update(V)
    
    # Cluster allocation
    S         <- S_update(y, as.numeric(strata)-1, p_RF, nu, theta_RF, eta_Fix)
    
    # Relabeling the strata according to S
    strata_DP         <- strata
    levels(strata_DP) <- S  # Levels are grouped.
    levels(strata_DP) <- c(levels(strata_DP), setdiff(levels(strata), levels(strata_DP)))
    strata_DP <- factor(strata_DP, levels = 1:H)
    
    # If missing, sample from the prior (zero mean and 1/tau variance)
    reg       <- tapply(y - 0.5 - omega*eta_Fix, strata_DP, sum); reg[is.na(reg)] <- 0 
    Sigma_RF  <- 1/(tapply(omega, strata_DP, sum) + tau); Sigma_RF[is.na(Sigma_RF)] <- 1/tau
    mu_RF     <- Sigma_RF * reg
    theta_RF  <- rnorm(H, mu_RF, sqrt(Sigma_RF))
    beta_RF   <- theta_RF[S]

    # Linear predictor of random effects
    eta_RF   <- beta_RF[as.numeric(strata)]
    
    # Updating Fixed Effects
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*eta_RF))
    
    A1       <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix <- mu_Fix + c(matrix(rnorm(1 * p_Fix), 1, p_Fix) %*% A1)
    
    # Linear predictor fixed effects
    eta_Fix  <- X_Fix %*% beta_Fix
    
    # Linear predictor - updating
    eta <- as.numeric(eta_RF + eta_Fix)
    
    # Updating alpha: cluster parameter
    alpha <- rgamma(1, a_alpha + H - 1, b_alpha - sum(log(1 - V[-H])))
    
    # Output
    if (r > burn_in) {
      rr                <- r - burn_in
      beta_RF_out[rr, ] <- beta_RF
      beta_Fix_out[rr, ]<- beta_Fix
      tau_out[rr]       <- tau
      alpha_out[rr]     <- alpha
      log_pdf           <-  y*eta - log(1 + exp(eta))
      log_pdf_hat       <- ((rr-1)*log_pdf_hat + log_pdf)/rr
      loglik_out[rr]    <-  sum(log_pdf)
      S_out[rr,]        <- S
      
      eta_hat           <-  ((rr-1)*eta_hat + eta)/rr
      prob_hat          <-  ((rr-1)*prob_hat + plogis(eta))/rr
      exp_lppd          <-  ((rr-1)*exp_lppd + exp(log_pdf))/rr
    }
    
    if (verbose) {
      if (r%%verbose_step == 0) 
        cat(paste("Sampling iteration: ", r, ".\n", 
                  sep = ""))
    }
  }
  
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  list(beta_Fix = beta_Fix_out, beta_RF = beta_RF_out, 
       tau = tau_out, 
       alpha=alpha_out,
       S =  S_out,
       eta_hat = eta_hat,
       prob_hat = prob_hat,
       log_pdf = log_pdf_hat,
       lppd    = log(exp_lppd), 
       loglik  = loglik_out,
       loglik_hat=loglik_hat,
       parameters=(p_Fix+p_RF)
  )
}

logit_ranefDP_spline <- function(formula, strata, x_spline, data,inner_knots, degree, dif,prior, R , burn_in, verbose) {
  
  verbose_step = ceiling(R/50)
  
  # Fixed quantities
  frame  <- model.frame(formula, data = data)
  X_Fix <- model.matrix(frame, data = data)[,-1]
  y      <- as.numeric(model.response(frame))-1
  strata <- as.factor(strata); 
  
  n     <- length(y)
  p_Fix <- NCOL(X_Fix)
  
  # Remove one element, since the first is used as baseline. The baseline will be regarded as "level 0"
  p_RF  <- length(levels(strata)) - 1; levels(strata) <- 0:p_RF
  
  # Hyperpriors
  P_Fix <- diag(prior$P_Fix_const,p_Fix)
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  a_alpha <- prior$a_alpha
  b_alpha <- prior$b_alpha
  H       <- prior$H
  a_lambda      <- prior$a_lambda
  b_lambda      <- prior$b_lambda
  
  # Knots placement
  xl    <- min(x_spline); xr <- max(x_spline); dx <- (xr - xl) / (inner_knots-1)
  knots <- seq(xl - degree * dx, xr + degree * dx, by = dx)
  
  # Fixed quantities
  B        <- spline.des(knots, x_spline, degree + 1, 0 * x_spline, outer.ok=TRUE)$design
  rankB    <- NCOL(B) - dif
  p_spline <- NCOL(B)
  
  if(dif==0) {D <- diag(p_spline)} else{ D <- diff(diag(p_spline),dif=dif)}
  DtD <- crossprod(D) 
  
  # Output
  beta_Fix_out <- matrix(0, R, p_Fix)
  beta_RF_out  <- matrix(0, R, p_RF)
  tau_out      <- numeric(R)
  alpha_out    <- numeric(R)
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
  S         <- sample(H,p_RF,replace=TRUE)
  theta_RF  <- numeric(H)
  alpha     <- 1
  
  # Gibbs sampling
  for (r in 1:(R + burn_in)) {
    
    # Beta quantities
    mahalanob     <- as.numeric(crossprod(D%*%beta_spline)) 
    
    # Smoothing component
    a_lambda_tilde  <- a_lambda + rankB/2
    b_lambda_tilde  <- b_lambda + mahalanob/2
    lambda          <- rgamma(1, a_lambda_tilde, b_lambda_tilde)
    
    # Updating random effects
    tau <- rgamma(1, a_tau + H/2, b_tau + sum(theta_RF^2)/2)
    
    # Updating Omega
    omega       <- rpg.devroye(num = n, n = 1, z = eta)
    
    # Updating the nu weights
    V     <- V_update(S, alpha, H); V[-H] <- pmin(V[-H],1-1e-16)
    nu    <- nu_update(V)
    
    # Cluster allocation. Notice that "strata" will start from "-1": the removed state.
    S     <- S_update(y, as.numeric(strata)-2, p_RF, nu, theta_RF, eta_Fix + eta_spline)
    
    # Relabeling the strata according to S
    strata_DP         <- strata
    levels(strata_DP) <- c(0,S)  # Levels are grouped.
    levels(strata_DP) <- c(levels(strata_DP), setdiff(levels(strata), levels(strata_DP)))
    strata_DP         <- factor(strata_DP, levels = 0:H)
    
    # If missing, sample from the prior (zero mean and 1/tau variance)
    reg       <- tapply(y - 0.5 - omega*(eta_Fix + eta_spline), strata_DP, sum)[-1]; # The first term (baseline), is  removed.
    reg[is.na(reg)] <- 0
    Sigma_RF  <- 1/(tapply(omega, strata_DP, sum) + tau)[-1]; # The first term (baseline), is removed.
    Sigma_RF[is.na(Sigma_RF)] <- 1/tau
    mu_RF     <- Sigma_RF * reg
    theta_RF  <- rnorm(H, mu_RF, sqrt(Sigma_RF))
    beta_RF   <- theta_RF[S]
    
    # Linear predictor of random effects. The "baseline" is fixed to 0.
    eta_RF   <- c(0,beta_RF)[as.numeric(strata)]
    
    # Updating Fixed Effects
    eig       <- eigen(crossprod(X_Fix * sqrt(omega)) + P_Fix, symmetric = TRUE)
    Sigma_Fix <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_Fix    <- Sigma_Fix %*% (crossprod(X_Fix, y - 0.5 - omega*(eta_RF + eta_spline)))
    
    A1       <- t(eig$vectors)/sqrt(eig$values)
    beta_Fix <- mu_Fix + c(matrix(rnorm(1 * p_Fix), 1, p_Fix) %*% A1)
    
    # Linear predictor fixed effects
    eta_Fix  <- X_Fix %*% beta_Fix
    
    # Updating Splines Effects
    eig          <- eigen(crossprod(B * sqrt(omega)) + lambda*DtD, symmetric = TRUE)
    Sigma_spline <- crossprod(t(eig$vectors)/sqrt(eig$values))
    mu_spline    <- Sigma_spline %*% (crossprod(B, y - 0.5- omega*(eta_RF + eta_Fix)))
    
    A1       <- t(eig$vectors)/sqrt(eig$values)
    beta_spline <- mu_spline + c(matrix(rnorm(1 * p_spline), 1, p_spline) %*% A1)
    
    # Linear predictor spline effect
    eta_spline    <- as.numeric(B%*%beta_spline) 
    
    # Linear predictor - updating
    eta <- as.numeric(eta_RF + eta_Fix + eta_spline)
    
    # Updating alpha: cluster parameter
    alpha <- rgamma(1, a_alpha + H - 1, b_alpha - sum(log(1 - V[-H])))
    
    # Output
    if (r > burn_in) {
      rr                 <- r - burn_in
      beta_RF_out[rr, ]  <- beta_RF
      beta_Fix_out[rr, ] <- beta_Fix
      beta_spline_out[rr,]<- beta_spline
      lambda_out[rr]     <- lambda
      tau_out[rr]        <- tau
      alpha_out[rr]      <- alpha
      
      log_pdf            <-  y*eta - log(1 + exp(eta))
      log_pdf_hat        <- ((rr-1)*log_pdf_hat + log_pdf)/rr
      loglik_out[rr]     <-  sum(log_pdf)
      S_out[rr,]         <- S
      
      prob_hat          <-  ((rr-1)*prob_hat + plogis(eta))/rr
      eta_hat           <-  ((rr-1)*eta_hat + eta)/rr
      exp_lppd          <-  ((rr-1)*exp_lppd + exp(log_pdf))/rr
    }
    
    if (verbose) {
      if (r%%verbose_step == 0) 
        cat(paste("Sampling iteration: ", r, ".\n", 
                  sep = ""))
    }
  }
  
  loglik_hat <- sum(y*eta_hat - log(1 + exp(eta_hat)))
  list(beta_Fix = beta_Fix_out, beta_RF = beta_RF_out, beta_spline=beta_spline_out,
       tau = tau_out, 
       alpha=alpha_out,
       lambda=lambda_out,
       S =  S_out,
       B = B,
       prob_hat = prob_hat,
       eta_hat = eta_hat,
       log_pdf = log_pdf_hat,
       lppd    = log(exp_lppd), 
       loglik  = loglik_out,
       loglik_hat=loglik_hat,
       parameters=(p_Fix+p_RF+p_spline)
  )
}

fit_logit <- function(formula, strata, x_spline, data, method, prior,  R = 20000, burn_in = 2000, inner_knots = min(round(length(strata)/4),40), degree=3, dif = 2, verbose = TRUE){
  
  if(method=="ranef") {
    return(logit_ranef(formula=formula, strata=strata, data=data, prior=prior, 
                       R=R, burn_in=burn_in,  verbose=verbose))
  } else if(method=="ranef_s")
  {
    return(logit_ranef_spline(formula=formula, x_spline=x_spline, strata=strata, data=data, prior=prior, 
                              inner_knots=inner_knots, degree=degree, dif=dif, R=R, burn_in=burn_in, verbose=verbose))
  } else if(method=="dp_ranef"){
    return(logit_ranefDP(formula=formula, strata=strata, data=data, 
                         prior=prior, R=R, burn_in=burn_in, verbose=verbose))
  }
    else if(method=="dp_ranef_s"){
      return(logit_ranefDP_spline(formula=formula, strata=strata, x_spline=x_spline, data=data, 
                                  inner_knots=inner_knots, degree=degree, dif=dif, 
                                  prior=prior, R =R, burn_in=burn_in, verbose=verbose))
    }
  else print("Provide a valid method.")
}


# Function for computing the Information Criteria
IC <- function(model){
  p      <- model$parameters
  p_DIC  <- 2*(model$loglik_hat - mean(model$loglik))
  DIC    <- -2*model$loglik_hat + 2*p_DIC

  p_WAIC <- 2*sum( model$lppd - model$log_pdf)
  WAIC   <- sum(model$lppd) - p_WAIC

  return(cbind(DIC,WAIC,p,p_DIC,p_WAIC))
}