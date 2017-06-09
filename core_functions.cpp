#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector nu_update(NumericVector V){    
  
  int H = V.length();
  NumericVector nu(H);
  NumericVector cumv(H);
  
  cumv[0] = 1-V[0];
  nu[0]= V[0];
  
  for(int h =1; h < H; h++) {
    cumv[h] = cumv[h-1]*(1-V[h]);
    nu[h]= V[h] * cumv[h-1];
  }
  return nu;
}


// [[Rcpp::export]]
NumericVector V_update(arma::vec S, double alpha, int H) {
  
  NumericVector V(H);
  arma::uvec nc;
  arma::uvec ncp;

  for(int h =1; h < H; h++) {
    nc  = find(S == h);
    ncp = find(S > h);
    V[h-1]  = R::rbeta(1.0 + nc.n_elem, alpha + ncp.n_elem);
  }
  V[H-1] = 1;
  return V;
}


// [[Rcpp::export]]
NumericVector S_update(arma::vec y, arma::vec strata, int p_RF, arma::vec nu, arma::vec theta, arma::vec eta_residual){
  
  // Initialization
  int H = nu.n_elem;
  
  NumericVector S(p_RF);               // Output
  IntegerVector clusters = Range(1,H); // Possible clusters
  NumericVector lprob(H);              // Log-probabilities
  NumericVector prob(H);               // Probabilities
  arma::uvec index;
  arma::vec y_i;
  double ltemp;
  
  for(int i=0; i < p_RF; i++) {    //Cycle the States
    index = find(strata == i);
    y_i = y(index);
    int n_i = index.n_elem;
    arma::vec eta_residual_i = eta_residual(index);
    
    for(int h =0; h < H; h++) {   // Cycle the possible intercept
      ltemp = 0;
      for(int j = 0; j < n_i; j++) { // Observations within the State
        double lin = theta(h) + eta_residual_i(j);
        double pi = 1/(1+ exp(-lin));
        ltemp = ltemp + R::dbinom(y_i(j),1,pi,TRUE);
        }
      lprob[h] = log(nu[h]) + ltemp;
      }
    lprob = lprob - max(lprob);
    prob  = exp(lprob);
    prob  = prob/sum(prob);
    S[i] =  RcppArmadillo::sample(clusters, 1, TRUE, prob)[0];
  }
  return S;
}


// [[Rcpp::export]]
arma::mat clusterdist(arma::mat cluster) {
  int R = cluster.n_rows;
  int n = cluster.n_cols;
  arma::mat out(n,n);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      out(i,j) = 0;
    }
  }
  
  for(int r = 0; r < R; r++) {
    for(int i = 0; i < n; i++) {
      for(int j = (i+1); j < n; j++) {
        if(cluster(r,i)!=cluster(r,j)){
          out(i,j) = out(i,j) + 1;
          out(j,i) = out(i,j);
        }
      }
    }
  }
  
  return out;  
}
