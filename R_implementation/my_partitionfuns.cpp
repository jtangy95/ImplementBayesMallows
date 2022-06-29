#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// Helper to compute the importance sampling smoothed fit
double compute_is_fit(double alpha, arma::vec fit){
  // Note : Here, input "fit" is a beta.values of polynomials estimating partition function log{Z_n(alpha)}
  // The partition function
  double logZ = 0;
  int n_items = fit.n_elem;

  for(int i = 0; i < n_items; ++i){
    logZ += std::pow(alpha, static_cast<double>(i)) * fit(i);
  }
  return(logZ);
}

// [[Rcpp::depends(RcppArmadillo)]]
int binomial_coefficient(int n, int k){

  // Special case:
  if( k > n ) return 0;

  int res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k )
    k = n - k;

  // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i)
  {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
double logz_cardinalities(const double& alpha, const int& n_items, const arma::vec& cardinalities, const std::string& metric){
  if(metric == "footrule"){
    arma::vec distances = arma::regspace(0, 2, std::floor(std::pow(static_cast<double>(n_items), 2.) / 2));
    return std::log(arma::sum(cardinalities % arma::exp(-alpha * distances / n_items)));
    // Note : `%` works as an element wise multiplication operator for Mat, Col, Row and Cube classes.
  } else if (metric == "spearman"){
    arma::vec distances = arma::regspace(0, 2, 2 * binomial_coefficient(n_items + 1, 3));
    // ! `binomial_coefficient` is defined in misc.cpp
    return std::log(arma::sum(cardinalities % arma::exp(-alpha * distances / n_items)));
  } else {
    Rcpp::stop("Cardinalities not implemented for the provided metric.");
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
double get_partition_function(int n_items, double alpha,
                              const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                              const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
                              std::string metric = "footrule"){


  if(cardinalities.isNotNull()){
    return logz_cardinalities(alpha, n_items, Rcpp::as<arma::vec>(cardinalities), metric);
  } else if(logz_estimate.isNotNull()) {
    return compute_is_fit(alpha, Rcpp::as<arma::vec>(logz_estimate));
    // Note : return value is estimatmed value of logZ_n(alpha) 
  } else {
      Rcpp::stop("Partition function cannot be calculated.") ;
  }
}