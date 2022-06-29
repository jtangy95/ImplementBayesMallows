#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
double footrule_distance(const arma::vec& r1, const arma::vec& r2){
  return arma::norm(r1 - r2, 1);
  // Note : `norm(X,p)` computes the p-norm of X where X can be a vector or matrix.
}

// [[Rcpp::depends(RcppArmadillo)]]
double spearman_distance(const arma::vec& r1, const arma::vec& r2){
    double a = arma::norm(r1 - r2, 2) ;
  return a * a ;
  
}

// [[Rcpp::depends(RcppArmadillo)]]
double  get_rank_distance(arma::vec r1, arma::vec r2, std::string metric){

  if (r1.n_elem != r2.n_elem){
    Rcpp::stop("r1 and r2 must have the same length");
  }

  if (metric == "footrule") {
    return footrule_distance(r1, r2);
  } else if (metric == "spearman") {
    return spearman_distance(r1, r2);
  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::vec rank_dist_vec(const arma::mat& rankings,
                        const arma::vec& rho,
                        const std::string& metric){

  int n = rankings.n_cols;
  // "n" is now the number of assessors
  arma::vec result = arma::zeros(n);

  for(int i = 0; i < n; ++i){
    result(i) = get_rank_distance(rankings.col(i), rho, metric);
  }
  // Note : `get_rank_distance` is defined above

  return(result);
}

// [[Rcpp::depends(RcppArmadillo)]]
double rank_dist_sum(const arma::mat& rankings, const arma::vec& rho,
                     const std::string& metric){
  return arma::sum(rank_dist_vec(rankings, rho, metric));
}

