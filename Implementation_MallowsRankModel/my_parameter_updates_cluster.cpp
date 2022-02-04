#include <RcppArmadillo.h>
#include "my_distances.h"
#include "my_leap_and_shift.h"
#include "my_partitionfuns.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void update_rho(arma::cube& rho, arma::uvec& rho_acceptance, arma::mat& rho_current,
                int& rho_index,  const int& cluster_index , const int& rho_thinning,
                const double& alpha_current, const int& leap_size, const arma::mat& rankings, const std::string& metric, 
                const int& n_items, const int& t) {
  
  arma::vec rho_cluster = rho_current.col(cluster_index) ;

  // Sample a rank proposal
  arma::vec rho_proposal;
  arma::uvec indices;
  double prob_backward, prob_forward;

  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho_cluster, leap_size);
  // ! `leap_and_shift` is defined in leapandshift.cpp
  // Note : then we get new "rho_proposal". Also, "indices" is the index where rho_old and rho_proposal are different. "prob_backward" and "prob_forward" are probabilities associated to the transition.

  // Compute the distances to current and proposed ranks
  double dist_new = rank_dist_sum(rankings.rows(indices), rho_proposal(indices), metric);
  double dist_old = rank_dist_sum(rankings.rows(indices), rho_cluster(indices), metric);
  // ! `rank_dist_sum` is defined in distances.cpp
  // `X.rows( vector_of_row_indices)` is a way of subsetting matrix
  // "dist_new" yields sum_{j=1}^N d(R_j, rho_proposal), where distance is caculated only on the indices where rho_proposal and rho_old are different.

  // Metropolis-Hastings ratio
  double ratio = - alpha_current / n_items * (dist_new - dist_old) +
    std::log(prob_backward) - std::log(prob_forward);
  // Note : Here a "ratio" is equal to log(r) in the paper.

  // Draw a uniform random number
  double u = std::log(arma::randu<double>());

  if(u < ratio){
    rho_current.col(cluster_index) = rho_proposal ;
    ++rho_acceptance(cluster_index) ;
    // Note : At first, rho_current was initially given. Now, rho_current has changed. 
    // Note : rho_accpetance holds the number of acceptance. Here, increment it by one.
  }

  // Save rho if appropriate
  if(t % rho_thinning == 0){
    if(cluster_index ==0) {
        ++rho_index;
    }
    rho.slice(rho_index).col(cluster_index) = rho_current.col(cluster_index) ;
  }

}

double update_alpha(arma::uvec& alpha_acceptance,
                  const double& alpha_current,
                  const arma::mat& rankings,
                  const int& cluster_index, 
                  const arma::vec& rho_current,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
                  double alpha_max = 1e6) {
  // Note : R_NilValue is equivalent to NULL value in R.

  // Set the number of assessors. Here we want the number of assessors in this cluster
  
  int n_items = rho_current.n_elem;
  int n_assessors = rankings.n_cols;

  double alpha_proposal = std::exp(arma::randn<double>() * alpha_prop_sd +
                              std::log(alpha_current));
   // Note : `randn` is a random number generated from standard normal distrn. By using simga_alpha and log(alpha_old), we can consider it as a random number generated from normal with parameter log(alpha_old) and sigma_alpha.
   // Note : Y=log(X)~N(u, s) --> X=exp(Y) where Y~N(u, s)

  double rank_dist = rank_dist_sum(rankings, rho_current, metric);


  // Difference between current and proposed alpha
  double alpha_diff = alpha_current - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    lambda * alpha_diff +
    n_assessors * (
        get_partition_function(n_items, alpha_current, cardinalities, logz_estimate, metric) -
          get_partition_function(n_items, alpha_proposal, cardinalities, logz_estimate, metric)
    ) + std::log(alpha_proposal) - std::log(alpha_current);
  // Note : Here `sum(obs_freq)` plays a rold of the number of assessors N. But I think it is better to replace it with n_assessors=ranking.n_cols because sub(obs_freq) is the total N while ranking.n_cols here is the number of assessors in this cluster.
  
  // Draw a uniform random number
  double u = std::log(arma::randu<double>());

  if(ratio > u && alpha_proposal < alpha_max){
    ++alpha_acceptance(cluster_index) ;
    return alpha_proposal;
  } else {
    return alpha_current;
  }
}