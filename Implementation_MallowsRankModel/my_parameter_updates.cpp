#include <RcppArmadillo.h>
#include "my_distances.h"
#include "my_leap_and_shift.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void update_rho(arma::mat& rho, int& rho_acceptance, arma::vec& rho_current,
                int& rho_index,  const int& rho_thinning,
                const double& alpha_current, const int& leap_size, const arma::mat& rankings, const std::string& metric, const int& n_items, const int& t) {

  // Sample a rank proposal
  arma::vec rho_proposal;
  arma::uvec indices;
  double prob_backward, prob_forward;

  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho_current, leap_size);
  // ! `leap_and_shift` is defined in leapandshift.cpp
  // Note : then we get new "rho_proposal". Also, "indices" is the index where rho_old and rho_proposal are different. "prob_backward" and "prob_forward" are probabilities associated to the transition.

  // Compute the distances to current and proposed ranks
  double dist_new = rank_dist_sum(rankings.rows(indices), rho_proposal(indices), metric);
  double dist_old = rank_dist_sum(rankings.rows(indices), rho_current(indices), metric);
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
    rho_current = rho_proposal ;
    ++rho_acceptance ;
    // Note : At first, rho_current was initially given. Now, rho_current has changed. 
    // Note : rho_accpetance holds the number of acceptance. Here, increment it by one.
  }

  // Save rho if appropriate
  if(t % rho_thinning == 0){
    ++rho_index;
    rho.col(rho_index) = rho_current ;
  }

}