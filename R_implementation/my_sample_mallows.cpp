#include "RcppArmadillo.h"
#include "my_leap_and_shift.h"
#include "my_distances.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rmallows(
    arma::vec rho0,
    double alpha0,
    int n_samples,
    int burnin,
    int thinning,
    int leap_size = 1,
    std::string metric = "footrule"
  ){

  // The number of items ranked
  int n_items = rho0.n_elem;

  // Declare the matrix to hold the sampled ranks
  arma::mat ranks(n_items, n_samples);

  // Other variables used
  int sample_index = 0;

  // Vector to hold the iteration value of rho
  // Initializing it to the modal ranking
  arma::vec rank_current = rho0;

  int t = 1;
  // This is the Metropolis-Hastings loop
  while(sample_index < n_samples){

    // Check if the user has tried to stop the running
    if (t % 1000 == 0) Rcpp::checkUserInterrupt();

    arma::vec rank_proposal;
    arma::uvec indices;
    double prob_backward, prob_forward;

    // Sample a proposal
    leap_and_shift(rank_proposal, indices, prob_backward, prob_forward,
                   rank_current, leap_size);

    // Compute the distances to current and proposed ranks
    double dist_new = get_rank_distance(rho0(indices), rank_proposal(indices), metric);
    double dist_old = get_rank_distance(rho0(indices), rank_current(indices), metric);

    // Metropolis-Hastings ratio
    double ratio = - alpha0 / n_items * (dist_new - dist_old) +
      std::log(prob_backward) - std::log(prob_forward);

    // Draw a uniform random number
    double u = std::log(arma::randu<double>());

    if(ratio > u){
      rank_current = rank_proposal;
    }

    // Save every thinning'th iteration after burnin
    if((t > burnin) & (t % thinning == 0)){
      ranks.col(sample_index) = rank_current;
      ++sample_index;
    }

    ++t;
  }
  return ranks;
}



