#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
void shift_step(arma::vec& rho_proposal, const arma::vec& rho_current,
                const int& u, double& delta_r, arma::uvec& indices){
  // Shift step:
  delta_r = rho_proposal(u - 1) - rho_current(u - 1);
  indices = arma::zeros<arma::uvec>(std::abs(delta_r) + 1);
  indices[0] = u-1;
  int index;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      index = arma::as_scalar(arma::find(rho_current == rho_current(u-1) + k));
      rho_proposal(index) -= 1;
      indices[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k = (-1); k >= delta_r; --k){
      index = arma::as_scalar(arma::find(rho_current == rho_current(u-1) + k));
      rho_proposal(index) += 1;
      indices[-(k)] = index;
    }
  }
  // Note : From the `shift_step`, "rho_proposal" changes from a rho_star to rho_star.
  // Note : "indices" were a NULL vector initially set in defining "update_rho" function in parameterupdates.cpp and after the `shift_step`, "indices" becomes a set of index among {1,2,...,n} on which "rho_old" and "rho_proposal" are distinct.
}

// [[Rcpp::depends(RcppArmadillo)]]
void leap_and_shift(arma::vec& rho_proposal, arma::uvec& indices,
                    double& prob_backward, double& prob_forward,
                    const arma::vec& rho_current, int leap_size){

  // Set proposal equal to current
  rho_proposal = rho_current;

  // Help vectors
  arma::vec support;

  // Number of items
  int n = rho_current.n_elem;

  // Other helper variables
  int u, index ;
  double delta_r, support_new_len ;

  // Leap 1
  // 1, sample u randomly between 1 and n
  u = arma::randi(arma::distr_param(1, n));

  // 2, compute the set S for sampling the new rank
  double dobL = static_cast<double>(leap_size);
  double dobn = static_cast<double>(n);

  // Defining linspace lengths here to avoid duplication in code
  double length1 = std::min(rho_current(u - 1) - 1, dobL);
  double length2 = std::min(n - rho_current(u - 1), dobL);

  if((rho_current(u - 1) > 1) & (rho_current(u - 1) < n)){
    support = arma::join_cols(arma::linspace(
      std::max(1.0, rho_current(u - 1) - leap_size), rho_current(u - 1) - 1, length1) ,
      arma::linspace(rho_current(u - 1) + 1, std::min(dobn, rho_current(u - 1) + leap_size), length2));
      // Note : `join_cols(A,B)` performs a vertical concatenation
  } else if(rho_current(u - 1) == 1){
    support = arma::linspace(rho_current(u - 1) + 1, std::min(dobn, rho_current(u - 1) + leap_size), length2);
  } else if(rho_current(u - 1) == n){
    support = arma::linspace(std::max(1.0, rho_current(u - 1) - leap_size), rho_current(u - 1) - 1, length1);
  }
  // Note : Now "support" becomes the set S in the paper. Briefly speaking, |S|= length1+length2 (where length1 = 0 if rho_u = 1 and length2 = 0 if rho_u = n)

  // 3. assign a random element of the support set, this completes the leap step
  index = arma::randi(arma::distr_param(0, support.n_elem-1));
  // Picked element index-1 from the support set
  rho_proposal(u-1) = support(index);
  // Note : Now "rho_proposal" becomes rho_star in the paper.

  // Compute the associated transition probabilities
  if(std::abs(rho_proposal(u - 1) - rho_current(u - 1)) == 1){
    // in this case the transition probabilities coincide! (and in fact for leap_size = 1 the L&S is symmetric)
    support_new_len = std::min(rho_proposal(u - 1) - 1, dobL) + std::min(n - rho_proposal(u - 1), dobL);
    // Note : This is similar to length1+length2 for S above
    prob_forward = 1.0 / (n * support.n_elem) + 1.0 / (n * support_new_len);
    prob_backward = prob_forward;
  } else {
    // P(proposed|current)
    prob_forward = 1.0 / (n * support.n_elem);
    // P(current|proposed)
    support_new_len = std::min(rho_proposal(u - 1) - 1, dobL) + std::min(n - rho_proposal(u-1), dobL);
    prob_backward = 1.0 / (n * support_new_len);
  }

  shift_step(rho_proposal, rho_current, u, delta_r, indices);

}



