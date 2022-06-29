#include <RcppArmadillo.h>
#include "my_leap_and_shift.h"
#include "my_distances.h"

// [[Rcpp::depends(RcppArmadillo)]]

void find_pairwise_limits(int& left_limit, int& right_limit, const int& item,
                          const Rcpp::List& assessor_constraints,
                          const arma::vec& current_ranking){

  // Find the items which are preferred to the given item
  // Items go from 1, ..., n_items, so must use [item - 1]
  arma::uvec items_above = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[item - 1]);
  arma::uvec items_below = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[item - 1]);

  // If there are any items above, we must find the possible rankings
  if(items_above.n_elem > 0){
    // Again subtracting 1 because of zero-first indexing
    // Find all the rankings of the items that are preferred to *item*
    arma::vec rankings_above = current_ranking.elem(items_above - 1);
    left_limit = arma::max(rankings_above);
  }

  // If there are any items below, we must find the possible rankings
  if(items_below.n_elem > 0){
    // Find all the rankings of the items that are disfavored to *item*
    arma::vec rankings_below = current_ranking.elem(items_below - 1);
    right_limit = arma::min(rankings_below);
  }

}

arma::vec propose_pairwise_augmentation(const arma::vec& ranking, const Rcpp::List& assessor_constraints){

  int n_items = ranking.n_elem;

  // Extract the constraints for this particular assessor
  arma::uvec constrained_items = Rcpp::as<arma::uvec>(assessor_constraints[0]);

  // Sample an integer between 1 and n_items
  int item = arma::randi<int>(arma::distr_param(1, n_items)) ;

  // Left and right limits of the interval we draw ranks from
  // Correspond to l_j and r_j, respectively, in Vitelli et al. (2018), JMLR, Sec. 4.2.
  int left_limit = 0, right_limit = n_items + 1;
  if(arma::any(constrained_items == item)){
    find_pairwise_limits(left_limit, right_limit, item, assessor_constraints, ranking) ;
  }

  // Now complete the leap step by sampling a new proposal uniformly between
  // left_limit + 1 and right_limit - 1
  int proposed_rank = arma::randi<int>(arma::distr_param(left_limit + 1, right_limit - 1)) ;

  // Assign the proposal to the (item-1)th item
  arma::vec proposal = ranking;
  proposal(item - 1) = proposed_rank;

  double delta_r;
  arma::uvec indices;

  // Do the shift step
  shift_step(proposal, ranking, item, delta_r, indices);

  return proposal;
}


void augment_pairwise(
    arma::mat& rankings,
    const arma::uvec& current_cluster_assignment,
    const arma::vec& alpha,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    arma::vec& aug_acceptance
){

  int n_assessors = rankings.n_cols;
  int n_items = rankings.n_rows;

  for(int i = 0; i < n_assessors; ++i){

    arma::vec proposal;

    // Sample a proposal
   
    proposal = propose_pairwise_augmentation(rankings.col(i), Rcpp::as<Rcpp::List>(constraints[i]));
    // Notes : constraints[i] is a list of constraints for i-th assessor
    // Notes : it consists of three lists : the first one is the constrained items
    // Notes : the second one is the list of "items above" and the third one is the list of "items below"
    // Notes : "items above" have n_item number of elements where each one is a vector. The same also holds for "items below" .
    
    // Finally, decide whether to accept the proposal or not
    // Draw a uniform random number
    double u = std::log(arma::randu<double>());

    // Find which cluster the assessor belongs to
    int cluster = current_cluster_assignment(i);

    double ratio = -alpha(cluster) / n_items *
      (get_rank_distance(proposal, rho.col(cluster), metric) -
      get_rank_distance(rankings.col(i), rho.col(cluster), metric));

    if(ratio > u){
      rankings.col(i) = proposal;
      ++aug_acceptance(i);
    }
  }

}

