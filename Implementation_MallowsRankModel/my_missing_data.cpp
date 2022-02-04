#include <RcppArmadillo.h>
#include "my_distances.h"

// [[Rcpp::depends(RcppArmadillo)]]
arma::uvec arma_setdiff(arma::uvec x, arma::uvec y){

  x = arma::unique(x);
  y = arma::unique(y);
  // Note :  `unique(x)` returns the unique elements of x, sorted in ascending order.

  for (size_t j = 0; j < y.n_elem; j++) {
    // `size_t` is unsigned integer type that is big enough to represent the size of the largest possible object on the target platform.
    arma::uvec q1 = arma::find(x == y[j]);
    // Note : since we took 'unique' function on x and y , q1 is a vector with only one element.
    if (!q1.empty()) {
      x.shed_row(q1(0));
      // Note : `shed_row(row_number)` removes the specified row.
    }
  }
  return x;
}



// [[Rcpp::depends(RcppArmadillo)]]
arma::vec propose_augmentation(const arma::vec& ranks, const arma::uvec& indicator){
  arma::vec proposal = ranks;
  proposal(arma::find(indicator == 1)) = arma::shuffle(ranks(arma::find(indicator == 1)));
  // Note : this function is used in "update_missing_ranks" , which is used in the situation where vector ranks does not have any missing value although those originally missing values (NA) are replaced by augmented integer values
  return(proposal);
}

void initialize_missing_ranks(arma::mat& rankings, const arma::umat& missing_indicator,
                              const arma::uvec& assessor_missing) {

  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0) {
      continue;
    } else {
      arma::vec rank_vector = rankings.col(i);
      arma::uvec present_inds = arma::find(missing_indicator.col(i) == 0);
      arma::uvec missing_inds = arma::find(missing_indicator.col(i) == 1);
      // Find the available ranks and permute them
      arma::uvec new_ranks = arma::shuffle(arma_setdiff(
        // `arma_setdiff` defined above is an analogue of R's setdiff() function.
        arma::regspace<arma::uvec>(1, rank_vector.size()),
        arma::conv_to<arma::uvec>::from(rank_vector(present_inds))
        // Note : angle brackets < > is used to specify a type
      ));
      rank_vector(missing_inds) = arma::conv_to<arma::vec>::from(new_ranks) ;
      rankings.col(i) = rank_vector;
    }
  }
}

void update_missing_ranks(arma::mat& rankings, 
                          arma::vec& aug_acceptance,
                          const arma::umat& missing_indicator,
                          const arma::uvec& assessor_missing,
                          const double& alpha, const arma::vec& rho,
                          const std::string& metric){

  int n_items = rankings.n_rows;
  int n_assessors = rankings.n_cols;
  
  // Note : we shall augment rank of each assessor one by one
  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0){
      ++aug_acceptance(i);
      continue;
    }

    // Sample an augmentation proposal
    arma::vec proposal = propose_augmentation(rankings.col(i), missing_indicator.col(i));

    // Draw a uniform random number
    double u = std::log(arma::randu<double>());

    // Find which cluster the assessor belongs to

    double ratio = -alpha / n_items *
      (get_rank_distance(proposal, rho, metric) -
      get_rank_distance(rankings.col(i), rho, metric));

    if(ratio > u){
      rankings.col(i) = proposal;
      ++aug_acceptance(i);
    }
  }
}
