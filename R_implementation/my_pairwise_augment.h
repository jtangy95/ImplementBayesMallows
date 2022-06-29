#ifndef PAIRWISE_AUGMENT_H
#define PAIRWISE_AUGMENT_H

#include <RcppArmadillo.h>
#include "my_leap_and_shift.h"
#include "my_distances.h"

void find_pairwise_limits(int& left_limit, int& right_limit, const int& item,
                          const Rcpp::List& assessor_constraints,
                          const arma::vec& current_ranking) ;

arma::vec propose_pairwise_augmentation(const arma::vec& ranking, const Rcpp::List& assessor_constraints) ;

void augment_pairwise(
    arma::mat& rankings,
    const arma::uvec& current_cluster_assignment,
    const arma::vec& alpha,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    arma::vec& aug_acceptance
) ;

#endif