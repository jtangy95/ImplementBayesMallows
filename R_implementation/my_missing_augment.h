#ifndef MISSING_AUGMENT_H
#define MISSING_AUGMENT_H

#include <RcppArmadillo.h>
#include "my_distances.h"

arma::uvec arma_setdiff(arma::uvec x, arma::uvec y) ;

arma::vec propose_augmentation(const arma::vec& ranks, const arma::uvec& indicator) ;

void initialize_missing_ranks(arma::mat& rankings, const arma::umat& missing_indicator,
                              const arma::uvec& assessor_missing) ;

void update_missing_ranks(arma::mat& rankings, const arma::uvec& current_cluster_assignment,
                          arma::vec& aug_acceptance,
                          const arma::umat& missing_indicator,
                          const arma::uvec& assessor_missing,
                          const arma::vec& alpha, const arma::mat& rho,
                          const std::string& metric) ;

#endif
