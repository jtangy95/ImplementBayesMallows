#ifndef MISSINGDATA_H
#define MISSINGDATA_H

#include <RcppArmadillo.h>
#include "my_distances.h"


arma::uvec arma_setdiff(arma::uvec x, arma::uvec y) ;

arma::vec propose_augmentation(const arma::vec& ranks, const arma::uvec& indicator) ;

void initialize_missing_ranks(arma::mat& rankings, const arma::umat& missing_indicator,
                              const arma::uvec& assessor_missing) ;

void update_missing_ranks(arma::mat& rankings, 
                          arma::vec& aug_acceptance,
                          const arma::umat& missing_indicator,
                          const arma::uvec& assessor_missing,
                          const double& alpha, const arma::vec& rho,
                          const std::string& metric) ;

#endif
