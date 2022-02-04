#ifndef PARAMETERUPDATES_CLUSTER_H
#define PARAMETERUPDATES_CLUSTER_H

#include <RcppArmadillo.h>
#include "my_distances.h"
#include "my_leap_and_shift.h"
#include "my_partitionfuns.h"

void update_rho(arma::cube& rho, arma::uvec& rho_acceptance, arma::mat& rho_current,
                int& rho_index,  const int& cluster_index , const int& rho_thinning,
                const double& alpha_current, const int& leap_size, const arma::mat& rankings, const std::string& metric, 
                const int& n_items, const int& t) ;

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
                  double alpha_max = 1e6) ;


#endif