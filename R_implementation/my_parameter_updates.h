#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include <RcppArmadillo.h>
#include "my_distances.h"
#include "my_leap_and_shift.h"
#include "my_partitionfuns.h"


void update_rho(arma::mat& rho, int& rho_acceptance, arma::vec& rho_current,
                int& rho_index,  const int& rho_thinning,
                const double& alpha_current, const int& leap_size, const arma::mat& rankings, const std::string& metric, const int& n_items, const int& t) ;

double update_alpha(int& alpha_acceptance,
                  const double& alpha_current,
                  const arma::mat& rankings,
                  const arma::vec& rho_current,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
                  double alpha_max = 1e6) ;

#endif