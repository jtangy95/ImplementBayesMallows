#ifndef MIXTURES_H
#define MIXTURES_H

#include "RcppArmadillo.h"
#include "my_partitionfuns.h"
#include "my_distances.h"

void update_dist_mat(arma::mat& dist_mat, const arma::mat& rankings,
                     const arma::mat& rho_current, const std::string& metric
                     ) ;

arma::vec update_cluster_probs(
    const arma::uvec& current_cluster_assignment,
    const int& n_clusters,
    const int& psi
) ;

arma::uvec update_cluster_labels(
    const arma::mat& dist_mat,
    const arma::vec& cluster_probs,
    const arma::vec& alpha_current,
    const int& n_items,
    const int& t,
    const std::string& metric,
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue
) ;



arma::vec update_wcd(const arma::uvec& current_cluster_assignment,
                     const arma::mat& dist_mat) ;

#endif