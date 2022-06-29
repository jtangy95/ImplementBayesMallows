#ifndef PARTITIONFUNS_H
#define PARTITIONFUNS_H

#include <RcppArmadillo.h>

double compute_is_fit(double alpha, arma::vec fit) ;

int binomial_coefficient(int n, int k) ;

double logz_cardinalities(const double& alpha, const int& n_items, const arma::vec& cardinalities, const std::string& metric) ;


double get_partition_function(int n_items, double alpha,
                              const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                              const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
                              std::string metric = "footrule") ;

#endif

