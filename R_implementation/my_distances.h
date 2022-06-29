#ifndef DISTANCES_H
#define DISTANCES_H

#include <RcppArmadillo.h>


double footrule_distance(const arma::vec& r1, const arma::vec& r2) ;

double spearman_distance(const arma::vec& r1, const arma::vec& r2) ;

double  get_rank_distance(arma::vec r1, arma::vec r2, std::string metric) ;

arma::vec rank_dist_vec(const arma::mat& rankings,
                        const arma::vec& rho,
                        const std::string& metric) ;

double rank_dist_sum(const arma::mat& rankings, const arma::vec& rho,
                     const std::string& metric) ;

#endif