#ifndef LEAPANDSHIFT_H
#define LEAPANDSHIFT_H

#include <RcppArmadillo.h>


void shift_step(arma::vec& rho_proposal, const arma::vec& rho_current,
                const int& u, double& delta_r, arma::uvec& indices) ;

void leap_and_shift(arma::vec& rho_proposal, arma::uvec& indices,
                    double& prob_backward, double& prob_forward,
                    const arma::vec& rho_current, int leap_size) ;

#endif