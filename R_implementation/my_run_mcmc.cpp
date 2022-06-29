#include <RcppArmadillo.h>
#include "my_parameter_updates.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat rankings,  int nmc,
                    Rcpp::Nullable<arma::vec> cardinalities,
                    Rcpp::Nullable<arma::vec> logz_estimate,
                    Rcpp::Nullable<arma::vec> rho_init,
                    std::string metric = "footrule",
                    int leap_size = 1,
                    double alpha_prop_sd = 0.5,
                    double alpha_init = 1,
                    int alpha_jump = 1,
                    double lambda = 0.1,
                    double alpha_max = 1e6,
                    int rho_thinning = 1,
                    bool verbose = false
                    ){

  // The number of items ranked
  int n_items = rankings.n_rows;

  // The number of assessors
  int n_assessors = rankings.n_cols;

  // Declare the matrix to hold the latent ranks
  arma::mat rho(n_items, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
  // Note : `static_cast<new_type>(target)` convert the type of target variable to the new type
  if(rho_init.isNotNull()){
      rho.col(0) = Rcpp::as<arma::vec>(rho_init) ;
  } else {
    rho.col(0) = arma::shuffle(arma::regspace(1, n_items)) ;
  }
  arma::vec rho_current = rho.col(0) ;
  

  // Declare the vector to hold the scaling parameter alpha
  arma::vec alpha(std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
  alpha(0) = alpha_init ;
  double alpha_current = alpha(0);
 

  // Declare indicator varaibles to hold acceptance or not
  int alpha_acceptance = 0 ;
  int rho_acceptance = 0 ;

  
  // Other variables used
  int alpha_index = 0, rho_index = 0 ;

  arma::uvec element_indices = arma::regspace<arma::uvec>(0, rankings.n_rows - 1);
  // Note : "element_indices" is same as 1:n in R where n is the number of items


  // This is the Metropolis-Hastings loop

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0,
  // and this has been done above
  for(int t = 1; t < nmc; ++t){

    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << t << " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }

    }

    
    update_rho(rho, rho_acceptance, rho_current, rho_index, 
                rho_thinning, alpha_current, leap_size,
                rankings,  metric, n_items, t) ;
    // ! `update_rho` is defined in my_parameter_updates.cpp
    
    // if(verbose){
    //   Rcpp::Rcout << "First " << t << " update of rho completed." << std::endl;
    // }
    

    if(t % alpha_jump == 0) {
      ++alpha_index;
      alpha(alpha_index) = update_alpha(alpha_acceptance, alpha_current ,
            rankings, rho_current, alpha_prop_sd, metric, lambda, cardinalities, logz_estimate, alpha_max);
      
    // if(verbose){
    //   Rcpp::Rcout << "First " << t << " update of alpha completed." << std::endl;
    // }
    // if(verbose){
    //   Rcpp::Rcout << "First " << t<< "accpeted alpha number is "<< alpha_acceptance << std::endl;
    // }
      // ! `update_alpha` is defined in my_parameterupdates.cpp
      // Update alpha_old
      alpha_current = alpha(alpha_index) ;
    }

  }



  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance,
    Rcpp::Named("n_items") = n_items,
    Rcpp::Named("n_assessors") = n_assessors
  );


}

