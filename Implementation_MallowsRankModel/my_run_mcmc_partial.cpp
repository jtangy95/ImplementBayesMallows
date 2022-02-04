#include <RcppArmadillo.h>
#include "my_parameter_updates.h"
#include "my_missing_data.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List run_mcmc_partial(arma::mat rankings,  int nmc,
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
                    int aug_thinning = 1,
                    bool save_aug = false,
                    bool verbose = false
                    ){

  // The number of items ranked
  int n_items = rankings.n_rows;

  // The number of assessors
  int n_assessors = rankings.n_cols;

  bool any_missing = !arma::is_finite(rankings); 
  // Note : `is_finite` checks whether all elements are finite.

  arma::umat missing_indicator;
  arma::uvec assessor_missing;

  if(any_missing){
    // Converting to umat will convert NA to 0, but might cause clang-UBSAN error, so converting explicitly.
    rankings.replace(arma::datum::nan, 0);
    // Note : `.replace(old_value, new_value)`    `nan` stands for not a number
    missing_indicator = arma::conv_to<arma::umat>::from(rankings);
    // Note : `conv_to<type>::from(X)`  converts from one matrix(cube) type of X to another. In this case, from `mat` to `umat`
    missing_indicator.transform( [](int val) { return (val == 0) ? 1 : 0; } );
    // Note : `.transform(lambda_function)` transforms each element using lambda function [](){}
    // Note : `?` is conditional operator -> `condition ? expression1 : expression 2` works like `ifelse` in R
    // Note : "missing indicator" is a matrix with same dimension with rankings where 1 is marked on NA and 0 is marked otherwise.
    assessor_missing = arma::conv_to<arma::uvec>::from(arma::sum(missing_indicator, 0));
    // Note : For matrix M, `sum(M, dim)` returns the sums of elements in each column (dim=0) or each row (dim=1)
    // Note : "assessor_missing" becomes N dim vector representing number of not answered item for each assessor
    initialize_missing_ranks(rankings, missing_indicator, assessor_missing);
    // ! `initialize_missing_ranks` is defined in missing_data.cpp 
    // Note : "rankings" has just become augmented
  } else {
    missing_indicator.reset();
    assessor_missing.reset();
    // Note : `.reset()` reset the size to zero so that the object will have no elements.
  }

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

  // If the user wants to save augmented data, we need a cube
  arma::cube augmented_data;
  if(save_aug){
    augmented_data.set_size(n_items, n_assessors, std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
    augmented_data.slice(0) = rankings;
  }
 

  // Declare indicator varaibles to hold acceptance or not
  int alpha_acceptance = 0 ;
  int rho_acceptance = 0 ;

  arma::vec aug_acceptance ;
  if(any_missing){
      aug_acceptance = arma::ones(n_assessors) ;
  } else{
      aug_acceptance.reset() ;
  }

  
  // Other variables used
  int alpha_index = 0, rho_index = 0 ;
  int aug_index = 0 ; 

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
    
    if(t % alpha_jump == 0) {
      ++alpha_index;
      alpha(alpha_index) = update_alpha(alpha_acceptance, alpha_current ,
            rankings, rho_current, alpha_prop_sd, metric, lambda, cardinalities, logz_estimate, alpha_max);
      
      // ! `update_alpha` is defined in my_parameterupdates.cpp
      // Update alpha_old
      alpha_current = alpha(alpha_index) ;
    }

    // Perform data augmentation of missing ranks, if needed
    if(any_missing){
        update_missing_ranks(rankings, aug_acceptance, missing_indicator,
                            assessor_missing, alpha_current, rho_current, metric);
    }
    // Save augmented data if the user wants this. Uses the same index as rho.
    if(save_aug & (t % aug_thinning == 0)){
        ++aug_index;
        augmented_data.slice(aug_index) = rankings;
    }

  }

  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance,
    Rcpp::Named("n_items") = n_items,
    Rcpp::Named("n_assessors") = n_assessors,
    Rcpp::Named("augmented_data") = augmented_data,
    Rcpp::Named("any_missing") = any_missing,
    Rcpp::Named("aug_acceptance") = aug_acceptance
  );


}

