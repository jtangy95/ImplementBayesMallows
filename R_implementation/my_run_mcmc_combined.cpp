#include <RcppArmadillo.h>
#include "my_parameter_updates_cluster.h"
#include "my_mixtures.h"
#include "my_missing_augment.h"
#include "my_pairwise_augment.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat rankings,  Rcpp::List constraints,
                    int nmc,
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
                    int n_clusters = 1,
                    double psi = 10,
                    bool include_wcd = false,
                    int clus_thin = 1,
                    int aug_thinning = 1,
                    bool save_aug = false,
                    bool verbose = false
                    ){

  // The number of items ranked
  int n_items = rankings.n_rows;

  // The number of assessors
  int n_assessors = rankings.n_cols;
  
  // Partial ranking ? -> any missing == true  /  Preference comparison ? -> augpair == true 
  bool augpair = (constraints.length() > 0);
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

  // Declare the cube to hold the latent ranks
  arma::cube rho(n_items, n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
  // Note : `static_cast<new_type>(target)` convert the type of target variable to the new type
  if(rho_init.isNotNull()){
      rho.slice(0) = arma::repmat(Rcpp::as<arma::mat>(rho_init), 1, n_clusters) ;
  } else {
    rho.slice(0) = arma::shuffle(arma::repmat(arma::regspace(1, n_items), 1, n_clusters)) ;
  }
  arma::mat rho_current = rho.slice(0) ;
  if(verbose){
        Rcpp::Rcout << "rho is initialized" << std::endl;
        Rcpp::Rcout << "size of rho is " << size(rho) << std::endl;
    }

  // Declare the matrix to hold the scaling parameter alpha
  arma::mat alpha(n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
  alpha.col(0).fill(alpha_init) ;
  arma::vec alpha_current = alpha.col(0) ;
  if(verbose){
        Rcpp::Rcout << "alpha is initialized" << std::endl;
    }

  // If the user wants to save augmented data, we need a cube
  arma::cube augmented_data;
  if(save_aug){
    augmented_data.set_size(n_items, n_assessors, std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
    augmented_data.slice(0) = rankings ;
  }

  if(verbose){
        Rcpp::Rcout << "augmented_data is initialized" << std::endl;
        Rcpp::Rcout << "size of augmented_data is " << size(augmented_data) << std::endl;
    }
 
  // Clustering
  bool clustering = n_clusters > 1;
  int n_cluster_assignments = clustering ? std::ceil(static_cast<double>(nmc * 1.0 / clus_thin)) : 1;
  if(verbose){
        Rcpp::Rcout << "number of montecarlo saved for clustering is " << n_cluster_assignments << std::endl;
    }

  arma::mat cluster_probs(n_clusters, n_cluster_assignments);
  cluster_probs.col(0).fill(1.0 / n_clusters);
  arma::vec current_cluster_probs = cluster_probs.col(0);
    if(verbose){
        Rcpp::Rcout << "cluster_probs is initialized" << std::endl;
    }
  
  arma::umat cluster_assignment(n_assessors, n_cluster_assignments);
  cluster_assignment.col(0) = arma::randi<arma::uvec>(n_assessors, arma::distr_param(0, n_clusters - 1));
  arma::uvec current_cluster_assignment = cluster_assignment.col(0);
   if(verbose){
        Rcpp::Rcout << "cluster_assignment is initialized" << std::endl;
    }
  // Note : `randi(n_elem, distr_param(a,b))` generates a vector with the elements set to random integer values in the [a,b] interval and the syntax is `vector_type v = randi<vector_type>(n_elem, distr_param(a,b))`
  // Note : cluster label is 0 , 1, ... , C - 1 instead of 1, 2, ... , C so that in the end we should add 1 elementwisely for the output
  
  // Matrix with precomputed distances d(R_j, \rho_j), used to avoid looping during cluster assignment`
  arma::mat dist_mat(n_assessors, n_clusters);
  update_dist_mat(dist_mat, rankings, rho_current, metric);
  if(verbose){
        Rcpp::Rcout << "dist_mat is initialized" << std::endl;
    }
  // ! `update_dist_mat` is defined in my_mixtures.cpp
  arma::mat within_cluster_distance(n_clusters, include_wcd ? nmc : 1);
  within_cluster_distance.col(0) = update_wcd(current_cluster_assignment, dist_mat);
   if(verbose){
        Rcpp::Rcout << "within_cluster_distance is initialized" << std::endl;
    }
  // ! `update_wcd` is defined in my_mixtures.cpp

  // Declare indicator varaibles to hold acceptance or not
  arma::uvec alpha_acceptance = arma::ones<arma::uvec>(n_clusters) ;
  arma::uvec rho_acceptance = arma::ones<arma::uvec>(n_clusters) ;

  arma::vec aug_acceptance ;
  if(any_missing | augpair){
      aug_acceptance = arma::ones(n_assessors) ;
  } else{
      aug_acceptance.reset() ;
  }

  
  // Other variables used
  int alpha_index = 0, rho_index = 0 , cluster_assignment_index = 0 , aug_index = 0 ;

  arma::uvec element_indices = arma::regspace<arma::uvec>(0, rankings.n_rows - 1) ;
  // Note : "element_indices" is same as 1:n in R where n is the number of items
  
  
  // This is the Metropolis-Hastings loop
  if(verbose){
        Rcpp::Rcout << "Metropolis_Hastings loop is started" << std::endl;
    }
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

    for(int i = 0 ; i < n_clusters ; ++i){
        update_rho(rho, rho_acceptance, rho_current, rho_index, 
                    i, rho_thinning, alpha_current(i), leap_size,
                    clustering ? arma::conv_to<arma::mat>::from(rankings.cols(arma::find(current_cluster_assignment == i))) : rankings , metric, n_items, t) ;
        // ! `update_rho` is defined in my_parameter_updates.cpp
    }

    
    if(t % alpha_jump == 0) {
      ++alpha_index;
      for(int i = 0 ; i < n_clusters ; ++i){
        alpha(i, alpha_index) = update_alpha(alpha_acceptance, alpha_current(i) ,
        clustering ? arma::conv_to<arma::mat>::from(rankings.cols(arma::find(current_cluster_assignment == i))) : rankings , i, rho_current.col(i) , alpha_prop_sd, metric, lambda, cardinalities, logz_estimate, alpha_max);
      }
      // ! `update_alpha` is defined in my_parameterupdates.cpp
      // Update alpha_current
      alpha_current = alpha.col(alpha_index) ;
    }

     // Note : I think that before updating cluster labels, we should update distance matrix since current rho matrix is updated
    if(clustering | include_wcd){
        update_dist_mat(dist_mat, rankings, rho_current, metric);
    }
    

    if(clustering){
        current_cluster_probs = update_cluster_probs(current_cluster_assignment, n_clusters, psi);
        // ! `update_cluster_probs" is defined in my_mixtures.cpp

        current_cluster_assignment = update_cluster_labels(dist_mat, current_cluster_probs,
                                                        alpha_current, n_items, t, metric, cardinalities,
                                                        logz_estimate);

        if(t % clus_thin == 0){
        ++cluster_assignment_index;
        cluster_assignment.col(cluster_assignment_index) = current_cluster_assignment;
        cluster_probs.col(cluster_assignment_index) = current_cluster_probs;
        }
    }

    if(include_wcd){
        // Update within_cluster_distance
        within_cluster_distance.col(t) = update_wcd(current_cluster_assignment, dist_mat);
    }

    // Perform data augmentation of missing ranks, if needed
    if(any_missing){
        update_missing_ranks(rankings, current_cluster_assignment, aug_acceptance, missing_indicator,
                            assessor_missing, alpha_current, rho_current, metric);
    }

    // Perform data augmentation of pairwise comparisons, if needed
    if(augpair){
        augment_pairwise(rankings, current_cluster_assignment, alpha_current, rho_current,
                        metric, constraints, aug_acceptance);

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
    Rcpp::Named("cluster_assignment") = cluster_assignment + 1,
    Rcpp::Named("cluster_probs") = cluster_probs ,
    Rcpp::Named("within_cluster_distance") = within_cluster_distance,
    Rcpp::Named("augmented_data") = augmented_data,
    Rcpp::Named("partial_rank") = any_missing,
    Rcpp::Named("preference_compare") = augpair,
    Rcpp::Named("aug_acceptance") = aug_acceptance
  );


}

