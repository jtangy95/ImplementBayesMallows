#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat rankings,  int nmc,
                    Rcpp::Nullable<arma::vec> cardinalities,
                    Rcpp::Nullable<arma::vec> logz_estimate,
                    Rcpp::Nullable<arma::vec> rho_init,
                    std::string metric = "footrule",
                    int leap_size = 1,
                    double alpha_prop_sd = 0.5,
                    double alpha_init = 5,
                    int alpha_jump = 1,
                    double lambda = 0.1,
                    double alpha_max = 1e6,
                    int rho_thinning = 1,
                    bool verbose = false,
                    ){

  // The number of items ranked
  int n_items = rankings.n_rows;

  // The number of assessors
  int n_assessors = rankings.n_cols;

  // Declare the matrix to hold the latent ranks
  arma::mat rho(n_items, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
  // Note : `static_cast<new_type>(target)` convert the type of target variable to the new type
  if(rho_init.isNotNull()){
      rho.col(0) = rho_init ;
  } else {
    rho.col(0) = arma::shuffle(arma::regspace(1, n_items)) ;
  }
  arma::vec rho_current = rho.col(0) ;
  

  // Declare the vector to hold the scaling parameter alpha
  arma::vec alpha(std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
  alpha(0) = alpha_init ;
  double alpha_current = alpha(0);
 

  // Declare indicator vectors to hold acceptance or not
  int alpha_acceptance = 0 ;
  int rho_acceptance = 0 ;

  
  // Other variables used
  int alpha_index = 0, rho_index = 0

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
    // ! `update_rho` is defined in parameterupdates.cpp
    

    if(t % alpha_jump == 0) {
      ++alpha_index;
      for(int i = 0; i < n_clusters; ++i){
        alpha(i, alpha_index) = update_alpha(alpha_acceptance, alpha_old(i),
              clustering ? rankings.submat(element_indices, arma::find(current_cluster_assignment == i)) : rankings,
              obs_freq, i, rho_old.col(i), alpha_prop_sd, metric, lambda, cardinalities, logz_estimate, alpha_max);
      }
      // ! `update_alpha` is defined in parameterupdates.cpp
      // Update alpha_old
      alpha_old = alpha.col(alpha_index);
    }

  if(clustering){
    current_cluster_probs = update_cluster_probs(current_cluster_assignment, n_clusters, psi);

    current_cluster_assignment = update_cluster_labels(dist_mat, current_cluster_probs,
                                                       alpha_old, n_items, t, metric, cardinalities,
                                                       logz_estimate, save_ind_clus);

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
                         assessor_missing, alpha_old, rho_old, metric);
  }

  // Perform data augmentation of pairwise comparisons, if needed
  if(augpair){
    augment_pairwise(rankings, current_cluster_assignment, alpha_old, 0.1, rho_old,
                     metric, constraints, aug_acceptance, clustering, error_model, Lswap);

  }

  // Save augmented data if the user wants this. Uses the same index as rho.
  if(save_aug & (t % aug_thinning == 0)){
    ++aug_index;
    augmented_data.slice(aug_index) = rankings;
  }

  if(clustering | include_wcd){
    update_dist_mat(dist_mat, rankings, rho_old, metric, obs_freq);
    }
  }



  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance / nmc,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance / nmc,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("shape1") = shape_1,
    Rcpp::Named("shape2") = shape_2,
    Rcpp::Named("cluster_assignment") = cluster_assignment + 1,
    Rcpp::Named("cluster_probs") = cluster_probs,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("within_cluster_distance") = within_cluster_distance,
    Rcpp::Named("augmented_data") = augmented_data,
    Rcpp::Named("any_missing") = any_missing,
    Rcpp::Named("augpair") = augpair,
    Rcpp::Named("aug_acceptance") = aug_acceptance / nmc,
    Rcpp::Named("n_assessors") = n_assessors,
    Rcpp::Named("obs_freq") = obs_freq
  );


}

