#include "RcppArmadillo.h"
#include "my_partitionfuns.h"
#include "my_distances.h"


// [[Rcpp::depends(RcppArmadillo)]]


void update_dist_mat(arma::mat& dist_mat, const arma::mat& rankings,
                     const arma::mat& rho_current, const std::string& metric
                     ){
  int n_clusters = dist_mat.n_cols;
  for(int i = 0; i < n_clusters; ++i)
    dist_mat.col(i) = rank_dist_vec(rankings, rho_current.col(i), metric);
    // ! `rank_dist_vec` is defined in distances.cpp
}

arma::vec update_cluster_probs(
    const arma::uvec& current_cluster_assignment,
    const int& n_clusters,
    const int& psi
){

  arma::vec cluster_probs(n_clusters);

  for(int i = 0; i < n_clusters; ++i){
    // Find the parameter for this cluster and provide it to the gamma distribution
    cluster_probs(i) = R::rgamma(arma::sum(current_cluster_assignment == i) + psi, 1.0);
  }
  // Finally, normalize cluster_probs with 1-norm.
  // result now comes from Dirichlet(tau_k(0), ..., tau_k(n_clusters))
  return arma::normalise(cluster_probs, 1);

}

arma::uvec update_cluster_labels(
    const arma::mat& dist_mat,
    const arma::vec& cluster_probs,
    const arma::vec& alpha_current,
    const int& n_items,
    const int& t,
    const std::string& metric,
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue
){
  int n_assessors = dist_mat.n_rows;
  int n_clusters = dist_mat.n_cols;
  arma::uvec new_cluster_assignment(n_assessors);


  arma::mat assignment_prob(n_assessors, n_clusters);
  for(int i = 0; i < n_clusters; ++i){
    // Compute the logarithm of the unnormalized probability
    assignment_prob.col(i) = std::log(cluster_probs(i)) -
      alpha_current(i) / n_items * dist_mat.col(i) -
      get_partition_function(n_items, alpha_current(i), cardinalities, logz_estimate, metric);
      // Note : For a scalar value k and a vector v , k + v is calculated as k * ones(v.n_elem) + v .
      // Note : Here, dist_mat.col(i) is a vector whereas the others are scalar values.
  }

  for(int i = 0; i < n_assessors; ++i){
    // Exponentiate to get unnormalized prob relative to max
    arma::rowvec probs = arma::exp(assignment_prob.row(i) -
      arma::max(assignment_prob.row(i)));

    // Normalize with 1-norm
    probs = arma::normalise(probs, 1);
    double u = arma::randu();
    arma::uvec sample_clus = arma::find(arma::cumsum(probs) > u , 1, "first") ;

    new_cluster_assignment(i) = arma::as_scalar(sample_clus) ;
  }
  
  return(new_cluster_assignment);
}



arma::vec update_wcd(const arma::uvec& current_cluster_assignment,
                     const arma::mat& dist_mat){
  int n_clusters = dist_mat.n_cols;
  arma::vec wcd(n_clusters);

  arma::uvec inds = arma::regspace<arma::uvec>(0, n_clusters - 1);
  for(int i = 0; i < n_clusters; ++i){
    arma::mat dist_vec = arma::conv_to<arma::mat>::from(dist_mat.rows(arma::find(current_cluster_assignment == i))).col(i) ;
    wcd(i) = arma::sum(arma::conv_to<arma::vec>::from(dist_vec));
    // Note : "wcd" is C-dim vector representing within cluster distnace for each cluster
  }

  return(wcd);
}
