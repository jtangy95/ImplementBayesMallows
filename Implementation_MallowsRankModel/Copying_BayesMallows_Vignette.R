# install.packages('BayesMallows')
library('BayesMallows')
data("potato_visual")
potato_visual
bmm_test<-compute_mallows(potato_visual)
bmm_test
attributes(bmm_test)
print(bmm_test$rho)
print(bmm_test$alpha)

assess_convergence(bmm_test)
assess_convergence(bmm_test, parameter='rho', items=1:5)

bmm_potato<-compute_mallows(potato_visual, nmc=501000, verbose = T)
attributes(bmm_potato)
bmm_potato$burnin<-1000
attributes(bmm_potato)
plot(bmm_potato)
compute_posterior_intervals(bmm_potato, decimals=2)
compute_consensus(bmm_potato, type = "CP")

compuplot(bmm_potato)
plot(bmm_potato, parameter='rho')
plot(bmm_potato, parameter='rho', items=1:20)

bmm_potato2<-compute_mallows(potato_visual, nmc=501000, alpha_jump = 10, verbose=T)

#install.packages('tictoc')
library(tictoc)
tic('mcmc')
bmm_potato<-compute_mallows(potato_visual, nmc=501000)
toc()
tic('mcmc2')
bmm_potato2<-compute_mallows(potato_visual, nmc=501000, alpha_jump = 10)
toc()

bmm_potato3<-compute_mallows(potato_visual, nmc=501000, metric='spearman')



data('beach_preferences')
head(beach_preferences)
str(beach_preferences)
beach_tc<-generate_transitive_closure(beach_preferences)
str(beach_tc)
head(beach_preferences,20)
head(beach_tc, 20)

beach_init_rank<-generate_initial_ranking(beach_tc)
colnames(beach_init_rank)<-paste('Beach', 1:ncol(beach_init_rank))
beach_init_rank

library(dplyr)
filter(beach_preferences, assessor==1)
filter(beach_preferences, assessor==1, bottom_item==2 | top_item==2)
filter(beach_tc, assessor==1, bottom_item==2 | top_item==2)

bmm_test<-compute_mallows(preferences=beach_tc, ranking=beach_init_rank, save_aug=TRUE)
assess_convergence(bmm_test)
assess_convergence(bmm_test, parameter='rho', items=c(1,3,5,7,9))
assess_convergence(bmm_test, parameter='Rtilde', items=c(1,3,5,7,9), assessors=c(1,2))

assess_convergence(bmm_test, parameter='Rtilde', items=c(2,6,15), assessors=1)
assess_convergence(bmm_test, parameter='Rtilde', items=c(1,15), assessors=2)

tic('mcmc')
bmm_beaches<-compute_mallows(preference=beach_tc, ranking=beach_init_rank, nmc=102000, save_aug=T, verbose=T)
toc()

attributes(bmm_beaches)
bmm_beaches$augmented_data
bmm_beaches$burnin<-2000

plot(bmm_beaches)
plot(bmm_beaches, parameter='rho', items=1:15)
compute_posterior_intervals(bmm_beaches)
compute_posterior_intervals(bmm_beaches, parameter='rho')

compute_consensus(bmm_beaches, type='CP')
compute_consensus(bmm_beaches, type='MAP')
compute_consensus(bmm_beaches, parameter = "Rtilde", type='CP', assessors = 2)
compute_consensus(bmm_beaches, parameter = "Rtilde", type='CP', assessors = c(3,5))


plot_top_k(bmm_beaches, k=3) # k=3 is default
predict_top_k(bmm_beaches)

predict_top_k(bmm_beaches) %>%
  filter(prob > 0.9, assessor %in% 1:10)



data('sushi_rankings')
head(sushi_rankings)
str(sushi_rankings)

library('parallel')
(nCores=detectCores()-1)
cl<-makeCluster(nCores)

tic('mcmc')
bmm<-compute_mallows_mixtures(n_clusters=c(1,4,7), cl=cl, rankings=sushi_rankings, nmc=5000, save_clus=F, include_wcd=F)
toc()
stopCluster(cl)

assess_convergence(bmm)
assess_convergence(bmm, parameter='rho')
assess_convergence(bmm, parameter='cluster_probs')

cl<-makeCluster(nCores)
tic('mcmc')
bmm<-compute_mallows_mixtures(n_clusters=1:10, cl=cl, rankings=sushi_rankings, nmc=100000, rho_thinning=10, save_clus=F, include_wcd=T)
toc()
stopCluster(cl)
# elapsed time is 21`43``

bmm
attributes(bmm)
plot_elbow(bmm, burnin=5000)

tic('mcmc')
bmm_sushi<-compute_mallows(rankings=sushi_rankings, n_clusters = 3,  nmc=100000, save_clus=T, clus_thin = 10, rho_thinning = 10, verbose =T )
toc()
# elapsed time is 5`28``

bmm_sushi$burnin<-5000
plot(bmm_sushi, parameter='cluster_probs')
plot(bmm_sushi, parameter='cluster_assignment')

## There is something wrong in implementing cluster model...
## All assessors happen to fall in only one cluster...

library(tidyr)
compute_consensus(bmm_sushi) %>%
  select(-cumprob) %>%
  spread(key=cluster, value=item)


## From estimate_partition_function_example.R file
## Importance sampling strategy for approximating partition function
alpha_vector <- seq(from = 0, to = 10, by = 0.5)
n_items <- 20
metric <- "spearman"
degree <- 10

# We start with 1e3 Monte Carlo samples
fit <- estimate_partition_function(method = "importance_sampling",
                                    alpha_vector = alpha_vector,
                                    n_items = n_items, metric = metric,
                                    nmc = 1e+4, degree = degree)
# A vector of polynomial regression coefficients is returned
fit


# We write a little function for storing the estimates in a dataframe
library(dplyr)
powers <- seq(from = 0, to = degree, by = 1)

compute_fit <- function(fit){
  tibble(alpha = alpha_vector) %>%
    rowwise() %>%
    mutate(logz_estimate = sum(alpha^powers * fit))
}

estimates <- compute_fit(fit)


# We can now plot the two estimates side-by-side
library(ggplot2)
ggplot(estimates, aes(x = alpha, y = logz_estimate)) +
  geom_line()
# We see that the two importance sampling estimates, which are unbiased,
# overlap. The asymptotic approximation seems a bit off. It can be worthwhile
# to try different values of n_iterations and K.

# When we are happy, we can provide the coefficient vector in the
# logz_estimate argument to compute_mallows
# Say we choose to use the importance sampling estimate with 1e4 Monte Carlo samples:
bmm_potato_is <- compute_mallows(potato_visual, metric = "spearman",
                             logz_estimate = fit, nmc = 501000, verbose = T)

bmm_potato_is$burnin <- 1000
plot(bmm_potato_is, parameter = 'rho', items = 1:20)
plot(bmm_potato_is, parameter = 'alpha')


























































