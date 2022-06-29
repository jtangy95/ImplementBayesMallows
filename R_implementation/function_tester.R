library(Rcpp)
library(RcppArmadillo)
library('BayesMallows')
data("potato_visual")
potato_visual

sourceCpp('/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_distances.cpp')

n=20
R1 = sample(1:n, n)
R2 = sample(1:n, n)

footrule_distance(R1, R2)
spearman_distance(R1, R2)

metric = "footrule"
get_rank_distance(R1, R2, metric)

rho_ref = 1:n
cbind(rho_ref, t(potato_visual))
rank_dist_vec(t(potato_visual), rho_ref, metric)
rank_dist_sum(t(potato_visual), rho_ref, metric)

sourceCpp('/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_leap_and_shift.cpp')

sourceCpp('/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_parameter_updates.cpp')

