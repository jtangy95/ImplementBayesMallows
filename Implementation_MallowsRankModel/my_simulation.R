source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_sample_mallows.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_mallows_combined.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_plot.BayesMallows_cluster.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_consensus.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_topk.R")
load(file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_potatoes.RData", verbose = T)
load(file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_sushi.RData", verbose = T)
load(file = "/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/dataset_beach.RData", verbose = T)
library(dplyr)


n_items = 30
n_samples = 200
n_clusters = 4
rho0 = matrix(0, ncol= n_clusters, nrow = n_items)
alpha0 = rep(0, n_clusters)
for(i in 1:n_clusters){
  rho0[, i] = sample(1:n_items, n_items)
  alpha0[i] = rexp(n=1 ,rate = 3)
}
samples = sample_mallows_mixed(rho0, alpha0, n_samples, n_clusters)

top5sample = create_topk(samples, 15)

fit = compute_mallows(rankings = top5sample , nmc = 1e5, verbose = T, n_clusters = 5,
                      leap_size = 2, alpha_prop_sd = 0.1, lambda = 0.1, alpha_jump = 100, 
                      save_aug = T, aug_thinning = 100, clus_thin = 10, rho_thinning  = 10)

fit$burnin = 5*1e3
compute_consensus(fit) %>% filter(cluster == "Cluster 5")
compute_consensus(fit, parameter= "Rtilde", assessors = 1)


assessor_index = 30
pred = compute_consensus(fit, parameter= "Rtilde", assessors = assessor_index) %>% filter(ranking >15 & ranking <=20) %>% pull(item)

real = colnames(samples)[which(samples[assessor_index,] %in% 16:20)]

sum(pred %in% real)

# comb = function(n, x) {
#   factorial(n) / factorial(n-x) / factorial(x)
# }


rightnum = rep(0, n_samples)
for(i in 1:n_samples){
  assessor_index = i
  pred = compute_consensus(fit, parameter= "Rtilde", assessors = assessor_index) %>% filter(ranking >15 & ranking <=20) %>% pull(item)
  real = colnames(samples)[which(samples[assessor_index,] %in% 16:20)]
  rightnum[i] = sum(pred %in% real)
}

counts = as.numeric(table(rightnum))
proportions = counts / sum(counts)

data = read.csv(file = "/Users/changtaeyeong/Downloads/ratings_small.csv")
length(as.numeric(unique(data %>% pull(movieId))))
summary(as.numeric(table(data %>% pull(movieId))))
sum(as.numeric(table(data %>% pull(movieId))) > 80)
pop_movies = dimnames(table(data %>% pull(movieId)))[[1]][which(as.numeric(table(data %>% pull(movieId))) > 80)]
length(pop_movies)
length(unique(data %>% filter(movieId %in% pop_movies) %>% pull(userId)))


numpairs = 15
combs = combn(1:30, 2)
ncombs = ncol(combs)
assessor = rep(1:n_samples, each = numpairs)
top_item = rep(0, numpairs * n_samples)
bottom_item = rep(0, numpairs * n_samples)
for(i in 1:n_samples){
  rank = samples[i,]
  pairs_indices = sample(1:ncombs, numpairs)
  pairs = combs[, pairs_indices]
  for(j in 1:numpairs){
    pair = pairs[,j]
    first = pair[1]
    second = pair[2]
    if(rank[first] < rank[second]){
      top_item[numpairs*(i-1)+j] = first
      bottom_item[numpairs*(i-1)+j] = second
    } else{
      top_item[numpairs*(i-1)+j] = second
      bottom_item[numpairs*(i-1)+j] = first
    }
  }
  
}
simul_pref = tibble(assessor= assessor, bottom_item = bottom_item, top_item = top_item)


fit = compute_mallows(preferences = simul_pref , nmc = 1e5, verbose = T, n_clusters = 4,
                      leap_size = 6, alpha_prop_sd = 0.1, lambda = 0.1, alpha_jump = 10, 
                      save_aug = T, aug_thinning = 10, clus_thin = 10, rho_thinning  = 10)
fit$burnin = 1e4

# comb = function(n, x) {
#   factorial(n) / factorial(n-x) / factorial(x)
# }

each_assessor <- split(simul_pref[, c("bottom_item", "top_item"), drop = FALSE],simul_pref$assessor)
constrained_items = lapply(each_assessor, function(x){unique(c(x[["bottom_item"]], x[["top_item"]]))} )

assessor_index = 8

con_items = paste("Item", constrained_items[[assessor_index]])
con_items

result = compute_consensus(fit, parameter= "Rtilde", assessors = assessor_index) 

pred = setdiff(result$item, con_items)[1:3]
real = names(sort(samples[assessor_index,][!(colnames(samples[assessor_index, , drop = F]) %in% con_items)]))[1:3]
sum(pred %in% real)

tic("summary")
rightnum = rep(0, n_samples)
for(i in 1:n_samples){
  assessor_index = i
  con_items = paste("Item", constrained_items[[assessor_index]])
  result = compute_consensus(fit, parameter= "Rtilde", assessors = assessor_index) 
  pred = setdiff(result$item, con_items)[1:3]
  real = names(sort(samples[assessor_index,][!(colnames(samples[assessor_index, , drop = F]) %in% con_items)]))[1:3]
  rightnum[i] = sum(pred %in% real)
  print(paste(i, "-th iteration completed"))
}
toc()

counts = as.numeric(table(rightnum))
proportions = counts / sum(counts)




