source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_sample_mallows.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_mallows_combined.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_plot.BayesMallows_cluster.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_compute_consensus.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/ImplementBayesMallows/Implementation_MallowsRankModel/my_topk.R")

source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/BayesMallows/R/compute_mallows_mixtures.R")
source("/Users/changtaeyeong/Desktop/BayesMallowsRankModel/BayesMallows/R/plot_elbow.R")

library(dplyr)
library(tictoc)

# Generating mallows rank samples 
set.seed(10)
n_items = 25
n_samples = 1000
n_clusters = 4
rho0 = matrix(0, ncol= n_clusters, nrow = n_items)
alpha0 = rep(0, n_clusters)
for(i in 1:n_clusters){
  rho0[, i] = sample(1:n_items, n_items)
  alpha0[i] = rexp(1 ,rate = 1)
}
samples = sample_mallows_mixed(rho0, alpha0, n_samples, n_clusters)

# Convert complete rank samples to top-10 rank samples
top10sample = create_topk(samples, 10)
head(top10sample)

# Fit mixture model for n_clusters=1:10
tic("mixture")
fit_mixed = compute_mallows_mixtures(n_clusters=1:10, rankings = top10sample , nmc=1e5, rho_thinning=10, save_aug = F, alpha_jump = 10,
                               leap_size = 5, lambda = 0.1, alpha_prop_sd = 0.1,   clus_thin = 10 , include_wcd=T)
toc()

# Find the elbow in the wcd plot to choose C
plot_elbow(fit_mixed, burnin = 1e4)

# Fit with determined C
tic("simul")
fit = compute_mallows(rankings = top10sample , nmc = 1e5, verbose = T, n_clusters = 4,
                      leap_size = 5, alpha_prop_sd = 0.1, lambda = 0.1, alpha_jump = 10, 
                      save_aug = T, aug_thinning = 20, clus_thin = 10, rho_thinning  = 10)
toc()

fit$burnin = 1e4
compute_consensus(fit) %>% filter(cluster == "Cluster 2")
compute_consensus(fit, parameter= "Rtilde", assessors = 1)


assessor_index = 10
pred = compute_consensus(fit, parameter= "Rtilde", assessors = assessor_index) %>% filter(ranking >10 & ranking <=15) %>% pull(item)

real = colnames(samples)[which(samples[assessor_index,] %in% 11:15)]

sum(pred %in% real)

## Calculate the number of correctly recommended items for each assessor by random draw
# comb = function(n, x) {
#   factorial(n) / factorial(n-x) / factorial(x)
# }
# proportions_true = rep(0, 6)
# for(j in 0:5){
#   proportions_true[j+1]  = comb(5, j) * comb(10, 5-j)
# }
# proportions_true = proportions_true / comb(15,5) 
# proportions_true = round(proportions_true, 5)

# Calculate the number of correctly recommended items for each assessor by our model
rightnum = rep(0, n_samples)
tic("summary")
for(i in 1:n_samples){
  assessor_index = i
  pred = compute_consensus(fit, parameter= "Rtilde", assessors = assessor_index) %>% filter(ranking >10 & ranking <=15) %>% pull(item)
  real = colnames(samples)[which(samples[assessor_index,] %in% 11:15)]
  rightnum[i] = sum(pred %in% real)
  print(paste(i, "-th iteration completed"))
}
toc()

table(rightnum)
counts = as.numeric(table(rightnum))
proportions = counts / sum(counts)
proportions

cumsum(rev(proportions))
cumsum(rev(proportions_true))

# temp_save$"zero" = table(rightnum)
# temp_save$"seed1" = table(rightnum)
# temp_save$"seed2" = table(rightnum)
# temp_save$"seed3" = table(rightnum)
# temp_save$"seed4" = table(rightnum)
# temp_save$"seed5" = table(rightnum)
# temp_save$"seed6" = table(rightnum)

#results
seed1 = c(71, 332, 415, 160, 22, 0)
seed2 = c(77, 303, 432, 171, 16, 1)
seed3 = c(47, 266, 391, 248, 47, 1)
seed4 = c(83, 315, 389, 191, 20, 2)
seed5 = c(36, 198, 370, 301, 90, 5)
seed6 = c(63, 326, 418, 163, 28, 2)
seed7 = c(64, 262, 419, 213, 41, 1)
seed8 = c(47, 225, 376, 276, 68, 8)
seed9 = c(54, 245, 337, 265, 93, 6)
seed10 = c(46, 288, 414, 204, 48, 0)

result = rbind(seed1, seed2, seed3, seed4, seed5, seed6, seed7, seed8, seed9, seed10)
result
apply(result, 1, sum)


comb = function(n, x) {
  factorial(n) / factorial(n-x) / factorial(x)
}
proportions_true = rep(0, 6)
for(j in 0:5){
  proportions_true[j+1]  = comb(5, j) * comb(10, 5-j)
}
proportions_true = proportions_true / comb(15,5)
proportions_true = round(proportions_true, 5)

true = cumsum(rev(proportions_true))

func<- function(x){
  res = x / sum(x)
  return(cumsum(rev(res)))
}

chart = apply(result, 1, func)
t(chart)

apply(chart, 1, mean)
apply(chart, 1, sd)
true


# Application to real data : movie ratings in TMDB

data = read.csv(file = "/Users/changtaeyeong/Downloads/ratings_small.csv")
data  = tibble(data)
data = data %>% select(-timestamp)
links = read.csv(file = "/Users/changtaeyeong/Downloads/links_small.csv")
links = tibble(links)
links = links %>% select(-imdbId)
data2 = inner_join(data, links, by = "movieId")
data2 = data2 %>% mutate(across(c(userId, movieId, tmdbId), as.character))

names = read.csv(file = "/Users/changtaeyeong/Downloads/movies_metadata.csv")
names = tibble(names)
names = names %>% select(c(id, title))

movies = inner_join(data2, names, by = c("tmdbId" = "id"))
movies %>% group_by(userId) %>% count() 
movies %>% group_by(movieId) %>% count() %>% filter(n >= 100)

pop_movies = movies %>% group_by(movieId) %>% count() %>% filter(n >= 100) %>% ungroup() %>% pull(movieId)
pop_movies_title = unique(movies %>% filter(movieId %in% pop_movies) %>% pull(title))

ratings = movies %>% filter(movieId %in% pop_movies) %>% select(-tmdbId)

df = ratings %>% group_by(movieId) %>% count() 
label = tibble(movieLabel = 1:nrow(df))
df = bind_cols(df, label)
movielabel = df %>% select(-n) %>% ungroup()
ratings = inner_join(ratings, movielabel, by = "movieId")

df = ratings %>% group_by(userId) %>% count()
label = tibble(userLabel = 1:nrow(df))
df = bind_cols(df, label)
userlabel = df %>% select(-n) %>% ungroup()
ratings = inner_join(ratings, userlabel, by = "userId")

n_movies = nrow(ratings %>% group_by(movieId) %>% count())
n_users =  nrow(ratings %>% group_by(userId) %>% count())
users = ratings %>% group_by(userLabel) %>% count() %>% pull(userLabel)
counts = ratings %>% group_by(userLabel) %>% count() %>% pull(n)


ratings_pref = tibble(assessor = double(), bottom_item = double(), top_item = double())
tic("creating preference dataframe")
for(j in 1:n_users){
  if(counts[j] < 2){
    next
  }
  combs = combn(1:counts[j], 2)
  dat = ratings %>% filter(userLabel == users[j])
  
  result = tibble(bottom_item = double(), top_item = double())
  
  for(i in 1:ncol(combs)){
    pair = combs[,i]
    pair_rating = dat[pair, ] %>% pull(rating)
    pair_movie = dat[pair, ] %>% pull(movieLabel)
    if(pair_rating[1] > pair_rating[2]){
      result = result %>% add_row(bottom_item = pair_movie[2], top_item = pair_movie[1])
    } else if(pair_rating[1] < pair_rating[2] ){
      result = result %>% add_row(bottom_item = pair_movie[1], top_item = pair_movie[2])
    }
  }
  
  result = result %>% mutate(assessor = rep(users[j], nrow(result))) %>% select(assessor, bottom_item, top_item)
  
  ratings_pref = bind_rows(ratings_pref, result)
}
toc()

ratings_pref_tc = generate_transitive_closure(ratings_pref)

tic("generating initial ranking")
ratings_pref_init_rank = generate_initial_ranking(ratings_pref_tc)
toc()

tic("generating constraints list")
ratings_pref_const = generate_constraints(ratings_pref_tc, n_movies)
toc()


alpha_vector <- seq(from = 0.01, to = 10, length = 100)
tic("estimating partition function beta parameters")
logz_estimate = estimate_partition_function(alpha_vector, n_movies, "footrule" , 1e5)
toc()

tic("mixture")
fit_mixed = compute_mallows_mixtures(n_clusters=1:8 , rankings = ratings_pref_init_rank , preferences = ratings_pref , nmc = 1e5, 
                                     verbose = T, leap_size = 5, alpha_prop_sd = 0.1, lambda = 0.1, alpha_jump = 10, 
                                     save_aug = T, aug_thinning = 20, clus_thin = 10, rho_thinning  = 10, 
                                     logz_estimate =  logz_estimate, constraints = ratings_pref_const, include_wcd = T)
toc()

plot_elbow(fit_mixed, burnin = 1e4)

tic("simul")
fit = compute_mallows(rankings = ratings_pref_init_rank ,preferences = ratings_pref , nmc = 1e5, verbose = T, n_clusters = 1,
                      leap_size = 5, alpha_prop_sd = 0.1, lambda = 0.1, alpha_jump = 10, 
                      save_aug = T, aug_thinning = 20, clus_thin = 10, rho_thinning  = 10, 
                      logz_estimate =  logz_estimate, constraints = ratings_pref_const)
toc()

labs = tibble(userLabel = unique(ratings_pref %>% pull(assessor)) , fit_userlabel = unique(fit$augmented_data %>% pull(assessor)))

ratings_fit = left_join(ratings, labs, by="userLabel")


fit$burnin = 1e4

movietitle = ratings %>% group_by(movieId, title) %>% count() %>% ungroup() %>% select(movieId, title)
movielabs = inner_join(movielabel, movietitle, by= "movieId")

movieitems = movielabs %>% mutate(item = paste("Item" , movieLabel) ) %>% select(item, title)

ass1 = compute_consensus(fit, parameter = "Rtilde" , assessor = 1)

ass1 = inner_join(ass1, movieitems, by= "item")
ratings_fit %>% filter(fit_userlabel == 1) %>% select(rating, title, userId)
ass1 %>% filter(ranking <20)

ass2 = compute_consensus(fit, parameter = "Rtilde" , assessor = 2)

ass2 = inner_join(ass2, movieitems, by= "item")
ratings_fit %>% filter(fit_userlabel == 2) %>% select(rating, title, userId)
ass2 %>% filter(ranking <20)

ass3 = compute_consensus(fit, parameter = "Rtilde" , assessor = 3)

ass3 = inner_join(ass3, movieitems, by= "item")
ratings_fit %>% filter(fit_userlabel == 3) %>% select(rating, title, userId) %>% print(n=100)
ass3 %>% filter(ranking <20)

tibble1 = ratings_fit %>% filter(fit_userlabel == 1) %>% select(rating, title, userId)
tibble2 = ass1 %>% filter(ranking <=20)
tibble1
tibble2
anti_join(tibble2, tibble1, by = "title")


tibble1 = ratings_fit %>% filter(fit_userlabel == 2) %>% select(rating, title, userId)
tibble2 = ass2 %>% filter(ranking <=20)
tibble1
tibble2
anti_join(tibble2, tibble1, by ="title")


ass = compute_consensus(fit, parameter = "Rtilde", assessor = 162)
ass = inner_join(ass, movieitems, by= "item")
ratings_fit %>% filter(fit_userlabel == 162) %>% select(rating, title, userId) 
ass %>% filter(ranking <20)

tibble1 = ratings_fit %>% filter(fit_userlabel == 162) %>% select(rating, title, userId) 
tibble2 = ass %>% filter(ranking <=20) 
tibble1 %>% print(n=Inf)
tibble2 %>% print(n = Inf)
anti_join(tibble2, tibble1, by = "title")