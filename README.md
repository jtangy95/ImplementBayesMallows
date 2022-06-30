# ImplementBayesMallows
I have studied about the article "[Probabilistic preference learning with the Mallows rank model](https://jmlr.org/papers/v18/15-481.html)" which suggests the Bayesian Mallows rank model. With understanding how this model works, I have simulated whether individual recommendation works well using the method.

In the `Notes_MallowsRankModel` folder, there are two tex files
- `Notes for Mallows Rank model Article` : My summary for the article
- `Seminar_MallowsRankModel` : Lab seminar presentation about the article  




In the `R_implementation` folder, there are R and C++ codes for implementation and simulation.
- Implementation : based on the `BayesMallows` package created by the authors of the article
    - `my_compute_mallows.R` & `my_run_mcmc.cpp` :  For the baseline setting - complemte rankings and homogeneous assessors
    - `my_compute_partial.R` & `my_run_mcmc_partial.cpp` : For partial rankings (ex. top-k ranks)
    - `my_compute_cluster.R` & `my_run_mcmc_cluster.cpp` : For the case of heterogenous assessors ; clustering involved
    - `my_compute_mallows.R` & `my_run_mcmc_combined.cpp` : Combining all the settings above
- Simulation and Application to real data : `my_simulation.R` 
    - Simulation to predict missing individual preference 
    - Applying our approach to a movie rating data of individual users