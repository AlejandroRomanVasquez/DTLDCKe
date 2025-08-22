This repository contains the source code for the paper "D-Vine Copula Based Knockoffs for Variable Selection in Gene Expression Studies"

The computational implementations of the proposed methods were carried out using R, employing versions 3.6.3.

The code for running the simulations and the application to a real dataset uses several packages from R. The packages (with versions) are the following:

ggplot2_3.3.5            survival_3.4-0         dplyr_1.0.9            doParallel_1.0.17     
iterators_1.0.14         knockoff_0.3.5         foreach_1.5.2          glmnet_4.1-4          
Matrix_1.5-3             TSP_1.2-4              VineCopula_2.4.5       rvinecopulib_0.6.2.1.3
seqknockoff_0.0.0.9000   latentcor_1.2.0 


The folder "Simulations" contains code related to three distinct data-generating processes (DGPs) for the predictors X considered in the paper: a multivariate normal distribution, a t-tailed Markov chain, and Survival Regression with Block Correlated AR1 Structure. 

The folder "Application" contains code and data for applying the proposed methodology to a real lung cancer dataset.
