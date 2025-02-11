#### libraries ####

library(rvinecopulib)
library(VineCopula)
library(TSP)
library(parallel)
library(glmnet)
library(foreach)
library(knockoff)
library(doParallel)
library(dplyr)
library(survival)
library(dplyr)
library(ggplot2)


##### Auxiliary functions #####


#A function that binns the continuous variables of a matrix X by 
#a specific quantile corresponding to the given probability prob

Continuous_to_binary_by_quantile <- function(X, column_type, prob=0.5){
  
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  
  if (is.null(column_type)) {
    stop("Argument column_type is missing")  
  }
  
  #Number of variables p 
  p <- dim(X)[2]
  
  X_mixed <- X
  for(i in 1:p) {   
    if (column_type[i]=="num"){
      X_mixed[,i] <- X[,i]
    }
    else {
      quantile_i <- as.numeric(quantile(X[,i],prob=prob)) 
      X_mixed[,i] <- ifelse(X[,i]<quantile_i,0,1)
    } 
  }
  return(X_mixed)
}

#Ordinal to a latent continuous uniform representation

ordinal_to_uniform <- function(column){
  
  if (is.null(column)) {
    stop("Argument column is missing")  
  }
  
  # An ordinal variable with more than 10 levels is considered a numerical variable
  if(length(unique(column)) < 11){
    
    column_lower <- column - 1
    
    #Empirical Cumulative Distribution Function
    ecdf_column <- ecdf(column)  
    
    u_upper <-   ecdf_column(column)
    u_lower <-  ecdf_column(column_lower)
    
    u_cond <- runif(n=length(column),min=u_lower, max= u_upper)
    
    return (u_cond)
  }
  
  else{
    return (column)
  }
}

#Function to identify if the variable is numerical(num) or ordinal (ord)
column_type_identification <- function(column){ 
  
  if (is.null(column)) {
    stop("Argument column is missing")  
  }
  
  # An ordinal variable with more than 10 levels is considered a numerical variable
  if(length(unique(column)) < 11) 
  {type<-"ord"}
  else
  {type<-"num"} 
  return(type)
}    

#Uniform to original (transformation)
uniform_to_original <- function(X, u_Xk, column_type){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  
  if (is.null(u_Xk)) {
    stop("Argument u_Xk is missing")  
  }
  
  if (is.null(column_type)) {
    stop("Argument column_type is missing")  
  }
  
  #Number of variables p 
  p <- dim(X)[2]
  
  Xk <- X
  for(i in 1:p) {   
    if (column_type[i]=="num"){
      Xk[,i] <- as.vector(quantile(X[,i], probs = punif(u_Xk[,i],min=0, max=1), type=8))
    }
    else {
      Xk[,i]<- round(as.vector(quantile(X[,i],  probs = u_Xk[,i], type=1)),0)
    } 
  }
  return(Xk)
  
}


#### Utility functions for the e-values procedure ####

#These functions are taken from
#https://github.com/zhimeir/derandomized_knockoffs_fdr


# The eBH procedure
### Input: 
###   E: e-values
###   alpha: target FDR level
### Output:
###   Variables selected by the e-BH procedure

ebh <- function(E, alpha){
  
  p <- length(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (p / alpha / (1:p))
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(list(rej = rej))
}

## Computing the early stopping time ##
### Input:
###   W: vector of knockoff feature importance statistics 
###   gamma: alpha_kn 
###   offset: value between 0 and 1
### Output: 
###   The modified knockoff stopping time defined in (14)

stop_early <- function(W, gamma, offset){
  
  tau <- alphakn_threshold(W, fdr =  gamma, offset = offset) 
  ord_W <- order(abs(W), decreasing = TRUE)
  sorted_W <- W[ord_W]
  
  if(sum(W>0) >= 1 / gamma){
    pos_ind <- which(sorted_W > 0)
    tau1 <- sorted_W[pos_ind[ceiling(1/gamma)-1]]
  }else{
    tau1 <- 0
  }
  tau <- min(tau,tau1) 
  
  return(tau)
}

## Compute stopping time w/ diff alpha_kn and offset ##
### Input:
###   W: a length p vector of knockoff feature importance statistics
###   fdr: the target FDR level
###   offset: 0 or 1 
### Output: 
###   the knockoff selection threshold

alphakn_threshold <- function(W, fdr, offset) {
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}


##### R functions related to Dvines #####

#get_dvine_order() -->   This function runs an heuristic procedure to determine the 
#order for the first tree in a D-vine structure using the TSP R package to solve 
#the traveling salesman problem. To solve it, we need to identify the shortest 
#Hamiltonian path by assigning weights based on the pairwise Kendall’s τ

#Arguments:
#X: matrix of predictors
#random_order: logic value for a nonrandom assignation order of variables using the 
#            argument "identity" in the solve_TSP() function.

get_dvine_order <- function(X, random_order=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is null")  
  }
  
  #Matrix transformation
  X <- as.matrix(X)
  
  #Matrix of 1 - tau_ij
  M_tau <- 1 - abs(TauMatrix(X))
  
  #Hamiltonian path and solution (functions of package TSP)
  hamilton <- insert_dummy(TSP(M_tau), label="cut")
  
  if (random_order==FALSE){
    sol <- solve_TSP(hamilton,method="identity")
  }
  else{
    sol <- solve_TSP(hamilton,method="repetitive_nn")
  }
  
  #Reordering
  TSP_order <- cut_tour(sol,"cut")
  names(TSP_order) <- NULL
  
  return(TSP_order)
  
}

#X_Xk_dvine_distributions()-->Function to fit the dvine distribution for X and X_X matrices
#Arguments:
#X_cont: matrix of continuous predictors (after transformation of ordinal to conditional uniforms)
#vinecop_family : String related to the family of pair copulas used in the dvine fitting. 
#                 Common options are "parametric", "nonparametric", "onepar". More details can be found
#                 in the documentation of R package rvinecopulib
#                  https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf
# n_cores: int -> number of cores for parallel processing

#Value: This function returns a list that contains objects of class vinecop_dist for X and X_X
#Note: more information about objects of class vinecop_dist can be found in 
#https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf

X_Xk_dvine_distributions <- function(X_cont, column_type, vinecop_family="parametric", n_cores=1){
  
  if (is.null(X_cont)) {
    stop("Argument X_cont is null")  
  }
  
  if (is.null(column_type)) {
    stop("Argument column_type is null")  
  }
  
  
  #Number of variables p and sample size n
  n <- dim(X_cont)[1]
  p <- dim(X_cont)[2]
  
  
  #dstructures for dvines
  X_X_dstructure <- dvine_structure((2*p):1)
  X_dstructure <- dvine_structure(p:1)
  
  #Dataset column binding
  X_X <- cbind(X_cont,X_cont)
  
  #Seudo-Observations u_X
  u_X <- matrix(data=rep(0,n*p),nrow=n,ncol=p)
  for(j in 1:p){
    if(column_type[j]=="ord") {
      u_X[,j] <- X_cont[,j]
    }
    else{
      u_X[,j] <- pseudo_obs(X_cont[,j])
    }
  }
  #Seudo-Observations for the duplicated data u_X_X  
  u_X_X <- cbind(u_X,u_X)
  
  #Fitting dvine distribution for X_X
  dvine_fitting_time <- system.time(
    fit_dvine_trunc <- vinecop(u_X_X, family_set=c(vinecop_family), structure= X_X_dstructure, presel=TRUE,
                               selcrit='mbicv', par_method='mle', psi0=0.95, show_trace=FALSE, cores=n_cores, trunc_lvl=p-1)
  )
  
  #Printing dvine X_X fitting time
  print("dvine fitting time in seconds:")
  print(dvine_fitting_time)
  
  #Pair-copula list for X_X
  X_X_dvine_pclist <- fit_dvine_trunc$pair_copulas
  
  #dvine distribution for X_Xk 
  X_X_dvine_dist <- vinecop_dist(X_X_dvine_pclist, X_X_dstructure)
  
  #Pair-copula list for X
  X_dvine_pclist <- list(rep(list(""),p-1))
  
  #Iniziating with Independent copula
  for (i in 1:(p-1)){
    bicop <- bicop_dist("indep",)
    X_dvine_pclist[i] <- list(rep(list(bicop),p-i))
  }
  
  #Pair copula list just for X dependencies
  for (i in 1:(p-1)){
    J <- p-i
    
    for (j in 1:J){
      X_dvine_pclist[[i]][j] <- X_X_dvine_pclist[[i]][j] 
      
    } 
  }
  
  # dvine distribution for X
  X_dvine_dist <- vinecop_dist(X_dvine_pclist, X_dstructure)
  
  #List with dvine distributions
  dvine_distributions <- list(X_dvine_dist=X_dvine_dist, 
                              X_X_dvine_dist=X_X_dvine_dist,
                              u_X=u_X)
  
  return(dvine_distributions)
}

# create_dvine_Knockoffs()--> Function to sample dvine knockoffs

#Arguments:
#X: matrix of original predictors
#X_cont: matrix of continuous predictors (after transformation of ordinal to conditional uniforms)
#column_type: vector with the type of each variable of X ("num" numeric, "ord" ordinal)
#dvine_distributions: a list containing: 
#   1) object of class vinecop_dist for X_X,
#   2) Object of class vinecop_dist for X_cont,
#   3) Pseudo observations of X_cont saved in u_X
#n_cores: int -> number of cores for parallel processing
#seed_Xk: integer value to reproduce the sampled Knockoffs. The default is NULL to ensure a random knockoffs simulation 
#Note: more information about objects of class vinecop_dist can be found in 
#https://cran.r-project.org/web/packages/rvinecopulib/rvinecopulib.pdf

#Value: This function returns a matrix Xk of knockoffs 

create_dvine_Knockoffs <- function(X, X_cont, column_type, dvine_distributions, n_cores=1, seed_Xk=NULL){
  
  if (is.null(X)) {
    stop("Argument X is null")  
  }
  if (is.null(X_cont)) {
    stop("Argument X_cont is null")  
  }
  if (is.null(column_type)) {
    stop("Argument colum_type is null")  
  }
  if (is.null( dvine_distributions)) {
    stop("Argument dvine_distributions is null")  
  }
  
  set.seed(seed_Xk) #For reproducibility
  
  #Number of variables p and sample size n
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Pseudo observations
  u_X <- dvine_distributions$u_X
  
  #Independent uniforms w
  w_X <- rosenblatt(x=u_X, model=dvine_distributions$X_dvine_dist, cores = n_cores)
  w_Xk <- matrix(runif(n=p*n,min=0,max=1),nrow=n,ncol=p)
  w_X_Xk <- cbind(w_X,w_Xk)
  
  #Knockoff sampling Xk
  u_X_Xk <- inverse_rosenblatt(u=w_X_Xk, model= dvine_distributions$X_X_dvine_dist, cores = n_cores)
  u_Xk <- u_X_Xk[,(p+1):(2*p)]
  
  #Transformation to original predictors
  Xk <- uniform_to_original(X, u_Xk, column_type ) 
  
  #Xk column names
  Xk_column_names <- sapply(1:p, function(number) paste0("Xk", number))
  colnames(Xk ) <- Xk_column_names
  
  return(Xk)
}

# stable_lasso_glmnet_parallel()--> Function to fit a regularized lasso regresion model using 
#some functions of the R package glmnet.
#It implements the stabilizing procedure of Roberts and Nowak (2014) to diminish sensitivity
#to the fold assignment used in cross-validation to select the hyperparameter lambda
# The running is in parallel using the function foreach()

#Arguments:
#X: matrix of predictors
#y: vector of response
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial"
#M_lasso: integer related to the number of runs for the stabilization against CV (Roberts and Nowak, 2014)
#n_folds: integer indicating the number of cross validations
#random_cv: logic value for a random assignation of the folds in the Cross-validation
#        TRUE is the default. FALSE is for reproducibility purposes. 
#Note: more information about the R package glmnet can be found in 
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
#Note 2: this function runs in parallel for the stabilization against CV (Roberts and Nowak, 2014)

#Value: This function returns a vector of the estimated coeficientes (without the intercept)


stable_lasso_glmnet_parallel <- function(X, y, lasso_family, M_lasso = 10, n_folds = 5, random_cv=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  if(lasso_family=="gaussian"){
    y_glmnet <- as.vector(y)}
  else if(lasso_family=="binomial"){
    y_glmnet  <- as.factor(y)}
  else {
    y_glmnet <- y}
  
  X_matrix <- as.matrix(X)
  
  #Stabilizing the lasso against CV (Roberts and Nowak, 2014)
  lambdas <- rep(0,M_lasso)
  
  #Random assignation of the folds in the Cross-validation
if(random_cv==TRUE){
    start_time <- Sys.time()  
    lambdas <- foreach(i = 1:M_lasso, .combine=c,.packages=c("glmnet")) %dopar% {
      set.seed(NULL)
      cvfit <- cv.glmnet(X_matrix, y_glmnet, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
      cvfit$lambda.min
    }
    end_time <- Sys.time()
    print("Time for running parallel Stabilizing the lasso:")
    print(end_time-start_time)  }
  
  #Non-random assignation of the folds in the Cross-validation
  else{
    start_time <- Sys.time() 
    lambdas <- foreach(i = 1:M_lasso, .combine=c,.packages=c("glmnet")) %dopar% {
      set.seed(i)
      cvfit <- cv.glmnet(X_matrix, y_glmnet, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
      cvfit$lambda.min
    }
    end_time <- Sys.time()
    print("Time for running parallel Stabilizing the lasso:")
    print(end_time-start_time)
  }

  #Selecting the median of the lambdas distribution
  lambda50 <- as.numeric(quantile(lambdas,probs=0.5))

  fit_coef <- coef(glmnet(X_matrix, y_glmnet, alpha = 1, lambda = lambda50, family = lasso_family, standardize = TRUE))
  fit_coef_vec <- as.vector(fit_coef)
  
  if(lasso_family=="cox")
    fit_coef_vec <- fit_coef_vec
  else {
    fit_coef_vec <- fit_coef_vec[-1]   
  }
  return(fit_coef_vec)
}


# stable_lasso_glmnet()--> Function to fit a regularized lasso regresion model using 
#some functions of the R package glmnet.
#It implements the stabilizing procedure of Roberts and Nowak (2014) to diminish sensitivity
#to the fold assignment used in cross-validation to select the hyperparameter lambda

#Arguments:
#X: matrix of predictors
#y: vector of response
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial"
#M_lasso: integer related to the number of runs for the stabilization against CV (Roberts and Nowak, 2014)
#n_folds: integer indicating the number of cross validations
#random_cv: logic value for a random assignation of the folds in the Cross-validation
#        TRUE is the default. FALSE is for reproducibility purposes. 
#Note: more information about the R package glmnet can be found in 
#https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
#Note 2: this function runs in parallel for the stabilization against CV (Roberts and Nowak, 2014)

#Value: This function returns a vector of the estimated coeficientes (without the intercept)


stable_lasso_glmnet <- function(X, y, lasso_family, M_lasso = 10, n_folds = 5, random_cv=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  if(lasso_family=="gaussian"){
    y_glmnet <- as.vector(y)}
  else if(lasso_family=="binomial"){
    y_glmnet  <- as.factor(y)}
  else {
    y_glmnet <- y}
  
  X_matrix <- as.matrix(X)
  
  #Stabilizing the lasso against CV (Roberts and Nowak, 2014)
  lambdas <- rep(0,M_lasso)
  
  #Random assignation of the folds in the Cross-validation
  if(random_cv==TRUE){
    start_time <- Sys.time()  
    for (i in 1:M_lasso){
      set.seed(NULL)
      cvfit <- cv.glmnet(X_matrix, y_glmnet, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
      lambdas[i] <- cvfit$lambda.min
    }
    end_time <- Sys.time()
    print("Time for running Stabilizing the lasso:")
    print(end_time-start_time)
  }
  #Non-random assignation of the folds in the Cross-validation
  else{
    start_time <- Sys.time() 
    for( i in 1:M_lasso){
      set.seed(i)
      cvfit <- cv.glmnet(X_matrix, y_glmnet, alpha=1, family = lasso_family, nfolds = n_folds, standardize = TRUE)
      lambdas[i] <- cvfit$lambda.min
    }
    end_time <- Sys.time()
    print("Time for running Stabilizing the lasso:")
    print(end_time-start_time)
  }
  
  #Selecting the median of the lambdas distribution
  lambda50 <- as.numeric(quantile(lambdas,probs=0.5))
  
  fit_coef <- coef(glmnet(X_matrix, y_glmnet, alpha = 1, lambda = lambda50, family = lasso_family, standardize = TRUE))
  fit_coef_vec <- as.vector(fit_coef)
  
  if(lasso_family=="cox")
    fit_coef_vec <- fit_coef_vec
  else {
    fit_coef_vec <- fit_coef_vec[-1]   
  }
  return(fit_coef_vec)
}

# ekn_dvines()--> Function to derandomized knockoffs using e-values for FDR control. This function
# considers the dvine knockoff procedure.
#The code to implement this function is adapted from 
#https://github.com/zhimeir/derandomized_knockoffs_fdr

#Arguments:
#X: matrix of predictors
#X_cont: matrix of continuous predictors (after transformation of ordinal to conditional uniforms)
#column_type: vector with the type of each variable of X ("num" numeric, "ord" ordinal)
#y: vector or matrix of response
#dvine_distributions: a list containing: 
#   1) object of class vinecop_dist for X_X,
#   2) Object of class vinecop_dist for X_cont,
#   3) Pseudo observations of X_cont saved in u_X
#M: integer denoting the number of generated copies of the knockff matrix Xk.
#M_lasso: integer related to the number of runs for the stabilzation against CV
#alpha: integer indicating FDR target level
#gamma: integer denoting target level for the knockoff threshold. According to Ren & Barber (2023),
#       experimentally, gamma=alpha/2 works well.           
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial" 
#n_cores: int -> number of cores for parallel processing
#random_cv: logic value for a random assignation of the folds in the Cross-validation used in LASSO regression
#        TRUE is the default. FALSE is for reproducibility purposes. 
#random_evalues: logic value for a random sampling of the knockoffs in the derandomized knockoff procedure
#        TRUE is the default. FALSE is for reproducibility purposes. 


#Note: the knockoff.threshold() function from the R knockoff package is used for 
#setting the Knockoff rejection threshold (https://cran.r-project.org/web/packages/knockoff/knockoff.pdf)

#Value: This function returns a list with the selected variables of the procedure

ekn_dvines <- function(X, X_cont, column_type, y, dvine_distributions, M=50, M_lasso=10, alpha=0.2, 
                       gamma=0.1, lasso_family, n_cores=1, n_folds = 5, random_cv=TRUE, random_evalues=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null( X_cont )) {
    stop("Argument X_cont is missing")  
  }
  if (is.null( column_type )) {
    stop("Argument column_type is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  
  if (is.null( dvine_distributions )) {
    stop("Argument dvine_distributions is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  #Number of variables p and sample size n  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Initial matrix of E-values 
  E <- matrix(0, M, p)
  
  # Creating empty lists to store the X_Xk matrices
  ls_X_Xk <- list()
  
  for(m in 1:M){
    
    if(random_evalues==TRUE){
      #dvine Knockoffs sampling
      Xk <- create_dvine_Knockoffs(X=X, X_cont=X_cont, column_type=column_type, 
                                   dvine_distributions=dvine_distributions, 
                                   n_cores=n_cores,seed_Xk=NULL)
    }
    else {
      #dvine Knockoffs sampling
      Xk <- create_dvine_Knockoffs(X=X, X_cont=X_cont, column_type=column_type, 
                                   dvine_distributions=dvine_distributions, 
                                   n_cores=n_cores, seed_Xk=m)    
    }
    #X and Xk column binding
    X_Xk <- cbind(X, Xk)
    
    # Storing the data in lists
    ls_X_Xk <- append(ls_X_Xk, list(X_Xk))
    
  }
  
  start_time <- Sys.time()  
  Z <- foreach(m = 1:M,.packages=c("glmnet"),.export = "stable_lasso_glmnet") %dopar% {
    
    stable_lasso_glmnet(X=ls_X_Xk[[m]], y=y, lasso_family=lasso_family,  
                                   M_lasso=M_lasso, n_folds=n_folds, 
                                   random_cv=random_cv)
  }
  end_time <- Sys.time()
  print("Time for getting Z:")
  print(end_time-start_time)
  
  for(m in 1:M){
    #Importance statistic
    W <- abs(Z[[m]][1:p]) - abs(Z[[m]][(p+1):length(Z[[m]])])
    
    #Knockoff rejection threshold - conservative procedure ("knockoffs+" offset = 1)
    tau <- stop_early(W, gamma, offset=1) 
    
    #E-vales for all the variables (columns) for m run
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))
    
  }
  
  #Averaging the e-values to select set of discoveries
  E <- p*colMeans(E)
  rej <- ebh(E, alpha)$rej
  
  return(list(rej = rej, E = E)) 
  
}

# ekn_second_order()--> Function to derandomized knockoffs using e-values for FDR control. This function
# considers the second-order approximation of the knockoff R package
#The code to implement this function is adapted from 
#https://github.com/zhimeir/derandomized_knockoffs_fdr

#Arguments:
#X: matrix of predictors
#y: vector or matrix of response
#M: integer denoting the number of generated copies of the knockff matrix Xk.
#M_lasso: integer related to the number of runs for the stabilization against Cross Validation  (Roberts and Nowak, 2014)
#alpha: integer indicating FDR target level
#gamma: integer denoting target level for the knockoff threshold. According to Ren & Barber (2023),
#       experimentally, gamma=alpha/2 works well.           
#lasso_family: a string to select linear regression "gaussian" or logistic regression "binomial" 
#n_cores: int -> number of cores for parallel processing
#random_cv: logic value for a random assignation of the folds in the Cross-validation used in LASSO regression
#        TRUE is the default. FALSE is for reproducibility purposes. 
#random_evalues: logic value for a random sampling of the knockoffs in the derandomized knockoff procedure
#        TRUE is the default. FALSE is for reproducibility purposes. 


#Note: the knockoff.threshold() function from the R knockoff package is used for 
#setting the Knockoff rejection threshold (https://cran.r-project.org/web/packages/knockoff/knockoff.pdf)

#Value: This function returns a list with the selected variables of the procedure


ekn_second_order <- function(X, y, M=50, M_lasso=10, alpha=0.2, gamma=0.1, 
                             lasso_family, n_cores=1,n_folds = 5, random_cv=TRUE, random_evalues=TRUE){
  
  if (is.null(X)) {
    stop("Argument X is missing")  
  }
  if (is.null(y)) {
    stop("Argument y is missing")  
  }
  if (is.null( lasso_family )) {
    stop("Argument lasso_family is missing")  
  }
  
  
  #Number of variables p and sample size n  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Matrix of E-values 
  E <- matrix(0, M, p)
  
  # Creating empty lists to store the X_Xk matrices
  ls_X_Xk <- list()
  
  for(m in 1:M){
    
    if(random_evalues==TRUE){
      set.seed(NULL) 
    }
    else {  
      set.seed(m) #The seed is set for reproducibility.
    }
    #Gaussian Knockoffs copy selection from the list object
    Xk <- create.second_order(X)
    
    #X and Xk column binding
    X_Xk <- cbind(X, Xk)
    
    # Storing the data in lists
    ls_X_Xk <- append(ls_X_Xk, list(X_Xk))
    
  }
  
  start_time <- Sys.time()  
  Z <- foreach(m = 1:M,.packages=c("glmnet"),.export = "stable_lasso_glmnet") %dopar% {
    
    stable_lasso_glmnet(X=ls_X_Xk[[m]], y=y, lasso_family=lasso_family,  
                        M_lasso=M_lasso, n_folds=n_folds, 
                        random_cv=random_cv)
  }
  end_time <- Sys.time()
  print("Time for getting Z:")
  print(end_time-start_time)
  
  for(m in 1:M){
    #Importance statistic
    W <- abs(Z[[m]][1:p]) - abs(Z[[m]][(p+1):length(Z[[m]])])
    
    #Knockoff rejection threshold - conservative procedure ("knockoffs+" offset = 1)
    tau <- stop_early(W, gamma, offset=1) 
    
    #E-vales for all the variables (columns) for m run
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))
    
  }
  
  #Averaging the e-values to select set of discoveries
  E <- p*colMeans(E)
  rej <- ebh(E, alpha)$rej
  
  return(list(rej = rej, E = E)) 
  
}




##### Data loading and preprocessing  #####


#Importing the CSV file
Lung_data_initial <- read.csv("lung cancer data.csv")

#Encoding for ordinal variables
Lung_data_initial$Pat_Gender <- ifelse(Lung_data_initial$Pat_Gender == "M", 0, 1)
Lung_data_initial$Pat_Stage <- ifelse(Lung_data_initial$Pat_Stage== "I_II", 0, 1)
Lung_data_initial$Pat_Grade <- ifelse(Lung_data_initial$Pat_Grade== 1, 0, ifelse(Lung_data_initial$Pat_Grade== 2,1,2))
dim(Lung_data_initial)
#[1] 442 506

Lung_data_initial <- Lung_data_initial  %>% dplyr::select(-c(Pat_Grade))

#Removing NAs from the dataframe (4.88% of missinges)
Lung_data_complete <- Lung_data_initial[complete.cases(Lung_data_initial), ]
dim(Lung_data_complete)
#[1] 442 506

#Few observations have a survival time equal to zero, which causes a problem when fitting the Cox-lasso model.
#Therefore, these values are modified to have a survival time of 1 day (1/30)
Lung_data_complete$Pat_Overall_Survival_Months[Lung_data_complete$Pat_Overall_Survival_Months == 0] <- 1 / 30

#Filtering the most variable genes in terms of their variability (variance)
#Notes: 
#1) In the Lung_data dataframe, the gene-expression column order is from the most variable to the less variable, 
#in terms of their variance. 
#2) The Lung_data contains the most 500 expressed genes of a total of 12 259 genes of the
# genomic dataset from the research of Shedden et al. (2008) located in the Lung Cancer Explorer (LCE) database 
# http://lce.biohpc.swmed.edu/.

# Number of most expressed genes selected
gen_p <- 100

#A final data frame with the most expressed genes, clinical variables 
#and survival information. 
#The first 2 columns are related to survival, the following 4 are clinical, 
#and the remaining are the most expressed genes. 

Lung_data <- Lung_data_complete[, 1:(3 + 2 + gen_p)]
dim(Lung_data)
#[1] 442 106

# Calculating the censoring rate for 'Pat_Died'
Censoring <- round(100 * as.numeric(table(Lung_data$Pat_Died)[1] / length(Lung_data$Pat_Died)), 2)
# Print censoring rate
print(Censoring)

# Percentage of dead patients
round(100 * as.numeric(table(Lung_data$Pat_Died)[2] / length(Lung_data$Pat_Died)), 2)
#[1] 53.39

# 50% percentile (median) of Overall Survival
median(Lung_data$Pat_Overall_Survival_Months)
#[1] 47

#Extracting X and y
X <- Lung_data %>% dplyr::select(-c("Pat_Died", "Pat_Overall_Survival_Months")) %>% as.matrix()
y <- Lung_data %>% dplyr::select(c("Pat_Died", "Pat_Overall_Survival_Months"))
y_surv <- Surv(y$Pat_Overall_Survival_Months,y$Pat_Died)
p <- dim(X)[2]

#### Applying the Model-X knockoff methodology ####


##### Non-parametric DTDCKe #####

# Initial setup
n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel

alpha_dvine <- 0.3 # Target FDR
nonpar_vinecop_family <- 'nonparametric' # Class of Dvine fitting
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 5 # Number of folds for Lasso CV


# D-vine order
TSP_order <- get_dvine_order(X, random_order = TRUE)
X_dvine_order <- X[, TSP_order]


# Column identification
column_type <- as.vector(apply(X_dvine_order, MARGIN = 2, FUN = column_type_identification))
  
# Ordinal latent transformation
X_cont_dvine_order <- apply(X = X_dvine_order, MARGIN = 2, FUN = ordinal_to_uniform)

# Noparametric D-vine fitting
nonpar_dvine_distributions <- X_Xk_dvine_distributions(X_cont = X_cont_dvine_order,column_type, nonpar_vinecop_family, n_cores)
  
# Matrix conversion
X_dvine_order_matrix <- as.matrix(X_dvine_order)
  
# Application of the derandomized procedure using e-values
start_time <- Sys.time()
res_nonpar_dvines <- ekn_dvines(X = X_dvine_order_matrix, X_cont = X_cont_dvine_order, 
                             column_type = column_type, y = y_surv, 
                         nonpar_dvine_distributions, M = M, M_lasso = M_lasso,
                             alpha = alpha_dvine, gamma = alpha_dvine / 2, 
                             lasso_family = lasso_family, n_cores = n_cores,
                             n_folds=n_folds, 
                             random_cv = TRUE, random_evalues = TRUE)
end_time <- Sys.time()
print("Time for running the non-DTDCKe procedure:")
print(end_time-start_time)
stopCluster(cl) #Stopping cluster

# Vector of integers that indicates the selected non-null positions 
just_rejections_nonpar_dvines <- sort(res_nonpar_dvines$rej)
print(paste0("The number of selected variables is: ", length(just_rejections_nonpar_dvines)))
  

# The vector that indicates the rejections considering all variables (0 null, 1 non-null)
rejections_nonpar_dvines <- rep(0, p)
rejections_nonpar_dvines[just_rejections_nonpar_dvines] <- 1
  

# Filtering column names based on rejections equal to 1
selected_columns_nonpar_dvines <- colnames(X_dvine_order)[rejections_nonpar_dvines == 1]
print(selected_columns_nonpar_dvines)


#### Parametric DTDCKe ####

# Initial setup
n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel

alpha_par_dvine <- 0.3 # Target FDR
par_vinecop_family <- 'parametric' # Class of Dvine fitting
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 5 # Number of folds for Lasso CV

# D-vine order
TSP_order <- get_dvine_order(X, random_order = TRUE)
X_dvine_order <- X[, TSP_order]

# Column identification
column_type <- as.vector(apply(X_dvine_order, MARGIN = 2, FUN = column_type_identification))

# Ordinal latent transformation
X_cont_dvine_order <- apply(X = X_dvine_order, MARGIN = 2, FUN = ordinal_to_uniform)

# Parametric D-vine fitting
par_dvine_distributions <- X_Xk_dvine_distributions(X_cont = X_cont_dvine_order,column_type, par_vinecop_family, n_cores)

# Matrix conversion
X_dvine_order_matrix <- as.matrix(X_dvine_order)

# Application of the derandomized procedure using e-values
start_time <- Sys.time()
res_par_dvines <- ekn_dvines(X = X_dvine_order_matrix, X_cont = X_cont_dvine_order, 
                         column_type = column_type, y = y_surv, 
                         par_dvine_distributions, M = M, M_lasso = M_lasso,
                         alpha = alpha_par_dvine, gamma = alpha_dvine / 2, 
                         lasso_family = lasso_family, n_cores = n_cores,
                         n_folds=n_folds, 
                         random_cv = TRUE, random_evalues = TRUE)
end_time <- Sys.time()
print("Time for running the par-DTDCKe procedure:")
print(end_time-start_time)
stopCluster(cl) #Stopping cluster

#Vector of integers that indicates the selected non-null positions 
just_rejections_par_dvines <- sort(res_par_dvines$rej)
print(paste0("The number of selected variables is: ", length(just_rejections_par_dvines)))


#Vector that indicates the rejections considering all variables (0 nulls, 1 non-null)
rejections_par_dvines <- rep(0, p)
rejections_par_dvines[just_rejections_par_dvines] <- 1


# Filtering column names based on rejections equal to 1
selected_columns_par_dvines <- colnames(X_dvine_order)[rejections_par_dvines == 1]
print(selected_columns_par_dvines)


#### Model-X knockoff second order ####
# Initial setup

n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel

alpha_second_order <- 0.5 # Target FDR
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 5 # Number of folds for Lasso CV

start_time <- Sys.time()
res_second_order <- ekn_second_order(X = X, y = y_surv, 
                                         M = M, M_lasso = M_lasso, alpha = alpha_second_order, 
                                         gamma = alpha_second_order / 2, 
                                         lasso_family = lasso_family,
                                         n_folds = n_folds,
                                         n_cores = n_cores, random_cv = FALSE, 
                                         random_evalues = FALSE)
end_time <- Sys.time()
print("Time for running the second-order Model-X:")
print(end_time-start_time)
stopCluster(cl) #Stopping cluster

# Vector of integers that indicates the selected non-null positions 
just_rejections_second_order <- sort(res_second_order$rej)
print(paste0("The number of selected variables is: ", length(just_rejections_second_order)))

#Vector that indicates the rejections considering all variables (0 null, 1 non-null)
rejections_second_order <- rep(0, p)
rejections_second_order[just_rejections_second_order] <- 1

# Filtering column names based on rejections equal to 1
selected_columns_second_order <- colnames(X)[rejections_second_order == 1]
print(selected_columns_second_order)

 
