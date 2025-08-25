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
library(MASS)
library(sn)

##### Auxiliary functions #####

#Function to find the upper limit of a uniform(0, u_max) distribution that permits the 
# establish of an approximate censoring percentage given by target_censoring
find_u_max <- function(target_censoring=10.0, t, tol = 0.1, max_iter = 1000) {
  
  if (is.null(t)) {
    stop("Argument t is missing")  
  }
  
  # Initialize u_max range
  u_min <- 0
  u_max <- 2*max(t) 
  iter <- 1
  
  # Guess Max as the midpoint of u_min and u_max
  u_mid <- (u_min + u_max) / 2
  
  for (i in 1:max_iter){
    
    # Simulate censored times with the current u_mid
    t_cens <- runif(length(t), min = 0, max = u_mid)
    
    # Calculate censoring indicator
    I_cens <- ifelse(t <= t_cens, 1, 0)
    
    # Calculate current censoring percentage
    censoring_percentage <- (1 - mean(I_cens)) * 100
    
    # Check if the current censoring percentage is close to the target
    if ( abs(censoring_percentage - target_censoring) < tol) {
      break
    }
    
    # Adjust the range for u_max based on the censoring percentage
    if (censoring_percentage > target_censoring) {
      u_mid  <- u_mid + u_mid/2  # Increase u_mid to decrease censoring
    } else {
      u_mid <- u_mid/2  # Decrease u_mid to increase censoring
    }
    iter <- iter + 1
  }
  return(u_mid) 
}

#A function that binns the continuous variables of a matrix X by 
#a specific quantile corresponding to the given probability prob
Continuous_to_binary_by_quantile <- function(X, column_type, prob){
  
  
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

#Ordinal to uniform transfomation
ordinal_to_uniform <- function(column){
  
  if (is.null(column)) {
    stop("Argument column is missing")  
  }
  
  # A discrete variable with more than 10 levels is considered a numerical variable
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
  
  # A discrete variable with more than 10 levels is considered a numerical variable
  if(length(unique(column)) < 11) 
  {type<-"ord"}
  else
  {type<-"num"} 
  return(type)
}    

#Uniform to original transformation
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


#create_sparse_coefficients()-->Function to Generate a set of sparse coefficients for single index or sparse additive models.
#The code to implement this function is adapted from:
#https://amspector100.github.io/knockpy/_modules/knockpy/dgp.html#create_sparse_coefficients

#Arguments:
#p : int (Dimensionality of coefficients )
#sparsity : num (Proportion of non-null coefficients
#coeff_size: num (Size of non-zero coefficients)
#sign_prob : num (The probability that each nonzero coefficient will be positive)
#iid_signs : logical (If TRUE, the signs of the coeffs are assigned independently, 
#           Else, exactly sign_prob*sparsity*p coefficients will be positive)
#coeff_dist : char Specifies the distribution of nonzero coefficients. Three options:
#- "None" : all coefficients have absolute value coeff_size.
#- "normal": all nonzero coefficients are drawn from Normal(coeff_size, 1). 
#- "uniform": nonzero coeffs are drawn from Unif(coeff_size/2, coeff_size).
#corr_signals : logical (If true, all of the nonzero coefficients will lie in a consecutive block.)

create_sparse_coefficients <- function(p,
                                       sparsity = 0.5,
                                       coeff_size = 1,
                                       coeff_dist = NULL,
                                       sign_prob = 0.5,
                                       iid_signs = TRUE,
                                       corr_signals = FALSE) {
  
  
  # A certain percentage of coefficients are nonzero
  num_nonzero <- floor(sparsity * p)
  if (corr_signals) {
    nonnull_start <- sample(0:(p - num_nonzero), 1)
    beta <- rep(0, p)
    beta[nonnull_start:(nonnull_start + num_nonzero - 1)] <- coeff_size
  } 
  else {
    beta <- c(rep(coeff_size, num_nonzero), rep(0, p - num_nonzero))
    beta <- sample(beta)
  }
  
  beta_nonzeros <- beta != 0
  num_nonzero <- sum(beta_nonzeros)
  
  # Now draw random signs
  if (iid_signs) {
    signs <- 1 - 2 * rbinom(p, 1, sign_prob)
    beta <- beta * signs
  } else {
    num_pos <- floor(num_nonzero * sign_prob)
    signs <- c(rep(1, num_pos), rep(-1, num_nonzero - num_pos))
    signs <- sample(signs)
    beta[beta_nonzeros] <- beta[beta_nonzeros] * signs
  }
  
  # Possibly change the absolute values of beta
  if (!is.null(coeff_dist)) {
    if (tolower(coeff_dist) == "normal") {
      beta <- (sqrt(coeff_size) * rnorm(p)) * beta_nonzeros
    } else if (tolower(coeff_dist) == "uniform") {
      beta <- beta * runif(p) / 2 + beta / 2
    }  else if (tolower(coeff_dist) == "none") {
      # Do nothing
    } else {
      stop("coeff_dist (", coeff_dist, ") must be 'none', 'normal', or 'uniform'")
    }
  }
  
  return(beta)
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

#Value: an integer vector with the new indices

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
      
      # Print the current iteration and the worker id
      worker_id <- Sys.getpid()  # Get the current process ID
      print(paste("Worker", worker_id, "is processing iteration", i))
      
      cvfit$lambda.min
    }
    end_time <- Sys.time()
    print("Time for running parallel Stabilizing the lasso:")
    print(end_time-start_time)
    print(lambdas)
  }

  #Selecting the median of the lambdas distribution
  lambda50 <- as.numeric(quantile(lambdas,probs=0.5))
  print(paste("lambda50 :", lambda50))
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
#y: vector of response (In the context of survival regression, y needs to be a survival object from the Survival package)
#lasso_family: a string to select linear regression "gaussian", logistic regression "binomial" 
#               or "cox" for survival cox regression model
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
  else { #This is the case for the survival regression model
         # y should be a survival object
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







##### Simulations #####

#Number of simulations
n_sim <- 100

#Censored limits
vec_censoring <- c(10.0, 20.0, 30.0, 40.0, 50.0)

#Number of covariables and sample size
p <- 100
n <- 200

#Skew-t distribution
alpha_level <- 2
nu_level <- 5

#Number of elements by block
b_by_block <- 10 

#Correlation AR1(by blocks)
rho <- 0.7

#Block Correlation matrix
n_blocks <- p%/%b_by_block
corMat <- diag(n_blocks) %x% matrix(rho^abs(matrix(1:b_by_block,b_by_block, b_by_block, byrow = TRUE) - matrix(1:b_by_block, b_by_block, b_by_block)), b_by_block, b_by_block)
diag(corMat) <- 1

#Sparsity of the non-null variables (proportion of non-nulls)
sp <- 0.2 

#Proportion of ordinals
prop_ord <- 0.2

#Proportion of zeros
prob <- 0.5

#Beta coefficients magnitude
amplitude <- 10
beta_factor <- amplitude / sqrt(n)

#Survival Weibull distribution parameters (scale --> sigma, shape --> nu, lambda=1/(scale^shape))
lambda_T = 0.05 # lamda=1/(scale^shape)=1/(sigma^nu)
nu_T = 1.5 # With nu_T=1(shape_T=1) we have the exponential distribution 

#X column names
X_column_names <- sapply(1:p, function(number) paste0("X", number))

#Beta column names
beta_column_names <- sapply(1:p, function(number) paste0("Beta", number))

#Column names of the data frame with all the simulations
df_column_names <- c("Censoring","Iteration", "y", X_column_names)
df_column_names_betas <- c("Censoring","Iteration", beta_column_names)

# Creating empty lists to store the simulations
ls_simulations <- list()
ls_simulations_betas <- list()

# Creating empty data.frame to store the data censoring percentage
df_simulations_censoring  <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df_simulations_censoring) <- c("Censoring","Iteration", "Data_censoring")

#A seed is set for reproducibility purposes
set.seed(123)
for ( censoring in vec_censoring ){
  for ( i in 1:n_sim){
    
    #Number of ordinal and continuous variables
    p_ord <- round(p*prop_ord)
    p_num <- p - p_ord
    
    #variable types for simulations
    ord <- rep("bin", p_ord)
    num <- rep("num", p_num)
    types <- sample(c(ord, num))
    
    #Multivariate normal distribution realizations
    X <- mvrnorm(n, mu=rep(0,p), Sigma = corMat)
    
    for(j in 1:p) {   if(types[j]=="num")
    {tryCatch({X[,j] <-as.vector(qst(pnorm(X[,j]),xi=0, omega=1,alpha=alpha_level,nu=nu_level))}, 
              error=function(e){X[,j] <-as.vector(qst(pnorm(X[,j]),xi=0, omega=1,alpha=alpha_level, nu=nu_level, solver="RFB"))}) }}
    
    #Binning the ordinal variables
    X_mixed <- Continuous_to_binary_by_quantile(X, column_type=types, prob=prob)
    
    colnames(X_mixed ) <- X_column_names
    
    # Creating random sparse coefficients with uniform distribution Unif(coeff_size/2, coeff_size).
    beta <- create_sparse_coefficients(p=p, sparsity=sp, 
                                       sign_prob=0.5, # probability of positive nonzero coefficient
                                       coeff_size=beta_factor, 
                                       coeff_dist='uniform')
    
    matrix_beta <- matrix(beta, nrow=1,ncol=p)
    colnames(matrix_beta) <- beta_column_names
    
    # Survival time simulations (Bender et al. 2006)
    t <- (-log(runif(n, min = 0, max = 1)) / (lambda_T * exp(X_mixed %*% beta)))^(1 / nu_T)
    
    u_max <- find_u_max(target_censoring=censoring, t, tol = 1, max_iter = 1000) 
    
    # Censored time and censored indicator
    t_cens <- runif(n, min = 0, max = u_max)
    I_cens <- ifelse(t <= t_cens, 1, 0)
    
    # Observed time
    t_obs <- pmin(t, t_cens)
    
    #Censoring's percentage
    data_censoring <- (1- sum(I_cens)/n)*100
    cat(sprintf("Censoring's percentange %.2f%%\n", data_censoring))
    
    # Response variable
    y_surv <- Surv(t_obs, I_cens)
    
    #Data frame with the information of this iteration
    data <- data.frame(Censoring=rep(censoring ,n), Iteration = rep(i,n),y_surv=y_surv, X_mixed)
    data_betas <- data.frame(Censoring=censoring , Iteration = i, matrix_beta )
    data_censoring <- data.frame(Censoring=censoring,Iteration = i, Data_censoring=data_censoring)
    
    # Storing the data in lists
    ls_simulations <- append(ls_simulations, list(data))
    ls_simulations_betas <- append(ls_simulations_betas, list(data_betas))
    
    # Merge the dataframes by stacking rows
    df_simulations_censoring <- rbind(df_simulations_censoring , data_censoring )
    
    
  }  
}

#df_simulations_censoring %>%  filter(Censoring==30.0) %>% summarise(average=mean(Data_censoring))

#kmfit <- survfit(y_surv~1, data= ls_simulations[[5]] )
#plot(kmfit, lty = c("solid", "dashed"), col = c("black", "grey"), xlab = "Survival Time In Days", ylab = "Survival Probabilities")

#### Non-parametric DTDCKe ####

# Initial setup
n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel
alpha_dvine <- 0.2 # Target FDR
vinecop_family <- 'nonparametric' # Class of Dvine fitting
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 3 # Number of folds for Lasso CV

# Number of simulations to iterate through
n_simulations <- length(ls_simulations)

# Data frame to store the power and FDP of the procedure
df_nondvine_results <- data.frame(matrix(ncol = 4, nrow = 0))
df_column_names_dvine_results <- c("Censoring", "Iteration", "Dvine_power", "Dvine_FDP")
colnames(df_nondvine_results) <- df_column_names_dvine_results

start_time <- Sys.time()

for (idx in 1:n_simulations) {
  # Extracting the data frame from the list
  df_sim <- ls_simulations[[idx]]
  df_beta_sim <- ls_simulations_betas[[idx]]
  
  # Censoring and iteration
  value <- unique(df_sim$Censoring)
  iteration <- unique(df_sim$Iteration)
  
  # Getting X, y, and beta for each simulation
  X <- as.matrix(df_sim[, 4:ncol(df_sim)])  # From column 4 to the final column
  y <- df_sim[, 3]              # Third column is y
  beta <- as.numeric(df_beta_sim[, 3:ncol(df_beta_sim)])  # From column 3 to final column for betas
  
  # D-vine order
  TSP_order <- get_dvine_order(X, random_order = FALSE)
  X_dvine_order <- X[, TSP_order]
  
  # Column identification
  column_type <- as.vector(apply(X_dvine_order, MARGIN = 2, FUN = column_type_identification))
  
  # Ordinal transformation
  X_cont_dvine_order <- apply(X = X_dvine_order, MARGIN = 2, FUN = ordinal_to_uniform)
  
  # D-vine fitting
  dvine_distributions <- X_Xk_dvine_distributions(X_cont = X_cont_dvine_order,column_type, vinecop_family, n_cores)
  
  # Matrix conversion
  X_dvine_order_matrix <- as.matrix(X_dvine_order)
  
  # Application of the derandomized procedure using e-values
  time_ekn_dvines <- system.time(
    res_dvines <- ekn_dvines(X = X_dvine_order_matrix, X_cont = X_cont_dvine_order, 
                             column_type = column_type, y = y, 
                             dvine_distributions, M = M, M_lasso = M_lasso,
                             alpha = alpha_dvine, gamma = alpha_dvine / 2, 
                             lasso_family = lasso_family, n_cores = n_cores,
                             n_folds=n_folds, 
                             random_cv = FALSE, random_evalues = FALSE)
  )
  print("Time for running the dvine e-kn procedure:")
  print(time_ekn_dvines)
  
  # Vector of integers that indicates the selected non-nulls position 
  just_rejections_dvines <- sort(res_dvines$rej)
  print(paste0("The number of selected variables is: ", length(just_rejections_dvines)))
  
  # Vector that indicates the rejections considering all variables (0 null, 1 non-null)
  rejections_dvines <- rep(0, p)
  rejections_dvines[just_rejections_dvines] <- 1
  
  # Adjusting the order of the beta vector to match the dvine order
  beta_dvine_order <- beta[TSP_order]
  
  # Power and FDP
  dvine_Power <- round(100 * sum(rejections_dvines * (beta_dvine_order != 0)) / sum(beta_dvine_order != 0), 2)
  dvine_FDP <- round(100 * sum(rejections_dvines * (beta_dvine_order == 0)) / max(1, sum(rejections_dvines)), 2)
  
  cat(sprintf("The knockoff DTDCke filter POWER %.2f%% with an FDP of %.2f%%\n", dvine_Power, dvine_FDP))
  
  # Data frame with the information of this iteration
  data <- data.frame(Censoring = value, Iteration = iteration, Dvine_power = dvine_Power, Dvine_FDP = dvine_FDP)
  
  # Merge the dataframes by stacking rows
  df_nondvine_results <- rbind(df_nondvine_results, data)
  
  cat(sprintf("The iteration is %d with censoring rate of %.2f \n", iteration, value))
}

end_time <- Sys.time()

stopCluster(cl) #Stopping cluster
print("Time for running the non-DTDCKe procedure for the simulations:")
nondvine_time <- end_time-start_time
print(nondvine_time)

#Exporting the results to csv file
write.csv(df_nondvine_results, file = "df_nondvine_results_Survival_nonGaussian_X.csv", row.names = FALSE)

results_nondvine <- df_nondvine_results%>% group_by(Censoring) %>%
  summarise(meanPower = mean(Dvine_power),meanFDP = mean(Dvine_FDP) )
View(results_nondvine)

#### Parametric DTDCKe ####

# Initial setup
n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel
alpha_dvine <- 0.2 # Target FDR
vinecop_family <- 'parametric' # Class of Dvine fitting
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 3 # Number of folds for Lasso CV

# Number of simulations to iterate through
n_simulations <- length(ls_simulations)

# Data frame to store the power and FDP of the procedure
df_dvine_results <- data.frame(matrix(ncol = 4, nrow = 0))
df_column_names_dvine_results <- c("Censoring", "Iteration", "Dvine_power", "Dvine_FDP")
colnames(df_dvine_results) <- df_column_names_dvine_results


start_time <- Sys.time()

for (idx in 1:n_simulations) {
  # Extracting the data frame from the list
  df_sim <- ls_simulations[[idx]]
  df_beta_sim <- ls_simulations_betas[[idx]]
  
  # Censoring and iteration
  value <- unique(df_sim$Censoring)
  iteration <- unique(df_sim$Iteration)
  
  # Getting X, y, and beta for each simulation
  X <- as.matrix(df_sim[, 4:ncol(df_sim)])  # From column 4th to the final column
  y <- df_sim[, 3]             # Third column is y
  beta <- as.numeric(df_beta_sim[, 3:ncol(df_beta_sim)])  # From column 3rd to final column for betas
  
  # D-vine order
  TSP_order <- get_dvine_order(X, random_order = FALSE)
  X_dvine_order <- X[, TSP_order]
  
  # Column identification
  column_type <- as.vector(apply(X_dvine_order, MARGIN = 2, FUN = column_type_identification))
  
  # Ordinal transformation
  X_cont_dvine_order <- apply(X = X_dvine_order, MARGIN = 2, FUN = ordinal_to_uniform)
  
  # D-vine fitting
  dvine_distributions <- X_Xk_dvine_distributions(X_cont = X_cont_dvine_order,column_type, vinecop_family, n_cores)
  
  # Matrix conversion
  X_dvine_order_matrix <- as.matrix(X_dvine_order)

  
  # Application of the derandomized procedure using e-values
  time_ekn_dvines <- system.time(
    res_dvines <- ekn_dvines(X = X_dvine_order_matrix, X_cont = X_cont_dvine_order, 
                             column_type = column_type, y = y, 
                             dvine_distributions, M = M, M_lasso = M_lasso,
                             alpha = alpha_dvine, gamma = alpha_dvine / 2, 
                             lasso_family = lasso_family, n_cores = n_cores, 
                             random_cv = FALSE, random_evalues = FALSE)
  )
  print("Time for running the dvine e-kn procedure:")
  print(time_ekn_dvines)
  
  # Vector of integers that indicates the selected non-nulls position 
  just_rejections_dvines <- sort(res_dvines$rej)
  print(paste0("The number of selected variables is: ", length(just_rejections_dvines)))
  
  # Vector that indicates the rejections considering all variables (0 null, 1 non-null)
  rejections_dvines <- rep(0, p)
  rejections_dvines[just_rejections_dvines] <- 1
  
  # Adjusting the order of the beta vector to match the dvine order
  beta_dvine_order <- beta[TSP_order]
  
  # Power and FDP
  dvine_Power <- round(100 * sum(rejections_dvines * (beta_dvine_order != 0)) / sum(beta_dvine_order != 0), 2)
  dvine_FDP <- round(100 * sum(rejections_dvines * (beta_dvine_order == 0)) / max(1, sum(rejections_dvines)), 2)
  
  cat(sprintf("The knockoff DTDCke filter POWER %.2f%% with an FDP of %.2f%%\n", dvine_Power, dvine_FDP))
  
  # Data frame with the information of this iteration
  data <- data.frame(Censoring = value, Iteration = iteration, Dvine_power = dvine_Power, Dvine_FDP = dvine_FDP)
  
  # Merge the dataframes by stacking rows
  df_dvine_results <- rbind(df_dvine_results, data)
  
  cat(sprintf("The iteration is %d with censoring rate of %.2f \n", iteration, value))
}

end_time <- Sys.time()

stopCluster(cl) #Stopping cluster
print("Time for running the DTDCKe procedure for the simulations:")
dvine_time <- end_time-start_time
print(dvine_time)

#Exporting the results to csv file
write.csv(df_dvine_results, file = "df_dvine_results_Survival_nonGaussian_X.csv", row.names = FALSE)

dvine_results <- df_dvine_results %>% group_by(Censoring) %>%
  summarise(meanPower = mean(Dvine_power),meanFDP = mean(Dvine_FDP) )
View(dvine_results)


#### Model-X knockoff second order ####
# Initial setup

n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel
alpha_second_order <- 0.2 # Target FDR
M <- 50
M_lasso <- 10
lasso_family <- 'cox'
n_folds <- 3 # Number of folds for Lasso CV

# Number of simulations to iterate through
n_simulations <- length(ls_simulations)

# Data frame to store the power and FDP of the procedure
df_results <- data.frame(matrix(ncol = 4, nrow = 0))
df_column_names_results <- c("Censoring", "Iteration", "Power", "FDP")
colnames(df_results) <- df_column_names_results

start_time <- Sys.time()

for (idx in 1:n_simulations) {
  
  # Extracting the data frame from the list
  df_sim <- ls_simulations[[idx]]
  df_beta_sim <- ls_simulations_betas[[idx]]
  
  # Censoring and iteration
  value <- unique(df_sim$Censoring)
  iteration <- unique(df_sim$Iteration)
  
  # Getting X, y, and beta for each simulation
  X <- as.matrix(df_sim[, 4:ncol(df_sim)])  # From column 4 to the final column
  y <- df_sim[, 3]              # Third column is y
  beta <- as.numeric(df_beta_sim[, 3:ncol(df_beta_sim)])  # From column 3 to final column for betas
  
  # Application of the derandomized procedure using e-values
  time_ekn_second_order <- system.time(
    res_second_order <- ekn_second_order(X = X, y = y, 
                                         M = M, M_lasso = M_lasso, alpha = alpha_second_order, 
                                         gamma = alpha_second_order / 2, lasso_family = lasso_family, 
                                         n_cores = n_cores, random_cv = FALSE, 
                                         random_evalues = FALSE)
  )
  print("Time for running the e-value procedure:")
  print(time_ekn_second_order)
  
  # Vector of integers that indicates the selected non-nulls position 
  just_rejections <- sort(res_second_order$rej)
  print(paste0("The number of selected variables is: ", length(just_rejections)))
  
  # Vector that indicates the rejections considering all the variables (0 null, 1 non-null)
  rejections <- rep(0, p)
  rejections[just_rejections] <- 1
  
  # Power and FDP
  Power <- round(100 * sum(rejections * (beta != 0)) / sum(beta != 0), 2)
  FDP <- round(100 * sum(rejections * (beta == 0)) / max(1, sum(rejections)), 2)
  cat(sprintf("The procedure POWER is %.2f%% with an FDP of %.2f%%\n", Power, FDP))
  
  # Data frame with the information of this iteration
  data <- data.frame(Censoring = value, Iteration = iteration, Power = Power, FDP = FDP)
  
  # Merging the dataframes by stacking rows
  df_results <- rbind(df_results, data)
  
  cat(sprintf("The iteration is %d with censoring rate of %.2f \n", iteration, value))
}

end_time <- Sys.time()
stopCluster(cl) # Stopping cluster

print("Time for running all the procedure for the simulations:")
second_order_time <- end_time - start_time
print(second_order_time)

df_second_order_results <- df_results

#Exporting the results to csv file
write.csv(df_second_order_results, file = "df_second_order_results_Survival_nonGaussian_X.csv", row.names = FALSE)

second_order_results <- df_second_order_results %>% group_by(Censoring) %>%
  summarise(meanPower = mean(Power),meanFDP = mean(FDP) )
View(second_order_results)


#### Lasso ####

n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel
lasso_family <- 'cox'
M_lasso <- 10
random_cv <- FALSE
n_folds <- 5 # Number of folds for Lasso CV

# Number of simulations to iterate through
n_simulations <- length(ls_simulations)

# Data frame to store the power and FDP of the procedure
df_results <- data.frame(matrix(ncol = 4, nrow = 0))
df_column_names_results <- c("Ordinals_Proportion", "Iteration", "Power", "FDP")
colnames(df_results) <- df_column_names_results

start_time <- Sys.time()

for (idx in 1:n_simulations) {
  
  # Extracting the data frame from the list
  df_sim <- ls_simulations[[idx]]
  df_beta_sim <- ls_simulations_betas[[idx]]
  
  # Censoring and iteration
  value <- unique(df_sim$Censoring)
  iteration <- unique(df_sim$Iteration)
  
  # Getting X, y, and beta for each simulation
  X <- as.matrix(df_sim[, 4:ncol(df_sim)])  # From column 4th to the final column
  y <- df_sim[, 3]             # Third column is y
  beta <- as.numeric(df_beta_sim[, 3:ncol(df_beta_sim)])  # From column 3rd to the final column for betas
  
  #Estimated lasso coefficients for X
  Z <- stable_lasso_glmnet(X=X, y=y, lasso_family=lasso_family,  
                           M_lasso=M_lasso, n_folds=n_folds, 
                           random_cv=random_cv)
  
  # Vector of integers that indicates the selected non-nulls position 
  just_rejections_lasso <- which(abs(Z)>0)
  just_rejections <- sort(just_rejections_lasso)
  
  print(paste0("The number of selected variables is: ", length(just_rejections)))
  
  # Vector that indicates the rejections considering all the variables (0 null, 1 non-null)
  rejections <- rep(0, ncol(X))
  rejections[just_rejections] <- 1
  
  # Power and FDP
  Power <- round(100 * sum(rejections * (beta != 0)) / sum(beta != 0), 2)
  FDP <- round(100 * sum(rejections * (beta == 0)) / max(1, sum(rejections)), 2)
  cat(sprintf("The procedure POWER is %.2f%% with an FDP of %.2f%%\n", Power, FDP))
  
  
  # Data frame with the information of this iteration
  data <- data.frame(Censoring = value, Iteration = iteration, Power = Power, FDP = FDP)
 
  # Merging the dataframes by stacking rows
  df_results <- rbind(df_results, data)
  
  cat(sprintf("The iteration is %d with censoring rate of %.2f \n", iteration, value))
}


end_time <- Sys.time()

stopCluster(cl) # Stopping cluster

print("Time for running all the procedure for the simulations:")
lasso_time <- end_time - start_time
print(lasso_time)

df_lasso_results <- df_results

#Exporting the results to csv file
write.csv(df_lasso_results, file = "df_lasso_results_Survival_nonGaussian_X.csv", row.names = FALSE)

results_lasso <- df_lasso_results %>% group_by(Censoring) %>%
  summarise(meanPower = mean(Power),meanFDP = mean(FDP) )
View(results_lasso)




#### SCAD ####

n_cores <- detectCores() - 1 # Parallel computing
cl <- makeCluster(n_cores)
registerDoParallel(cl) # Registering a parallel backend for doParallel
penalty_type <- "SCAD"
n_folds <- 10 # Number of folds for CV

# Number of simulations to iterate through
n_simulations <- length(ls_simulations)

# Data frame to store the power and FDP of the procedure
df_results <- data.frame(matrix(ncol = 4, nrow = 0))
df_column_names_results <- c("Ordinals_Proportion", "Iteration", "Power", "FDP")
colnames(df_results) <- df_column_names_results

start_time <- Sys.time()

for (idx in 1:n_simulations) {
  
  # Extracting the data frame from the list
  df_sim <- ls_simulations[[idx]]
  df_beta_sim <- ls_simulations_betas[[idx]]
  
  # Censoring and iteration
  value <- unique(df_sim$Censoring)
  iteration <- unique(df_sim$Iteration)
  
  # Getting X, y, and beta for each simulation
  X <- as.matrix(df_sim[, 4:ncol(df_sim)])  # From column 4th to the final column
  y <- df_sim[, 3]             # Third column is y
  beta <- as.numeric(df_beta_sim[, 3:ncol(df_beta_sim)])  # From column 3rd to the final column for betas
  
  #Estimated coefficients
  cvfit <- cv.ncvsurv(X=X, y=y, penalty = penalty_type, cluster=cl, nfolds=n_folds)
  fit <- ncvsurv(X=X, y=y, penalty = penalty_type, cluster=cl, lambda=cvfit$lambda.min)
  Z <- coef(fit) #Thereis not intercept to exclude

  # Vector of integers that indicates the selected non-nulls position 
  just_rejections_SCAD <- which(abs(Z)>0)
  just_rejections <- sort(just_rejections_SCAD)
  
  print(paste0("The number of selected variables is: ", length(just_rejections)))
  
  # Vector that indicates the rejections considering all the variables (0 null, 1 non-null)
  rejections <- rep(0, ncol(X))
  rejections[just_rejections] <- 1
  
  # Power and FDP
  Power <- round(100 * sum(rejections * (beta != 0)) / sum(beta != 0), 2)
  FDP <- round(100 * sum(rejections * (beta == 0)) / max(1, sum(rejections)), 2)
  cat(sprintf("The procedure POWER is %.2f%% with an FDP of %.2f%%\n", Power, FDP))
  
  
  # Data frame with the information of this iteration
  data <- data.frame(Censoring = value, Iteration = iteration, Power = Power, FDP = FDP)
  
  # Merging the dataframes by stacking rows
  df_results <- rbind(df_results, data)
  
  cat(sprintf("The iteration is %d with censoring rate of %.2f \n", iteration, value))
}


end_time <- Sys.time()

stopCluster(cl) # Stopping cluster

print("Time for running all the procedure for the simulations:")
SCAD_time <- end_time - start_time
print(SCAD_time)

df_SCAD_results <- df_results

#Exporting the results to csv file
write.csv(df_SCAD_results, file = "df_SCAD_results_Survival_nonGaussian_X.csv", row.names = FALSE)

results_SCAD <- df_SCAD_results %>% group_by(Censoring) %>%
  summarise(meanPower = mean(Power),meanFDP = mean(FDP) )
View(results_SCAD)


