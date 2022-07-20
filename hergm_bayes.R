# hergm_bayes
# Author: Johan Koskinen
# Email: jkoskinen@unimelb.edu.au
# Description:
# This script, that can be sourced using
#     source("https://raw.githubusercontent.com/johankoskinen/BayesERGM/main/hergm_bayes.R")
# contains functions for fitting a Bayesian hierarchical ERGMs to repeated obeservations on networks
# The main workhorses are 
# A. h.ergm
# B. h.ergm.slow
# Both functions take a standard stanet formula object, but
# the networks need to be formatted as a list of networks, and
# covariate networks need to be formatted using get.covnets

# EXAMPLE USAGE
# In Agneessens, Trincado-Munoz, Koskinen, the model with the main effect of expertise was specified as
###=== define model terms ====
# you define the formula for one of the networks
# form <- net ~
#   #1. Structure
#   edges+
#   mutual+
#   gwodegree(decay=.2,fixed = T)+
#   gwidegree(decay=.2,fixed = T)+
#   twopath+
#   dgwesp(type="OTP",decay =.69,fixed = T)+
#   dgwesp(type="ITP",decay =.69,fixed = T)+
#   #2. Nodal attributes
#   nodeicov('AgeS')+
#   nodeocov('AgeS')+
#   nodeicov('Gend1FS')+
#   nodeocov('Gend1FS')+
#   nodematch('Gend1FS')+
#   nodeicov('ExpS')+
#   nodeocov('ExpS')
# 
# ###=== define fixed effects ====
# fixed <- c(
#   TRUE,#edges+
#   FALSE,# mutual+
#   FALSE,# gwodegree(decay=.2,fixed = T)+
#   FALSE,# gwidegree(decay=.2,fixed = T)+
#   FALSE,# twopath+
#   FALSE,# dgwesp(type="OTP",decay =.69,fixed = T)+
#   FALSE,# dgwesp(type="ITP",decay =.69,fixed = T)+
#   TRUE,# nodeicov('AgeS')+
#   TRUE,# nodeocov('AgeS')+
#   TRUE,# nodeicov('Gend1FS')+
#   TRUE,# nodeocov('Gend1FS')+
#   TRUE,# nodematch('Gend1FS')+
#   TRUE,# nodeicov('ExpS')+
#   TRUE)# nodeicov('ExpS')
# 
# ###=== define random effects ====
# isalpha <- fixed==FALSE
# # a random effect theta_k for group k is drawn from the distribution:
# # theta_k = alpha_k ~ MVN(mu, Sigma)
# ###=== define network effects ====
# isbeta <- c(
#   FALSE,#edges+
#   FALSE,# mutual+
#   FALSE,# gwodegree(decay=.2,fixed = T)+
#   FALSE,# gwidegree(decay=.2,fixed = T)+
#   FALSE,# twopath+
#   FALSE,# dgwesp(type="OTP",decay =.69,fixed = T)+
#   FALSE,# dgwesp(type="ITP",decay =.69,fixed = T)+
#   FALSE,# nodeicov('AgeS')+
#   FALSE,# nodeocov('AgeS')+
#   FALSE,# nodeicov('Gend1FS')+
#   FALSE,# nodeocov('Gend1FS')+
#   FALSE,# nodematch('Gend1FS')+
#   FALSE,# nodeicov('ExpS')+
#   FALSE)# nodeicov('ExpS')
# # a network effect theta_k for group k is
# theta_k = alpha_k + beta*W_k,
# where
# alpha_k ~ MVN(mu, Sigma)
# is ALWAYS a random effect


get.covnets <- function(DyadCovs,directed=TRUE)
{
	require('sna')
	require('network')
	require('ergm')


NumCovNets <- length(DyadCovs)
n <- dim(DyadCovs[[1]])
covnet <- network(matrix(1,n,n), directed=directed)
for (i in c(1:NumCovNets))
{
	set.edge.value(covnet, names(DyadCovs)[i], DyadCovs[[i]])
}

covnet

}

get.covnet.list <- function(DyadCovs)
{
	Num.Nets <- length(DyadCovs)
	covnets <- vector('list',Num.Nets)
	for (i in c(1:Num.Nets))
	{
		covnets[[i]] <- get.covnets(DyadCovs[[i]])
	}
	covnets
}

get.clist.main <- function(net=NULL,covnet=NULL,form=NULL)
{
	formula <- as.formula(paste("net ~", as.character(form)[2]))
	y <- ergm.getnetwork(formula)
	model <- ergm_model(formula, y)
	sy <- summary(formula)
	p <- length( sy )
	Clist <- ergm.Cprepare(y, model)
	control <- control.ergm(MCMC.burnin = aux.iters, MCMC.interval = 1, 
        MCMC.samplesize = 1)
    proposal <- ergm_proposal(object = ~., constraints = ~., 
        arguments = control$MCMC.prop.args, nw = y)  
	args <- list(Clist=Clist,control=control,proposal=proposal)
	args
}

### The function to carry out work 
sim_net_list <- function(form, nets, thetas,covnets,
                         statsonly = FALSE,
                         sample_size = 1, 
                         burnin = 1e+4,
                         interval = 1000,
                         num_cores = 1) { 
suppressMessages(library(ergm))
suppressMessages(library(parallel))
# note to self: needs testing
  # Create list for simulation  
  sim_list <- rep(list(NULL), length(nets))
  for (i in 1:length(sim_list)) { 
    sim_list[[i]] <- list(form = NULL, net = NULL)
    sim_list[[i]]$form = form 
    sim_list[[i]]$net = nets[[i]]
    sim_list[[i]]$theta = thetas[[i]]
    sim_list[[i]]$covnets = covnets[[i]]
  }

  
  sim <- mclapply(sim_list, 
                  sim_net_slave, 
                  statsonly = statsonly, 
                  sample_size = sample_size, 
                  burnin = burnin, 
                  interval = interval, 
                  mc.cores = num_cores)
  return(sim)
}


sim_net_slave <- function(cur_sim,
							sample_size = 1,
                          burnin = 1e+4,
                          interval = 1000) { 
# note to self: not tested with mcapply
  
  Clist <- cur_sim$Clist
  proposal <- cur_sim$proposal
  theta <- cur_sim$theta
  control <- cur_sim$control
delta <- ergm_MCMC_slave(	Clist = Clist,
							proposal =proposal,
							samplesize = 1,
							eta = theta,
							control = control,
							verbose = FALSE,
							burnin = burnin)
  return(delta)
}


h.ergm <- function(formula, tot.iterations = NULL , aux.var.iter = NULL , DyadCovs =NULL,nets=NULL)
{
Num.Nets <- length(nets)# have to be list of network objects
covnets <- 	get.covnet.list(DyadCovs)
# nets <- vector('list',Num.Nets) # have to be list of network objects
thetas <- vector('list',Num.Nets)
args <- vector('list',Num.Nets)
for (i in c(1:Num.Nets)){
args[[i]] <- get.clist.main(net=nets[[i]],covnet=covnets[[i]],form=formula)
}
# --- what's needed to sumulate using slave:
#delta <- ergm_MCMC_slave(Clist = ClistMain, proposal =proposal, samplesize = 1,
#               eta = theta1, control = control, verbose = FALSE,
 #              burnin = aux.var.iter)
	
}

#source('hergm_flippers.R')

h.ergm.slow <- function(formula, 
                        tot.iterations = NULL , 
                        aux.var.iter = 1e+4 , 
                        covnets =NULL,
                        nets=NULL, 
                        W=NULL,
                        fixed=NULL,
                        isalpha=NULL, 
                        isbeta=NULL,
                        mu.alpha = NULL,
                        sigma.alpha =NULL,
                        sigma.epsilon.alpha = NULL,
                        sigma.epsilon.eta =NULL,
                        sigma.epsilon.beta = NULL,
                        mu.eta =NULL,
                        sigma.eta =NULL,
                        mu.beta =NULL,
                        sigma.beta =NULL ,
                        priorDf =NULL,
                        priorKappa = NULL,
                        priorSigma = NULL,
                        priorMu = NULL,
                        parallel = FALSE,
                        printFreq = NULL,
                        verbose =TRUE,
                        burnin = 1,
                        thinning = 1,
                        startEta = NULL,
                        startBeta = NULL,
                        startAlpha = NULL,
                        obs.stats = NULL,
                        test.beta =TRUE)
{
	require('ergm')
  require('sna')
  require('network')
	 require('mvtnorm')
	 require('expm')
	 
	 if (is.null(printFreq))
	 {
	 	printFreq <- round(tot.iterations/10)
	 }
Num.Nets <- length(nets)# have to be list of network objects

# nets <- vector('list',Num.Nets) # have to be list of network objects
thetas <- vector('list',Num.Nets)
net <- nets[[1]] 
covnet <- covnets[[1]] 
modelSpec <- as.formula(paste("net ~", as.character(formula)[3]))

# mod <- ergm.getmodel(modelSpec , net) #depreciated
mod <- ergm_model(modelSpec , net)
# PosStat <- ergm.getglobalstats(net, mod) #depreciated
PosStat <- summary( mod, net)

#n <- dim(y[,])[1]
# figure out total number of iterations
p <- length(PosStat)

theta.save.ind <- seq( from = burnin , to = tot.iterations, by =  thinning )

globalStats <- matrix(0,p,Num.Nets)


for (k in c(1:Num.Nets) ){
 
  net <- nets[[k]] 
  covnet <- covnets[[k]] 
  modelSpec <- as.formula(paste("net ~", as.character(formula)[3]))
  mod <- ergm_model(modelSpec , net)
  globalStats[,k] <- summary( mod, net)
cat('initialised net ',k,', ')
}


# then we can define p1 and p2 and q
p1 <- sum(isalpha)  
p2 <- sum(isbeta)
q <- sum(fixed)
if (is.null(startEta ))
{
	
	Eta <- matrix(0,q,1)
}
if (!is.null(startEta ))
{
	
	Eta <- startEta
	
	if (length(Eta)>1)
	{
	  Eta <- matrix( Eta,q,1)
	}
}
if (is.null(startBeta  ))
{
	
	Beta <- matrix(0,p2,1)
}
if (!is.null(startBeta  ))
{
	
	Beta <- startBeta
	if (length(Beta)>1)
	{
	  Beta <- matrix( Beta,p2,1)
	}
	
}
if (is.null(startAlpha  ))
{
	
	Alpha <- matrix(0,p1,Num.Nets)
}
if (!is.null(startAlpha  ))
{
	
	Alpha <- startAlpha
}

if (!is.null(obs.stats))
{
  do.MPLE <- TRUE
  like <- 0
  like.sav <- matrix(0, length(theta.save.ind))
}
                        
if (is.null(obs.stats))
{
  do.MPLE <- FALSE
  like <- NULL
}



Alpha1 <- Alpha
Beta1 <- Beta
Eta1 <- Eta
theta <- matrix(0,p,Num.Nets)
theta1 <- theta

Post.Eta <- matrix(0, length(theta.save.ind) , q )
colnames(Post.Eta) <- names(PosStat)[fixed]
Post.Beta <- matrix(0, length(theta.save.ind) , p2 )
colnames(Post.Beta) <- names(PosStat)[isbeta]
Post.Alpha <- array(0, dim = c(length(theta.save.ind) , p1,Num.Nets) )

Post.Mu <- matrix(0, length(theta.save.ind) , p1 )
colnames(Post.Mu) <- names(PosStat)[isalpha]
Post.Sigma <- array(0, dim=c(length(theta.save.ind) , p1,p1) )
Post.Theta <- array(0, dim=c(length(theta.save.ind) , p,Num.Nets) )
t <- 1


for (k in c(1:Num.Nets) ){
  
  theta[fixed==TRUE,k] <-  Eta
  theta[isalpha,k] <- Alpha[,k]
  if (p2>0){
    theta[isbeta,k] <-  theta[isbeta,k] + Beta*W[k,]# isbeta always is a subset of isalpha
  }
  
  thetas[[k]] <- theta[,k]
}


for ( inter in c(1: tot.iterations) )
{

# update varying alphas:
# only fixed==FALSE
 
for (k in c(1:Num.Nets) ){
  Alpha1[,k] <-  rmvnorm(1, Alpha[,k] , sigma = sigma.epsilon.alpha)[1,]
  theta1[,k] <- theta[,k]
  theta1[isalpha,k] <- Alpha1[,k]
  if (p2>0){
  theta1[isbeta,k] <-  theta1[isbeta,k] + Beta*W[k,]# isbeta always is a subset of isalpha
  }
	
thetas[[k]] <- theta1[,k]
}


if (parallel == TRUE & do.MPLE==FALSE)
{
sim2 <- sim_net_list.slow(form, nets, thetas,covnets,
                     statsonly = TRUE,
                     sample_size = 1,
                     burnin =aux.var.iter,
                     num_cores = 2)
 sim2 <- matrix(unlist(sim2),p,Num.Nets)                     
 }  
 
 if (parallel == FALSE & do.MPLE==FALSE)
{
	sim2  <- matrix( 0 , p, Num.Nets )
	for (k in c(1:Num.Nets) ){
	net <- nets[[k]] 
  covnet <- covnets[[k]] 
  modelSpec <- as.formula(paste("net ~", as.character(formula)[3]))
  #net
  mod <- ergm_model(modelSpec , net)
  sim2[,k] <- simulate( modelSpec, coef=thetas[[k]], output = 'stats',control=control.simulate(
    MCMC.burnin=aux.var.iter))
  }
   
	
}                  
                    
                    

delta <- matrix(0,p,1)
for (k in c(1:Num.Nets) ){
	# for TRUE random use hierarchial MVN
	# for regression parameters use mu=beta*z and sigma
	 pr <- dmvnorm(rbind(Alpha1[,k] , Alpha[,k] ), 
                          mean = mu.alpha, 
                          sigma = sigma.alpha,
                          log=TRUE)
	 
	if (do.MPLE==FALSE){ 
	delta <- sim2[,k] - globalStats[,k]
	accratio <- (theta[,k] - theta1[,k]) %*% delta + pr[1] - pr[2] 
	if (accratio >= log(runif(1))) {
              #  theta[isalpha,k] <- theta1[isalpha,k]# we might as well update the entire vector
                theta[isalpha,k] <- theta1[isalpha,k]
                Alpha[,k] <- Alpha1[,k]
	}
	
	
	}
	if (do.MPLE){
	  like.old <- get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta[,k])
	  like.new <- get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta1[,k])
	  accratio <- (like.new - like.old) + pr[1] - pr[2] 
	  if (accratio >= log(runif(1))) {
	    theta[isalpha,k] <- theta1[isalpha,k]
	    Alpha[,k] <- Alpha1[,k]
	  }
	} 
	thetas[[k]] <- theta[,k]
}
# update constant theta:
# only fixed==TRUE

 if (dim(Eta)[1]==1)
 {
 	Eta1[1,1] <- rnorm(1, Eta, sd = sqrt(sigma.epsilon.eta) )
 }
 if (dim(Eta)[1]>1)
 {
Eta1[,1] <- rmvnorm(1, Eta, sigma = sigma.epsilon.eta)[1,]
 }



if (p2>0){
  if (test.beta){
if (dim(Beta)[1]==1)
 {
 Beta1[1,1] <-  rnorm(1, Beta , sd = sqrt(sigma.epsilon.beta ) )
 
 }
 if (dim(Beta)[1]>1)
 {
 
Beta1[,1] <-  rmvnorm(1, Beta , sigma = sigma.epsilon.beta)[1,]
 }
  }
  if (test.beta==FALSE)
  {
    Beta1 <-  Beta
  }
}
for (k in c(1:Num.Nets) ){
  theta1[,k] <- theta[,k]
  theta1[fixed==TRUE,k] <-  Eta1
  theta1[isalpha,k] <- Alpha[,k]
  if (p2>0){
  theta1[isbeta,k] <-  theta1[isbeta,k] + Beta1*W[k,]# isbeta always is a subset of isalpha
  }
thetas[[k]] <- theta1
# need to get dimesions right

}


### ---- determine prior for eta
prior1 <- 0
prior2 <- 0

if (dim(Eta)[1]>1)
 {

pr <- dmvnorm(rbind(t(Eta1) , t(Eta )), 
              mean = mu.eta, 
              sigma = sigma.eta,
              log=TRUE)

# dmvnorm(t(Eta1),
#         mean = mu.eta,
#         sigma = sigma.eta,
#         log=TRUE)
# 
# dmvnorm(t(Eta),
#         mean = mu.eta,
#         sigma = sigma.eta,
#         log=TRUE)

 }    
 if (dim(Eta)[1]==1)
 {
pr <- dnorm(rbind(Eta1 , Eta ), 
              mean = mu.eta, 
              sd = sqrt(sigma.eta),
              log=TRUE)
 }

### ----- determine prior for beta
prior1 <-pr[1]
prior2 <- pr[2]
if (p2>0){
if (dim(Beta)[1]>1)
 {

pr <- dmvnorm(rbind(t(Beta1) , t(Beta) ), 
              mean = mu.beta, 
              sigma = sigma.beta,
              log=TRUE)
 }
 if (dim(Beta)[1]==1)
 {
pr <- dnorm(rbind(Beta1 , Beta ), 
              mean = mu.beta, 
              sd = sqrt(sigma.beta),
              log=TRUE)
 	
 }
}
### ---- multiply priors
prior1 <-prior1+pr[1]
prior2 <- prior2+pr[2]

### ---- draw networks and evaluate acceptance ratio

delta <- matrix(0,p,1)
accratio <- 0

if (parallel == TRUE & do.MPLE==FALSE)
{
  sim2 <- sim_net_list.slow(form, nets, thetas,covnets,
                            statsonly = TRUE,
                            sample_size = 1,
                            burnin =aux.var.iter,
                            num_cores = 2)
  sim2 <- matrix(unlist(sim2),p,Num.Nets)      
  
  Success <- TRUE
  for (k in c(1:Num.Nets) ){
    delta <- sim2[,k] - globalStats[,k]
    #	accratio <- accratio+(theta[,k] - theta1[,k]) %*% delta
    # needs to be done group by group: accratio >= log(runif(1))
    accratio <- (theta[,k] - theta1[,k]) %*% delta
    accratio <- accratio + prior1  - prior2
    Success <- (accratio >= log(runif(1)))
  }
  
  
}  


#### === check for acceptance ratio
do.min = FALSE
if (parallel == FALSE & do.MPLE==FALSE & do.min == TRUE)
{
  sim2  <- matrix( 0 , p, Num.Nets )
  Success <- TRUE
  for (k in c(1:Num.Nets) ){
    if (Success){
    net <- nets[[k]] 
    covnet <- covnets[[k]] 
    modelSpec <- as.formula(paste("net ~", as.character(formula)[3]))
    #net
    mod <- ergm_model(modelSpec , net)
    sim2[,k] <- simulate( modelSpec, coef=theta1[,k], output = 'stats',control=control.simulate(
      MCMC.burnin=aux.var.iter))
    delta <- sim2[,k] - globalStats[,k]
    accratio <- (theta[,k] - theta1[,k]) %*% delta
    accratio <- accratio + prior1  - prior2
    Success <- (accratio >= log(runif(1)))
    }
  }
  
  
} 

if (parallel == FALSE & do.MPLE==FALSE & do.min == FALSE)
{
  sim2  <- matrix( 0 , p, Num.Nets )
  Success <- TRUE
  accratio <- 0
  for (k in c(1:Num.Nets) ){
    
      net <- nets[[k]] 
      covnet <- covnets[[k]] 
      modelSpec <- as.formula(paste("net ~", as.character(formula)[3]))
      #net
      mod <- ergm_model(modelSpec , net)
      sim2[,k] <- simulate( modelSpec, coef=theta1[,k], output = 'stats',control=control.simulate(
        MCMC.burnin=aux.var.iter))
      delta <- sim2[,k] - globalStats[,k]
      accratio <-accratio+ (theta[,k] - theta1[,k]) %*% delta
      
     
    
  }
  accratio <- accratio + prior1  - prior2
  Success <- (accratio >= log(runif(1)))
} 

###=== go back to main algo

if (do.MPLE)
{
  
  like.old <- 0
  like.new <- 0
  for (k in c(1:Num.Nets) ){
  
    like.old <- like.old + get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta[,k])
    like.new <- like.new + get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta1[,k])
    
  }
  accratio <- (like.new - like.old) +pr[1] - pr[2] 
  Success <- (accratio >= log(runif(1)))
  if (Success)
  {
    like <- like.new 
  }
  if (Success==FALSE)
  {
    like <- like.old 
  }
}

#browser()






#accratio <- (theta - theta1) %*% delta + prior1  - prior2
#accratio <- accratio + prior1  - prior2
if (Success) {
  
                Eta <- Eta1
                Beta <- Beta1
            }

# make sure that theta as a whole is updated:
for (k in c(1:Num.Nets) ){
 
  theta[fixed==TRUE,k] <-  Eta
  theta[isalpha,k] <- Alpha[,k]
  if (p2>0){
  theta[isbeta,k] <-  theta[isbeta,k] + Beta*W[k,]# isbeta always is a subset of isalpha
  }
  thetas[[k]] <- theta
  # need to get dimesions right
}


### now we can update  mu.alpha and sigma.alpha
updates <- sampleMuSigma(Thetas =Alpha, p1=p1,nGroup=Num.Nets,priorDf=priorDf,priorSigma=priorSigma, priorKappa=priorKappa, priorMu = priorMu)
# out <-list( SigmaTemp =  SigmaTemp ,  muTemp =  muTemp)
mu.alpha <- updates$muTemp
sigma.alpha <- updates$SigmaTemp

if (inter %in% theta.save.ind){
 				 			Post.Eta[t, ] <- Eta
 				 			
 				 			Post.Beta[t,] <- Beta
Post.Alpha[t, ,] <- Alpha
Post.Mu[t,] <- mu.alpha
Post.Sigma[t,,] <-sigma.alpha
rownames(theta) <- names(PosStat)
Post.Theta[t,,] <- theta

if (do.MPLE){
like.sav[t] <- like
}
 				 			t <- t +1
 			}

if ((verbose==TRUE) && (inter %% printFreq == 0) )
      {
  
      	cat('at iteration ', inter,' out of ',tot.iterations, ' parameters:\n')
      	
      	cat('mu_alpha',mu.alpha ,'\n')
      	cat('eta ',Eta ,'\n')
      	cat('beta ',Beta ,'\n')
      	}
      	
}


if (do.MPLE){
  like <- like.sav
}

Results <- list(Eta = Post.Eta,Beta=Post.Beta,Alpha= 
Post.Alpha, mu.alpha =Post.Mu, sigma.alpha = Post.Sigma, thetas = Post.Theta, likelihood = like)

Results
}


### The function to carry out work 
sim_net_list.slow <- function(form, nets, thetas,covnets,
                         statsonly = FALSE,
                         sample_size = 1, 
                         burnin = 1e+4,
                         interval = 1000,
                         num_cores = 1) { 
suppressMessages(library(ergm))
suppressMessages(library(parallel))

  # Create list for simulation  
  sim_list <- rep(list(NULL), length(nets))
  for (i in 1:length(sim_list)) { 
    sim_list[[i]] <- list(form = NULL, net = NULL)
    sim_list[[i]]$form = form 
    sim_list[[i]]$net = nets[[i]]
    sim_list[[i]]$theta = thetas[[i]]
    sim_list[[i]]$covnets = covnets[[i]]
  }

  
  sim <- mclapply(sim_list, 
                  sim_net_slave.slow, 
                  statsonly = statsonly, 
                  sample_size = sample_size, 
                  burnin = burnin, 
                  interval = interval, 
                  mc.cores = num_cores)
  return(sim)
}

sim_net_slave.slow <- function(cur_sim,
                          statsonly = FALSE, 
                          sample_size = 1,
                          burnin = 1e+4,
                          interval = 1000) { 

  cur_net <- cur_sim$net 
  covnet <- cur_sim$covnets
  theta <- cur_sim$theta
  form <- cur_sim$form
  form <- as.formula(paste("cur_net ~", as.character(form)[3]))
  sim <- simulate(form, 
                  coef = theta, 
                  control = control.simulate(MCMC.burnin = burnin,
                                             MCMC.interval = interval),
                  nsim = sample_size,
                  statsonly = statsonly)

  return(sim)
}

##@sampleMuSigma internal sienaBayes;
## Gibbs algorithm to sample new Mu and Sigma parameters
sampleMuSigma <- function(Thetas =NULL,
p1=NULL,
nGroup=NULL,
priorDf=NULL,
priorSigma=NULL,
priorKappa=NULL,
priorMu=NULL)
{		
  ##@protectedInverse internal sampleMuSigma inverse of p.s.d matrix
  protectedInverse <- function(x)
  {
    if (inherits(try(xinv <- chol2inv(chol(x)),
                     silent=TRUE), "try-error"))
    {
      # Now make this x positive definite, if it is not.
      # See above for a more extensive treatment of the same.
      # Adapted from function make.positive.definite in package corpcor
      # which uses a method by Higham (Linear Algebra Appl 1988)
      # but changed to make the matrix positive definite (psd)
      # instead of nonnegative definite.
      # The idea is to left-truncate all eigenvalues to delta0.
      # The construction with tol, not used now,
      # is to ensure positive definiteness given numerical inaccuracy.
      es <- eigen(x)
      esv <- es$values
      delta0 <- 1e-6
      cat("Eigenvalues Sigma = ", sort(esv), "\n")
      if (min(esv) < delta0)
      {
        delta <- delta0
        tau <- pmax(delta, esv)
        #						cat("Smallest eigenvalue of Sigma now is ",
        #									min(esv),"; make posdef.\n")
        xinv <- es$vectors %*% 
          diag(1/tau, dim(x)[1]) %*% t(es$vectors)
      }
    }
    xinv
  }
  invWish <- function(v,S){
    # Draw from the inverse Wishart distribution with df v
    # and scale matrix S
    # inefficient if drawn multiply for the same S
    protectedInverse(rWishart(1,v,protectedInverse(S))[,,1])
  }
  
  muhat <- matrix( c(0), p1, 1)
 
						
  muhat <- rowMeans( Thetas )# the average across groups
  matQ <- (nGroup - 1)*cov(t(Thetas))
  # prior Lambda = z$priorDf*z$priorSigma
  SigmaTemp <<- invWish(priorDf + nGroup,
                          priorDf*priorSigma + matQ +
                            (priorKappa*nGroup/(priorKappa + nGroup))*
                            tcrossprod( muhat - priorMu, muhat - priorMu ) )
  muTemp <<- t(chol( ( priorKappa + nGroup )^(-1)*SigmaTemp )) %*%
    rnorm( p1 , 0 , 1 ) +
    (nGroup/( priorKappa + nGroup ))*muhat +
    (priorKappa/( priorKappa + nGroup ))*priorMu
  out <-list( SigmaTemp =  SigmaTemp ,  muTemp =  muTemp)
  out
}

get.test.dyad <- function(adj=NULL, covnet=NULL,sex=NULL)
{
  n <- dim(adj)[1]
  cov.mat <- array(data=NA,dim=c(n*(n-1)/2,4,4))
  put.here <- 1
  for (i in c(1:(n-1)))
  {
    for (j in c((i+1):n ) )
    {
      cov.mat[put.here,1,] <- c(0,1,1,2)#edges
      cov.mat[put.here,2,] <- c(0,0,0,1)# mutual
      cov.mat[put.here,3,] <- c(0,sex[i],sex[j],sex[i]+sex[j])# sex
      cov.mat[put.here,4,] <- c(0,covnet[i,j],covnet[j,i],covnet[i,j]+covnet[j,i] )#covnet
      put.here <- put.here +1
    }
  }
  cov.mat
}

get.resp.mat <- function(adj=NULL)
{
  
  n <- dim(adj)[1]
  resp.mat <- matrix(0,n*(n-1)/2,4)
  put.here <- 1
  for (i in c(1:(n-1)))
  {
    for (j in c((i+1):n ) )
    {
      if (adj[i,j]==1 & adj[j,i]==1)
          {
      resp.mat[put.here,4] <-  1
      }
      if (adj[i,j]==1 & adj[j,i]==0)
      {
        resp.mat[put.here,2] <-  1
      }
      if (adj[i,j]==0 & adj[j,i]==1)
      {
        resp.mat[put.here,3] <-  1
      }
      if (adj[i,j]==0 & adj[j,i]==0)
      {
        resp.mat[put.here,1] <-  1
      }
      put.here <- put.here+1
    }
  }
  
  resp.mat
}

get.mut.like <- function(resp.mat=NULL,cov.mat=NULL,theta=NULL)
{
  log.like <- 0
  M <- dim(resp.mat)[1]
  for (k in c(1:M))
  {
    zobs <- cov.mat[k,,resp.mat[k,]==1]
   
    this.like <- theta %*% zobs
    like.denom <- exp(theta %*% cov.mat[k,,1])
    like.denom <-like.denom +exp(theta %*% cov.mat[k,,2])
    like.denom <-like.denom +exp(theta %*% cov.mat[k,,3])
    like.denom <-like.denom +exp(theta %*% cov.mat[k,,4])
    this.like <-  this.like - log(like.denom )
    log.like <-log.like +this.like
  }
 
  log.like    
}

get.mut.stats <- function(resp.mat=NULL,cov.mat=NULL)
{
  stats <- matrix(0,1,4)
  
  M <- dim(resp.mat)[1]
  for (k in c(1:M))
  {
    stats <- stats + cov.mat[k,,resp.mat[k,]==1]
    
   
  }
  
  stats   
}


r.mvnmat <- function(n=NULL, sex=NULL,coef=NULL,covnet=NULL)
  {
  # generate n by 4 matrix
  net <- matrix(0,n,4)
  net[,1] <- rnorm(n, mean = coef[1],sd=1)# mimic edges
  net[,2] <- rnorm(n, mean = coef[2],sd=1)# mimic mutual
  net[,3] <- rnorm(n, mean = coef[3]*sex,sd=1)# mimic sex
  act.cov <- colSums(as.sociomatrix(covnet, 'colaborate'))
  net[,4] <- rnorm(n, mean = coef[4]*act.cov ,sd=1)# mimic sex
  
  
  net
}

d.mvnmat <- function(net=NULL, sex=NULL,coef=NULL,covnet=NULL)
{
  # generate n by 4 matrix
  like <- matrix(0,1,4)
  like[,1] <- sum(dnorm(net[,1], mean = coef[1],sd=1,log=TRUE))# mimic edges
  like[,2] <- sum(dnorm(net[,2], mean = coef[2],sd=1,log=TRUE))# mimic mutual
  like[,3] <- sum(dnorm(net[,3], mean = coef[3]*sex,sd=1,log=TRUE))# mimic sex
  act.cov <- colSums(as.sociomatrix(covnet, 'colaborate'))
  like[,4] <- sum(dnorm(net[,4], mean = coef[4]*act.cov ,sd=1,log=TRUE))# mimic sex
  
  
  like
}



h.ergm.slow.fake <- function(formula, 
                        tot.iterations = NULL , 
                        aux.var.iter = 1e+4 , 
                        covnets =NULL,
                        nets=NULL, 
                        W=NULL,
                        fixed=NULL,
                        isalpha=NULL, 
                        isbeta=NULL,
                        mu.alpha = NULL,
                        sigma.alpha =NULL,
                        sigma.epsilon.alpha = NULL,
                        sigma.epsilon.eta =NULL,
                        sigma.epsilon.beta = NULL,
                        mu.eta =NULL,
                        sigma.eta =NULL,
                        mu.beta =NULL,
                        sigma.beta =NULL ,
                        priorDf =NULL,
                        priorKappa = NULL,
                        priorSigma = NULL,
                        priorMu = NULL,
                        parallel = FALSE,
                        printFreq = NULL,
                        verbose =TRUE,
                        burnin = 1,
                        thinning = 1,
                        startEta = NULL,
                        startBeta = NULL,
                        startAlpha = NULL,
                        obs.stats = NULL,
                        test.beta =TRUE)
{
  require('ergm')
  require('sna')
  require('network')
  require('mvtnorm')
  require('expm')
  
  if (is.null(printFreq))
  {
    printFreq <- round(tot.iterations/10)
  }
  Num.Nets <- length(nets)# have to be list of network objects
  
  # nets <- vector('list',Num.Nets) # have to be list of network objects
  thetas <- vector('list',Num.Nets)
  net <- nets[[1]] 
  
  PosStat <- colSums(net)
  names( PosStat) <- c('edge','mut','sex','cov')
  #n <- dim(y[,])[1]
  # figure out total number of iterations
  p <- dim(net)[2]
  
  theta.save.ind <- seq( from = burnin , to = tot.iterations, by =  thinning )
  
  globalStats <- matrix(0,p,Num.Nets)
  
  
  for (k in c(1:Num.Nets) ){
    
    net <- nets[[k]] 
    covnet <- covnets[[k]] 
   
    globalStats[,k] <- colSums(net)
    cat('initialised net ',k,', ')
  }
  
  
  # then we can define p1 and p2 and q
  p1 <- sum(isalpha)  
  p2 <- sum(isbeta)
  q <- sum(fixed)
  if (is.null(startEta ))
  {
    
    Eta <- matrix(0,q,1)
  }
  if (!is.null(startEta ))
  {
    
    Eta <- startEta
  }
  if (is.null(startBeta  ))
  {
    
    Beta <- matrix(0,p2,1)
  }
  if (!is.null(startBeta  ))
  {
    
    Beta <- startBeta
  }
  if (is.null(startAlpha  ))
  {
    
    Alpha <- matrix(0,p1,Num.Nets)
  }
  if (!is.null(startAlpha  ))
  {
    
    Alpha <- startAlpha
  }
  
  
    do.MPLE <- TRUE
    like <- 0
    like.sav <- matrix(0, length(theta.save.ind))
  
  
  
  
  
  Alpha1 <- Alpha
  Beta1 <- Beta
  Eta1 <- Eta
  theta <- matrix(0,p,Num.Nets)
  theta1 <- theta
  
  Post.Eta <- matrix(0, length(theta.save.ind) , q )
  colnames(Post.Eta) <- names(PosStat)[fixed]
  Post.Beta <- matrix(0, length(theta.save.ind) , p2 )
  colnames(Post.Beta) <- names(PosStat)[isbeta]
  Post.Alpha <- array(0, dim = c(length(theta.save.ind) , p1,Num.Nets) )
  
  Post.Mu <- matrix(0, length(theta.save.ind) , p1 )
  colnames(Post.Mu) <- names(PosStat)[isalpha]
  Post.Sigma <- array(0, dim=c(length(theta.save.ind) , p1,p1) )
  Post.Theta <- array(0, dim=c(length(theta.save.ind) , p,Num.Nets) )
  t <- 1
  
  
  for (k in c(1:Num.Nets) ){
    
    theta[fixed==TRUE,k] <-  Eta
    theta[isalpha,k] <- Alpha[,k]
    if (p2>0){
      theta[isbeta,k] <-  theta[isbeta,k] + Beta*W[k,]# isbeta always is a subset of isalpha
    }
    
    thetas[[k]] <- theta[,k]
  }
  
  
  for ( inter in c(1: tot.iterations) )
  {
    
    # update varying alphas:
    # only fixed==FALSE
    
    for (k in c(1:Num.Nets) ){
      Alpha1[,k] <-  rmvnorm(1, Alpha[,k] , sigma = sigma.epsilon.alpha)[1,]
      theta1[,k] <- theta[,k]
      theta1[isalpha,k] <- Alpha1[,k]
      if (p2>0){
        theta1[isbeta,k] <-  theta1[isbeta,k] + Beta*W[k,]# isbeta always is a subset of isalpha
      }
      
      thetas[[k]] <- theta1[,k]
    }
    
    
    
    
          
    
    
    
    delta <- matrix(0,p,1)
    for (k in c(1:Num.Nets) ){
      # for TRUE random use hierarchial MVN
      # for regression parameters use mu=beta*z and sigma
      pr <- dmvnorm(rbind(Alpha1[,k] , Alpha[,k] ), 
                    mean = mu.alpha, 
                    sigma = sigma.alpha,
                    log=TRUE)
      
     
      if (do.MPLE){
        net <- nets[[k]] 
        covnet <- covnets[[k]] 
        sex <- obs.stats[[k]]$sex 
        like.old <- sum(d.mvnmat(net=net, sex=sex,coef=theta[,k],covnet=covnet))
        like.new <- sum(d.mvnmat(net=net, sex=sex,coef=theta1[,k],covnet=covnet))
       # like.old <- get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta[,k])
        #like.new <- get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta1[,k])
        accratio <- (like.new - like.old) + pr[1] - pr[2] 
        if (accratio >= log(runif(1))) {
          theta[isalpha,k] <- theta1[isalpha,k]
          Alpha[,k] <- Alpha1[,k]
        }
      } 
      thetas[[k]] <- theta[,k]
    }
    # update constant theta:
    # only fixed==TRUE
    
    if (dim(Eta)[1]==1)
    {
      Eta1[1,1] <- rnorm(1, Eta, sd = sqrt(sigma.epsilon.eta) )
    }
    if (dim(Eta)[1]>1)
    {
      Eta1 <- rmvnorm(1, Eta, sigma = sigma.epsilon.eta)[1,]
    }
    
    
    if (p2>0){
      if (test.beta){
        if (dim(Beta)[1]==1)
        {
          Beta1[1,1] <-  rnorm(1, Beta , sd = sqrt(sigma.epsilon.beta ) )
          
        }
        if (dim(Beta)[1]>1)
        {
          
          Beta1 <-  rmvnorm(1, Beta , sigma = sigma.epsilon.beta)[1,]
        }
      }
      if (test.beta==FALSE)
      {
        Beta1 <-  Beta
      }
    }
    for (k in c(1:Num.Nets) ){
      theta1[,k] <- theta[,k]
      theta1[fixed==TRUE,k] <-  Eta1
      theta1[isalpha,k] <- Alpha[,k]
      if (p2>0){
        theta1[isbeta,k] <-  theta1[isbeta,k] + Beta1*W[k,]# isbeta always is a subset of isalpha
      }
      thetas[[k]] <- theta1
      # need to get dimesions right
      
    }
    
    
    ### ---- determine prior for eta
    prior1 <- 0
    prior2 <- 0
    if (dim(Eta)[1]>1)
    {
      pr <- dmvnorm(rbind(Eta1 , Eta ), 
                    mean = mu.eta, 
                    sigma = sigma.eta,
                    log=TRUE)
    }    
    if (dim(Eta)[1]==1)
    {
      pr <- dnorm(rbind(Eta1 , Eta ), 
                  mean = mu.eta, 
                  sd = sqrt(sigma.eta),
                  log=TRUE)
    }
    
    ### ----- determine prior for beta
    prior1 <-pr[1]
    prior2 <- pr[2]
    if (p2>0){
      if (dim(Beta)[1]>1)
      {
        
        pr <- dmvnorm(rbind(Beta1 , Beta ), 
                      mean = mu.beta, 
                      sigma = sigma.beta,
                      log=TRUE)
      }
      if (dim(Beta)[1]==1)
      {
        pr <- dnorm(rbind(Beta1 , Beta ), 
                    mean = mu.beta, 
                    sd = sqrt(sigma.beta),
                    log=TRUE)
        
      }
    }
    ### ---- multiply priors
    prior1 <-prior1+pr[1]
    prior2 <- prior2+pr[2]
    
    ### ---- draw networks and evaluate acceptance ratio
    
    delta <- matrix(0,p,1)
    accratio <- 0
    
   
    
   
    if (do.MPLE)
    {
      
      like.old <- 0
      like.new <- 0
      for (k in c(1:Num.Nets) ){
        net <- nets[[k]] 
        covnet <- covnets[[k]] 
        sex <- obs.stats[[k]]$sex 
        like.old <- like.old + sum(d.mvnmat(net=net, sex=sex,coef=theta[,k],covnet=covnet))
        like.new <- like.new +sum(d.mvnmat(net=net, sex=sex,coef=theta1[,k],covnet=covnet))
        #like.old <- like.old + get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta[,k])
        #like.new <- like.new + get.mut.like(resp.mat=obs.stats[[k ]]$resp.mat,cov.mat=obs.stats[[k ]]$cov.mat,theta=theta1[,k])
        
      }
      accratio <- (like.new - like.old) +pr[1] - pr[2] 
      Success <- (accratio >= log(runif(1)))
      if (Success)
      {
        like <- like.new 
      }
      if (Success==FALSE)
      {
        like <- like.old 
      }
    }
    
    #browser()
    
    
    
    
    
    
    #accratio <- (theta - theta1) %*% delta + prior1  - prior2
    #accratio <- accratio + prior1  - prior2
    if (Success) {
      
      Eta <- Eta1
      Beta <- Beta1
    }
    
    # make sure that theta as a whole is updated:
    for (k in c(1:Num.Nets) ){
      
      theta[fixed==TRUE,k] <-  Eta
      theta[isalpha,k] <- Alpha[,k]
      if (p2>0){
        theta[isbeta,k] <-  theta[isbeta,k] + Beta*W[k,]# isbeta always is a subset of isalpha
      }
      thetas[[k]] <- theta
      # need to get dimesions right
    }
    
    
    ### now we can update  mu.alpha and sigma.alpha
    updates <- sampleMuSigma(Thetas =Alpha, p1=p1,nGroup=Num.Nets,priorDf=priorDf,priorSigma=priorSigma, priorKappa=priorKappa, priorMu = priorMu)
    # out <-list( SigmaTemp =  SigmaTemp ,  muTemp =  muTemp)
    mu.alpha <- updates$muTemp
    sigma.alpha <- updates$SigmaTemp
    
    if (inter %in% theta.save.ind){
      Post.Eta[t, ] <- Eta
      
      Post.Beta[t,] <- Beta
      Post.Alpha[t, ,] <- Alpha
      Post.Mu[t,] <- mu.alpha
      Post.Sigma[t,,] <-sigma.alpha
      rownames(theta) <- names(PosStat)
      Post.Theta[t,,] <- theta
      
      if (do.MPLE){
        like.sav[t] <- like
      }
      t <- t +1
    }
    
    if ((verbose==TRUE) && (inter %% printFreq == 0) )
    {
      
      cat('at iteration ', inter,' out of ',tot.iterations, ' parameters:\n')
      
      cat('mu_alpha',mu.alpha ,'\n')
      cat('eta ',Eta ,'\n')
      cat('beta ',Beta ,'\n')
    }
    
  }
  
  
  if (do.MPLE){
    like <- like.sav
  }
  
  Results <- list(Eta = Post.Eta,Beta=Post.Beta,Alpha= 
                    Post.Alpha, mu.alpha =Post.Mu, sigma.alpha = Post.Sigma, thetas = Post.Theta, likelihood = like)
  
  Results
}

gof.bayes.hergm <- function(form = NULL, nets = NULL, covnets = NULL,posteriors=NULL,MCMC.burnin=1000,nsim=1000)
{
  Num.Nets <- length(nets)
  Numpost <- dim(posteriors$thetas)[1]
  do.these <- seq(1,Numpost,length.out = nsim)
  gofnets <- vector('list',Num.Nets)
  for (k in c(1:Num.Nets)){
  net <- nets[[k]] 
  covnet <- covnets[[k]] 
  modelSpec <- as.formula(paste("net ~", as.character(form)[3]))
  gofnets[[k]] <- vector('list',nsim)
  for (i in c(1:nsim))
  {
    gofnets[[k]][[i]] <- simulate(modelSpec,  seed = NULL, coef=posteriors$thetas[do.these[i],,k],nsim =1,control=control.simulate(
    MCMC.burnin=MCMC.burnin))
  }
  }
  gofnets
}

result.table <- function(posteriors=NULL,d=3,fixed=NULL,
                         isalpha=NULL, 
                         isbeta=NULL,burnin =1,thinning=1 ,name='HERGM',perc=0.05)
{
  require('xtable')
  N <- dim(posteriors$mu.alpha)[1]
  pick.these <- seq(burnin,N,by=thinning)
  p <- dim(posteriors$thetas)[2]
  alphanames <- colnames(posteriors$mu.alpha)
  etanames <- colnames(posteriors$Eta)
  name.vec <- matrix(NA,p,1)
  name.vec[fixed] <- etanames
  name.vec[isalpha] <- alphanames
  p1 <- sum(isalpha)  
  p2 <- sum(isbeta)
  q <- sum(fixed)
  # radnom effects
  alpha.mean <- matrix(NA,p,1)
  alpha.mean[isalpha] <-  colMeans(posteriors$mu.alpha[pick.these ,])
  alpha.var.tmp <- matrix(NA,p1,1)
  alpha.var <- matrix(NA,p,1)
  alpha.var.tmp.sd <- matrix(NA,p1,1)
  alpha.var.sd <- matrix(NA,p,1)
  
  alpha.ci.tmp <- matrix(NA,p1,2)
  alpha.ci <- matrix(NA,p,2)
  for (k in c(1:p1))
  {
    alpha.var.tmp[k] <- mean(posteriors$sigma.alpha[pick.these ,k,k] )
    alpha.var.tmp.sd[k] <- sd(posteriors$sigma.alpha[pick.these ,k,k] )
    alpha.ci.tmp[k,] <- get.ci(x=posteriors$mu.alpha[pick.these ,k],perc=perc)
  }
  alpha.var[isalpha] <-  alpha.var.tmp
  alpha.ci[isalpha,] <- alpha.ci.tmp
  alpha.var.sd[isalpha] <-alpha.var.tmp.sd
  
  
  alpha.sd.tmp <- matrix(NA,p1,1)
  alpha.sd <- matrix(NA,p,1)
  for (k in c(1:p1))
  {
    alpha.sd.tmp[k] <- sd(posteriors$mu.alpha[pick.these ,k] )
  }
  alpha.sd[isalpha] <-  alpha.sd.tmp
  ### fixed effects
  eta.mean <- matrix(NA,p,1)
  eta.mean[fixed] <- colMeans(posteriors$Eta[pick.these ,])
  eta.sd <- matrix(NA,p,1)
  eta.sd.tmp <- matrix(NA,q,1)
  
  eta.ci.tmp <- matrix(NA,q,2)
  eta.ci <- matrix(NA,p,2)
  
  for (k in c(1:q))
  {
    eta.sd.tmp[k] <- sd(posteriors$Eta[pick.these ,k] )
    eta.ci.tmp[k,] <- get.ci(x=posteriors$Eta[pick.these ,k],perc=perc)
  }
  eta.sd[fixed] <-  eta.sd.tmp
  eta.ci[fixed,] <- eta.ci.tmp
  ## beta
  
  if (sum(isbeta)>0){
  name.beta<-name.vec[isbeta] 
  beta.mean <- matrix(NA,p2,1)
  beta.sd <- matrix(NA,p2,1)
  beta.ci <- matrix(NA,p2,2)
  for (k in c(1:p2))
  {
    beta.mean[k] <- mean( posteriors$Beta[pick.these ,k])
    beta.sd[k] <- sd(posteriors$Beta[pick.these ,k])
    beta.ci[k,] <- get.ci(x=posteriors$Beta[pick.these ,k],perc=perc)
  }
  my.tab <- data.frame(Effects = c(name.vec,name.beta) ,
                       Alpha.mean = c(alpha.mean,matrix(NA,p2,1)),
                       Alpha.sd = c(alpha.sd,matrix(NA,p2,1)),
                       Alpha.low = c(alpha.ci[,1],matrix(NA,p2,1)),
                       Alpha.high = c(alpha.ci[,2],matrix(NA,p2,1)),
                       Alpha.var = c(alpha.var,matrix(NA,p2,1)),
                       Alpha.var.sd = c(alpha.var.sd,matrix(NA,p2,1)),
                       Eta.mean = c(eta.mean,matrix(NA,p2,1)),
                       Eta.sd = c(eta.sd,matrix(NA,p2,1)),
                       Eta.low = c(eta.ci[,1],matrix(NA,p2,1)),
                       Eta.high = c(eta.ci[,2],matrix(NA,p2,1)),
                       Beta.mean = c(matrix(NA,p,1),beta.mean ),
                       Beta.sd = c(matrix(NA,p,1),beta.sd ),
                       Beta.low = c(matrix(NA,p,1),beta.ci[,1]) ,
                       Beta.high = c(matrix(NA,p,1),beta.ci[,2]))
  names(my.tab) <- c("Effect","Mean alpha","SD alpha",
                     "2.5 perc","97.5 perc",
                     "RE var","RE var sd",
                     "Mean eta","SD eta",
                     "2.5 perc","97.5 perc",
                     "Mean beta","SD beta",
                     "2.5 perc","97.5 perc")
  }
  
  if (sum(isbeta)==0){
    my.tab <- data.frame(Effects = name.vec,
                         Alpha.mean = alpha.mean,
                         Alpha.sd = alpha.sd,
                         Alpha.low = alpha.ci[,1],
                         Alpha.high = alpha.ci[,2],
                         Alpha.var = alpha.var,
                         Alpha.var.sd = alpha.var.sd,
                         Eta.mean = eta.mean,
                         Eta.sd = eta.sd,
                         Eta.low = eta.ci[,1],
                         Eta.high = eta.ci[,2])
    names(my.tab) <- c("Effect","Mean alpha","SD alpha",
                       "2.5 perc","97.5 perc",
                       "RE var","RE var sd",
                       "Mean eta","SD eta",
                       "2.5 perc","97.5 perc")
  }
  
  for (k in c(2:dim(my.tab)[2] )){
  my.tab[,k] <- round(my.tab[,k] ,d)
  }
 
    
  my.tab <- xtable(my.tab,auto=TRUE)
  sink(paste(name,'.html',sep='') )
  print(my.tab,type='html')
  sink()
}


plot.alpha <- function(posteriors=NULL,panel=c(2,3)){
par(mfrow=panel)
Num.Nets <- dim(posteriors$thetas)[3]
plot.name <- colnames( posteriors$mu.alpha)
p1 <- dim(posteriors$mu.alpha)[2]
cat('iterations ',dim(posteriors$mu.alpha)[1],'\n')
for (k in c(1:p1))
{
  plot( posteriors$Alpha[,k,1] ,type='l',ylab=plot.name[k],ylim=range(posteriors$Alpha[,k,]),xlab='iteration')# first group
  for (i in c(1:Num.Nets))
  {
    lines( posteriors$Alpha[,k,i])# rest of groups
  }
  lines( posteriors$mu.alpha[,k], col='red')# the posterior for the population-level mean
}
}


plot.eta <- function(posteriors=NULL,panel=c(2,5))
{  
  par(mfrow=panel)
plot.name <- colnames( posteriors$Eta)
q <- dim(posteriors$Eta)[2]
Num.Nets <- dim(posteriors$thetas)[3]
cat('iterations ',dim(posteriors$mu.alpha)[1],'\n')
for (k in c(1:q))
{
  plot(ts(posteriors$Eta[,k]),ylab=plot.name[k],xlab='iteration')
}
### === network-level effects  === #####
p2 <- dim(posteriors$Beta)[2]
if (p2>0){
for (k in c(1:p2)){
plot(ts(posteriors$Beta[,k]),ylab=expression(beta),xlab='iteration')
}
}

}

get.ci <- function(x=NULL,perc=0.05)
{
  x <- sort(x)
  N <- length(x)
  ci.lowhigh <- c( x[ceiling(N*(perc/2))]  , x[floor(N*(1-perc/2))]  )
  ci.lowhigh
  
}

triad.gof <- function(gofnets=NULL,nets=NULL,reorder=TRUE)
{
  require(vioplot)
  Num.Nets <- length(nets) 
  n.sims <- length(gofnets[[k]])
  Big.T <- array(0,dim=c(n.sims,16,Num.Nets))
  for (i in c(1:Num.Nets))
  {
    Big.T[,,i] <- triad.census(gofnets[[i]])
    obs.tri <- triad.census(nets[[i]])
    for (triad.ind in c(1:16))
    {
      Big.T[,triad.ind,i] <- Big.T[,triad.ind,i] -obs.tri[1,triad.ind] 
    }
  }
  par(mfrow=c(4,4))
  for (triad.ind in c(1:16))
  {
    if (reorder)
    {
      
      m.order <- colMeans(Big.T[,triad.ind,])
      vioplot(Big.T[,triad.ind,order(m.order )],xlab=colnames(obs.tri)[triad.ind])
      abline(h=0,col='red')
    }
    
    if (reorder==FALSE)
    {
    vioplot(Big.T[,triad.ind,],xlab=colnames(obs.tri)[triad.ind])
    abline(h=0,col='red')
    }
  }
  
}

degree.gof <- function(gofnets=NULL,nets=NULL,cmode='outdegree',numdegs=6,reorder=TRUE)
{
  require(vioplot)
  Num.Nets <- length(nets) 
  n.sims <- length(gofnets[[k]])
  # outdegree distribution
  n.sims <- length(gofnets[[k]])
 
  Big.Out <- array(0,dim=c(n.sims,numdegs,Num.Nets))
  for (i in c(1:Num.Nets))
  {
    obsdeg <- degree(nets[[i]], cmode = cmode)
    simdeg <- matrix(0,n.sims,length(obsdeg))
    
    for (k in c(1:n.sims))
    {
      simdeg[k,] <- degree(gofnets[[i]][[k]], cmode = cmode)
    }
    
    for (d in c(1:numdegs))
    {
      Big.Out[,d,i] <- rowSums(simdeg==(d-1) )-sum(obsdeg ==(d-1) )
    }
    
  }
  
  par(mfrow=c(2,3))
  for (d in c(1:numdegs))
  {
    if (reorder)
    {
      
      m.order <- colMeans(Big.Out[,d,])
      vioplot( Big.Out[,d,order(m.order) ] ,xlab= paste('degree: ',(d-1)))
      abline(h=0,col='red')
    }
    
    if (reorder==FALSE)
    {
    vioplot( Big.Out[,d, ] ,xlab= paste('degree: ',(d-1)))
    abline(h=0,col='red')
    }
    
  }
  
  
  
}

paste.posterior <- function(pos1=NULL, pos2=NULL, pos3=NULL, pos4=NULL, pos5=NULL, pos6=NULL, pos7=NULL)
{
  posteriors <- pos1
  Num.Nets <- dim(pos1$thetas)[3]
  p <- dim(pos1$thetas)[2]
  p1 <- dim(pos1$mu.alpha)[2]
  q <- dim(posteriors$Eta)[2]
  p2 <- dim(posteriors$Beta)[2]
  
  
  for (s in c(2:7) )
  {
    #paste(assign('pos', paste('pos',s,sep=''))
    pos <- get( paste('pos',s,sep='') )
    if (!is.null(pos)){
    N <- dim(posteriors$mu.alpha)[1]
    Nn <- dim(pos$mu.alpha)[1]
    posteriors$mu.alpha <- rbind(posteriors$mu.alpha,pos$mu.alpha)
    posteriors$Eta <- rbind(posteriors$Eta,pos$Eta)
    
    thetatmp <- array( NA, dim = c(N+Nn,p,Num.Nets) )
    thetatmp[1:N,,] <- posteriors$thetas
    thetatmp[(N+1):(N+Nn),,] <- pos$thetas
    posteriors$thetas <- thetatmp
    posteriors$Beta <- rbind(posteriors$Beta,pos$Beta)
    sigmatmp <- array( NA, dim = c(N+Nn,p1,p1) )
    sigmatmp[1:N,,] <- posteriors$sigma.alpha
    sigmatmp[(N+1):(N+Nn),,] <- pos$sigma.alpha
    posteriors$sigma.alpha <- sigmatmp
    # Alpha
    sigmatmp <- array( NA, dim = c(N+Nn,p1,Num.Nets) )
    sigmatmp[1:N,,] <- posteriors$Alpha
    sigmatmp[(N+1):(N+Nn),,] <- pos$Alpha
    posteriors$Alpha <- sigmatmp
    
    
    }
    
  }

  posteriors  
}



simulate.data.hergm <- function(Num.Nets = 40)
{
  require('ergm')
  require('sna')
  require('mvtnorm')
  require('expm')
  #### ----- now with Beta different from 0
  Num.Nets <-Num.Nets
  nnet <- sample(c(7,8), Num.Nets, replace=TRUE)
  nets <- vector('list',Num.Nets)
  covnets <- vector('list',Num.Nets)
  
  # nets <- vector('list',Num.Nets) # have to be list of network objects
  thetas <- vector('list',Num.Nets)
  W <- matrix(runif(Num.Nets)>0.5,Num.Nets,1)+0
  fixed <- c(TRUE,FALSE,FALSE,FALSE,TRUE)
  isalpha <- c(FALSE,TRUE,TRUE,TRUE,FALSE)
  isbeta <- c(FALSE,TRUE,FALSE,FALSE,FALSE)
  
  p1 <- sum(isalpha)
  ##
  
  
  
  mu.alpha <- matrix( c(.5,-1,0) ,p1,1)
  priorMu <- mu.alpha
  sigma.alpha <- diag(.5,p1)
  priorSigma <- sigma.alpha 
  
  priorSigma <- matrix(0,p1,p1)
  diag(priorSigma) <- 1
  Eta <- matrix(c(-1.5,0.9),2,1)
  Beta <- matrix(-1,1,1)
  #Beta <- matrix(0,1,1)
  priorKappa <- 1
  p<- 5
  Alpha <- matrix(0,p1,Num.Nets)
  my.thetas <- matrix(0,p,Num.Nets)
  my.stats <- matrix(0,p,Num.Nets)
  
  
  obs.stats <- vector('list',Num.Nets)
  your.stats <- matrix(0,p,Num.Nets)
  for (k in c(1:Num.Nets))
  {
    n <- nnet[k]
    net <- rgraph(n,m=1, tprob=0.2)
    net <- as.network(net)
    sex <- matrix(runif(n)>0.5,n,1)+0
    net %v% 'sex' <- sex[,1]
    colaborate <- rgraph(n,m=1, tprob=0.2)
    covnet <- network(matrix(1,n,n), directed=TRUE)
    
    set.edge.value(covnet, 'colaborate', colaborate)
    
    form <- net ~ edges+mutual+nodeocov('sex')+edgecov(covnet,'colaborate' )+gwesp(decay =.69,fixed = T)
    theta <- matrix(0, p, 1)
    theta[fixed ,1] <- Eta
    Alpha[,k] <-  matrix(rmvnorm(1, mu.alpha , sigma = sigma.alpha)[1,],p1,1)
    theta[isalpha,1] <- Alpha[,k]
    theta[isbeta,1] <-  theta[isbeta,1] + Beta*W[k,]
    my.thetas[,k] <- theta
    net <- simulate(form, nsim = 1, seed = NULL, coef=theta)
    nets[[k]] <- net
    covnets[[k]] <- covnet
    mod <- ergm_model(form , net)
    my.stats[,k] <- summary( mod, net)
    # get.mut.like(resp.mat=resp.mat,cov.mat=cov.mat,theta=my.thetas[,k])
  }
  
  my.data <- list(
    form=form, 
    covnets = covnets,
    nets=nets, 
    W=W,
    fixed=fixed,
    isalpha=isalpha, 
    isbeta=isbeta,
    mu.alpha = mu.alpha,
    sigma.alpha =sigma.alpha,
    mu.eta =matrix(0,2,1),
    sigma.eta = diag(100,2),
    mu.beta =0,
    sigma.beta =100 ,
    priorDf =1,
    priorSigma=priorSigma,
    priorKappa = priorKappa,
    priorMu = priorMu,
    startEta = Eta,
    startBeta = Beta,
    startAlpha = Alpha)
  my.data 
}


plot.eta.hist <- function(posteriors=NULL,panel=c(2,5))
{  
  par(mfrow=panel)
  plot.name <- colnames( posteriors$Eta)
  q <- dim(posteriors$Eta)[2]
  Num.Nets <- dim(posteriors$thetas)[3]
  cat('iterations ',dim(posteriors$mu.alpha)[1],'\n')
  for (k in c(1:q))
  {
   hist(posteriors$Eta[,k],xlab=plot.name[k],yaxt='n',ylab=NA,main=NA)
    abline(v=mean(posteriors$Eta[,k]),col='red')
    abline(v=0,col='grey',lty=2)
  }
  ### === network-level effects  === #####
  p2 <- dim(posteriors$Beta)[2]
  if (p2>0){
    for (k in c(1:p2)){
      hist(posteriors$Beta[,k],ylab=expression(beta),xlab='iteration')
    }
  }
  
}



plot.alpha.hist <- function(posteriors=NULL,panel=c(2,3),burnin=1,thin=1){
  par(mfrow=panel)
  N <- dim(posteriors$mu.alpha)[1]
  pick.these <- seq(burnin,N,by=thin)
  Num.Nets <- dim(posteriors$thetas)[3]
  plot.name <- colnames( posteriors$mu.alpha)
  p1 <- dim(posteriors$mu.alpha)[2]
  cat('iterations ',length(pick.these),'\n')
  for (k in c(1:p1))
  {
    #plot( posteriors$Alpha[,k,1] ,type='l',ylab=plot.name[k],ylim=range(posteriors$Alpha[,k,]),xlab='iteration')# first group
    HK <- density( posteriors$mu.alpha[pick.these,k])
    xlim <- range(HK$x)
    ylim <- range(HK$y)
    densities <- vector('list',Num.Nets)
    for (i in c(1:Num.Nets))
    {
      densities[[i]] <- density( posteriors$Alpha[pick.these,k,i] )
      xlim <- range(xlim,densities[[i]]$x)
      ylim <- range(ylim,densities[[i]]$y)
     
    }
    plot(HK$x,HK$y,xlim=xlim,ylim=ylim ,col='red',type='l',xlab=plot.name[k],yaxt='n',ylab=NA)
    for (i in c(1:Num.Nets))
    {
      lines(densities[[i]]$x,densities[[i]]$y ,col='grey')
     # abline( v = mean(posteriors$Alpha[pick.these,k,i] ), col='grey' )
    }
   abline( v=mean(posteriors$mu.alpha[pick.these,k]), col='red')# the posterior for the population-level mean
  }
  
  
 
}

plot.alpha.catepillar <- function(posteriors=NULL,panel=c(2,3),burnin=1,thin=1){
  require(vioplot)
  par(mfrow=panel)
  N <- dim(posteriors$mu.alpha)[1]
  pick.these <- seq(burnin,N,by=thin)
  Num.Nets <- dim(posteriors$thetas)[3]
  plot.name <- colnames( posteriors$mu.alpha)
  p1 <- dim(posteriors$mu.alpha)[2]
  cat('iterations ',length(pick.these),'\n')
  for (k in c(1:p1))
  {
    #plot( posteriors$Alpha[,k,1] ,type='l',ylab=plot.name[k],ylim=range(posteriors$Alpha[,k,]),xlab='iteration')# first group
    ci.lowhigh <- get.ci(x=posteriors$mu.alpha[pick.these,k],perc=0.05)
    mu.mean <- mean(posteriors$mu.alpha[pick.these,k])
    ylim <- range(c(posteriors$Alpha[pick.these,k,],posteriors$mu.alpha[pick.these,k]) )
    m.order <- colMeans(posteriors$Alpha[pick.these,k,])
    vioplot(posteriors$Alpha[pick.these,k,order(m.order )],xlab='group',ylab=plot.name[k],ylim=ylim)
   abline(h=mu.mean,col='red')
    abline(h=ci.lowhigh[1])
    abline(h=ci.lowhigh[2])
   
   }
  
  
  
}


plot.eta.acf <- function(posteriors=NULL,panel=c(5,6),thinning=1,burnin=1,thisperc=0.1)
{  
  par(mfrow=panel,cex.lab=.75,
  oma = c(1,0,0,0) + 0.1,
  mar = c(1,2,4,2) + 0.1)
  plot.name <- colnames( posteriors$Eta)
  N <- dim(posteriors$mu.alpha)[1]
  pick.these <- seq(burnin,N,by=thinning)
  q <- dim(posteriors$Eta)[2]
  Num.Nets <- dim(posteriors$thetas)[3]
  cat('iterations ',dim(posteriors$mu.alpha)[1],'\n')
  for (k in c(1:q))
  {
    ESS <- effectiveSize(posteriors$Eta[,k])
    Hk <- acf(posteriors$Eta[,k],lag.max=1000,plot=FALSE)
    lagminper <- min(which(Hk$acf<thisperc))-1
    acf(posteriors$Eta[,k],lag.max=500,main=paste('N:',N,';ESS:',round(ESS),';c_10:', lagminper ) ) 
    Hk <- acf(posteriors$Eta[pick.these,k],lag.max=1000,plot=FALSE)
    plot(ts(posteriors$Eta[pick.these,k]),ylab=plot.name[k],xlab='iteration',
    main=paste('N:',length(posteriors$Eta[pick.these,k]),';thin:', thinning,';r_30:',round(Hk$acf[30+1],2) ) ) 
    
    abline(h = mean(posteriors$Eta[pick.these,k]),col='red')
    abline( h = 0, col='grey')
    HK <- density(posteriors$Eta[pick.these,k] )
  plot(HK$x,HK$y,type='l',main=plot.name[k],yaxt='n')
    ci1<- get.ci(x=posteriors$Eta[pick.these,k] ,perc=0.05)
    ci2<- get.ci(x=posteriors$Eta[,k] ,perc=0.05)
    lines(ci2,c(0,0),lwd=4,col='grey')
    lines(ci1,c(0,0),lwd=2,col='red')
   
  }
  ### === network-level effects  === #####
  p2 <- dim(posteriors$Beta)[2]
  if (p2>0){
    for (k in c(1:p2)){
      plot(ts(posteriors$Beta[,k]),ylab=expression(beta),xlab='iteration')
      
      
      
    }
  }
  
}

