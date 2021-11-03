ergm.gcd <- function(formula,ergm.est,nsim=NULL,burnin=NULL,thinning =NULL,returnsims=FALSE, silent=FALSE)
{
	require('ergm')
	
	observedNetwork <-  ergm.getnetwork(formula)
	therearemissing <- any(get.edge.attribute(observedNetwork,'na'))
	notergmobject <- class(ergm.est)!='ergm'
	if ( c('m.prior') %in% names(ergm.est) )
	{
		notergmobject <- FALSE
		isbergm <- TRUE
	}
	else {
		isbergm <- FALSE
	}
	
	if (isbergm)
	{
		cat('ergm.est is bergm output, Bayesian GCD will be tried\n')
		 if (is.null(nsim))
  {
  	nsim <- dim(ergm.est$Theta)[1]
  }
  if (is.null(burnin))
  {
  	
  	burnin <- ergm.est$control$MCMC.burnin
  }
		if (length(dim(ergm.est$Theta))==3 )
		{
			Theta <- ergm.est$Theta[seq(from=1,to=dim(ergm.est$Theta)[1],length.out = min(nsim,dim(ergm.est$Theta)[1]) ),,1]
		}
		if (length(dim(ergm.est$Theta))==2 )
		{
			Theta <- ergm.est$Theta[seq(from=1,to=dim(ergm.est$Theta)[1],length.out = min(nsim,dim(ergm.est$Theta)[1]) ),]
		}
ergm.est$covar <- cov(Theta)
		ergm.est$coef <- colMeans(Theta)
		names(ergm.est$coef) <- est$specs
		
	}
	if(therearemissing | notergmobject){
		
		if (therearemissing)
		{
			err.mess <- "missing edges are not handled yet"
		}
		if (notergmobject)
		{
			err.mess <- "ergm.est is not an ergm estimation object"
		}
	 stop(err.mess)
	 }
	 else{
  mod <- ergm.getmodel(formula, observedNetwork)
  n <- network.size(observedNetwork)
  M <- network.edgecount(observedNetwork)
  stat <- ergm.getglobalstats(observedNetwork, mod)
  mrow <- length(stat)
  stat <- matrix(stat,mrow,1)
  GCD <- matrix(NA,n,mrow+1)
  stats <- matrix(NA,mrow,n)
  # GCD defined as ||E(z(Y)|Y_{-ij}=y_{-ij}) - z(y_{obs})||_{I(theta)}
  # ||A||_B = t(A) B^{-1} A

  if (is.null(nsim))
  {
  	nsim <- ergm.est$control$MCMC.samplesize
  }
  if (is.null(burnin))
  {
  	
  	burnin <- ergm.est$control$MCMC.burnin
  }
  if (is.null(thinning))
  {
  	thinning <- ergm.est$control$MCMC.interval
  }
  if (returnsims)
  {
  	SimulatedStats <- array(data = NA, dim = c(n,mrow,nsim), dimnames = list(c(1:n),names(ergm.est$coef),c(1:nsim)))
  }

if (silent==FALSE){  
  cat('Calculating GCD for actor:\n')}
  for (i in c(1:n))
  {
  	cat(i,' ')
  	# simulate with fixallbut(free.dyads)
  	free.dyads <- matrix(0,n,2)
  	free.dyads[,1] <- rep(i,n)
  	free.dyads[,2] <- c(1:n)
  	free.dyads<- free.dyads[free.dyads[,2]!=i,]
  	free.dyads<- as.edgelist(free.dyads,n=n, directed =is.directed(observedNetwork))
  	
 
  	if (isbergm==FALSE){
  	test.net <- simulate(formula, coef=ergm.est$coef , nsim=nsim, statsonly=TRUE, constraint =~fixallbut(free.dyads),control=control.simulate(MCMC.burnin=burnin, MCMC.interval=thinning) )#
  	}
  	if (isbergm)
  	{
  		test.net <- matrix(0,nsim,mrow)
  		for (k in c(1:nsim))
  		{
  			test.net[k,] <- simulate(formula, coef=Theta[k,] , nsim=1, statsonly=TRUE, constraint =~fixallbut(free.dyads),control=control.simulate(MCMC.burnin=burnin) )#
  		}
  		
  	}
  	
  	if (returnsims)
  {
  	SimulatedStats[i,,] <- test.net
  	}
  	  stats[,i] <- matrix(colMeans(test.net),mrow,1)
  	  GCD[i,mrow+1] <- t( stats[,i] - stat) %*% ergm.est$covar %*% (stats[,i] - stat)
  }	
  
  
  GCD[,1:mrow] <- t(stats)
  if (is.null(get.vertex.attribute(observedNetwork,'vertex.names'))!=TRUE)
  {
  	rownames(GCD) <- get.vertex.attribute(observedNetwork,'vertex.names')
  }
 colnames(GCD) <- c(names(ergm.est$coef),'Generalised Cooks')
 GCooksDistance <- list(observed.stats= stat,expected.stats = GCD[,1:mrow],GCDMVP=GCD[,mrow+1])
 if (returnsims)
 {
 	GCooksDistance$sim.stats <- SimulatedStats
 }
 GCooksDistance
 }

}


test.gcd <- function(formula,ergm.est,nsim=NULL,burnin=NULL,thinning =NULL,returnsims=FALSE,numreps=100,burnin.rep=10000,print.stats=FALSE)
{
	# calculate observed GCD
  cat('Calculating observed GCD for actor: \n')
	GCDobs <- ergm.gcd(formula,ergm.est,burnin=burnin,thinning =thinning,returnsims=FALSE,silent=TRUE)
	GCDReps <- matrix(NA, length(GCDobs$GCDMVP),numreps)
	observedNetwork <-  ergm.getnetwork(formula)
	 mod <- ergm.getmodel(formula, observedNetwork)
	 stat <- ergm.getglobalstats(observedNetwork, mod)
  n <- network.size(observedNetwork)
  M <- network.edgecount(observedNetwork)
 
  mrow <- length(stat)

	cat('\n Doing replicate number: \n')
	 if (is.null(nsim))
  {
  	nsim <- ergm.est$control$MCMC.samplesize
  }
  if (is.null(burnin))
  {
  	
  	burnin <- ergm.est$control$MCMC.burnin
  }
   if (is.null(thinning))
  {
  	thinning <- ergm.est$control$MCMC.interval
  }
	for (k in c(1:numreps))
	{
		if (k%%10 == 0)
		{
			cat(k,' ')
		}
		rep.net <- simulate(formula, coef=ergm.est$coef , nsim=1, statsonly=FALSE,control=control.simulate(MCMC.burnin=burnin.rep) )
		 modBen <- ergm.getmodel(formula, rep.net )
        	stat <- ergm.getglobalstats(rep.net , modBen)
    if (print.stats){    	
				print(stat)
    }
		for (i in c(1:n))
		{
			free.dyads <- matrix(0,n,2)
  	free.dyads[,1] <- rep(i,n)
  	free.dyads[,2] <- c(1:n)
  	free.dyads<- free.dyads[free.dyads[,2]!=i,]
  	free.dyads<- as.edgelist(free.dyads,n=n, directed =is.directed(observedNetwork))
  	
#### === PAVEL broke this
#  	control <- control.simulate.formula(MCMC.burnin = burnin, MCMC.interval= thinning)
#  	control$MCMC.samplesize <- nsim
#MHproposal <- MHproposal.ergm(object= modBen, 
 #                                 constraints = ~fixallbut(free.dyads), arguments = control$MCMC.prop.args, 
  #                                nw = rep.net, weights = control$MCMC.prop.weights, class = "c", 
   #                               reference = ~Bernoulli, response = NULL) 
 
  	
#  	test.net <-          as.matrix(ergm.getMCMCsample(nw=rep.net, model=modBen,MHproposal= MHproposal, eta0= matrix(as.numeric(ergm.est$coef),1,mrow), verbose=FALSE,control=control)$statsmatrices[[1]])
#### === So let us do it the hard way
  	form <- as.formula(paste("rep.net ~", as.character(formula)[3]))
  #	netone.stats <- summary(form)
  #	discrepancy <- sum( netone.stats != netone.stats.alt )	
  #	netone.stats.alt <- delta$s+	sy.obs.copy
  	
  	  	test.net <- simulate(form,
  	  	                     coef=ergm.est$coef , 
  	  	                     nsim=nsim, 
  	  	                     statsonly=TRUE, 
  	  	                     constraint =~fixallbut(free.dyads),
  	  	                     control=control.simulate(MCMC.burnin=burnin, MCMC.interval=thinning) )#
  	
  	
  	statsRep <- matrix(colMeans(test.net),mrow,1)
  	GCDReps[i,k] <- t( statsRep- stat) %*% ergm.est$covar %*% (statsRep - stat)
  	 # GCDReps[i,k] <- t( statsRep ) %*% ergm.est$covar %*% statsRep 
  	  #browser()
		}
		
	}

MCMC.test.GCD <- list(GCDobs=GCDobs,GCDReps=GCDReps)
}
	
	
example.gcd <- function()
{	require('ergm')
	data(flo)
	flonet <- as.network(flo,directed=FALSE)
	est <- ergm( flonet ~ edges+kstar(2:3)+triangles)
	convergence.test(flonet ~ edges+kstar(2:3)+triangles,est)
	GCD <- ergm.gcd(flonet ~ edges+kstar(2:3)+triangles,est,returnsims=FALSE)
	MCMC.GCD.text <- test.gcd(flonet ~ edges+kstar(2:3)+triangles,est,numreps=10,nsim=100,burnin=100,thinning =10)
	require('Bergm')
	data(florentine)
	est <- bergm(flomarriage ~ edges +kstar(2:3)+triangles,burn.in=10, aux.iters=500, main.iters=500,  gamma=1)

GCD <- ergm.gcd(flomarriage ~ edges+kstar(2:3)+triangles,est,returnsims=FALSE)

	# Lazega
	library(network)
library(sna)
library(igraph)
library(ergm)
library(data.table)
temp <- tempfile()
download.file("https://www.stats.ox.ac.uk/~snijders/siena/LazegaLawyers.zip", temp)
LazLawer <- fread(unzip(temp, files = "ELwork36.dat"))
LasAtt <- fread(unzip(temp, files = "ELattr.dat"))
rm(temp)
LazLawer <- as.matrix(LazLawer)
senior <- LasAtt[,1]
sex <- LasAtt[,3]-1
office <- LasAtt[,4]
LazNet <- network(LazLawer, 
		directed = FALSE)
LazNet %v% "practice" <- LasAtt[,7]-1
LazNet %v% "senior" <- senior/36
LazNet %v% "sex" <- sex
LazNet %v% "office" <- office

m3 <- ergm(LazNet  ~ edges 
		+ nodecov("practice")
		+ nodecov("senior")
		 + match("practice")
		 + match("sex") 
		 + match("office")
		+ gwesp(0.7781, fixed = TRUE),control=control.ergm(MCMC.samplesize=5000,MCMC.interval=5000))
summary(m3)
GCD <- ergm.gcd(LazNet  ~ edges 
		+ nodecov("practice")
		+ nodecov("senior")
		 + match("practice")
		 + match("sex") 
		 + match("office")
		+ gwesp(0.7781, fixed = TRUE),m3)

	
}

convergence.test <- function(formula,ergm.est,multiplication.factor=30,nsim=NULL,burnin=NULL)
{
	require('ergm')
	
	observedNetwork <-  ergm.getnetwork(formula)
	therearemissing <- any(get.edge.attribute(observedNetwork,'na'))
	notergmobject <- class(ergm.est)!='ergm'
	if (therearemissing | notergmobject){
		
		if (therearemissing)
		{
			err.mess <- "missing edges are not handled yet"
		}
		if (notergmobject)
		{
			err.mess <- "ergm.est is not an ergm estimation object"
		}
	 stop(err.mess)
	 }
	 else
	 {
	 	 mod <- ergm.getmodel(formula, observedNetwork)
  n <- network.size(observedNetwork)
  M <- network.edgecount(observedNetwork)
  stat <- ergm.getglobalstats(observedNetwork, mod)
  mrow <- length(stat)
  if (is.null(nsim))
  {
  	nsim <- ergm.est$control$MCMC.samplesize
  }
  if (is.null(burnin))
  {
  	
  	burnin <- ergm.est$control$MCMC.burnin
  }
  
  	thinning <- network.density(observedNetwork)*(1-network.density(observedNetwork))*(n^2)*multiplication.factor
  	test.net <- simulate(formula, coef=ergm.est$coef , nsim=nsim, statsonly=TRUE,control=control.simulate(MCMC.burnin=burnin, MCMC.interval=thinning) )
  	
 t.stat <- stat
 SACF <- matrix(0,mrow,1)
 for (k in c(1:mrow))
 {
 	 t.stat[k] <- (t.stat[k]-mean(test.net[,k]))/sd(test.net[,k])
 	  SACF[k] <-  acf(test.net[,k],plot=FALSE,lag.max=1)$acf[2]
 
 }
 
  
con.stats <- cbind(ergm.est$coef ,sqrt(diag(ergm.est$covar)),t.stat,SACF)
colnames(con.stats) <- c('coef','se','convergence t','SACF')

con.stats

	 }


}

test.gcd.plot <- function(MCMC.test.GCD,do.plot=TRUE)
{
	isnotlist <- class(MCMC.test.GCD)!='list'
	hasnotObs <- ('GCDobs' %in% names(MCMC.test.GCD))==FALSE
	hasnotReps <- ('GCDReps' %in% names(MCMC.test.GCD))==FALSE
	if (isnotlist | hasnotObs | hasnotReps )
	{
		stop('argument MCMC.test.GCD needs to be object produced by test.gcd()')

	}
	else
	{
		n <- length(MCMC.test.GCD$GCDobs$GCDMVP)
		extreme.act <- which(MCMC.test.GCD$GCDobs$GCDMVP==max(MCMC.test.GCD$GCDobs$GCDMVP))
		obs.max <- MCMC.test.GCD$GCDobs$GCDMVP[extreme.act]
		if (length(obs.max)>1)
		{
			obs.max <- obs.max[1]
		}
		mmax <-apply(MCMC.test.GCD$GCDReps,2,max)
		gcd.table <- as.data.frame( matrix(NA,n,3) )
		gcd.table[,1] <- names(MCMC.test.GCD$GCDobs$GCDMVP)
		gcd.table[,2] <- MCMC.test.GCD$GCDobs$GCDMVP
		for (i in c(1:n))
		{
			gcd.table[i,3] <- mean(MCMC.test.GCD$GCDReps[i,]<=obs.max)
		}
		names(gcd.table) <- c('Actor','GCDVMP','p-val')
		print(gcd.table)
		if (do.plot)
		{
			F12 <- ecdf(mmax)
			plot(F12,main='CDF for maximal GCD',xlab=expression(S^(n)))
			lines(c(obs.max,obs.max),c(0, F12(obs.max) ) , col='red')
			lines(c(0,obs.max),c(F12(obs.max), F12(obs.max) ) , col='green')
			text(0,F12(obs.max)+.03,F12(obs.max))
			text(obs.max+1,.1,expression({S^(n)} (y[obs]) ) )
			
		}
	
	}
	gcd.table
}