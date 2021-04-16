# miss_ERGM_Bayes
miss_ERGM_Bayes <- function(formula, 
                            tot.iterations = NULL , 
                            aux.var.iter = NULL , 
                            imput.var.iter = NULL, 
                            do.test.sim = FALSE, 
                            sigma.epsilon = NULL, 
                            tuning.const = 1, 
                            thinning = 1, 
                            burnin = 1 , 
                            save.nets = FALSE , 
                            num.nets = NULL, 
                            printFreq = NULL, 
                            MultFac.aux =30 ,
                            MultFac.miss =30, 
                            theta= NULL, 
                            diagonal.proposal=FALSE,
                            save.name=NULL,
                            warm.burnin=10000,
                            warm.thin = 5000,
                            start.net =NULL)
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
	do.calculate.propovar <- FALSE
	do.diag.covar <- FALSE
	if (is.null(sigma.epsilon))
	{
		do.calculate.propovar <- TRUE
		#print('sigma.epsilon: a proposal variance or variance covariance is needed')
	}
	
	aux.iters <- 1
	# extract model, network, and observed stats
	y <- ergm.getnetwork(formula)
	model <- ergm_model(formula, y)
	sy <- summary(formula)
	p <- length( sy )
	#n <- dim(y[,])[1]
	n <- network.size(y)
	#numMiss <- sum(is.na(y[,] ) )
	numMiss <- network.naedgecount(y)
	if (numMiss==0)
	{
		print('no missing ties - pls use other routine')
	}
  	#densObs <- sum(y[,],na.rm = TRUE)/(n*(n-1)-numMiss)
  	densObs <- network.edgecount(y)/(n*(n-1)-numMiss)
	
	if ( is.null(aux.var.iter) )
	{
		
		aux.var.iter <- round(max(MultFac.aux * densObs * (1 - densObs) * n * n,1000))
		cat('burnin auxiliary var draws: ',aux.var.iter,' with multiplication factor: ',MultFac.aux,'\n')
	}
	if ( is.null(imput.var.iter) )
	{
	imput.var.iter <- round(max(MultFac.miss * densObs * (1 - densObs) * numMiss,1000))
	cat('burnin missing var draws: ',imput.var.iter,' with multiplication factor: ',MultFac.miss,'\n')
	}
	
	theta.save.ind <- seq( from = burnin , to = tot.iterations, by =  thinning )
	Post.theta <- matrix(0, length(theta.save.ind) , p )
	
	if (save.nets)
	{
		if (is.null(num.nets) || num.nets >(tot.iterations- burnin))
		{
			print('the number of networks you pick up must be leq to number of paramters')
			
		}
		nets.save.ind <- round(seq( from = burnin , to = tot.iterations, length.out = num.nets ))
		aux.nets <- vector('list',length(nets.save.ind))
		mix.nets <- vector('list',length(nets.save.ind))
	}
	# ------- initialise main part for updating theta -------- 
	if (do.calculate.propovar==FALSE && length(sigma.epsilon)==1)
	{
		do.diag.covar <- TRUE
	}
	if (do.calculate.propovar==FALSE && length(sigma.epsilon)>1)
	{
		do.diag.covar <- FALSE
	}
	if ( do.calculate.propovar==FALSE && do.diag.covar )
	{
		sigma.epsilon <- diag( sigma.epsilon ,p)
	}
	
	if ( do.calculate.propovar==FALSE && do.diag.covar==FALSE )
	{
		sigma.epsilon <- sigma.epsilon*tuning.const
	}
	
	if (do.calculate.propovar)
	{
		
		if (is.null(theta)){
  		thetaL <- - log(1 / densObs - 1)
  		theta <- matrix(0,1,p)
  		theta[1,1] <- thetaL
		}
	
	  
	  
	  actual.burning <- ceiling(warm.burnin/ warm.thin)
	  totalSimNets <-  actual.burning+1000

	  
	  
	  
	  if (is.null(start.net))
	  {
	    start.net <- y
	  }
	  
	cat('simulating ',totalSimNets,' networks with thinning ',warm.thin,', total iterations: ',totalSimNets*warm.thin ,'\n')
		tempStats <- simulate( formula, coef = theta, 
                  output = "stats", basis = start.net ,nsim = totalSimNets, control = control.simulate(MCMC.burnin = 1,MCMC.prop.weights='default',MCMC.interval = warm.thin) )
    # plot(ts(tempStats[,c(1:10)]))  
		tempStats <- tempStats[actual.burning:totalSimNets,]
		
		sigma.epsilon <-  solve(cov(tempStats)) *tuning.const 
         #sigma.epsilon <-  sqrtm(cov(tempStats))*tuning.const
         if (diagonal.proposal)
         {
         	sigma.epsilon <- diag(diag(sigma.epsilon ))
         }
             cat('warming phase done\n')
	}
	
	# define missing based on extracted network
	missingTies <- matrix(0, n,n)
	#missingTies[is.na( y[,] )] <- 1
	missingTies[is.na(as.matrix.network(y))] <- 1
	Ass <- dyad.census(y)[2]
	if (Ass==0){
	missingTies <- as.edgelist(as.network(missingTies,directed=FALSE), n = n)
	}
	if (Ass>0){
	  missingTies <- as.edgelist(as.network(missingTies,directed=TRUE), n = n)
	}
	
	# if we are checking consistency of algorithm:
	if (do.test.sim)
	{
		# keep copy of original observed network
		y.obs.copy <- y
		# and stats
		sy.obs.copy <- sy
	}
	stats.curr <- sy
	# get control object from 'simulate'
	if (is.null(theta)){
	theta <- matrix(0,1,p)
	}
	
	cat('initialise c-list etc \n')
	#browser()
	inits <- simulate(formula,coef =  theta,
				basis = y , nsim = 1,control = control.simulate(MCMC.burnin = 1 ),do.sim=FALSE) # tried to add, MCMC.prop.weights = "TNT"), 
					
	# initialise C-side       
	Clist <- ergm.Cprepare(y, model)
	# initiallise proposal
	myprop <- ergm_proposal(~fixallbut(missingTies), arguments=NULL,
						nw= y , weights = "default", class = "c",
						reference = ~Bernoulli,
						response = NULL)
# ------- initialise main part for updating theta -------- 
	y.no.na <- y
	sy.no.na <- sy
	ClistMain <- ergm.Cprepare(y.no.na, model)
	control <- control.ergm(MCMC.burnin = aux.iters, MCMC.interval = 1, 
        MCMC.samplesize = 1)
    proposal <- ergm_proposal(object = ~., constraints = ~., 
        arguments = control$MCMC.prop.args, nw = y.no.na)  
 # ------- hot loop  --------  #        
 k <- 1
 k1 <- 1
 k2 <- 1
 curr.aux.net <- y.no.na
 #curr.miss.net <- y
 curr.miss.net <- as.matrix.network(y)
 
 curr.miss.net[is.na(curr.miss.net)] <- 0
 curr.miss.net <- as.network( curr.miss.net, directed=FALSE)
 save(curr.aux.net,curr.miss.net,inits,file=paste('missERGMinits',save.name,'par',p,'.RData',sep=''))
 gio.test <- FALSE
 cat('Start of estimation \n')
        for (iter in c(1:tot.iterations))
        {
        	# ------- propose new parameters  -------- 
        	theta1 <- rmvnorm(1,theta, sigma = sigma.epsilon)[1,]
 # ------- simulate auxiliary variable  --------     
        
	delta <- ergm_MCMC_slave(Clist = ClistMain, proposal =proposal, samplesize = 1,
               eta = theta1, control = control, verbose = FALSE,
               burnin = aux.var.iter)
	
	
# ------- the following test that code is consistent:   -------- 
# check that the statistics for the auxiliary network are calculated correctly
             
   
# ------- back to the main algorithm   -------- 

# ------- statistics for simulate auxiliary variable  --------     
	sim.stats <- delta$s+	sy.no.na
	diff.stats <- sim.stats - stats.curr
	
	if (gio.test ==TRUE)
	{
	 test.start.net <- update( y.no.na , cbind(delta$newnwtails[2:(delta$newnwtails[1]+1) ],
	                                          delta$newnwheads[2:(delta$newnwheads[1]+1)] ) , matrix.type = 'edgelist' ) 
	 form <- as.formula(paste("test.start.net ~", as.character(formula)[3]))
	 test.start.stats <- summary(form)
	 if (sim.stats[1] != test.start.stats[1] )
	 {
	   print('missmatch between simulated stats and stats of simulated network')
	 }
	}
	
	 beta <- (theta - theta1) %*% t(diff.stats) 
        if (beta >= log(runif(1))) {
          
          theta <- theta1
          if (save.nets )
          {
          	
          	if (sim.stats[1]==0)
          	{
          		y.no.na[,]<-0
          	}
            # update the auxilairy network
            if (sim.stats[1]>0)
          	{
          		if (sum(delta$newnwtails)>0)
 		{
 			#print('no ties to add to network auxiliary')
 			
 		#	browser()
 	

          	curr.aux.net <- update( y.no.na , cbind(delta$newnwtails[2:(delta$newnwtails[1]+1) ],
 		delta$newnwheads[2:(delta$newnwheads[1]+1)] ) , matrix.type = 'edgelist' ) 
          	# overwrite the current one
          	y.no.na <- update( y.no.na , cbind(delta$newnwtails[2:(delta$newnwtails[1]+1) ],
          	                                  delta$newnwheads[2:(delta$newnwheads[1]+1)] ) , matrix.type = 'edgelist' ) 
          	                                  }
          
          	}
          	
          	# these lines work for calulating statistics for new network
          	form <- as.formula(paste("curr.aux.net ~", as.character(formula)[3]))
          	curr.aux.stats <- summary(form)# this works
          	# let's see if we can just overwrite y
          
          ClistMain <- ergm.Cprepare(y.no.na, model)
          control <- control.ergm(MCMC.burnin = aux.iters, MCMC.interval = 1, 
                                  MCMC.samplesize = 1)
          proposal <- ergm_proposal(object = ~., constraints = ~., 
                                    arguments = control$MCMC.prop.args, nw = y.no.na) 
          sy.no.na <- curr.aux.stats
          # now the main algorithm starts in y with sy current aux
          # simulate and see where simulation starts from 
          # delta.test <- ergm_MCMC_slave(Clist = ClistMain, proposal =proposal, samplesize = 1,
          #                         eta = theta1, control = control, verbose = FALSE,
          #                         burnin = aux.var.iter)
         # calulculate two new networks:
          # curr.aux.net.1 <- update( y.no.na , cbind(delta.test$newnwtails[2:(delta.test$newnwtails[1]+1) ],
          #                                  delta.test$newnwheads[2:(delta.test$newnwheads[1]+1)] ) , matrix.type = 'edgelist' )
          # 
          # form <- as.formula(paste("curr.aux.net.1 ~", as.character(formula)[3]))
          # curr.aux.stats.1 <- summary(form)#
          #curr.aux.stats.1.alt <-  delta.test$s+	sy# should be different because sy has not been overwritten
          #curr.aux.stats.1.alt2 <-  delta.test$s+		curr.aux.stats # identical to  curr.aux.stats.1, so started in updated y
          ###
          	
          	
          	
          }
        }

# ------- simulate missing data  --------  		
 	delta <- ergm_MCMC_slave(Clist = Clist, proposal = myprop, 
 			samplesize = 1,
 			eta = theta, control = inits$control,
 			verbose = FALSE, burnin = imput.var.iter)
# ------- statistics for simulate missing data  --------   			
 	stats.curr <- delta$s+	sy
 	if (save.nets && (iter %in% nets.save.ind) ){
 		
 		if (sum(delta$newnwtails)==0)
 		{
 			print('no ties to add to network missing')
 		}
 	curr.miss.net <- update( y , cbind(delta$newnwtails[2:(delta$newnwtails[1]+1) ],
 		delta$newnwheads[2:(delta$newnwheads[1]+1)] ) , matrix.type = 'edgelist' ) 
 	
 	}
# ------- the following test that code is consistent:   -------- 
# check that the statistics for the imputed network are calculated correctly
# and that constraint is respected 	
 	if (do.test.sim)
 	{
 		netone <- update( y , cbind(delta$newnwtails[2:(delta$newnwtails[1]+1) ],
 		delta$newnwheads[2:(delta$newnwheads[1]+1)] ) , matrix.type = 'edgelist' ) 
 		discrepancy <- sum(  netone[!is.na(y[,])] !=  y[!is.na(y[,])]) 
 		print( summary( netone ) )
 		if (discrepancy != 0)
	{
		print(paste('missing data algorithm does not respect constraint, diff: ', discrepancy) )
	}
		form <- as.formula(paste("netone ~", as.character(formula)[3]))
		netone.stats <- summary(form)
		netone.stats.alt <- delta$s+	sy.obs.copy
		discrepancy <- sum( netone.stats != netone.stats.alt )	
		if (discrepancy != 0)
		{
			cat('stats differ by ',netone.stats[1]-netone.stats.alt[1],'\n')
			cat('and ',netone.stats[2]-netone.stats.alt[2],'\n')	
	
		}


 		
 	}
# ------- if the test show that inconsistencies the algorithm is wrong  --------  	
 			if (iter %in% theta.save.ind){
 				 			Post.theta[k, ] <- theta	
 				 			k <- k +1
 			}
 			
 			if (save.nets && (iter %in% nets.save.ind) ){
 				 			aux.nets[[k1]] <- curr.aux.net
 				 			#cat('perform test of update\n')
 				 			#browser()
 				 			
 				 			k1 <- k1+1
 				 		#	browser()
 				 			mix.nets[[k2]] <- curr.miss.net
 		k2 <- k2 +1 
 			}
 			
 	if ((verbose=TRUE) && (iter %% printFreq == 0) )
      {
  
      	cat('at iteration ', iter,' out of ',tot.iterations, ' parameters:\n')
      	cat(names(sy),'\n')
      	cat(round(theta,3))
      	cat('\n')
      	
      	save(Post.theta, aux.nets ,mix.nets,file=paste('missERGM',save.name,'par',p,'.RData',sep=''))
      	if (save.nets)
      	{
      		
      	
      		cat('num auxiliary ties: ',network.edgecount(curr.aux.net),' and imputed ties: ',  network.edgecount(curr.miss.net ),'\n' )
      	}
      	}
      	
	}

colnames(Post.theta) <- names(sy)
out.put <- list(Post.theta= Post.theta)
if (save.nets)
{
	if (network.naedgecount(aux.nets[[1]])>0)
	{
	  aux.nets[[1]] <- aux.nets[[2]]
	}
out.put$aux.nets <- aux.nets
	out.put$mix.nets <- mix.nets
}	

	out.put	
}


creatTab <- function(mat, thinning = 1, burnin=1)
{
  p <- dim(mat)[2]
  N <- dim(mat)[1]
  
  myTab <- as.data.frame(matrix(0,p,5))
  colnames(myTab) <- c('effect','mean','sd','CI 0.025','CI 0.975')
  useThese <- seq(burnin,N, by=  thinning )
  newN <- length(useThese )
  for (k in c(1:p))
  {
    datSorted <- sort(mat[ useThese,k])
    myTab[k,1] <- colnames(mat)[k]
    myTab[k,2] <- round(mean(datSorted),3)
    myTab[k,3] <- round(sd(datSorted),3)
    myTab[k,4] <- round(datSorted[ceiling(0.025*   newN)], 3)
    myTab[k,5] <- round(datSorted[floor( (1- 0.025) *   newN)], 3)
    
  }
  
  myTab
}


plot.gof.degree <- function(gofdegs,max.deg.plot.ind = 11)
{
  N <- dim( gofdegs$aux.indeg )[1]
  n <- dim( gofdegs$aux.indeg )[2]
  
  lowerp <- ceiling(0.025*   N)
  upperp <- floor( (1- 0.025) *   N)

  par( mfrow = c(1,2))
  plot( c(0:(max.deg.plot.ind-1) ), seq(0,1,length.out = max.deg.plot.ind), type='n'  ,
        bty ='n',yaxt='n', ylab=NA, xlab='in-degree')
  for (k in c(1:max.deg.plot.ind ))
  {
    sort.dat <- sort( gofdegs$aux.indeg[,k] )
    l_ci <- sort.dat [lowerp]
    u_ci <- sort.dat [upperp]
    # lines( c(k-1, k-1 ), range(gofdegs$aux.indeg[,k]) ,  lwd = 3 )
    lines( c(k-1, k-1 ), c(l_ci,u_ci) ,  lwd = 3 )
    gofdegs$mix.indeg[,k] <- sort( gofdegs$mix.indeg[,k] )
  }
  lines( c(0:(max.deg.plot.ind-1) ),  gofdegs$mix.indeg[lowerp, c(1: max.deg.plot.ind ) ] , col ='grey')
  lines( c(0:(max.deg.plot.ind-1) ),  gofdegs$mix.indeg[upperp, c(1: max.deg.plot.ind ) ] , col ='grey')
  
  plot( c(0:(max.deg.plot.ind-1) ), seq(0,1,length.out = max.deg.plot.ind),
        type='n'  , bty ='n',yaxt='n', ylab=NA, xlab='out-degree')
  for (k in c(1:max.deg.plot.ind ))
  {
    sort.dat <- sort( gofdegs$aux.oudeg[,k])
    l_ci <- sort.dat [lowerp]
    u_ci <- sort.dat [upperp]
    #lines( c(k-1, k-1 ), range(gofdegs$aux.oudeg[,k]) , lwd = 3 )
    lines( c(k-1, k-1 ), c(l_ci,u_ci) ,  lwd = 3 )
    gofdegs$mix.oudeg[,k] <- sort( gofdegs$mix.oudeg[,k] )
  }
  lines( c(0:(max.deg.plot.ind-1) ),  gofdegs$mix.oudeg[lowerp, c(1: max.deg.plot.ind ) ] , col ='grey' )
  lines( c(0:(max.deg.plot.ind-1) ),  gofdegs$mix.oudeg[upperp, c(1: max.deg.plot.ind ) ], col ='grey' )
  
}

gof.degree <- function(mat, directed=TRUE)
{
  N <- length(mat$aux.nets)
  n <-  network.size(mat$aux.nets[[1]])
  aux.indeg <- matrix(0,N,n)
  aux.oudeg <- matrix(0,N,n)
  
  mix.indeg <- matrix(0,N,n)
  mix.oudeg <- matrix(0,N,n)
  
  for ( k in c(1:N))
  {
    
    degtemp <- degree(dat=mat$aux.nets[[k  ]],cmode="indegree")
    aux.indeg[k,] <- degree.table(bin= degtemp,min.deg = 0, max.deg = n-1)/n
    
    degtemp <- degree(dat=mat$aux.nets[[k]],cmode="outdegree")
    aux.oudeg[k,] <- degree.table(bin= degtemp,min.deg = 0, max.deg = n-1)/n
    
    degtemp <- degree(dat=mat$mix.nets[[k]],cmode="indegree")
    mix.indeg[k,] <- degree.table(bin= degtemp,min.deg = 0, max.deg = n-1)/n
    
    degtemp <- degree(dat=mat$mix.nets[[k]],cmode="outdegree")
    mix.oudeg[k,] <- degree.table(bin= degtemp,min.deg = 0, max.deg = n-1)/n
    
    
  }
  
  gofdegs <- list(aux.indeg=aux.indeg,aux.oudeg=aux.oudeg, mix.indeg= mix.indeg,  mix.oudeg = mix.oudeg)
  
  gofdegs
}

degree.table <- function(bin,min.deg = 0, max.deg = NULL)
{
  if (is.null(max.deg))
  {
    max.deg = max(bin)
  }
  
  nbins <- max.deg - min.deg+1
  freq <- matrix(0,1,nbins)
  curr.deg <- min.deg
  for ( k in c(1:(nbins) ) )
  {
    freq[k] <- sum(bin==curr.deg)
    curr.deg <- curr.deg +1
    
  }
  freq
  
}


get.geodist.net <- function( net , max.path = NULL )
{
  myDist <- geodist(net)$gdist
  geoDistTable <- matrix(0, max.path +1 , 1 )
  for (k in c(1:max.path ))
  {
    geoDistTable[ k ] <- sum( ( myDist[ upper.tri( myDist  ) ] )==k   )
  }
  geoDistTable[ max.path + 1 ] <- sum(  is.infinite( myDist[ upper.tri( myDist  ) ] ))
  
  
  geoDistTable
}

gof.dist <- function(mat, directed=TRUE , max.path = 20)
{
  N <- length(mat$aux.nets)
  n <-  network.size(mat$aux.nets[[1]])
  aux.dist <- matrix(0,N,max.path +1 )
  mix.dist <- matrix(0,N, max.path +1)
  num.dyads <- n*(n-1)/2
  
  for ( k in c(1:N))
  {
    
    
    aux.dist[k,] <- get.geodist.net( net = mat$aux.nets[[k  ]] , max.path = max.path )/num.dyads
    
    mix.dist[k,] <- get.geodist.net( net = mat$mix.nets[[k]] , max.path = max.path )/num.dyads
    
    
    
    
  }
  
  gofdist <- list( aux.dist = aux.dist ,  mix.dist =  mix.dist)
  
  gofdist
}


plot.gof.dist <- function(gofdist,max.dist.plot.ind = 11,do.inf.dist=FALSE)
{
  N <- dim( gofdist$aux.dist )[1]
  n <- dim( gofdist$aux.dist )[2]
  lowerp <- ceiling(0.025*   N)
  upperp <- floor( (1- 0.025) *   N)
  
  max.freq <- max( max(gofdist$aux.dist[, 1:(max.dist.plot.ind-1)]), max(gofdist$mix.dist[, 1:(max.dist.plot.ind-1)]) )
  par(mfrow=c(1,1))
  plot( c(1:(max.dist.plot.ind) ), seq(0,max.freq,length.out = max.dist.plot.ind), type='n'  ,
        bty ='n',yaxt='n', ylab=NA, xlab='geodesic distance')
  for (k in c(1: (max.dist.plot.ind-1) ))
  {
    sort.dat <- sort( gofdist$aux.dist[,k] )
    l_ci <- sort.dat [lowerp]
    u_ci <- sort.dat [upperp]
    # lines( c(k-1, k-1 ), range(gofdegs$aux.indeg[,k]) ,  lwd = 3 )
    lines( c(k, k ), c(l_ci,u_ci) ,  lwd = 3 )
    
    #lines( c(k, k ), range(gofdist$aux.dist[,k]) ,  lwd = 3 )
    gofdist$mix.dist[,k] <- sort( gofdist$mix.dist[,k] )
  }
  lines( c(1:(max.dist.plot.ind-1) ), gofdist$mix.dist[lowerp, c(1: (max.dist.plot.ind-1) ) ] , col ='grey')
  lines( c(1:(max.dist.plot.ind-1) ), gofdist$mix.dist[upperp, c(1:(max.dist.plot.ind-1)  ) ] , col ='grey')
  #lines( c(max.dist.plot.ind , max.dist.plot.ind ), range(gofdist$aux.dist[,n]) ,  lwd = 5 )
  #lines( c(max.dist.plot.ind , max.dist.plot.ind), range(gofdist$mix.dist[,n]) ,  lwd = 3 , col = 'grey')
}


gof.obs.degree <- function(obs.net,mat,is.resopndent ,directed=TRUE)
{
  N <- length(mat$aux.nets)
  n <-  network.size(mat$aux.nets[[1]])
  obs.indeg <- matrix(0,1,n)
  obs.oudeg <- matrix(0,1,n)
  
  mix.indeg <- matrix(0,N,n)
  mix.oudeg <- matrix(0,N,n)
  
  n.obs <- sum(is.resopndent)
  
  degtemp <- degree(dat=obs.net,cmode="indegree")
  aux.indeg <- degree.table(bin= degtemp[is.resopndent],min.deg = 0, max.deg = n-1)/ n.obs
  
  degtemp <- degree(dat=obs.net,cmode="outdegree")
  aux.oudeg <- degree.table(bin= degtemp[is.resopndent],min.deg = 0, max.deg = n-1)/ n.obs
  
  for ( k in c(1:N))
  {
    
   
    
    degtemp <- degree(dat=mat$mix.nets[[k]],cmode="indegree")
    mix.indeg[k,] <- degree.table(bin= degtemp,min.deg = 0, max.deg = n-1)/n
    
    degtemp <- degree(dat=mat$mix.nets[[k]],cmode="outdegree")
    mix.oudeg[k,] <- degree.table(bin= degtemp,min.deg = 0, max.deg = n-1)/n
    
    
  }
  
  gofdegs <- list(aux.indeg=aux.indeg,aux.oudeg=aux.oudeg, mix.indeg= mix.indeg,  mix.oudeg = mix.oudeg)
  
  gofdegs
}

numESP <- function(networkObject=NULL, directed=TRUE)
{
  
  ADJ <- as.matrix.network(networkObject)
  n <- dim(ADJ)[1]
  Fnbd2 <- matrix(0,n,1)
  if (directed == FALSE){
  for (i in c(1:(n-1)))
  {
    for (j in c((i+1):n))
    {
      if (ADJ[i,j]==1)
      {
        nbd2 <- sum( ADJ[i,] * ADJ[j,] ) -1 
        Fnbd2[nbd2 +1 ] <- Fnbd2[nbd2 +1 ] +1
        
      }
    }
  }
  }
  
  if (directed == TRUE){
    for (i in c(1:n))
    {
      for (j in c(1:n))
      {
        if (ADJ[i,j]==1)
        {
          nbd2 <- sum( ADJ[i,] * ADJ[j,] ) -1 
          Fnbd2[nbd2 +1 ] <- Fnbd2[nbd2 +1 ] +1
          
        }
      }
    }
  }
  Fnbd2
  
}


getGOFESP <- function(networklist=NULL,directed=TRUE )
{
  # simDegs[,k]
  NNets <- length(networklist)
  n <- dim(as.matrix.network(networklist[[1]]))[1]
  simDegs <- matrix(0,NNets, n)
  for (k in c(1:NNets))
  {
    simDegs[k, ] <- numESP(networkObject=networklist[[k]])
    
    
  }
  simDegs
}

gof.esp <- function(ans=NULL,directed=TRUE)
{
  gof.esp.aux <- getGOFESP(networklist=ans$aux.nets,directed=directed )
  gof.esp.mix <- getGOFESP(networklist=ans$mix.nets,directed=drirected )
  gofesp <- list( aux.esp =  gof.esp.aux ,  mix.esp =  gof.esp.mix)
  gofesp
  
}

plot.gof.esp <- function(gofesp,max.dist.plot.ind = 11)
{
  N <- dim( gofesp$aux.esp )[1]
  n <- dim( gofesp$mix.esp )[2]
  lowerp <- ceiling(0.025*   N)
  upperp <- floor( (1- 0.025) *   N)
  
  max.freq <- max( max(gofesp$aux.esp[, 1:(max.dist.plot.ind-1)]), max(gofesp$mix.esp[, 1:(max.dist.plot.ind-1)]) )
  

  par(mfrow=c(1,1))
  plot( c(0:(max.dist.plot.ind-1) ), seq(0,max.freq,length.out = max.dist.plot.ind), type='n'  ,
        bty ='n',yaxt='n', ylab=NA, xlab='edgewise shared partner')
  for (k in c(1: (max.dist.plot.ind) ))
  {
    sort.dat <- sort( gofesp$aux.esp[,k] )
    l_ci <- sort.dat [lowerp]
    u_ci <- sort.dat [upperp]
    # lines( c(k-1, k-1 ), range(gofdegs$aux.indeg[,k]) ,  lwd = 3 )
    lines( c(k-1, k -1), c(l_ci,u_ci) ,  lwd = 3 )
    
    #lines( c(k, k ), range(gofdist$aux.dist[,k]) ,  lwd = 3 )
    gofesp$mix.esp[,k] <- sort( gofesp$mix.esp[,k] )
  }
  lines( c(0:(max.dist.plot.ind-1) ), gofesp$mix.esp[lowerp, c(1: (max.dist.plot.ind) ) ] , col ='grey')
  lines( c(0:(max.dist.plot.ind-1) ), gofesp$mix.esp[upperp, c(1:(max.dist.plot.ind)  ) ] , col ='grey')
  #lines( c(max.dist.plot.ind , max.dist.plot.ind ), range(gofdist$aux.dist[,n]) ,  lwd = 5 )
  #lines( c(max.dist.plot.ind , max.dist.plot.ind), range(gofdist$mix.dist[,n]) ,  lwd = 3 , col = 'grey')
}

