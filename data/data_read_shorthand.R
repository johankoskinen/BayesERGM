#read in networks demo
rebNet <-read.csv('StysRebels.csv',header=TRUE,stringsAsFactors=FALSE) # n
#rebNet[rebNet$relation %in% c(1,2,3,5,6,7,1,2),] <- NA
# rebNet <- rebNet[rebNet$relation %in% c(4,8),]
DRC.rebels <- as.network(rebNet,directed=FALSE,matrix.type="edgelist",ignore.eval=FALSE)
#DRC.rebels[DRC.rebels %in% c(1,2,3,5,6,7,1,2) ] <- 0
#DRC.rebels[DRC.rebels %e%'relation' %in% c(1,2,3,5,6,7,1,2) ]
DRC.mat <- as.sociomatrix(DRC.rebels,attrname='relation')
DRC.coop <- DRC.mat
DRC.coop[DRC.coop %in% c(2,3,5,6,7) ] <- 0
DRC.coop.net <- as.network(DRC.coop,directed=FALSE)

plot(DRC.coop.net,displaylabels =TRUE)
