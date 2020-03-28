calc_RNDmin <- function(bial, populations, outgroup){

if(outgroup[1]==FALSE || length(outgroup[1])==0){
stop("This statistic needs an outgroup ! (set.outgroup)")
}

if(length(populations)!=2){
stop("This statistic requires 2 populations")
}


# First pop vs outgroup
n1 	<- length(populations[[1]])
n2 	<- length(outgroup)

diff    <- rep(NaN, n1*n2) 

# Get the mean distance sequence between pop1 and outgroup 
zz <- 1 
for(xx in 1:n1){

	seq1 <- bial[populations[[1]][xx],]	
	
	for(yy in 1:n2){

		seq2     <- bial[outgroup[yy],]
		diff[zz] <- sum(seq1 != seq2, na.rm=TRUE)
		zz	 <- zz + 1
	}
	
}

d1O <- mean(diff, na.rm=TRUE)


# Second pop vs outgroup
n1 	<- length(populations[[2]])
n2 	<- length(outgroup)

diff    <- rep(NaN, n1*n2) 

# Get the mean distance sequence between pop2 and outgroup 
zz <- 1 
for(xx in 1:n1){

	seq1 <- bial[populations[[2]][xx],]	
	
	for(yy in 1:n2){

		seq2     <- bial[outgroup[yy],]
		diff[zz] <- sum(seq1 != seq2, na.rm=TRUE)
		zz	 <- zz + 1
	}
	
}

d2O <- mean(diff, na.rm=TRUE)


# pop 1 vs pop 2 
n1 	<- length(populations[[1]])
n2 	<- length(populations[[2]])

diff    <- rep(NaN, n1*n2) 
ppp1    <- rep(NaN, n1*n2)
ppp2    <- rep(NaN, n1*n2)

# Get the min distance sequence between pop1 and pop2 
zz <- 1 
for(xx in 1:n1){

	seq1 <- bial[populations[[1]][xx],]	
	
	for(yy in 1:n2){

		seq2     <- bial[populations[[2]][yy],]
		diff[zz] <- sum(seq1 != seq2, na.rm=TRUE)
                ppp1[zz] <- xx
		ppp2[zz] <- yy
		zz	 <- zz + 1
	}
	
}

RNDmin  <- min(diff, na.rm=TRUE)/((d1O+d2O)/2) 
Gmin    <- min(diff, na.rm=TRUE)/mean(diff, na.rm=TRUE)

return(list(RNDmin=RNDmin,Gmin=Gmin))

}
############################# new #############################################
#id       <- which(diff==min(diff,na.rm=TRUE))[1]
#seqmin1  <- bial[populations[[1]][ppp1[id]],] 
#seqmin2  <- bial[populations[[2]][ppp2[id]],]

# seqmin1 Calc distance to the outgroup 
#diff1    <- rep(NaN, length(outgroup))
#diff2    <- rep(NaN, length(outgroup))
#zz       <- 1	
#for(yy in 1:length(outgroup)){
#
#		out       <- bial[outgroup[yy],]
#		diff1[zz] <- sum(seqmin1 != out, na.rm=TRUE)
#		diff2[zz] <- sum(seqmin2 != out, na.rm=TRUE)
#                zz	  <- zz + 1
#}
#	
#ONE <- sum(diff1, na.rm=TRUE)
#TWO <- sum(diff2, na.rm=TRUE)

#ONEX <- ONE/(ONE+TWO)
#TWOX <- TWO/(ONE+TWO)

#RNDmin <- (ONEX*RNDmin - TWOX*RNDmin)
############################# end new #########################################


