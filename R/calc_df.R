calc_df <- function(bial, populations, outgroup, keep.site.info=FALSE, dxy.table=FALSE, l.smooth){


if(l.smooth){
#print(l.smooth)
## add one ! Laplace ####################################################
pop1 <- bial[populations[[1]],,drop=FALSE]
pop2 <- bial[populations[[2]],,drop=FALSE]

out  <- bial[outgroup, ,drop=FALSE]
# get the ancestral state 
# return(out)
res  <- apply(out,2,function(x){
        vec <- x[!is.na(x)]
	if(length(unique(vec))==1){return(!vec[1])}else{return(NaN)}	
	})
#check if pop1 all NaN sites
pp1  <- apply(pop1,2,function(x){
        return(all(is.na(x)))
	})
#check if pop2 all NaN sites
pp2  <- apply(pop2,2,function(x){
        return(all(is.na(x)))
	})

res[pp1] <- NaN
res[pp2] <- NaN

# add res to pop1 and pop2
res2   <- as.numeric(!res)
bial   <- rbind(bial,res, res2)
populations[[1]] <- c(populations[[1]],dim(bial)[1],dim(bial)[1]-1)  
populations[[2]] <- c(populations[[2]],dim(bial)[1],dim(bial)[1]-1)

#print(populations[[1]])
#print(populations[[2]])
#return(bial)
}
##########################################################################

# Patterson's D
# f

# pop 1: populations[[1]] :
# pop 2: populations[[2]] : 
# pop 3: populations[[3]] : archaic population

if(length(populations)!=3){
stop("This statistic requires 3 populations; the third population is the archaic population")
}

if(outgroup[1]==FALSE || length(outgroup[1])==0){
stop("This statistic needs an outgroup ! (set.outgroup)")
}

# calculate frequencies of the derived alleles.
freqs  <- jointfreqdist(bial,populations,outgroup,keep.all.sites=TRUE)
freqs  <- freqs$jfd
#print(freqs)


# calc the B_d
p      <- freqs[1,]
q      <- freqs[2,]

# dxy for each site and each pop vs archaic pop 3

dxy_pop13 <- apply(bial,2,function(x){
	pop1 <- x[populations[[1]]]
	pop3 <- x[populations[[3]]]
	ones <- sum(pop1==1,na.rm=TRUE) * sum(pop3==0,na.rm=TRUE)
	zero <- sum(pop1==0,na.rm=TRUE) * sum(pop3==1,na.rm=TRUE)
	dxy  <- (ones + zero)/(sum(!is.na(pop1))*sum(!is.na(pop3)))
	return(dxy)		
})

dxy_pop23 <- apply(bial,2,function(x){
	pop2 <- x[populations[[2]]]
	pop3 <- x[populations[[3]]]
	ones <- sum(pop2==1,na.rm=TRUE) * sum(pop3==0,na.rm=TRUE)
	zero <- sum(pop2==0,na.rm=TRUE) * sum(pop3==1,na.rm=TRUE)
	dxy  <- (ones + zero)/(sum(!is.na(pop2))*sum(!is.na(pop3)))
	return(dxy)		
})

dxy_pop12 <- apply(bial,2,function(x){
	pop1 <- x[populations[[1]]]
	pop2 <- x[populations[[2]]]
	ones <- sum(pop1==1,na.rm=TRUE) * sum(pop2==0,na.rm=TRUE)
	zero <- sum(pop1==0,na.rm=TRUE) * sum(pop2==1,na.rm=TRUE)
	dxy  <- (ones + zero)/(sum(!is.na(pop1))*sum(!is.na(pop2)))
	return(dxy)		
})


dxy_pop123 <- apply(bial,2,function(x){
	pop1 <- x[c(populations[[1]],populations[[2]])]
	pop2 <- x[populations[[3]]]
	ones <- sum(pop1==1,na.rm=TRUE) * sum(pop2==0,na.rm=TRUE)
	zero <- sum(pop1==0,na.rm=TRUE) * sum(pop2==1,na.rm=TRUE)
	dxy  <- (ones + zero)/(sum(!is.na(pop1))*sum(!is.na(pop2)))
	return(dxy)		
})


if(dxy.table[[1]]!=FALSE){
#print("Using table")
ids       <- match(dxy_pop12,as.numeric(names(dxy.table[[1]])))
dxy_pop12 <- dxy.table[[1]][ids]

ids       <- match(dxy_pop13,as.numeric(names(dxy.table[[2]])))
dxy_pop13 <- dxy.table[[2]][ids]

ids       <- match(dxy_pop23,as.numeric(names(dxy.table[[3]])))
dxy_pop23 <- dxy.table[[3]][ids]
}else{
#print("Not using table")
}

d13   <- dxy_pop13 #/(dxy_pop13+dxy_pop23)
d23   <- dxy_pop23 #/(dxy_pop13+dxy_pop23)	
d12   <- dxy_pop12
d123  <- dxy_pop123

#d12   <- abs(p-q)
#d13   <- abs(p-freqs[3,])
#d23   <- abs(q-freqs[3,])

alpha   <- (d12-d23)^2
beta    <- (d12-d13)^2
theta   <- (d13-d23)^2 

root <- freqs[3,]
root[root!=0] <- 1

# ---------------------------------
# google docs version
#BABA    <- ( (p     + alpha )     * ((1-q) + beta) ) * freqs[3,]
#ABBA    <- ( ((1-p) + alpha )     * (q     + beta) ) * freqs[3,] 


#Bd-fraction

BABA  <-  p*d23
ABBA  <-  q*d13  


DENOM   <- ABBA + BABA 

######################################################

#print(BABA)
sum_ABBA <- sum(ABBA,na.rm=TRUE)
sum_BABA <- sum(BABA,na.rm=TRUE) 

# 
D        <- (sum_ABBA - sum_BABA)/sum(DENOM, na.rm=TRUE)#(sum_ABBA + sum_BABA) #/valid.sites

### Bayes-----------------------------------
# sample sizes
n1 <- length(populations[[1]])
n2 <- length(populations[[2]])
n3 <- length(populations[[3]])

B_ABBA <- (1-p)*q*freqs[3,]
B_BABA <- p*(1-q)*freqs[3,]
B_BBAA <- p*q*(1-freqs[3,])

scale      <- (n1+n2+n3)/3

alpha_ABBA <- scale   + sum( q*d13  * scale, na.rm=TRUE) 
alpha_BABA <- scale   + sum( p*d23  * scale, na.rm=TRUE)
beta_BBAA  <- scale   + sum( B_BBAA * scale, na.rm=TRUE)

#test
#beta_ABBA  <- scale   + sum( B_ABBA * scale, na.rm=TRUE)
#beta_BABA  <- scale   + sum( B_BABA * scale, na.rm=TRUE)
#

L_abba      <- lbeta(alpha_ABBA,beta_BBAA) #(lgamma(alpha_ABBA)*lgamma(beta_BBAA))/lgamma(alpha_ABBA+beta_BBAA)
L_baba      <- lbeta(alpha_BABA,beta_BBAA) #(lgamma(alpha_BABA)*lgamma(beta_BBAA))/lgamma(alpha_BABA+beta_BBAA)


#D_bayes            <- L_abba/L_baba

D_bayes             <- D
D_bayes[D==0]       <- 1
D_bayes[D >0]       <- 1 + exp((L_abba/L_baba)) - exp(1)
D_bayes[D <0]       <- -(1 + exp((L_baba/L_abba)) - exp(1))

#D_bayes            <- sum_ABBA/sum_BABA

#-------------------------------------------

#if(keep.site.info){
 D_site <- (ABBA - BABA)/DENOM #(ABBA + BABA)
 ABBA_site <- ABBA
 BABA_site <- BABA

 x      <- D_site
 Bd_dir <- (sum((d12-d23), na.rm=TRUE) + sum((d12-d13), na.rm=TRUE))

#D <- mean(D_site, na.rm=TRUE)

#}else{
# D_site      <- NULL
# ABBA_site   <- NULL
# BABA_site   <- NULL
#}

# calc f_BD
freqs23    <- freqs[2:3,,drop=FALSE]
maxfreqs23 <- apply(freqs23,2,max)


# f_d denominator
maxBABA <- (d23 * p     + d13 * (1-maxfreqs23))  *maxfreqs23
maxABBA <- (d23 * (1-p) + d13 * maxfreqs23)      *maxfreqs23

sum_maxABBA <- sum(maxABBA, na.rm=TRUE)
sum_maxBABA <- sum(maxBABA, na.rm=TRUE)

f <- (sum_ABBA - sum_BABA)/(sum_maxABBA - sum_maxBABA)

#D3
D3 <- (d23-d13)/(d23+d13)

return(list(Bd_dir=Bd_dir, D=D, D_bayes=D_bayes, f=f, D_site=D_site, ABBA=ABBA_site, BABA=BABA_site, alpha_ABBA=alpha_ABBA, alpha_BABA=alpha_BABA, beta_BBAA=beta_BBAA,D3=D3))

}
