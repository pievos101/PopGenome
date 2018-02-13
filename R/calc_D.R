calc_D <- function(bial, populations, outgroup, keep.site.info=FALSE){

# Patterson's D
# f

# pop 1: populations[[1]] :
# pop 2: populations[[2]] : 
# pop 3: populations[[3]] : archaic population

ARCHAIC <- TRUE

if(length(populations)!=3){
warning("The original statistic requires 3 populations; the third population is the archaic population")
warning("Proceeding. Ignoring the arachaic population !")
ARCHAIC <- FALSE
}

if(outgroup[1]==FALSE || length(outgroup[1])==0){
stop("This statistic needs an outgroup ! (set.outgroup)")
}

# calculate frequencies of the derived alleles.
freqs  <- jointfreqdist(bial,populations,outgroup,keep.all.sites=TRUE)
freqs  <- freqs$jfd
#print(freqs)

###!!!!!!!!!!!!!!!!!!!!
#freqs[freqs==0 | freqs==1] <- NaN
###!!!!!!!!!!!!!!!!!!!!

if(ARCHAIC){
# calc D
ABBA   <- (1-freqs[1,])*freqs[2,]*freqs[3,]
#print(ABBA)
BABA   <- freqs[1,]*(1-freqs[2,])*freqs[3,]
#print(BABA)
}else{
# calc D
ABBA   <- (1-freqs[1,])*freqs[2,]
#print(ABBA)
BABA   <- freqs[1,]*(1-freqs[2,])
#print(BABA)
}


sum_ABBA <- sum(ABBA,na.rm=TRUE)
sum_BABA <- sum(BABA,na.rm=TRUE) 

D <- (sum_ABBA - sum_BABA)/(sum_ABBA + sum_BABA)

if(keep.site.info){
 D_site <- (ABBA - BABA)/(ABBA + BABA)
 ABBA_site <- ABBA
 BABA_site <- BABA
}else{
 D_site      <- NULL
 ABBA_site   <- NULL
 BABA_site   <- NULL
}

# calc f
if(ARCHAIC){
freqs23    <- freqs[2:3,,drop=FALSE]
maxfreqs23 <- apply(freqs23,2,max)

#print(maxfreqs23)

maxABBA <- (1-freqs[1,])*maxfreqs23*maxfreqs23
maxBABA <- freqs[1,]*(1-maxfreqs23)*maxfreqs23

#print(maxABBA)
#print(maxBABA)

sum_maxABBA <- sum(maxABBA, na.rm=TRUE)
sum_maxBABA <- sum(maxBABA, na.rm=TRUE)

f <- (sum_ABBA - sum_BABA)/(sum_maxABBA - sum_maxBABA)
}else{
f <- NaN
}

return(list(D=D, f=f, D_site=D_site, ABBA=ABBA_site, BABA=BABA_site))

}
