D_jacknife <- function(Dvalues, D.base, block.size=FALSE){

D_sites 	<- Dvalues[1,]
ABBA_sites 	<- Dvalues[2,]
BABA_sites	<- Dvalues[3,]

D_sites    <- D_sites[!is.na(D_sites)]
ABBA_sites <- ABBA_sites[!is.na(D_sites)]
BABA_sites <- BABA_sites[!is.na(D_sites)]
 
# If block.size is not specified calc the traditional Z-test
if(!block.size){
#after Martin (two tailed z-test)
D_sd  <- sd(D_sites,na.rm=TRUE)
D_err <- D_sd/sqrt(length(D_sites))
D_Z   <- D.base / D_err
D_p   <- 2*pnorm(-abs(D_Z))
return(list(z=D_Z, pval=D_p))
}
###############


n.sites         <- length(D_sites)
width           <- block.size
jump 		<- 1


if(n.sites<block.size){
return(list(z=NaN, pval=NaN))
}


repeatlength <- ceiling( (n.sites-width+1)/jump )
  
D_sim 	     <- rep(NaN,repeatlength)      

for(zz in 1:repeatlength){
 
        
        start      <- ((zz-1) * jump + 1)
        end 	   <- ((zz-1) * jump + width) 
        window	   <- start:end 

	ABBA       <- ABBA_sites[-window]
	BABA	   <- BABA_sites[-window]
	D_sim[zz]  <- ( sum(ABBA, na.rm=TRUE) - sum(BABA, na.rm=TRUE) ) / ( sum(ABBA, na.rm=TRUE) + sum(BABA, na.rm=TRUE) )

}

#after Martin (two tailed z-test)
D_sd  <- sd(D_sim,na.rm=TRUE)
D_err <- D_sd/sqrt(repeatlength)
D_Z   <- D.base / D_err
D_p   <- 2*pnorm(-abs(D_Z))

#after Eaton and Ree 2013 
#D_Z    <- abs(D.base/sd(D_sim,na.rm=TRUE))
#D_p    <- 2 * (1 - pnorm(D_Z)) 


return(list(z=D_Z, pval=D_p))

}
