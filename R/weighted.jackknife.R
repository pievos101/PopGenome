setGeneric("weighted.jackknife", function(object, do.D=TRUE, do.BDF=TRUE, per.region=FALSE, block.size=1) standardGeneric("weighted.jackknife"))
 setMethod("weighted.jackknife", "GENOME",

 function(object, do.D, do.BDF, per.region, block.size){
  

region.names                 <- object@region.names
n.region.names               <- length(region.names)

#n.region.names=10

BDF <- object@BDF
D   <- object@D


if(per.region){ #drop-one jackknife 
  
  # init 
  D.z      <- matrix(NaN,n.region.names,1) # jacknife
  D.pval   <- matrix(NaN,n.region.names,1) # jacknife
  BDF.z    <- matrix(NaN,n.region.names,1) # jacknife
  BDF.pval <- matrix(NaN,n.region.names,1) # jacknife
  BDF.SE   <- matrix(NaN,n.region.names,1) # jacknife

  # get site specific values 
  site.D    <- object@region.stats@D
  site.BDF  <- object@region.stats@BDF

## PROGRESS #########################
 progr <- progressBar()
#####################################  

  for(xx in 1:n.region.names){

	if(do.D){
	res        <- D_jacknife(site.D[[xx]], D[xx], block.size=block.size)
	D.z[xx]	   <- res$z
	D.pval[xx] <- res$pval	
	}
	if(do.BDF){	
	res          <- D_jacknife(site.BDF[[xx]], BDF[xx], block.size=block.size)	
	BDF.z[xx]    <- res$z
	BDF.pval[xx] <- res$pval	
	}

 ## PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
 ###################################################################

  }

object@BDF.z    <- BDF.z
object@BDF.pval <- BDF.pval
object@D.z      <- D.z
object@D.pval   <- D.pval

return(object)

}

if(!per.region){ # weighted jackknife 

# Calculate the global D and BDF 
D_ABBA <- sum( sapply(object@region.stats@D, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
D_BABA <- sum( sapply(object@region.stats@D, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

BDF_ABBA <- sum( sapply(object@region.stats@BDF, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
BDF_BABA <- sum( sapply(object@region.stats@BDF, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

D   <- (D_ABBA-D_BABA)/(D_ABBA+D_BABA)
BDF <- (BDF_ABBA-BDF_BABA)/(BDF_ABBA+BDF_BABA)

#print(D)
#print(BDF)

#get informative sites per region
BDF.valid <- sapply(object@region.stats@BDF, function(x){if(length(x)==0){return(NaN)}else{return(sum(!is.na(x[1,])))}})
D.valid   <- sapply(object@region.stats@D, function(x){if(length(x)==0){return(NaN)}else{return(sum(!is.na(x[1,])))}})

BDF_M     <- sum(BDF.valid, na.rm=TRUE)
D_M       <- sum(D.valid, na.rm=TRUE)

#print(BDF_M)
#print(D_M)

BDF.valid <- BDF.valid/BDF_M
D.valid   <- D.valid/D_M


reset_jackBDF   <- object@region.stats@BDF
reset_jackD     <- object@region.stats@D


jackBDF   <- object@region.stats@BDF
jackD     <- object@region.stats@D

Di      <- rep(NaN,n.region.names)
BDFi    <- rep(NaN,n.region.names)

## PROGRESS #########################
 progr <- progressBar()
#####################################

 for (xx in 1:n.region.names){
	
        # jack one region 
	jackBDF[[xx]] <- NULL
        jackD[[xx]]   <- NULL

	#recalculate BDF and D 
	D_ABBA <- sum( sapply(jackD, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
	D_BABA <- sum( sapply(jackD, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

	BDF_ABBA <- sum( sapply(jackBDF, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
	BDF_BABA <- sum( sapply(jackBDF, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

	Di[xx]   <- (D_ABBA-D_BABA)/(D_ABBA+D_BABA)
	BDFi[xx] <- (BDF_ABBA-BDF_BABA)/(BDF_ABBA+BDF_BABA)

	#print(Di[xx])
	#print(BDFi[xx])


	#reset
	jackBDF[[xx]] <- reset_jackBDF[[xx]]
	jackD[[xx]]   <- reset_jackD[[xx]]

## PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
###################################################################

 }

Du   <- sum(Di, na.rm=TRUE)/n.region.names
BDFu <- sum(BDFi, na.rm=TRUE)/n.region.names


#print("------")

#print(Du)
#print(BDFu)

SE_D   <- sqrt( n.region.names*sum(D.valid*(Di-Du)^2, na.rm=TRUE) )
SE_BDF <- sqrt( n.region.names*sum(BDF.valid*(BDFi-BDFu)^2, na.rm=TRUE) )

#print(SE_D)
#print(SE_BDF)

object@BDF.z    <- as.matrix(BDF/SE_BDF)
object@BDF.pval <- as.matrix(2*pnorm(-abs(object@BDF.z)))
object@D.z      <- as.matrix(D/SE_D)
object@D.pval   <- as.matrix(2*pnorm(-abs(object@D.z)))
object@BDF.SE   <- as.matrix(SE_BDF)
object@D.SE     <- as.matrix(SE_D)

return(object)

}



})
 
 
