setGeneric("weighted.jackknife", function(object, do.D=TRUE, do.df=TRUE, per.region=FALSE, block.size=1) standardGeneric("weighted.jackknife"))
 setMethod("weighted.jackknife", "GENOME",

 function(object, do.D, do.df, per.region, block.size){
  

region.names                 <- object@region.names
n.region.names               <- length(region.names)

#n.region.names=10

df  <- object@df
D   <- object@D


if(per.region){ #drop-one jackknife 
  
  # init 
  D.z      <- matrix(NaN,n.region.names,1) # jacknife
  D.pval   <- matrix(NaN,n.region.names,1) # jacknife
  df.z    <- matrix(NaN,n.region.names,1) # jacknife
  df.pval <- matrix(NaN,n.region.names,1) # jacknife
  df.SE   <- matrix(NaN,n.region.names,1) # jacknife

  # get site specific values 
  site.D    <- object@region.stats@D
  site.df  <- object@region.stats@df

## PROGRESS #########################
 progr <- progressBar()
#####################################  

  for(xx in 1:n.region.names){

	if(do.D){
	res        <- D_jacknife(site.D[[xx]], D[xx], block.size=block.size)
	D.z[xx]	   <- res$z
	D.pval[xx] <- res$pval	
	}
	if(do.df){	
	res          <- D_jacknife(site.df[[xx]], df[xx], block.size=block.size)	
	df.z[xx]    <- res$z
	df.pval[xx] <- res$pval	
	}

 ## PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
 ###################################################################

  }

object@df.z    <- df.z
object@df.pval <- df.pval
object@D.z      <- D.z
object@D.pval   <- D.pval

return(object)

}

if(!per.region){ # weighted jackknife 

# Calculate the global D and df 
D_ABBA <- sum( sapply(object@region.stats@D, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
D_BABA <- sum( sapply(object@region.stats@D, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

df_ABBA <- sum( sapply(object@region.stats@df, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
df_BABA <- sum( sapply(object@region.stats@df, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

D   <- (D_ABBA-D_BABA)/(D_ABBA+D_BABA)
df <- (df_ABBA-df_BABA)/(df_ABBA+df_BABA)

#print(D)
#print(df)

#get informative sites per region
df.valid <- sapply(object@region.stats@df, function(x){if(length(x)==0){return(NaN)}else{return(sum(!is.na(x[1,])))}})
D.valid   <- sapply(object@region.stats@D, function(x){if(length(x)==0){return(NaN)}else{return(sum(!is.na(x[1,])))}})

df_M     <- sum(df.valid, na.rm=TRUE)
D_M       <- sum(D.valid, na.rm=TRUE)

#print(df_M)
#print(D_M)

df.valid <- df.valid/df_M
D.valid   <- D.valid/D_M


reset_jackdf   <- object@region.stats@df
reset_jackD     <- object@region.stats@D


jackdf   <- object@region.stats@df
jackD     <- object@region.stats@D

Di      <- rep(NaN,n.region.names)
dfi    <- rep(NaN,n.region.names)

## PROGRESS #########################
 progr <- progressBar()
#####################################

 for (xx in 1:n.region.names){
	
        # jack one region 
	jackdf[[xx]] <- NULL
        jackD[[xx]]   <- NULL

	#recalculate df and D 
	D_ABBA <- sum( sapply(jackD, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
	D_BABA <- sum( sapply(jackD, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

	df_ABBA <- sum( sapply(jackdf, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[2,], na.rm=TRUE))}}), na.rm=TRUE)
	df_BABA <- sum( sapply(jackdf, function(x){if(length(x)==0){return(NaN)}else{return(sum(x[3,], na.rm=TRUE))}}), na.rm=TRUE)

	Di[xx]   <- (D_ABBA-D_BABA)/(D_ABBA+D_BABA)
	dfi[xx] <- (df_ABBA-df_BABA)/(df_ABBA+df_BABA)

	#print(Di[xx])
	#print(dfi[xx])


	#reset
	jackdf[[xx]] <- reset_jackdf[[xx]]
	jackD[[xx]]   <- reset_jackD[[xx]]

## PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
###################################################################

 }

Du   <- sum(Di, na.rm=TRUE)/n.region.names
dfu <- sum(dfi, na.rm=TRUE)/n.region.names


#print("------")

#print(Du)
#print(dfu)

SE_D   <- sqrt( n.region.names*sum(D.valid*(Di-Du)^2, na.rm=TRUE) )
SE_df <- sqrt( n.region.names*sum(df.valid*(dfi-dfu)^2, na.rm=TRUE) )

#print(SE_D)
#print(SE_df)

object@df.z    <- as.matrix(df/SE_df)
object@df.pval <- as.matrix(2*pnorm(-abs(object@df.z)))
object@D.z      <- as.matrix(D/SE_D)
object@D.pval   <- as.matrix(2*pnorm(-abs(object@D.z)))
object@df.SE   <- as.matrix(SE_df)
object@D.SE     <- as.matrix(SE_D)

return(object)

}



})
 
 
