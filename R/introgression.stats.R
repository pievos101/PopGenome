setGeneric("introgression.stats", function(object, subsites=FALSE, do.D=TRUE, do.df=TRUE, keep.site.info=TRUE, block.size=FALSE, do.RNDmin=FALSE, l.smooth=TRUE) standardGeneric("introgression.stats"))

setMethod("introgression.stats","GENOME",function(object, subsites, do.D, do.df, keep.site.info, block.size, do.RNDmin, l.smooth){

do.BD=FALSE
dxy.table=FALSE
D.global=FALSE
do.CLR=FALSE
dgt=FALSE
lambda=FALSE

  
if(do.CLR){
keep.site.info <- TRUE
do.BD          <- TRUE
if(length(D.global)<=1){
stop("To calculate the Bd-clr >>D.global<< has to be specified !")
}
}

if(do.D){ 
do.BD     <- FALSE
#do.df   <- FALSE
}

if(do.BD){ 
do.D      <- FALSE
do.df    <- FALSE
#do.RNDmin <- FALSE
}

if(do.df){ 
do.BD    <- FALSE
#do.D     <- FALSE
}


  region.names                 <- object@region.names
  n.region.names               <- length(region.names)
  if(object@big.data){region.names <- NULL} # because of memory space
 
    NEWPOP <- FALSE
    npops                       <- length(object@populations)     # total number of populations
     

# If Jacknife should be done we need keep.site.info to be set
if(block.size[1]!=FALSE){
keep.site.info=TRUE
}

   
# init
  D      <- matrix(NaN,n.region.names,1)
  BD     <- matrix(NaN,n.region.names,1)
  df    <- matrix(NaN,n.region.names,1)
  D3    <- matrix(NaN,n.region.names,1)
  df_bayes  <- matrix(NaN,n.region.names,1)
  alpha_ABBA  <- matrix(NaN,n.region.names,1)
  alpha_BABA  <- matrix(NaN,n.region.names,1)
  beta_BBAA  <- matrix(NaN,n.region.names,1)


  Bd_dir <- matrix(NaN,n.region.names,1)
  Bd_clr <- matrix(NaN,n.region.names,1)
P.Bd_clr <- matrix(NaN,n.region.names,1)
  f      <- matrix(NaN,n.region.names,1)
  RNDmin <- matrix(NaN,n.region.names,1)
  Gmin   <- matrix(NaN,n.region.names,1)
  D.z    <- matrix(NaN,n.region.names,1) # jacknife
  D.pval <- matrix(NaN,n.region.names,1) # jacknife
  df.z    <- matrix(NaN,n.region.names,1) # jacknife
  df.pval <- matrix(NaN,n.region.names,1) # jacknife
#--------------------------------------------------

  # Names ----------------------------------------
  rownames(D)       <- region.names
  colnames(D)       <- "D statistic"
  rownames(BD)      <- region.names
  colnames(BD)      <- "BD statistic"
  rownames(df)     <- region.names
  colnames(df)     <- "df statistic"
  rownames(D3)     <- region.names
  colnames(D3)     <- "D3 statistic"
  rownames(df_bayes)     <- region.names
  colnames(df_bayes)     <- "df bayes factor"
  rownames(Bd_dir)  <- region.names
  colnames(Bd_dir)  <- "Direction"
  rownames(Bd_clr)  <- region.names
  colnames(Bd_clr)  <- "Bd-clr"
  rownames(P.Bd_clr)  <- region.names
  colnames(P.Bd_clr)  <- "P.Bd_clr"
  rownames(f)       <- region.names
  colnames(f)       <- "f statistic"
  rownames(RNDmin)  <- region.names
  colnames(RNDmin)  <- "RNDmin"
  rownames(Gmin)  <- region.names
  colnames(Gmin)  <- "RNDmin"
  # ----------------------------------------------

  site.D         <- vector("list",n.region.names) # region stats
  site.df       <- vector("list",n.region.names)

  change    	 <- object@region.stats

 # Bd-CLR
 if(do.CLR==TRUE){
 print("Calculate Bd-Transition Table")
 D.global     <- D.global[!is.na(D.global)]
 #D.global[is.na(D.global)] <- 0	
 D.global     <- round(D.global, digits=dgt)
 Trans.matrix <- calc_Bd_CLR_table(D.global, lambda)
 }


## PROGRESS #########################
 progr <- progressBar()
#####################################
   
for(xx in 1:n.region.names){

### if Subsites ----------------------------------

bial <- popGetBial(object,xx)

if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
  # object@Pop_FSTH$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_FSTH$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
  # object@Pop_FSTH$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
  # object@Pop_FSTH$sites <- "nonsynonymous"
}

if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
  # object@Pop_Linkage$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
  # object@Pop_FSTH$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
  # object@Pop_FSTH$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]])==TRUE)
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
  # object@Pop_FSTH$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]]==TRUE)
   bial             <- bial[,gene,drop=FALSE]
  # object@Pop_FSTH$sites <- "gene"
}

if(subsites=="intergenic"){
  intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   if(length(intron)==0){
     intron <- !object@region.data@ExonSNPS[[xx]]	  
   }

  utr            <- object@region.data@UTRSNPS[[xx]]
  exon           <- object@region.data@ExonSNPS[[xx]]
  gene           <- object@region.data@GeneSNPS[[xx]]
  coding         <- !is.na(object@region.data@synonymous[[xx]])  

  inter          <- !(intron|utr|exon|gene|coding)
  bial           <- bial[,inter,drop=FALSE]
  #object@Pop_FSTH$sites <- "intergenic"
}

if(subsites=="included" & length(bial!=0)){
   included         <- which(object@region.data@included[[xx]]==TRUE)
   bial             <- bial[,included,drop=FALSE]
}

}# End if subsites
############### ---------------------------------


 if(length(bial)!=0){ # if a biallelic position exists
           
    populations <- object@region.data@populations[[xx]] 
    outgroup    <- object@region.data@outgroup[[xx]]  

    if(do.D){
	
	    res       <- calc_D(bial, populations, outgroup, keep.site.info) 
    	    D[xx]     <- res$D
    	    f[xx]     <- res$f	 
		
   # fill detailed Slots --------------------------------# 
     if(keep.site.info){	 	
     site.D[[xx]]  		<- rbind(res$D_site,res$ABBA,res$BABA)
     rownames(site.D[[xx]])     <- c("D","ABBA","BABA")
     }		  
   # ----------------------------------------------------# 
    }

    if(do.BD){
	    res        <- calc_BD(bial, populations, outgroup, keep.site.info, dxy.table) 
    	    BD[xx]     <- res$D
    	    #f[xx]      <- res$f
            Bd_dir[xx] <- res$Bd_dir
		
   # fill detailed Slots --------------------------------# 
     if(keep.site.info){	 	
     site.D[[xx]]  		<- rbind(res$D_site,res$ABBA,res$BABA)
     rownames(site.D[[xx]])     <- c("D","ABBA","BABA")
     }		  
   # ----------------------------------------------------# 
    }

    if(do.df){
	    res            <- calc_df(bial, populations, outgroup, keep.site.info, dxy.table, l.smooth) 
    	    df[xx]        <- res$D
            D3[xx]        <- res$D3
            df_bayes[xx]  <- res$D_bayes
            alpha_ABBA[xx] <- res$alpha_ABBA
	    alpha_BABA[xx] <- res$alpha_BABA
	    beta_BBAA[xx]  <- res$beta_BBAA
		
    	    #f[xx]      <- res$f
            Bd_dir[xx] <- res$Bd_dir
		
   # fill detailed Slots --------------------------------# 
     if(keep.site.info){	 	
     site.df[[xx]]  		<- rbind(res$D_site,res$ABBA,res$BABA)
     rownames(site.df[[xx]])   <- c("df","ABBA","BABA")
     }		  
   # ----------------------------------------------------# 
    }


   if(do.RNDmin){
        res        <- calc_RNDmin(bial, populations, outgroup)
  	RNDmin[xx] <- res[[1]]
	Gmin[xx]   <- res[[2]]	
   }		


  # Perform Bd_clr
  if(do.CLR==TRUE){
	 D.local      <- site.D[[xx]][1,]
	 D.local      <- round(site.D[[xx]][1,], digits=dgt)
	 rett         <- calc_Bd_clr(D.local, Trans.matrix, lambda)
  	 Bd_clr[xx]   <- rett[[1]]       
         P.Bd_clr[xx] <- rett[[2]]	
  }		
   	
  # Perform jacknife if requested
  if(block.size[1]!=FALSE){
        
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
 }	



  ## PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}

 change@D      			<- site.D
 change@df			<- site.df	
 object@region.stats            <- change
 rm(change)
 gc()
 object@D        <- D
 object@BD       <- BD
 object@df       <- df
 object@D3       <- D3
 object@df_bayes      <- df_bayes
 object@alpha_ABBA     <- alpha_ABBA
 object@alpha_BABA     <- alpha_BABA
 object@beta_BBAA      <- beta_BBAA
 
 object@Bd_clr   <- Bd_clr
 object@Bd_dir   <- Bd_dir
 object@P.Bd_clr <- P.Bd_clr
 object@f        <- f
 object@D.z      <- D.z
 object@D.pval   <- D.pval
 object@df.z      <- df.z
 object@df.pval   <- df.pval

 object@RNDmin   <- RNDmin
 object@Gmin     <- Gmin

 return(object)
 })


