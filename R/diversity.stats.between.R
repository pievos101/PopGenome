setGeneric("diversity.stats.between", function(object,new.populations=FALSE,subsites=FALSE, keep.site.info=FALSE, haplotype.mode=FALSE, nucleotide.mode=TRUE) standardGeneric("diversity.stats.between"))

setMethod("diversity.stats.between","GENOME",function(object,new.populations,subsites,keep.site.info,haplotype.mode,nucleotide.mode){
  
  region.names                 <- object@region.names
  n.region.names               <- length(region.names)
  if(object@big.data){region.names <- NULL} # because of memory space
 
  #object@Pop_FSTH$sites        <- "ALL"
  #object@Pop_FSTH$calculated   <- TRUE
   
  if(!missing(new.populations)){
    NEWPOP <- TRUE
    populations <- vector("list",length(new.populations))
    npops       <- length(populations)            # Wenn mehr Pops definiert werden
    object@Pop_FSTH$Populations  <- new.populations
  }else{
    NEWPOP <- FALSE
    npops                       <- length(object@populations)     # alte Anzahl der Populationen
    object@Pop_FSTH$Populations <- object@populations
  }
 
 #########################################
 # INIT
 #########################################
  
   # Get the names ----------for pairwaise comparison --------------------------------------------------
  #############################################################################
  
if(npops>1){

  poppairs <- choose(npops,2)   # Outgroup is not included !!
  pairs    <- combn(1:(npops),2) 
 
#### --- Names of population pairs --- #### 
 nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
 if(dim(pairs)[2]>1){ # more than 2 Populations
  for(xx in 2:dim(pairs)[2]){
    m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
  } 
 }#END if
}# End npops > 1
else{poppairs <- 1;nn <- "pop1"} 
##### ------------------------------ ####------------------------------------------------ 
#########################################################################################  
  

  #nam    <- paste("pop",1:npops)
  init1  	 <- matrix(0,n.region.names, poppairs)
  
  nucb   	 <- init1
  hapb		 <- init1
  site.nucb      <- vector("list",n.region.names) # region stats
 
  change    <- object@region.stats
  #Pop_FSTH  <- vector("list",n.region.names)
  
#--------------------------------------------------

  # Names ----------------------------------------
  rownames(nucb)     <- region.names
  colnames(nucb)     <- nn
  rownames(hapb)     <- region.names
  colnames(hapb)     <- nn
  # ----------------------------------------------


## PROGRESS #########################
 progr <- progressBar()
#####################################



   
for(xx in 1:n.region.names){

### if Subsites ----------------------------------

bial <- popGetBial(object,xx)
if(length(bial)==0){next} 

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
     
     if(NEWPOP){ # Wenn eine andere Population definiert !
          
       for(yy in 1:npops){       
           if(is.character(new.populations[[yy]])){            
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]             
           }else{ # numeric values
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}
           }                      
       }     
       
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   

       if(length(populations)==0){next} # Keine Population vorhanden
       
        
    }else{populations <- object@region.data@populations[[xx]]} # Wenn keine neue Pop definiert !
    

	if(NEWPOP) {temp       <- checkpoppairs(npops,popmissing,pairs,nn)} # welche populationen wurden ueberhaupt berechnet
        if(!NEWPOP){temp       <- checkpoppairs(npops,object@region.data@popmissing[[xx]],pairs,nn)} 
   
        respop     <- temp$respop
        respairpop <- temp$respairpop    

   
      res                    <- site_diversity_between(bial,populations,keep.site.info,haplotype.mode,nucleotide.mode)

    
    # fill detailed Slots --------------------------------#
     if(nucleotide.mode & keep.site.info){	   	
     	site.nucb[[xx]]  		<- res$site_div  
     	if(npops>2){
     	rownames(site.nucb[[xx]])	<- nn
     	}	 	
     }
    # ----------------------------------------------------# 
     
     if(nucleotide.mode){
     nucb[xx,respairpop]   <- res$window_div
     }

     if(haplotype.mode){	
     hapb[xx,respairpop]   <- res$hapavek	   
     }

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}

 
 change@nuc.diversity.between      <- site.nucb
 #change@Pop_FSTH                  <- Pop_FSTH
 object@region.stats               <- change
 rm(change)
 gc()
 object@nuc.diversity.between  <- nucb
 object@hap.diversity.between  <- hapb
 
 return(object)
 })


