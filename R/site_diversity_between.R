site_diversity_between <- function(matrix, populations, keep.site.info=TRUE, haplotype.mode=FALSE, nucleotide.mode=TRUE){
  
window_div <- NULL
site_div   <- NULL
hapavek    <- NULL

# nucleotide diversity ----------------------------------------------------------------------------- 
if(nucleotide.mode){

  n.individuals   <- dim(matrix)[1]
  if(n.individuals==1){
	  return(list(window_div=NaN,site_div=NaN,hapavek=NaN))
	  #return(0)
  }

  if(length(populations)<2){
	return(list(window_div=NaN,site_div=NaN,hapavek=NaN))
  	#return(NaN)
  }
  
  npops       <- length(populations)
  pop.pairs   <- combn(npops,2)
  n.pop.pairs <- choose(npops,2)  

  erg <- apply(matrix,2,function(x){
	
	res <- apply(pop.pairs,2,function(y){
	      pop1vec  <- x[populations[[y[1]]]]
              pop2vec  <- x[populations[[y[2]]]]
              pop1zero <- sum(pop1vec==0, na.rm=TRUE)
              pop2zero <- sum(pop2vec==0, na.rm=TRUE)
              pop1one  <- sum(pop1vec==1, na.rm=TRUE)
              pop2one  <- sum(pop2vec==1, na.rm=TRUE)
	      pop1size <- pop1zero + pop1one
              pop2size <- pop2zero + pop2one
	      vergl    <- pop1size * pop2size	
              divers   <- (pop1zero*pop2one + pop1one*pop2zero)/vergl
	      return(divers)	  	
	})
 #return(mean(res))
  }) 
 
 if(npops>2){
 	window_div <- rowSums(erg, na.rm=TRUE)
 }else{
	window_div <- sum(erg, na.rm=TRUE)
 }
 
 #print(window_div)

 if(keep.site.info){
 	site_div   <- erg
 }else{
	erg 	 <- NULL
	site_div <- NULL
 }	

 #return(list(window_div=window_div,site_div=site_div))
} # END OF IF NUCLEOTIDE MODE

# haplotype diversity --------------------------------------------------------------------------

if(haplotype.mode){

npops       <- length(populations)

matrix_hap <- matrix
rownames(matrix_hap) <- rownames(matrix)
matrix_hap_sub       <- matrix_hap[unique(unlist(populations)),,drop=FALSE]#### 

duplids               <- .Call("my_unique_C", matrix_hap_sub)
uniquematrix          <- matrix_hap_sub[!duplids,,drop=FALSE]

nhgesamt             <- dim(uniquematrix)[1]
sfreqh               <- matrix(0,npops,nhgesamt)
 
nam       	     <- paste("pop", 1:npops)
rownames(sfreqh)     <- nam
colnames(sfreqh)     <- rownames(uniquematrix)
rownames(matrix_hap) <- NULL



 for(xx in 1:npops){
  
  sfreqh[xx,]  <- .Call("C_get_sfreqh_C",uniquematrix,matrix_hap[populations[[xx]],,drop=FALSE]) 
  #nh[xx]       <- length(which(sfreqh[xx,]!=0)) # number of haplotypes for each population
 
 }


 pairs <- combn(npops,2)
 #--------------

 # Apply
 hapavek     <- apply(pairs,2,function(x){
     
      hapa <- NaN
     
      m1_size <- length(populations[[ x[1] ]]) # size of population 1
      m2_size <- length(populations[[ x[2] ]]) # size of population 2     
     
     freqpop1 <- sfreqh[x[1],,drop=FALSE]
     freqpop2 <- sfreqh[x[2],,drop=FALSE]
     
     hapa  <- .Call("combnsum2_C",freqpop1,freqpop2)
         
     hapa  <- hapa/(m1_size*m2_size)
     
     return(hapa)

}) 

#hapavek           <- as.matrix(hapavek)
} # END OF HAPLOTYPE MODE

return(list(window_div=window_div,site_div=site_div,hapavek=hapavek))
}# End of function
