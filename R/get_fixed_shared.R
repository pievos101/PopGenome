get_fixed_shared <- function(MAF, fixed.threshold, fixed.threshold.fst, bial, populations){


  npops       <- dim(MAF)[1]
  pop.pairs   <- combn(npops,2)
  #n.pop.pairs <- choose(npops,2)  
	
	res <- apply(pop.pairs,2,function(y){


	      pop1vec  <- MAF[y[1],]
              pop2vec  <- MAF[y[2],]
	      #n.snps   <- sum(!is.na(pop1vec) & !is.na(pop2vec))

              if(fixed.threshold.fst){
	      fst        <- site_FST(bial,list(populations[[y[1]]],populations[[y[2]]]))
              fixedTF    <- fst>=fixed.threshold.fst
	      nonfixedTF <- !fixedTF	
	      }else{
              fixedTF    <- abs(pop1vec - pop2vec)>=fixed.threshold 
              nonfixedTF <- !fixedTF
	      }	
	
              fixed    <- sum(fixedTF, na.rm=TRUE)
	      mono     <- (pop1vec==0) & (pop2vec==0)		
	      mono     <- sum(mono,na.rm=TRUE)
	      shared   <- (pop1vec>0) & (pop2vec>0) & nonfixedTF
              shared   <- sum(shared,na.rm=TRUE)
	      #cat("pop",y[1],"vs","pop",y[2],"\n")
	      #cat("fixed:",fixed,"\n")
	      #cat("shared:",shared,"\n")
	      #cat("mono:",mono,"\n")
	      return(c(fixed,shared,mono))	  	
	})

return(res)

}
