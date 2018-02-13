setGeneric("set.filter", function(object, missing.freqs=TRUE, minor.freqs=FALSE, maf.lower.bound=0, maf.upper.bound=1, miss.lower.bound=0, miss.upper.bound=1) standardGeneric("set.filter"))
 setMethod("set.filter", "GENOME",
 function(object, missing.freqs, minor.freqs, maf.lower.bound, maf.upper.bound, miss.lower.bound, miss.upper.bound){


if((missing.freqs==TRUE) & (minor.freqs==FALSE)){

	print("Using missing frequencies")
	object   <- count.unknowns(object)
	change   <- object@region.data

	MAF      <- object@region.stats@missing.freqs

	included <- vector("list",length(object@region.names))

	## PROGRESS #########################
	progr <- progressBar()
	#####################################

	# Iterate over the genomic regions
	for ( xx in 1:length(object@region.names) ){
		
		if(length(object@region.data@biallelic.sites[[xx]])==0){next}

		truefalse <- apply(MAF[[xx]],2,function(y){
			  check <- (y >= miss.lower.bound) & (y <= miss.upper.bound)
			  if(all(check)){return(TRUE)}else{return(FALSE)}	
		})

	included[[xx]] <- truefalse

	# PROGRESS #######################################################
        progr <- progressBar(xx,length(object@region.names), progr)
  	##################################################################

	}

change@included    <- included 
object@region.data <- change

return(object)
}

if((minor.freqs==TRUE) & (missing.freqs==FALSE) ){

	print("Using minor allele frequencies")
	object   <- detail.stats(object)
	change   <- object@region.data

	MAF      <- object@region.stats@minor.allele.freqs

	included <- vector("list",length(object@region.names))

	## PROGRESS #########################
	progr <- progressBar()
	#####################################

	# Iterate over the genomic regions
	for ( xx in 1:length(object@region.names) ){
		
		if(length(object@region.data@biallelic.sites[[xx]])==0){next}

		truefalse <- apply(MAF[[xx]],2,function(y){
			  check <- (y >= maf.lower.bound) & (y <= maf.upper.bound)
			  if(all(check)){return(TRUE)}else{return(FALSE)}	
		})

	included[[xx]] <- truefalse

	# PROGRESS #######################################################
        progr <- progressBar(xx,length(object@region.names), progr)
  	##################################################################

	}

change@included    <- included 
object@region.data <- change

return(object)
}

if((minor.freqs==TRUE) & (missing.freqs==TRUE) ){

	print("Using minor & missing allele frequencies")
	object   <- detail.stats(object)
	change   <- object@region.data

	MAF      <- object@region.stats@minor.allele.freqs

	included <- vector("list",length(object@region.names))

	## PROGRESS #########################
	progr <- progressBar()
	#####################################

	# Iterate over the genomic regions
	for ( xx in 1:length(object@region.names) ){
		
		if(length(object@region.data@biallelic.sites[[xx]])==0){next}

		truefalse <- apply(MAF[[xx]],2,function(y){
			  check <- (y >= maf.lower.bound) & (y <= maf.upper.bound) & (y >= miss.lower.bound) & (y <= miss.upper.bound)
			  if(all(check)){return(TRUE)}else{return(FALSE)}	
		})

	included[[xx]] <- truefalse

	# PROGRESS #######################################################
        progr <- progressBar(xx,length(object@region.names), progr)
  	##################################################################

	}

change@included    <- included 
object@region.data <- change

return(object)
}

print("Nothing filtered")

return(object)

})
