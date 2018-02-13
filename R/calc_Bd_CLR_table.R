calc_Bd_CLR_table <- function(D, lambda){


D       <- D[!is.na(D)]
#D       <- abs(D)

#D[is.na(D)] <- 0
#------------------------------
#newD <- rep(0,length(D)-1)
#for  ( xx in 1:(length(D)-1) ){
#newD[xx] <- abs(D[xx]-D[xx+1])
#}
#------------------------------
D.table <- table(D)
#return(D.table)

D.names <- as.numeric(names(D.table))

Trans   <- matrix(0,length(D.names),length(D.names))
rownames(Trans) <- D.names
colnames(Trans) <- D.names

for  ( xx in 1:(length(D)-1) ){

	#runs	
	#if(xx>1){
	#	if((D[xx]==D[xx-1]) && (D[xx]==D[xx+1])){next}
	#}

	#if(D[xx]!=D[xx+1]){next}

	row <- match(D[xx],D.names)
	col <- match(D[xx+1],D.names)
	#a <- 1 - abs(D[xx])
	#b <- 1 - abs(D[xx+1]) 

	Trans[row,col]  <- Trans[row,col]     +  1/exp(lambda*abs(D[xx]+D[xx+1])) #exp( (2-abs(D[xx]+D[xx+1]) ) )


	#Trans[col,row] <- Trans[col,row]     +  1 #+ abs(row-col) #2-(abs(D[xx]+D[xx+1])) #+ abs(D[xx]-D[xx+1])/2 
}


return(Trans)

}
