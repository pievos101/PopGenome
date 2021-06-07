calc_pairwise_Tv <- function(bial, populations, outgroup, transitions, shared.only){

indis        <- unique(unlist(populations))
n.ind        <- length(indis)
if(n.ind==0){

	n.ind <- dim(bial)[1]

}else{

	bial <- bial[indis,] 

}	

n.bial.sites <- dim(bial)[2]

TVmatrix <- matrix(0, n.ind, n.ind)
rownames(TVmatrix) <- rownames(bial)
colnames(TVmatrix) <- rownames(bial)

PAIRS <- combn(n.ind,2)

for(xx in 1:dim(PAIRS)[2]){

	ind1_id <- PAIRS[1,xx]
	ind2_id <- PAIRS[2,xx]

	bial_ind1  <- bial[ind1_id,]
	bial_ind2  <- bial[ind2_id,]

	diff = which(bial_ind1 != bial_ind2)
	n.tv = sum(!transitions[diff])
	TVmatrix[ind1_id, ind2_id] <- n.tv/n.bial.sites
	TVmatrix[ind2_id, ind1_id] <- n.tv/n.bial.sites

}

return(TVmatrix)

}# End of function