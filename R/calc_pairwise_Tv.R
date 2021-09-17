calc_pairwise_Tv <- function(bial, populations, outgroup, transitions, 
								shared.only, by.genotype){

if(!by.genotype){

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

}

### BY GENOTYPE

if(by.genotype){


indis        <- unique(unlist(populations))
n.ind        <- length(indis)

if(n.ind==0){

	n.ind <- dim(bial)[1]

}else{

	bial <- bial[indis,] 

}	

n.bial.sites <- dim(bial)[2]

TVmatrix <- matrix(0, n.ind/2, n.ind/2)
iii      <- seq(1,dim(bial)[1],2)
rownames(TVmatrix) <- rownames(bial)[iii]
colnames(TVmatrix) <- rownames(bial)[iii]


idx  <- sort(rep(1:(dim(bial)[1]/2),2))

geno <- split(rownames(bial), idx)


PAIRS <- combn(length(geno),2)

for(xx in 1:dim(PAIRS)[2]){

	ind1_id <- PAIRS[1,xx]
	ind2_id <- PAIRS[2,xx]

	geno1   <- geno[[ind1_id]] 
	geno2   <- geno[[ind2_id]]

	geno1bial <- bial[geno1,]
	geno2bial <- bial[geno2,]
	
	check1 <- geno1bial[1,]+geno1bial[2,]
	check2 <- geno2bial[1,]+geno2bial[2,]

	SCA  <- rep(0.5, dim(bial)[2])
	ONE  <- which((check1==0 & check2==2) || (check1==2 & check2==0))
	ZERO <- which((check1==0 & check2==0) || (check1==2 & check2==2))
	SCA[ONE]  <- 1
	SCA[ZERO] <- 0

	transv <- as.numeric(!transitions)
	n.tv   <- sum(SCA*transv)
	n.tv   <- sum(SCA)


	TVmatrix[ind1_id, ind2_id] <- n.tv/dim(bial)[2]
	TVmatrix[ind2_id, ind1_id] <- n.tv/dim(bial)[2]

}

return(TVmatrix)

}

}# End of function