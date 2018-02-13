setGeneric("PG_plot.biallelic.matrix", function(object, region, ind.names = FALSE , cex.axis = 0.5, title="") standardGeneric("PG_plot.biallelic.matrix"))
 setMethod("PG_plot.biallelic.matrix", "GENOME",
 function(object, region, ind.names, cex.axis,title){


bial   <- get.biallelic.matrix(object,region)
names  <- rownames(bial)

if(ind.names[1]!=FALSE){
 ids   <- match(ind.names, names)
 names <- names[ids]
 bial  <- bial[ids,]
}


bial.sites <- object@region.data@biallelic.sites[[region]]

# Plot first individual
CO   <- colors()
ind1 <- bial[1,]
ind1[ind1==0] <- NA
plot(bial.sites,ind1,ylim=c(0,length(names)),yaxt="n",ylab="", xlab="polymorphic sites",cex=0.2,pch=19, main=title) 
axis(2,1:length(names), labels=names, lwd=0.5, cex.axis=cex.axis, las=1)
abline(h=1, col="grey")

for (xx in 2:length(names) ){

	ind  <- bial[xx,]
        ind2 <- ind
	ind2[ind==0] <- xx
	ind2[ind!=0] <- NA
	points(bial.sites,ind2,cex=0.2,col="gray",pch=19) #col=CO[xx])
	ind2 <- ind
	ind2[ind==1] <- xx
	ind2[ind!=1] <- NA
	points(bial.sites,ind2,cex=0.2,col="black",pch=19)
	abline(h=xx, col="grey")

}

})
