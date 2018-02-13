#object, gff.file, chr, attribut
split_data_into_GFF_features <- function(object, gff.file, chr, feature){

chr <- as.character(chr)

 region      <- .Call("find_lines_GFF_Human2",gff.file,chr)
 start       <- region[1]
 end         <- region[2]
 gff.table   <- read.table(gff.file,sep="\t",colClasses=c(rep("NULL",2),"character",rep("integer",2),rep("NULL",4)),
                           skip = start - 1, nrows = end - start + 1)
 featids     <- which(gff.table[,1]==feature)
 if(length(featids)==0){stop("This feature was not found")}
 poslist     <- apply(gff.table[featids,2:3],1,function(x){return(x[1]:x[2])})

print("Split the data into Features")
object               <- splitting.data(object, positions=poslist, type=2)  
return(object)


}
