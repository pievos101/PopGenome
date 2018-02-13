\name{split_data_into_GFF_attributes}

\docType{methods}
\alias{split_data_into_GFF_attributes}

\title{Split the data into GFF attributes}


\description{
 Splits the data into GFF attributes defined by the user.
}

\usage{
split_data_into_GFF_attributes(object, gff.file, chr, attribute)
}

\arguments{ 
  \item{object}{An object of class GENOME}
  \item{gff.file}{The corresponding GFF file}
  \item{chr}{The chromosome/scaffold identifier}
  \item{attribute}{The attribute to use for splitting}
}

\details{
The algorithm splits the data into attributes. \cr 
An attribute can be \code{"gene_name"}, \code{"Parent"} or just a single \cr
gene name like \code{"geneXYZ"}. 
}

\value{
The returned value is an object of class \code{"GENOME"} \cr
See \code{GENOME.class.split@region.names} and \code{GENOME.class.split@region.names} \cr
after splitting the data.
}


\examples{
# GENOME.class <- readVCF("chr1.vcf.gz",1000,"1",1,100000)
# GENOME.class.split <- split_data_into_GFF_attributes(GENOME.class,"Homo_sapiens.GRCh37.73.gtf", 
# "1", "gene_name")
# GENOME.class.split@region.names
# GENOME.class.split@feature.names
}

