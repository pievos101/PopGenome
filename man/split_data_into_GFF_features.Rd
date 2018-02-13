\name{split_data_into_GFF_features}

\docType{methods}
\alias{split_data_into_GFF_features}

\title{Split the data into GFF features}


\description{
 Splits the data into GFF feautures defined by the user.
}

\usage{
split_data_into_GFF_features(object, gff.file, chr, feature)
}

\arguments{ 
  \item{object}{An object of class GENOME}
  \item{gff.file}{The corresponding GFF file}
  \item{chr}{The chromosome/scaffold identifier}
  \item{feature}{The feature used for splitting}
}

\details{
The algorithm splits the data into features. \cr 
A feature can be \code{"gene"}, \code{"exon"} etc. \cr 
depending on what is specified in the GFF3 file.
}

\value{
The returned value is an object of class \code{"GENOME"} \cr
See \code{GENOME.class.split@region.names} and \code{GENOME.class.split@region.names} \cr
after splitting the data.
}


\examples{
# GENOME.class <- readVCF("chr1.vcf.gz",1000,"1",1,100000)
# GENOME.class.split <- split_data_into_GFF_features(GENOME.class,"Homo_sapiens.GRCh37.73.gtf", 
# "1", "gene")
# GENOME.class.split@region.names
}

