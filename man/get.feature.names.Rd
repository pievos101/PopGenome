\name{get.feature.names}

\docType{methods}
\alias{get.feature.names}

\title{Feature informations and GFF-attributes}


\description{
 Returns feature names and additional attributes for a given region
}

\usage{
get.feature.names(object, gff.file, chr)
}

\arguments{ 
  \item{object}{An object of class GENOME}
  \item{gff.file}{The corresponding GFF file}
  \item{chr}{The chromosome/scaffold identifier}
  
}

\details{
The algorithm uses the information stored in \code{GENOME.class.split@region.names} to iterate \cr
over the GFF file and returns attribute plus feature informations for each given region. \cr
Note, the functions splitting.data, split_data_into_GFF_attributes or sliding.window.transform \cr
should be performed prior to that.\cr\cr

The slot \code{region.names} must have the following form: "pos1 - pos2". 

}

\value{
The returned value is a character vector of length \code{length(GENOME.class.split@region.names)}
}


\examples{
# GENOME.class <- readVCF("chr1.vcf.gz",1000,"1",1,100000)
# GENOME.class.split <- split_data_into_GFF_attributes(GENOME.class,"Homo_sapiens.GRCh37.73.gtf", 
# "1", "gene_name")
# GENOME.class.split@region.names
# info <- get.feature.names(GENOME.class.split, gff.file="Homo_sapiens.GRCh37.73.gtf", chr="1")
# GENOME.class.split <- splitting.data(GENOME.class, subsites="gene")
# GENOME.class.split@region.names
# info <- get.feature.names(GENOME.class.split, gff.file="Homo_sapiens.GRCh37.73.gtf", chr="1")
}

