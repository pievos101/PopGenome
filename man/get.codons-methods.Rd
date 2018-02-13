\name{get.codons-methods}
\docType{methods}
\alias{get.codons,GENOME-method}
\alias{get.codons-methods}


\title{Detailed information about the nature of codon changes}


\description{
This generic function returns some information about the codon changes resulting 
from the observed SNPs.
}


\usage{
\S4method{get.codons}{GENOME}(object, regionID)
}

\arguments{
  \item{object}{an object of class \code{"GENOME"}}
  \item{regionID}{what region/alignment should be analyzed ?}
}
\details{
The slot GENOME.class@region.data@synonymous and GENOME.class@region.data@codons \cr
have to be set.\cr The data have to be read in with the correponding GFF file. \cr
The function set.synnonsyn(..., save.codons=TRUE) sets the syn/nonsny \cr 
sites in case of SNP data and stores the corresponding codon changes. 
}

\value{
The function \code{get.codons} returns a data.frame with the following information \cr

\item{1}{Position of the SNPs} 
\item{2}{Major Codon} 
\item{3}{Minor Codon}   
\item{4}{Major amino acid}    
\item{5}{Minor amino acid} 
\item{6}{synonymous (TRUE/FALSE)} 
\item{7}{Polarity of the major amino acid} 
\item{8}{Polarity of the minor amino acid}   

}

\examples{

# Alignments
# GENOME.class <- readData("FASTA",gffpath="GFF")
# get.codons(GENOME.class,1)
# SNP data
# GENOME.class <- readData("VCF",gffpath="GFF")
# GENOME.class <- set.synnonsyn(GENOME.class, ref.chr="ref.fas",save.codons=TRUE)
# get.codons(GENOME.class,1)

}



\keyword{methods}
