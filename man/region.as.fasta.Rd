\name{region.as.fasta-methods}
\docType{methods}

\alias{region.as.fasta,GENOME-method}
\alias{region.as.fasta-methods}


\title{Extract a region and write it to a FASTA file}


\description{
This generic function writes a FASTA file of the observed biallelic positions to the current workspace.
}

\usage{

\S4method{region.as.fasta}{GENOME}(object,region.id=FALSE,filename=FALSE,type=1,ref.chr=FALSE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{region.id}{region of the genome}
  \item{filename}{name of the FASTA file}
  \item{type}{1: extract SNPs; 2: extract all nucleotides}
  \item{ref.chr}{reference sequence}
}


\details{
In case of \code{type=2} we recommend to use the function 
splitting.data(,positions=list( ... ), type=2) before and apply 
the region.as.fasta() to this splitted object afterwards. The \code{type=1}
method will write a FASTA file including only the biallelic.sites. \cr
\code{region.id} is the the region number specified in the PopGenome class \code{GENOME}.
}


\examples{

#GENOME.class       <- readSNP("Arabidopsis",CHR=1)
# split the data into the genmic positions 100 to 2000
#GENOME.class.split <- splitting.data(GENOME.class, positions=list(100:2000),type=2)
#GENOME.class.split@region.names
#region.as.fasta(GENOME.class.split,1,"my_fasta_file.fas",type=2, ref.chr="chrom1.fas")

}



\keyword{methods}

