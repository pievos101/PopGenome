\name{MKT-methods}
\docType{methods}
\alias{MKT,GENOME-method}
\alias{MKT-methods}
\alias{get.MKT,GENOME-method}
\alias{get.MKT-methods}

\title{McDonald-Kreitman Test (McDonald \& Kreitman 1991)}


\description{
 This generic function calculates an approximate version of the McDonald-Kreitman Test.
}

\usage{

\S4method{MKT}{GENOME}(object,
new.populations=FALSE,
do.fisher.test=FALSE,
fixed.threshold.fst=FALSE,
subsites=FALSE)

\S4method{get.MKT}{GENOME}(object)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{do.fisher.test}{P-value calculation out of the Dn,Ds,Pn,Ps table}
  \item{fixed.threshold.fst}{Fixed threshold}
  \item{subsites}{Subsites}
  
}

\details{

This approximate version of the McDonald-Kreitman test assumes that the probability that \cr 
two single nucleotide polymorphisms (SNPs) occur in the same codon is very small. Thus, \cr
only codons with a single SNP are examined. \cr 
If no gff-file was specified when the data was read in, it is assumed that the \cr
alignment is in the correct reading frame (starting at a first codon position). \cr
The outgroup has to be defined as a population !

}



\value{


  The return value is a modified object of class "GENOME" \cr \cr
  --------------------------------------------------------- \cr
  The following slots will be modified in the "GENOME" object \cr
  --------------------------------------------------------- \cr

  \item{MKT}{a matrix which includes the following values:}

   \tabular{rll}{
            \tab Columns               	                 \tab Description\cr
          1.\tab \code{P_nonsyn}                         \tab nonsynonymous sites\cr
          2.\tab \code{P_syn}         	                 \tab synonymous sites\cr
          3.\tab \code{D_nonsyn}                         \tab fixed nonsynonymous sites\cr      
          4.\tab \code{D_syn}         	                 \tab fixed synonymous sites\cr
	  5.\tab \code{neutrality.index}           	 \tab $(P_{nonsyn}/P_{syn})/(D_{nonsyn}/D_{syn})$\cr
          6.\tab \code{alpha}                            \tab 1-neutrality.index\cr        
        }
  
  
}

\references{ 

McDonald, J. H.; Kreitman, M. (1991). \emph{Adaptive protein evolution at the Adh locus in Drosophila.} Nature 351 (6328): 652-654 \cr \cr 

}


\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- MKT(GENOME.class)
# GENOME.class <- MKT(GENOME.class,list(1:7,8:12))
# show the result:
# get.MKT(GENOME.class)
}



\keyword{methods}
