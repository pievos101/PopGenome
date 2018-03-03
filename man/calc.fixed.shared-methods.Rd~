\name{calc.fixed.shared-methods}

\docType{methods}

\alias{calc.fixed.shared,GENOME-method}
\alias{calc.fixed.shared-methods}


\title{Fixed and shared polymorphisms}

\description{
 A generic function to calculate the number of fixed and shared polymorphisms.
}

\usage{

\S4method{calc.fixed.shared}{GENOME}(object,
subsites=FALSE,
new.populations=FALSE,
fixed.threshold=1,
fixed.threshold.fst=1)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{subsites}{ 

                  \code{"transitions"}: SNPs that are transitions. \cr
                  \code{"transversions"}: SNPs that are transversions. \cr
                  \code{"syn"}: synonymous sites. \cr
                  \code{"nonsyn"}:  nonsynonymous sites. \cr
                  \code{"exon"}: SNPs in exon regions. \cr
                  \code{"intron"}: SNPs in intron regions. \cr
                  \code{"coding"}: SNPs in coding regions (CDS). \cr
                  \code{"utr"}:  SNPs in UTR regions. \cr
                  \code{"gene"}: SNPs in genes. \cr
		  \code{"intergenic"} : SNPs in intergenic regions. \cr
  } 
  \item{new.populations}{list of populations. default=\code{FALSE}}
  \item{fixed.threshold}{Polymorphisms are considered as fixed >= threshold value}
  \item{fixed.threshold.fst}{Polymorphisms are considered as fixed >= threshold value}
  

}


\value{

Returned value is a modified object of class \code{"GENOME"} \cr  \cr
--------------------------------------------------------- \cr
The following slots will be modified in the \code{"GENOME"}  object \cr
--------------------------------------------------------- \cr

\tabular{rlll}{
            \tab Slot                	    \tab Reference   \tab Description\cr
          1.\tab \code{n.fixed.sites}	\tab [x]       \tab Number of fixed sites \cr
          2.\tab \code{n.shared.sites}	                \tab [x]         \tab Number of shared sites    \cr 
          3.\tab \code{n.monomorphic.sites}	\tab [x]         \tab Number of monomorphic sites \cr
          }

}

\details{
The nucleotide diversities have to be devided by \code{GENOME.class@n.sites} to give diversities per site.
}

\references{
[x]
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# set.populations
# GENOME.class <- calc.fixed.shared(GENOME.class)
}


\keyword{methods}
