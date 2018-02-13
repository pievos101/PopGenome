\name{recomb.stats-methods}
\docType{methods}
\alias{recomb.stats,GENOME-method}
\alias{recomb.stats-methods}
\alias{get.recomb,GENOME-method}
\alias{get.recomb-methods}


\title{Recombination statistics}


\description{
 This generic function calculates the Four-Gamete test (Hudson 1985).
}

\usage{

\S4method{recomb.stats}{GENOME}(object,new.populations=FALSE,subsites=FALSE)
\S4method{get.recomb}{GENOME}(object)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
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
          default:\code{FALSE}}
  
}



\value{


The return value is a modified object of class \code{"GENOME"} \cr \cr
-------------------------------------------------------------- \cr
The following slots will be modified in the \code{"GENOME"} object \cr
-------------------------------------------------------------- \cr

   \item{Hudson.RM}{Four-gamete test}
  
}

\references{
Hudson, R. K. (1985). \emph{Statistical Properties of the Number of Recombination Events in the History of a Sample of DNA Sequences}
Genetics 111 (1): 147-164.
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- recomb.stats(GENOME.class)
# GENOME.class <- recomb.stats(GENOME.class,list(1:7,8:12))
# show the result:
# recomb.values <- get.recomb(GENOME.class)
# recomb.values[[1]] # first population !
# GENOME.class@region.stats@Hudson.RM
 

}



\keyword{methods}
