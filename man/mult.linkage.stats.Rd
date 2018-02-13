\name{mult.linkage.stats-methods}
\docType{methods}
\alias{mult.linkage.stats,GENOME-method}
\alias{mult.linkage.stats-methods}


\title{Multilocus linkage statistics}


\description{
 This generic function calculates the linkage disequilibrium between regions.
}

\usage{

\S4method{mult.linkage.stats}{GENOME}(object,lower.bound=0,upper.bound=1,pairs=FALSE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{lower.bound}{sites with minor-allele-frequency>=lower.bound are considered}
  \item{upper.bound}{sites with minor-allele-frequency<=upper.bound are considered}
  \item{pairs}{permutation matrix of pairwise comparisons}
}

\details{

\code{pairs} is a matrix. Each column contains the pairwise comparison region IDs.\cr

1 1 \cr
2 3 \cr

compares region 1 with 2, and region 1 with 3.



}

\value{


The return value is a modified object of class \code{"GENOME"} \cr \cr
-------------------------------------------------------------- \cr
The following slots will be modified in the \code{"GENOME"} object \cr
-------------------------------------------------------------- \cr

   \item{mult.Linkage}{Some linkage statistics for each pair of regions (R2, P-value, Distance)}
  
The Fisher-Exact-Test is used to calculate the P-values.

}

\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class
# GENOME.class <- mult.linkage.stats(GENOME.class)
# show the result:
# GENOME.class@mult.Linkage
 

}



\keyword{methods}
