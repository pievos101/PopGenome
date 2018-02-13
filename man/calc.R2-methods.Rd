\name{calc.R2-methods}
\docType{methods}
\alias{calc.R2,GENOME-method}
\alias{calc.R2-methods}


\title{Linkage statistics (R2, P-value, Distance)}


\description{
 This generic function calculates some linkage disequilibrium statistics.
}

\usage{

\S4method{calc.R2}{GENOME}(object,subsites=FALSE,lower.bound=0,upper.bound=1)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{subsites}{same as in the other modules}
  \item{lower.bound}{sites with minor-allele-frequency>=lower.bound are considered}
  \item{upper.bound}{sites with minor-allele-frequency<=upper.bound are considered}
}

\details{
Note, the pairwise comparisons are computed via \code{combn(n.snps,2)}.
}

\value{
The slot GENOME.class@region.stats@linkage.disequilibrium will be filled.\cr 
(R2,P-value,Distance) \cr  
Fisher's Exact Test is used for the P-values.

}

\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class
# GENOME.class <- calc.R2(GENOME.class)
# show the result:
# GENOME.class@region.stats@linkage.disequilibrium
# [[x]][[y]] x:region, y:population 

}



\keyword{methods}
