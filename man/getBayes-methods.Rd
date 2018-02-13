 \name{getBayes-methods}

\docType{methods}

\alias{getBayes,GENOME-method}
\alias{getBayes-methods}


\title{Get values for BayeScanR}

\description{
This function returns the values that are necessary to run BayeScanR.
}

\usage{

\S4method{getBayes}{GENOME}(object,snps=FALSE)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{snps}{SNPs are considered seperately}

}


\value{

coming soon !

}

\references{

Foll M and OE Gaggiotti (2008). \emph{A genome scan method to identify selected loci appropriate for both dominant and codominant markers: A Bayesian perspective.} Genetics 180: 977-993 \cr \cr

}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class <- F_ST.stats(GENOME.class,list(1:4,5:10))
# Bayes.input  <- getBayes(GENOME.class)
# Bayes.class  <- BayeScanR(Bayes.input)

}


\keyword{methods}
