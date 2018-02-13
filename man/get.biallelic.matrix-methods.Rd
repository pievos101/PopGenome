\name{get.biallelic.matrix-methods}

\docType{methods}

\alias{get.biallelic.matrix,GENOME-method}
\alias{get.biallelic.matrix-methods}


\title{Get the biallelic matrix}

\description{
This function returns the biallelic matrix of a specific region.
}

\usage{

\S4method{get.biallelic.matrix}{GENOME}(object,region)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{region}{ID of the region}

}


\value{

Biallelic matrix \cr
rows:    names of individuals \cr
columns: biallelic sites \cr

}


\examples{

# GENOME.class <- readData("\home\Alignments")
# get.biallelic.matrix(GENOME.class,7) # biallelic matrix of the 7th alignment
}


\keyword{methods}
