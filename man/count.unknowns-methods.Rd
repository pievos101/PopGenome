\name{count.unknowns-methods}

\docType{methods}

\alias{count.unknowns,GENOME-method}
\alias{count.unknowns-methods}


\title{Calculate missing nucleotide frequencies}

\description{
 A generic function to calculate the missing nucleotide frequencies.
}

\usage{

\S4method{count.unknowns}{GENOME}(object)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}

}


\value{

Returned value is a modified object of class \code{"GENOME"} \cr  \cr

The slot \code{GENOME.class@missing.freqs} for the missing frequencies for the whole region. \cr
The slot \code{GENOME.class@region.stats@missing.freqs} for the missing frequencies for each SNP in a given region \cr

}

\examples{

# GENOME.class <- readData("VCF", format="VCF", include.unknown=TRUE)
# GENOME.class@region.stats
# GENOME.class <- count.unknowns(GENOME.class)
# GENOME.class@missing.freqs
# GENOME.class@region.stats@missing.freqs

}


\keyword{methods}
