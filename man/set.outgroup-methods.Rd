\name{set.outgroup-methods}
\docType{methods}

\alias{set.outgroup,GENOME-method}
\alias{set.outgroup-methods}


\title{Define an outgroup}


\description{
This generic function defines the outgroup by matching the specified vector against
each region.  \cr 
}

\usage{

\S4method{set.outgroup}{GENOME}(object,new.outgroup=FALSE, diploid=FALSE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.outgroup}{a vector of outgroup individuals}
  \item{diploid}{if diploid data is present}
}



\examples{

# GENOME.class <- readData("\home\Alignments")
# outgroup <- c("seq1","seq2")
# GENOME.class <- set.outgroup(GENOME.class,new.outgroup=outgroup)
# GENOME.class <- neutrality.stats(GENOME.class)
}



\keyword{methods}

