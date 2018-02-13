\name{concatenate.regions}

\docType{methods}
\alias{concatenate.regions}

\title{Concatenate regions}


\description{
 This function concatenates the regions/chunks contained in one \code{GENOME} object. 
}

\usage{

concatenate.regions(object)

}

\arguments{
 
  \item{object}{object of class \code{GENOME}}
   
}


\value{
  The function creates an object of class "GENOME".
}


\examples{

# GENOME.class <- readData("Three_Alignments/")
# WHOLE        <- concatenate.regions(GENOME.class)
# WHOLE        <- neutrality.stats(WHOLE)
# WHOLE@Tajima.D

}

\keyword{methods}
