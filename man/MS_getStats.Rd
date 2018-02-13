\name{MS_getStats}

\docType{methods}
\alias{MS_getStats}


\title{Get the simulated MS/MSMS statistics}


\description{
This function extracts the simulated values from the class \code{cs.stats}
}

\usage{

MS_getStats(object,locus=1,population=1)

}

\arguments{ 
 \item{object}{object of class \code{"cs.stats"}}
 \item{locus}{the locus ID}
 \item{population}{the population ID}
}  


\value{
The return value is a matrix containing 
the simulation results of different statistical tests.\cr
(see MS())
}


\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class <- neutrality.stats(GENOME.class)
# ms <- MS(GENOME.class,thetaID="Tajima",neutrality=TRUE)
# MS_getStats(ms)

}



\keyword{methods}
