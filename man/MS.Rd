\name{MS}
\docType{methods}
\alias{MS}

\title{Coalescent simulation with or without selection}


\description{

 This function uses Hudson's MS and Ewing's MSMS to compare simulated data with the observed data. 
 
}

\usage{

MS(GENO,niter=10,thetaID="user",params=FALSE,detail=FALSE,
neutrality=FALSE,linkage=FALSE,F_ST=FALSE,MSMS=FALSE,big.data=FALSE)

}

\arguments{

  \item{GENO}{an object of class \code{"GENOME"}}
  \item{niter}{number of samples per locus}
  \item{thetaID}{\code{"Tajima"},\code{"Watterson"} or \code{"user"}. default:\code{"user"}}
  \item{neutrality}{Calculate neutrality tests. default=\code{FALSE}}
  \item{linkage}{Calculate linkage disequilibrium. default=\code{FALSE}}
  \item{F_ST}{Calculate fixation index. default=\code{FALSE}}
  \item{params}{an object of class \code{"test.params"}. see \code{?test.params}}
  \item{detail}{detailed statistics. Note: slower computations! default=\code{FALSE}}
  \item{MSMS}{specify parameter for MSMS simulation with selection (has to be specified as a string)}
  \item{big.data}{if \code{TRUE} the ff-package is used}

}


\details{

You can choose different mutation rate estimators to generate simulation data.
When \code{thetaID="user"}, you have to define the theta values in an object of class \code{"test.params"}. 
The \code{"test.params"} class can also be used to specify some additional parameter like migration 
and/or recombination rates... (\code{?test.params}). \cr \cr 

Please read the MSMS documentation for the correct use of coalescent simulations to assess statistical significance.
 
}

\note{The executable file ms has to be stored in the current workspace. \cr If you want to use the MSMS application, put the msms folder including the corresponding executable files in the current workspace.\cr
Both programs can be obtained from their websites (see references).}

\value{

  The function creates an object of class \code{"cs.stats"} 

}


\references{ 

Hudson, R. R. (2002). \emph{Generating samples under a Wright-Fisher neutral model of genetic variation.}
Bioinformatics 18: 337-338 \cr \cr 

Gregory Ewing and Joachim Hermisson, \emph{MSMS: A Coalescent Simulation Program Including Recombination, Demographic Structure, and Selection at a Single Locus.} Bioinformatics 2010, doi: 10.1093/bioinformatics/btq322 

}




\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class <- neutrality.stats(GENOME.class,list(1:6))
# MS.class <- MS(GENOME.class,thetaID="Tajima",neutrality=TRUE)
# MS.class <- MS(GENOME.class,thetaID="Tajima",neutrality=TRUE,
#                MSMS="-N 1000 -SAA 200 -SaA 100 -SF 1e-2")
# MS.class
# MS.class@obs.val
# MS.class@locus[[1]] 

}

\keyword{methods}
