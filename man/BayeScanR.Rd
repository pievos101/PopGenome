\name{BayeScanR}

\docType{methods}
\alias{BayeScanR}


\title{An R implementation of BayeScan (Foll \& Gagiotti 2008)}


\description{

BayeScanR is an R implementation of BayeScan for analysis of codominant markers.

}

\usage{

BayeScanR(input,nb.pilot=10,pilot.runtime=2500,main.runtime=100000, discard=50000)

}

\arguments{

  \item{input}{textfile or an R-object returned by getBayes()}
  \item{nb.pilot}{number of pilot runs}
  \item{pilot.runtime}{length of pilot runs}
  \item{main.runtime}{length of main runs} 
  \item{discard}{how many runs in the main.loop should be discarded?} 

} 

\value{

returned value is an object of class \code{"BAYESRETURN"} \cr \cr
-------------------------------------------------------------- \cr
The following Slots will be filled \cr
-------------------------------------------------------------- \cr

   \item{alpha}{alpha effects}
   \item{beta}{beta effects}
   \item{var_alpha}{variance of alpha values}
   \item{a_inc}{which alpha is included in the model}
   \item{fst}{FST values}
   \item{P}{P-value}


}

\references{ 
[1] Foll M and OE Gaggiotti (2008). \emph{A genome scan method to identify selected loci appropriate for both dominant and codominant markers: A Bayesian perspective.} Genetics 180: 977-993 \cr \cr
}




\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class <- F_ST.stats(GENOME.class,list(1:5,6:10))
# Bayes.input  <- getBayes(GENOME.class)
# BAYES.class  <- BayeScanR(Bayes.input)
# BAYES.class

}



\keyword{methods}
