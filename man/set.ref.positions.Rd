\name{set.ref.positions-methods}
\docType{methods}
\alias{set.ref.positions,GENOME-method}
\alias{set.ref.positions-methods}


\title{Set reference positions for SNP data}


\description{
 This generic function sets the positions of the SNP data.
 Should be used if you use alignment formats to store SNP data (i.e., data restricted to the polymorphic positions). 
}

\usage{

\S4method{set.ref.positions}{GENOME}(object, positions)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{positions}{a list of reference positions}
}



\value{


returned value is a modified object of class \code{"GENOME"} \cr 
  
  
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- set.ref.positions(GENOME.class,list(c(1000,2001,3000),
#                 c(3200,12000)))
 

}



\keyword{methods}
