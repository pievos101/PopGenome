\name{load.session}

\docType{methods}
\alias{load.session}


\title{Loading a PopGenome session}


\description{
This function loads a PopGenome session (more precisely: the corresponding \code{"GENOME"} object) from the current workspace.
}

\usage{

load.session(folder)

}

\arguments{
 \item{folder}{name of the folder/object}
}  

\details{
This function has to be used in the same workspace (folder) where \cr
the object of class \code{"GENOME"} was saved.
}

\value{
An object of class \code{"GENOME"}.
}


\examples{

# GENOME.class <- readData("...\Alignments")
# save.session(GENOME.class,folder="GENOME.class")
# q()
# R
# library(PopGenome)
# load.session("GENOME.class")

}



\keyword{methods}
