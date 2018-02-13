\name{save.session}

\docType{methods}
\alias{save.session}


\title{Save the \code{"GENOME"} object of a PopGenome session}


\description{
This function saves the \code{"GENOME"} object of a PopGenome session to the current workspace.
The object can be loaded again with load.session().
}

\usage{

save.session(object,folder)

}

\arguments{

  \item{object}{object of class \code{"GENOME"}}
  \item{folder}{name of the folder/object}

}  

\details{
Saving R and ff-objects created by the ff-package in a folder.
}


\examples{

# GENOME.class <- readData("...\Alignments")
# save.session(GENOME.class,"GENOME.class")
# load.session("GENOME.class")
}



\keyword{methods}
