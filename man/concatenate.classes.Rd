\name{concatenate.classes}

\docType{methods}
\alias{concatenate.classes}

\title{Concatenate GENOME classes}


\description{
 This function concatenates objects of class GENOME, allowing to stitch together larger datasets from smaller objects.
}

\usage{

concatenate.classes(classlist)

}

\arguments{
 
  \item{classlist}{a list of GENOME objects}
  
}


\value{
  The function creates an object of class "GENOME".
}


\examples{

# a   <- readData("Three_Alignments/")
# b   <- readData("Two_Alignments/")
# ab  <- concatenate.classes(list(a,b))
# ab  <- neutrality.stats(ab)
# ab@Tajima.D
# ab@region.names

}

\keyword{methods}
