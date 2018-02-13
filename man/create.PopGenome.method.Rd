\name{create.PopGenome.method}

\docType{methods}
\alias{create.PopGenome.method}


\title{Integration of own functions into the PopGenome-framework}


\description{
This function generates a skeleton for a PopGenome function.
It thereby facilitates the effortless integration of new methods into the PopGenome framework.  
}

\usage{

create.PopGenome.method(function.name,population.specific=TRUE)

}

\arguments{

  \item{function.name}{name of your function}
  \item{population.specific}{\code{TRUE}:function returns one value per population.\code{FALSE}:function returns one value calculated across all populations (as in the case of FST measurements)}

}  

\details{
This mechanism enables you to use your own functions in the PopGenome environment. The functions can also be applied to sliding windows or subsites.\cr
Please look at the generated function, which documents where to place your own function in detail. 
}


\examples{

# GENOME.class <- readData(".../Alignments")
# create.PopGenome.method("myFunction")
# edit myFunction.R
# source("myFunction")
# value <- myFunction(test)
# value

}



\keyword{methods}
