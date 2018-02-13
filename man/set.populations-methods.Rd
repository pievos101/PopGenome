\name{set.populations-methods}
\docType{methods}

\alias{set.populations,GENOME-method}
\alias{set.populations-methods}


\title{Define populations}


\description{
This generic function defines the populations.\cr
Using this function obviates the need to specify the populations for each calculation separately.\cr
The populations can be set differently for different PopGenome modules by applying the function between module calls.
}

\usage{

\S4method{set.populations}{GENOME}(object,new.populations=FALSE, diploid=FALSE,
triploid=FALSE,tetraploid=FALSE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{diploid}{if diploid data is present}
  \item{triploid}{if triploid data is present}
  \item{tetraploid}{if tetraploid data is present}
}



\examples{

# GENOME.class <- readData("\home\Alignments")
# pop.1 <- c("seq1","seq2")
# pop.2 <- c("seq3","seq4","seq1")
# GENOME.class <- set.populations(GENOME.class,list(pop.1,pop.2))
# GENOME.class@region.data@populations2
# GENOME.class <- neutrality.stats(GENOME.class)
}



\keyword{methods}

