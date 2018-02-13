\name{set.filter-methods}

\docType{methods}

\alias{set.filter,GENOME-method}
\alias{set.filter-methods}


\title{Setting filter to the analysis}

\description{
 A generic function to set filter regarding e.g missing data.
}

\usage{

\S4method{set.filter}{GENOME}(object,
missing.freqs=TRUE,
minor.freqs=FALSE,
maf.lower.bound=0,
maf.upper.bound=1,
miss.lower.bound=0,
miss.upper.bound=1)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{missing.freqs}{Set filter for missing data}
  \item{minor.freqs}{Set filter for the MAF}
  \item{maf.lower.bound}{frequency of the MAF}
  \item{maf.upper.bound}{...}
  \item{miss.lower.bound}{frequency of the missing freq.}
  \item{miss.upper.bound}{...}
  

}


\value{

Returned value is a modified object of class \code{"GENOME"} \cr  \cr
--------------------------------------------------------- \cr
The following slots will be modified in the \code{"GENOME"}  object \cr
--------------------------------------------------------- \cr

\tabular{rlll}{
            \tab Slot                	    \tab Reference   \tab Description\cr
          1.\tab \code{region.data@included}	\tab [x]       \tab ... \cr
          }

}

\details{
This function sets the slot \code{region.data@included}.
}

\references{
[x]
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class <- set.filter(GENOME.class, missing.freqs=TRUE,
# miss.lower.bound=0, miss.upper.bound=0.2)
# now apply any statistic to the filtered data set.
# GENOME.class <- calc.fixed.shared(GENOME.class, subsites="included")
}


\keyword{methods}
