\name{PG_plot.biallelic.matrix-methods}
\docType{methods}
\alias{PG_plot.biallelic.matrix,GENOME-method}
\alias{PG_plot.biallelic.matrix-methods}


\title{Plot the biallelic matrix}


\description{
This function plots the biallelic matrix of a specific region. 
}

\usage{

\S4method{PG_plot.biallelic.matrix}{GENOME}(object,region, ind.names = FALSE , cex.axis = 0.5,
title="")

}

\arguments{

  \item{object}{object of class \code{"GENOME"}}
  \item{region}{the region ID}
  \item{ind.names}{individual names/IDs. default:ALL}
  \item{cex.axis}{size of text (y-axis)}
  \item{title}{Title of the plot}


}  


\examples{

# GENOME.class <- readData("...\Alignments")
# PG_plot.biallelic.matrix(GENOME.class, region = 1)

}



\keyword{methods}
