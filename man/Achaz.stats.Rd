\name{Achaz.stats-methods}
\docType{methods}
\alias{Achaz.stats,GENOME-method}
\alias{Achaz.stats-methods}

\title{Achaz statistic}


\description{
Achaz statistic
}

\usage{

\S4method{Achaz.stats}{GENOME}(object,new.populations=FALSE,new.outgroup=FALSE,subsites=FALSE)


}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{new.outgroup}{outgroup vector. default:\code{FALSE}}
  \item{subsites}{ 

                  \code{"transitions"}: SNPs that are transitions. \cr
                  \code{"transversions"}: SNPs that are transversions. \cr
                  \code{"syn"}: synonymous sites. \cr
                  \code{"nonsyn"}:  nonsynonymous sites. \cr
                  \code{"exon"}: SNPs in exon regions. \cr
                  \code{"intron"}: SNPs in intron regions. \cr
                  \code{"coding"}: SNPs in coding regions (CDS). \cr
                  \code{"utr"}:  SNPs in UTR regions. \cr
                  \code{"gene"}: SNPs in genes. \cr
          default:\code{FALSE}}
  
}



\value{


returned value is a modified object of class \code{"GENOME"} \cr \cr
-------------------------------------------------------------- \cr
The following Slots will be modified in the \code{"GENOME"} object \cr
-------------------------------------------------------------- \cr

   \item{Yach}{Achaz Y statistic}
   
}

\references{

Achaz G.,2008 \emph{Testing for neutrality in samples with sequencing errors.} Genetics 179: 1409.

}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- Achaz.stats(GENOME.class)
# GENOME.class <- Achaz.stats(GENOME.class,list(1:7,8:12))
# show the result:
# GENOME.class@Yach
 

}



\keyword{methods}
