 \name{F_ST.stats.2-methods}

\docType{methods}

\alias{F_ST.stats.2,GENOME-method}
\alias{F_ST.stats.2-methods}


\title{Fixation Index (2)}

\description{
 A generic function to calculate some FST measurenments.
}

\usage{

\S4method{F_ST.stats.2}{GENOME}(object,new.populations="list",subsites=FALSE,snn=TRUE,Phi_ST=FALSE)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default=\code{FALSE}}
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
		  \code{"intergenic"} : SNPs in intergenic regions. \cr
  } 

  \item{snn}{Snn statistic from Hudson}
  \item{Phi_ST}{Statistic from Excoffier et al.}

}


\value{

Returned value is an modified object of class \code{"GENOME"} \cr  \cr
--------------------------------------------------------- \cr
Following slots will be modified in the \code{"GENOME"}  object \cr
--------------------------------------------------------- \cr

\tabular{rlll}{
            \tab Slot                	\tab Reference   \tab Description\cr
          1.\tab \code{Hudson.Snn}	\tab [1]         \tab Snn statistic from Hudson (2000) \cr
          2.\tab \code{Phi_ST}	        \tab [2]         \tab Phi_ST from Excoffier (1992)      \cr 
          }

}

\references{
[1] Hudson, R. R. (2000).\emph{A new statistic for detecting genetic differentiation.} Genetics 155: 2011-2014.\cr

[2] Excoffier, L., Smouse, P., Quattro, J. (1992),\emph{Analysis of molecular variance inferred from met- \cr
    ric distances among DNA haplotypes: application to human mitochondrial DNA restriction data.} \cr
    Genetics 131: 479-91
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- F_ST.stats.2(GENOME.class)
# GENOME.class <- F_ST.stats.2(GENOME.class,list(1:4,5:10))
# GENOME.class <- F_ST.stats.2(GENOME.class,
# list(c("seq1","seq5","seq3"),c("seq2","seq8")))
# show the result:
# GENOME.class@Hudson.Snn

}


\keyword{methods}
