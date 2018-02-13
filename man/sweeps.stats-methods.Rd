\name{sweeps.stats-methods}
\docType{methods}
\alias{sweeps.stats,GENOME-method}
\alias{sweeps.stats-methods}
\alias{get.sweeps,GENOME-method}
\alias{get.sweeps-methods}


\title{Selective Sweeps}


\description{
 This module calculates some statistics to detect selective sweeps.
}

\usage{

\S4method{sweeps.stats}{GENOME}(object,new.populations=FALSE,subsites=FALSE,
 freq.table=FALSE, FST=FALSE)
\S4method{get.sweeps}{GENOME}(object)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{subsites}{ 

                  \code{"transitions"}: SNPs that are transitions. \cr
                  \code{"transversions"}: SNPs that are transversions. \cr
                  \code{"syn"}: synonymous sites. \cr
                  \code{"nonsyn"}:  non-synonymous sites. \cr
                  \code{"exon"}: SNPs in exon regions. \cr
                  \code{"intron"}: SNPs in intron regions. \cr
                  \code{"coding"}: SNPs in coding regions (CDS). \cr
                  \code{"utr"}:  SNPs in UTR regions. \cr
                  \code{"gene"}: SNPs in genes. \cr
          default:\code{FALSE}}
\item{freq.table}{the frequency counts for the CLR test. \code{"list"}}
\item{FST}{use FST values instead of the minor allele frequencies}
  
}



\value{


The return value is a modified object of class \code{"GENOME"} \cr \cr
-------------------------------------------------------------- \cr
The following slots will be modified in the \code{"GENOME"} object \cr
-------------------------------------------------------------- \cr

   \item{CL}{Composite Likelihood of SNPs}
   \item{CLR}{Nielsen's CLR test}

}

\details{
The \code{freq.table} contains the global sets of frequency counts.
It can be produced with the module \code{detail.stats}.
The values in the slot \code{GENOME.class@region.stats@minor.allele.frequencies}
can be used to create this global set. (use the R function \code{table})
\code{freq.table} is a list of length n.pops.
}


\references{
Cai JJ (2008) \emph{PGEToolbox: A Matlab toolbox for population genetics and evolution}
Journal of Heredity Jul-Aug;99(4):438-40.doi:10.1093/jhered/esm127 \cr \cr 
Nielson, R. (2005). \emph{Genomic scans for selective sweeps using SNP data}
Genome Res. 2005 15: 1566-1575 \cr \cr
}

\examples{

# Reading one alignment stored in the folder Aln
# GENOME.class <- readData("\home\Aln")
#
# CL
# GENOME.class <- sweeps.stats(GENOME.class)
# GENOME.class@CL
#
# CLR
# create global set
# GENOME.class  <- detail.stats(GENOME.class)
# freq <- GENOME.class@region.stats@minor.allele.freqs[[1]]
# freq.table <- list()
# freq.table[[1]] <- table(freq)
# define the region of interest
# GENOME.class.split <- splitting.data(GENOME.class, positions= ...)
# calculate CLR
# GENOME.class.split <- sweeps.stats(GENOME.class.split, freq.table=freq.table)
# GENOME.class@CLR
 

}



\keyword{methods}
