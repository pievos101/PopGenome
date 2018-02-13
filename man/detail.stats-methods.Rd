\name{detail.stats-methods}
\docType{methods}
\alias{detail.stats,GENOME-method}
\alias{detail.stats-methods}
\alias{get.detail,GENOME-method}
\alias{get.detail-methods}


\title{Several statistics}


\description{
 This generic function calculates some mixed statistics.
}



\usage{

\S4method{detail.stats}{GENOME}(
	object,
	new.populations=FALSE,
	new.outgroup=FALSE,
	subsites=FALSE,
	biallelic.structure=FALSE,
	mismatch.distribution=FALSE,
	site.spectrum=TRUE,
        site.FST=FALSE
        )
\S4method{get.detail}{GENOME}(object, biallelic.structure=FALSE)



}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations.}
  \item{new.outgroup}{outgroup sequences.}
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
  }

  \item{biallelic.structure}{fixed and shared polymorphisms (stored in \code{GENOME.class@region.stats}).}
  \item{mismatch.distribution}{statistics based on mismatch distribution}
  \item{site.spectrum}{minor allele frequency of each SNP}
  \item{site.FST}{computes FST for each SNP}
}


\value{


The return value is a modified object of class \code{"GENOME"} \cr \cr
------------------------------------------------------------------ \cr
The following Slots will be modified in the \code{"GENOME"} object \cr
------------------------------------------------------------------ \cr

   \item{MDSD}{...}
   \item{MDG1}{...}
   \item{MDG2}{...}
   \item{region.stats}{the slot  \code{biallelic.structure} and \code{minor.allele.freqs} will be filled} \cr

The function \code{get.detail(GENOME.class, biallelic.structure=TRUE)} \cr
returns a matrix for each region, where \cr

\item{0}{population is polymorphic, the remaining individuals are polymorphic} \cr
\item{1}{population is polymorphic, the remaining individuals are monomorphic} \cr
\item{2}{population is monomorphic, the remaining individuals are polymorphic} \cr  
\item{3}{population is monomorphic, the remaining individuals are monomorphic with the same value}  \cr  
\item{4}{population is monomorphic, the remaining individuals are monomorphic with different values} \cr

  
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- set.populations(GENOME.class,list(1:10))
# GENOME.class <- detail.stats(GENOME.class)
# show the result:
# mismatch.values   <- get.detail(GENOME.class)
# bial.struc.values <- get.detail(GENOME.class, biallelic.structure=TRUE)
# GENOME.class@region.stats@biallelic.structure
# GENOME.class@region.stats@biallelic.structure[[1]]
 

}



\keyword{methods}
