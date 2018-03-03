\name{readSNP}

\docType{methods}
\alias{readSNP}

\title{Read data in .SNP format}


\description{
 This function reads data in .SNP (quality_variant) format, as distributed by the 1001 Genomes project (Arabidopsis). 
}

\usage{

readSNP(folder,populations=FALSE,outgroup=FALSE,gffpath=FALSE,
	CHR=FALSE,ref.chr=FALSE,snp.window.size=FALSE,
	parallized=FALSE,ffpackagebool=TRUE,
	include.unknown=FALSE
	)

}

\arguments{
 
  \item{folder}{the basepath of the variant_calls}
  \item{outgroup}{vector of outgroup sequences}
  \item{populations}{list of populations}
  \item{gffpath}{the corresponding GFF file}
  \item{CHR}{which chromosome ?, default: all chromosomes}
  \item{ref.chr}{reference chromosome (to classify synonymous/non-synonymous positions)}
  \item{snp.window.size}{scan SNP chunks}
  \item{parallized}{multicore computation}
  \item{ffpackagebool}{use the ff-package to save memory space. (slower)}
  \item{include.unknown}{include positions with unknown nucleotides}
}

\details{
The ff-package we use to store the SNP information limits the data size to\cr
individuals * (number of SNPs) <= .Machine$integer.max \cr
The text files containing the SNP information of each individual have to be stored in one folder. \cr
The slots \code{transitions}, \code{biallelic.sites}, and \code{biallelic.substitutions} of the class \cr
\code{"regions.data"} will be filled. \cr
At this time, if a GFF/GTF is used the data should be organized in a way that \cr
the "CHR" is a numerical value. The prefix "Chr" or "chr" is also supported.
}

\value{


  The function creates an object of class \code{"GENOME"} \cr \cr
  --------------------------------------------------------- \cr
  Following Slots will be filled in the \code{"GENOME"} object \cr
  --------------------------------------------------------- \cr

            \tabular{rll}{
            \tab Slot                	        \tab Description\cr
          1.\tab \code{n.sites}                 \tab total number of sites\cr
          2.\tab \code{n.biallelic.sites}       \tab number of biallelic sites \cr
          3.\tab \code{region.data}            	\tab some detailed information about the data read\cr
          4.\tab \code{region.names}            \tab names of regions \cr 
        }


}


\examples{

# GENOME.class <- readSNP("...\SNPData")
# GENOME.class <- readSNP("...\SNPData", CHR=1)
# GENOME.class <- readSNP("...\SNPData", CHR=1, gffpath="Gff_file.gff")
# GENOME.class
# GENOME.class <- neutrality.stats(GENOME.class,FAST=TRUE)
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
