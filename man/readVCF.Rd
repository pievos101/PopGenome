\name{readVCF}

\docType{methods}
\alias{readVCF}

\title{Read SNP data in tabixed VCF format}


\description{
 This function reads tabixed VCF-files, as distributed from the 1000 Genomes project (human). 
}

\usage{

readVCF(filename, numcols, tid, frompos, topos,
        samplenames=NA, gffpath = FALSE, include.unknown=FALSE, approx=FALSE,
	out="", parallel=FALSE)

}

\arguments{
 
  \item{filename}{the corresponding tabixed VCF-file}
  \item{numcols}{number of SNPs that should be read in as a chunk}
  \item{tid}{which chromosome ? (character)}
  \item{frompos}{start of the region}
  \item{topos}{end of the region}
  \item{samplenames}{a vector of individuals}
  \item{gffpath}{the corresponding GFF file}
  \item{include.unknown}{includ positions with unknown/missing nucleotides}
  \item{approx}{see details !}
  \item{out}{a folder suffix where the temporary files should be saved} 
  \item{parallel}{parallel computation using mclapply} 
}

\details{
The readVCF function expects a tabixed VCF file with a diploid GT field. \cr
In case of haploid data, the GT field has to be transformed to a pseudo-diploid \cr
field (such as 0 -> 0|0). An alternative is to use readData(..., format="VCF"),\cr
which can read non-tabixed haploid and any kind of polyploid VCFs directly. \cr 
When \code{approx=TRUE}, the algorithm will apply a logical OR to the GT-field:\cr
(0|0=0,1|0=1,0|1=1,1|1=1). Note, this is an approximation for diploid data, which will \cr 
speed up calculations. In case of haploid data, \code{approx} should be switched to \code{TRUE}. \cr  
If \code{approx=FALSE}, the full diploid information will be considered.\cr
The ff-package PopGenome uses to store the SNP information limits total data size to \cr
individuals * (number of SNPs) <= .Machine$integer.max \cr
In case of very large data sets, the bigmemory package will be used; \cr
this will slow down calculations (e.g. this package have to be installed first !!!).\cr
Use the function \code{vcf_handle <-.Call("VCF_open", filename)} \cr 
to open a VCF-file and \code{.Call("VCF_getSampleNames",vcf_handle)} \cr 
to get and define the individuals which should be considered in the analysis. \cr
See also readData(..., format="VCF") ! \cr
}

\value{


  The function creates an object of class \code{"GENOME"} \cr \cr
  --------------------------------------------------------- \cr
  The following slots will be filled in the \code{"GENOME"} object \cr
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

# GENOME.class <- readVCF("...\chr1.vcf.gz", 1000, "1", 1, 100000)
# GENOME.class
# GENOME.class@region.names
# GENOME.class <- neutrality.stats(GENOME.class,FAST=TRUE)
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
