\name{readData}

\docType{methods}
\alias{readData}


\alias{get.sum.data,GENOME-method}
\alias{get.sum.data-methods}

\title{Read alignments and calculate summary data}


\description{
 This function reads alignments/SNP data in several formats and calculates some summary data. 
}


\usage{

readData(path,populations=FALSE,outgroup=FALSE,include.unknown=FALSE,
         gffpath=FALSE,format="fasta",parallized=FALSE,
         progress_bar_switch=TRUE, FAST=FALSE,big.data=FALSE,
         SNP.DATA=FALSE
        )

\S4method{get.sum.data}{GENOME}(object)

}

\arguments{
 
  \item{object}{object of class \code{"GENOME"}}
  \item{path}{the basepath (folder) of the alignments}
  \item{outgroup}{vector of outgroup sequences}
  \item{include.unknown}{if positions with unknown nucleotides should be considered.}
  \item{populations}{list of populations. default:\code{FALSE}}
  \item{gffpath}{the basepath (folder) of the corresponding GFF-files. default:\code{FALSE}}
  \item{format}{data formats. \code{"fasta"} is default. See details !}
  \item{parallized}{parallel processing to accelerate the reading process. See details !}
  \item{progress_bar_switch}{progress_bar}
  \item{FAST}{fast computation. See details !}
  \item{big.data}{use the ff-package}
  \item{SNP.DATA}{important for reference positions; should be TRUE if you use SNP-data in alignment format}

}

\details{
All data (alignments or SNP-files) have to be stored in one folder. The folder is the input of this \cr 
function. If no GFF file (which also have to be stored in a folder) is specified, an alignment in \cr
the correct reading frame (starting at a first codon position) is expected. \cr
Otherwise synonymous and non-synonymous positions are not identified correctly. \cr \cr
Note: \cr 
The GFF-files have to be EXACTLY the same names (without any extensions like .fas or .gff) \cr 
as the files storing the nucleotide data to ensure correct matching \cr \cr

format: \cr 
\code{"fasta"},\code{"nexus"},\code{"phylip"}, \cr
\code{"MAF"},\code{"MEGA"} \cr
\code{"HapMap"},\code{"VCF"} \cr
\code{"RData"} \cr 
Valid nucleotides are T,t,U,u,G,g,A,a,C,c,N,n,- \cr \cr

parallized: \cr 
- will speed up calculations if you use a very large amount of alignments \cr \cr

FAST:   \cr 
- will not classify synonymous/non-synonymous SNPs directly \cr
- fast computation (via compiled C code) of biallelic matrix, biallelic sites, transversions/transitions \cr
  and biallelic substitutions \cr
- can be switched to \code{TRUE} in case of SNP data without loss of information \cr \cr                  

big.data: \cr 
- use the ff-package \cr 
- ff mechanism is used for biallelic.matrix and GFF/GTF information \cr 
- is automatically activated for readVCF or readSNP \cr
- Note! you should set this to TRUE if you use big chunks of data \cr 
  and you want to later concatenate them in the PopGenome framework \cr
  (for example: sliding windows of the whole dataset).\cr \cr

SNP.DATA: \cr 
- should be switched to \code{TRUE} if you use SNP-data in alignment format. \cr 
- the corresponding SNP positions can be set via \code{set.ref.positions}
}

\value{

  The function creates an object of class \code{"GENOME"} \cr \cr
  --------------------------------------------------------- \cr
  The following slots will be filled in the \code{"GENOME"} object \cr
  --------------------------------------------------------- \cr

            \tabular{rll}{
            \tab Slot                	                \tab Description\cr
          1.\tab \code{n.sites}                         \tab total number of sites\cr
          2.\tab \code{n.biallelic.sites}          	\tab number of biallelic sites \cr
          3.\tab \code{n.gaps}                          \tab number of sites with gaps\cr      
          4.\tab \code{n.unknowns}          	        \tab number of sites with unknown nucleotides\cr
	  5.\tab \code{n.valid.sites}            	\tab number of valid sites\cr
          6.\tab \code{n.polyallelic.sites}             \tab number of sites with >2 nucleotides \cr      
          7.\tab \code{trans.transv.ratio}          	\tab transition/transversion ratio of biallelic sites\cr
          8.\tab \code{region.names}          	        \tab names of regions\cr
	  9.\tab \code{region.data}            	        \tab some detailed information about the data read\cr
          
        }


}


\examples{

# GENOME.class <- readData("...\Alignments", FAST=TRUE)
# GENOME.class <- readData("VCF", format="VCF")
# Note, "Alignments" and "VCF" are folders !
# GENOME.class@region.names
# GENOME.class <- readData("...\Alignments", big.data=TRUE)
# object.size(GENOME.class)
# GENOME.class <- readData("...\Alignments",gffpath="...\Alignments_GFF")
# GENOME.class
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
