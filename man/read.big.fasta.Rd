\name{read.big.fasta}

\docType{methods}
\alias{read.big.fasta}

\title{Reading large FASTA alignments}


\description{
 This function splits FASTA alignments that are too large to fit into the computer memory into chunks. 
}

\usage{

read.big.fasta(filename,populations=FALSE,outgroup=FALSE,window=2000,
               SNP.DATA=FALSE,include.unknown=FALSE,
               parallized=FALSE,FAST=FALSE,big.data=TRUE)

}

\arguments{
 
  \item{filename}{the basepath of the FASTA alignment}
  \item{outgroup}{vector of outgroup sequences}
  \item{populations}{list of populations}
  \item{window}{chunk size: number of columns/nucleotide sites}
  \item{SNP.DATA}{should be switched to TRUE if you use SNP data in alignment format}
  \item{include.unknown}{include unknown positions in the biallelic.matrix}
  \item{parallized}{Use parallel computations to speed up the reading - works only on UNIX systems! }
  \item{FAST}{Fast computation. see readData()}
  \item{big.data}{use the ff-package}
}

\details{
The algorithm reads the data for each individual and stores the information \cr
on disk. The data can be analyzed as regions of the defined window size, or can \cr
be concatenated in the PopGenome framework via the function concatenate.regions. \cr
This function should only be used when the FASTA file does not fit into the RAM; \cr
else, use the function \code{readData}. \cr
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
          3.\tab \code{region.names}            \tab names of regions \cr
          4.\tab \code{region.data}            	\tab some detailed information about the data\cr
        }


}


\examples{

# GENOME.class <- read.big.fasta("Alignment.fas", big.data=TRUE)
# GENOME.class
# GENOME.class@region.names
# CON <- concatenate.regions(GENOME.class)
# CON@region.data@biallelic.sites
# GENOME.class.slide <- sliding.window.transform(GENOME.class,100,100)
# GENOME.class <- neutrality.stats(GENOME.class,FAST=TRUE)
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}

