\name{readHapMap}

\docType{methods}
\alias{readHapMap}

\title{Read SNP data from the HapMap consortium}


\description{
 This function reads HapMap data. 
}

\usage{

readHapMap(folder,hap_gffpath,populations=FALSE,outgroup=FALSE)

}

\arguments{
 
  \item{folder}{the basepath of the variant_calls}
  \item{hap_gffpath}{the basepath of the corresponding GFF files. Note! The HapMap GFF file does not contain information about subsites. see details!}
  \item{populations}{list of populations}
  \item{outgroup}{vector of outgroup sequences}
  
}

\details{
PopGenome reads the GFF file distributed on the HapMap plattform only to verify the reference positions of the chromosomes.
In the next release, this function will also handle GFF/GTF files to get information about subsites (exons, introns, ...).
The input folder should include the files of different individuals for one chromosome. This facilitates FST calculations of the HapMap data. \cr
See also readData("...", format="HapMap") which can read the files of single populations directly.
}

\value{


  The function creates an object of class \code{"GENOME"} \cr \cr
  --------------------------------------------------------- \cr
  The following slots will be filled in the \code{"GENOME"} object \cr
  --------------------------------------------------------- \cr

            \tabular{rll}{
            \tab Slot                	        \tab Description\cr
          1.\tab \code{n.sites}                        \tab total number of sites\cr
          2.\tab \code{n.biallelic.sites}          	\tab number of biallelic sites \cr
          3.\tab \code{region.data}            	\tab some detailed information about the data read\cr
        }


}


\examples{

# GENOME.class <- readHapMap("...\HapMapData")
# GENOME.class
# show the result:
# get.sum.data(GENOME.class)
# GENOME.class@region.data

}



\keyword{methods}
