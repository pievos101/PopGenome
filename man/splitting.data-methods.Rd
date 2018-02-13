\name{splitting.data-methods}
\docType{methods}
\alias{splitting.data,GENOME-method}
\alias{splitting.data-methods}

\title{Split data into subsites}


\description{
 This generic function splits the data into subsites, if GFF/GTF information is present or if
 positions are defined accordingly.
 
}

\usage{

\S4method{splitting.data}{GENOME}(object,subsites=FALSE,positions=FALSE,type=1,
          whole.data=TRUE)


}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{positions}{list of positions}
 
 
  \item{subsites}{ 
                  \code{"exon"}: SNPs in exon regions. \cr
                  \code{"intron"}: SNPs in intron regions. \cr
                  \code{"coding"}: SNPs in coding regions (CDS). \cr
                  \code{"utr"}:  SNPs in UTR regions. \cr
                  \code{"gene"}: SNPs in genes. \cr
  }  
 \item{type}{\code{1}: SNP positions \cr
              \code{2}: Genome positions   
  }
  
 \item{whole.data}{Scan the whole data by concatenating the regions. 
                   If FALSE, the regions are scanned seperately}
}

\details{
Note, if \code{whole.data=FALSE} data with \code{n.biallelic.sites==0} should be removed.
}


\value{


The return value is a modified object of class \code{"GENOME"}. \cr

  
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class.split <- splitting.data(GENOME.class,subsites="exon")
# GENOME.class.split@region.names
# GENOME.class.split <- splitting.data(GENOME.class,positions=list(1:7,8:12))
# GENOME.class.split <- splitting.data(GENOME.class,
# positions=list(2000:3000,12000:13000),type=2)
# GENOME.class.split 

}



\keyword{methods}
