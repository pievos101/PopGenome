\name{get_gff_info}

\docType{methods}
\alias{get_gff_info}


\title{Annotation info}


\description{
This function extracts annotation information from a GTF/GFF file.
}

\usage{

get_gff_info(object=FALSE,gff.file,chr,position,feature=FALSE,extract.gene.names=FALSE)

}

\arguments{

  \item{object}{object of class \code{GENOME}} 
  \item{gff.file}{basepath of the GTF/GFF file}
  \item{chr}{the chromosome}
  \item{position}{reference positions or region id (when \code{object} is specified)}
  \item{feature}{feature to search for in the gff-file. returns a list of positions}
  \item{extract.gene.names}{returns the gene names of the chromosome}
}  

\details{
This function extracts annotation information from a GTF/GFF file.
}


\examples{
# get_gff_info("Arabidopsis.gff",chr=1,200202)
# get_gff_info(GENOME.class,"Arabidopsis.gff",chr=1,position=3)
}



\keyword{methods}
