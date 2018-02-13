\name{GFF_split_into_scaffolds}

\docType{methods}
\alias{GFF_split_into_scaffolds}

\title{Split a GFF file into multiple scaffold-GFFs}


\description{
 This function splits a GFF file into multiple GFFs including data for exactly one scaffold each.
}

\usage{
GFF_split_into_scaffolds(GFF.file, output.folder)
}

\arguments{ 
  \item{GFF.file}{the basepath of the GFF file}
  \item{output.folder}{name of the folder where the GFFs should be stored}
}

\details{
The algorithm splits the GFF into multiple scaffold based GFFs and \cr 
stores the files in a given folder. This folder can be used as an input for readData(,gffpath="") 
}

\value{
\code{TRUE}
}


\examples{
# GFF_split_into_scaffolds("GFFfile.gff","scaffoldGFFs")
# test <- readData("scaffoldVCFs", format="VCF", gffpath="scaffoldGFFs")
}

