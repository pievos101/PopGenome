\name{VCF_split_into_scaffolds}

\docType{methods}
\alias{VCF_split_into_scaffolds}

\title{Split a VCF file into multiple scaffold-VCFs}


\description{
 This function splits a VCF file into multiple VCFs including data for exactly one scaffold each.
}

\usage{
VCF_split_into_scaffolds(VCF.file, output.folder)
}

\arguments{ 
  \item{VCF.file}{the basepath of the VCF file}
  \item{output.folder}{name of the folder where the VCFs should be stored}
}

\details{
The algorithm splits the VCF into multiple scaffold based VCFs and \cr 
stores the files in a given folder. This folder can be used as an input for readData(,format="VCF") 
}

\value{
\code{TRUE}
}


\examples{
# VCF_split_into_scaffolds("VCFfile.vcf","scaffoldVCFs")
# test <- readData("scaffoldVCFs", format="VCF")
}

