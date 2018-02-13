\name{vcf_file}
\docType{data}
\alias{vcf_file}
\title{VCF file (subdirectory "data")}
\description{
A VCF file stored in the directory "data" of the PopGenome package.
The file (unpacked) has to be stored in a folder (for example "VCF").
Note that many VCF-files can be stored in this folder and are read consecutively. If the 
VCF file is too large to fit into the computer's main memory, split it into chunks (by position) !
PopGenome is able to concatenate these chunks afterwards.
\code{readData("VCF",format="VCF",FAST=TRUE)}
}
\keyword{datasets}
