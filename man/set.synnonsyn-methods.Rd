\name{set.synnonsyn-methods}
\docType{methods}
\alias{set.synnonsyn,GENOME-method}
\alias{set.synnonsyn-methods}


\title{Set synonymous positions for SNP data}


\description{
 This generic function classifies the observed biallelic positions read from SNP data files into synonymous and non-synonymous SNPs.\cr
}

\usage{

\S4method{set.synnonsyn}{GENOME}(object,ref.chr,save.codons=FALSE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{ref.chr}{the reference chromosome in FASTA format}
  \item{save.codons}{save codon changes}
}

\value{
The return value is a modified object of class \code{"GENOME"} storing \cr
syn/nonsyn informations in the slot \code{GENOME.class@region.data@synonymous}\cr
for each SNP. (1=synonymous,0=non-synonymous) \cr 
When \code{save.codons} is \code{TRUE} the SNP related codon changes \cr
are saved in the corresponding slot \code{GENOME.class@region.data@codons}. \cr 
(see also get.codons(), codontable() and codonise64())  
}

\note{
The data has to be read in with a corresponding GFF/GTF file (CDS fields must be specified); 
otherwise a correct classification is not possible. The set.synnonsyn() function does 
not work for splitted objects e.g produced via sliding.window.transform() or splitting.data(). 
Note, transcripts which are in the same CDS region but have different reading frames are \cr
are not specified correctly. PopGenome can also handle coding regions on reverse strands.
We have used the program SNPeff to validate our results.
}

\examples{
# GENOME.class <- readData("VCF",format="VCF",gffpath="GFF.Folder")
# GENOME.class <- set.synnonsyn(GENOME.class,ref.chr="ref.fas")
# GENOME.class@region.data@synonymous 

}



\keyword{methods}
