\name{PopGenome}
\alias{PopGenome}
\alias{PopGenome}
\docType{package}
\title{

 PopGenome



}
\description{

R-package for Population genetic & genomic analyses

}
\details{

Index:
\preformatted{
F_ST.stats               Fixation index
diversity.stats          Diversities
MKT                      McDonald & Kreitman test
MS                       Coalescent simulations
detail.stats             Several misc. statistics
linkage.stats            Linkage disequilibrium
neutrality.stats         Neutrality statistics
readData                 Reading alignments and calculating summary data
readSNP                  Read data in .SNP format (e.g., from the 1001 Arabidopsis Genomes project)
readVCF                  Read data in VCF format (e.g., from the 1000 human Genomes project)
readHapMap               Read data in HapMap format
sliding.window.transform Sliding window transformation
splitting.data           Split data into subsites
test.params              Set parameters for coalescent simulations.

}

}

\author{
Bastian Pfeifer
Maintainer: Bastian Pfeifer <Bastian.Pfeifer@uni-duesseldorf.de>

}

\keyword{ package }
\seealso{

?readData
\link{readData}

}
\examples{

# GENOME.class <- readData("...\Alignments")
# GENOME.class <- neutrality.stats(GENOME.class)
# values       <- get.neutrality(GENOME.class)
# GENOME.class <- F_ST.stats(GENOME.class,list(1:5,6:10))
# values       <- get.F_ST(GENOME.class)

}
