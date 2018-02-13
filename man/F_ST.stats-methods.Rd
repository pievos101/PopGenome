\name{F_ST.stats-methods}

\docType{methods}
\alias{F_ST.stats,GENOME-method}
\alias{F_ST.stats-methods}

\alias{get.diversity,GENOME-method}
\alias{get.diversity-methods}

\alias{get.F_ST,GENOME-method}
\alias{get.F_ST-methods}

\title{Fixation Index}

\description{
 A generic function to calculate some F-statistics and nucleotide/haplotype diversities.
}

\usage{

\S4method{F_ST.stats}{GENOME}(
object,
new.populations=FALSE,
subsites=FALSE,
detail=TRUE,
mode="ALL",
only.haplotype.counts=FALSE,
FAST=FALSE
)

\S4method{get.diversity}{GENOME}(object,between=FALSE)
\S4method{get.F_ST}{GENOME}(object,mode=FALSE,pairwise=FALSE)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{subsites}{ 

                  \code{"transitions"}: SNPs that are transitions. \cr
                  \code{"transversions"}: SNPs that are transversions. \cr
                  \code{"syn"}: synonymous sites. \cr
                  \code{"nonsyn"}:  nonsynonymous sites. \cr
                  \code{"exon"}: SNPs in exon regions. \cr
                  \code{"intron"}: SNPs in intron regions. \cr
                  \code{"coding"}: SNPs in coding regions (CDS). \cr
                  \code{"utr"}:  SNPs in UTR regions. \cr
                  \code{"gene"}: SNPs in genes. \cr
                  \code{"intergenic"} : SNPs in intergenic regions. \cr
  } 
  \item{detail}{detail statistics. Note: slower!}
  \item{between}{\code{TRUE}: show between-diversities. \code{FALSE}: show within-diversities}
  \item{mode}{mode=\code{"haplotype"} or mode=\code{"nucleotide"}}
  \item{only.haplotype.counts}{only calculate the haplotype counts}
  \item{FAST}{if \code{TRUE} only calculate a subset of statistics. see details!}
  \item{pairwise}{show paiwise comparisons. default:\code{FALSE}}
}


\value{

\tabular{rlll}{
            \tab Slot                	        \tab Reference   \tab Description\cr
          1.\tab \code{haplotype.F_ST} 	        \tab [1]         \tab Fixation Index based on haplotype      frequencies\cr
          2.\tab \code{nucleotide.F_ST}         \tab [1]         \tab Fixation Index based on minor.allele   frequencies\cr
          3.\tab \code{Nei.G_ST}             	\tab [2]         \tab Nei's Fixation Index\cr
          4.\tab \code{Hudson.G_ST}             \tab [3]         \tab see reference ... \cr
          5.\tab \code{Hudson.H_ST}             \tab [3]         \tab see reference ... \cr
          6.\tab \code{Hudson.K_ST}             \tab [3]         \tab see reference ... \cr
          7.\tab \code{nuc.diversity.within}    \tab [1,5]       \tab Nucleotide diversity (within the population)\cr
          8.\tab \code{hap.diversity.within}    \tab [1]         \tab Haplotype diversity (within the population) \cr
          9.\tab \code{Pi}          	        \tab [4]         \tab Nei's diversity (within the population)  \cr
          10.\tab \code{hap.F_ST.vs.all}        \tab [1]         \tab Fixation Index for each population against all other individuals (haplotype)\cr
          11.\tab \code{nuc.F_ST.vs.all}        \tab [1]         \tab Fixation Index for each population against tall other individuals (nucleotide)\cr
          12.\tab \code{hap.diversity.between}  \tab [1]         \tab Haplotype diversities between populations\cr
          13.\tab \code{nuc.diversity.between}  \tab [1,5]       \tab Nucleotide diversities between populations\cr
          14.\tab \code{nuc.F_ST.pairwise}      \tab [1]         \tab Fixation Index for every pair of populations (nucleotide)\cr
          15.\tab \code{hap.F_ST.pairwise}      \tab [1]         \tab Fixation Index for every pair of populations (haplotype)\cr
          16.\tab \code{Nei.G_ST.pairwise}      \tab [2]         \tab Fixation Index for every pair of populations (Nei)\cr
          17.\tab \code{region.stats}           \tab             \tab an object of class "region.stats" for detailed statistics\cr
       }

}

\seealso{
# methods?F_ST.stats.2  
#\link{F_ST.stats.2}
}

\details{

If \code{FAST} is switched on, this module only calculates \code{nuc.diversity.within}, \code{hap.diversity.within}, \code{haplotype.F_ST}, \code{nucleotide.F_ST} and \code{pi}.\cr
Note: \cr 
1) The nucleotide diversities have to be devided by the size of region considered (e.g. \code{GENOME@n.sites}) to give diversities per site. \cr
2) When missing or unknown nucleotides are included (include.unknown=TRUE) those sites are completely deleted in case of haplotype based statistics. \cr
3) The function \code{detail.stats(...,site.FST=TRUE)} will compute SNP specific FST values which are then 
stored in the slot \code{GENOME.class@region.stats@site.FST}.\cr
4) We recommend to use \code{mode="nucleotide"} in case you have many unknowns included in your dataset.
}

\references{

[1] Hudson, R. R., M. Slatkin, and W.P. Maddison (1992). \emph{Estimating levels of gene flow from DNA sequence data.} Gentics 13(2),583-589 \cr \cr
[2] Nei, M. (1973). \emph{Analysis of gene diversity in subdivided populations.} Proc.Natl. Acad. Sci. USA 70: 3321-3323 \cr \cr
[3] Hudson, R. R., Boos, D.D. and N. L. Kaplan (1992). \emph{A statistical test for detecting population subdivison.} Mol. Biol. Evol. 9: 138-151. \cr \cr
[4] Nei, M. (1987). \emph{Molecular Evolutionary Genetics.} Columbia Univ. Press, New York. \cr \cr
[5] Wakeley, J. (1996).\emph{The Variance of Pairwise Nucleotide Differences in Two Populations with Migration.}
THEORETICAL POPULATION BIOLOGY. 49, 39-57. \cr 

}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- F_ST.stats(GENOME.class)
# GENOME.class <- F_ST.stats(GENOME.class,list(1:4,5:10),subsites="syn")
# GENOME.class <- F_ST.stats(GENOME.class,list(c("seq1","seq5","seq3"),c("seq2","seq8")))
# show the result:
# get.F_ST(GENOME.class)
# get.F_ST(GENOME.class, pairwise=TRUE)
# get.diversity(GENOME.class, between=TRUE)
# GENOME.class@Pi --> population specific view
# GENOME.class@region.stats

}


\keyword{methods}

