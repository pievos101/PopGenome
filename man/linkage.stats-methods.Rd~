\name{linkage.stats-methods}

\docType{methods}

\alias{linkage.stats,GENOME-method}
\alias{linkage.stats-methods}

\alias{get.linkage,GENOME-method}
\alias{get.linkage-methods}

\title{Linkage Disequilibrium}


\description{
 A generic function to calculate some linkage disequilibrium statistics.
}

\usage{

\S4method{linkage.stats}{GENOME}(object,new.populations=FALSE,subsites=FALSE,detail=FALSE,
			do.ZnS,do.WALL=TRUE)
\S4method{get.linkage}{GENOME}(object)

}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default=\code{FALSE}}
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
          default:\code{FALSE}}
 
  \item{detail}{if you want to calculate some detailed statistics. This can be considerably slower! default:\code{FALSE}}  
  \item{do.ZnS}{calculate ZnS, ZA and ZZ} 
  \item{do.WALL}{calculate Wall's B/Q}
}

\details{
Note, the pairwise comparisons are computed via \code{combn(n.snps,2)}.
}

\value{

The return value is a modified object of class \code{"GENOME"} \cr  \cr
--------------------------------------------------------- \cr 
The following slots will be modified in the \code{"GENOME"}  object \cr 
--------------------------------------------------------- \cr 

\tabular{rlll}{
            \tab Slot                	\tab Reference   \tab Description\cr
          1.\tab \code{Wall.B}	                \tab [2]         \tab Wall $B$ statistic (only adjacent positions are considered)\cr
          2.\tab \code{Wall.Q}            	\tab [2]         \tab Wall $Q$ statistic (only adjacent positions are considered)\cr
          3.\tab \code{Kelly.Z_nS}              \tab [3]         \tab Kelly $Z_{n}S$ statistic (if detail==\code{TRUE})\cr      
          4.\tab \code{Rozas.ZA}            	\tab [1]         \tab Rozas $ZA$ statistic (adjacent positions, if detail==\code{TRUE})\cr
	  5.\tab \code{Rozas.ZZ}             	\tab [1]         \tab Rozas $ZZ$ statistic ($ZZ=ZA-Z_{n}S$, if detail==\code{TRUE})\cr
        }
  
}

\references{ 
[1] Rozas, J., M.Gullaud, G.Blandin, and M.Aguade(2001). \emph{DNA variation at the rp49 gene region of Drosophila simulans: evolutionary inferences from an unusual haplotype structure.} Genetics 158(3),1147-1155 \cr \cr
[2] Wall, J.(1999). \emph{Recombination and the power of statistical tests of neutrality.} Genet Res 74, 65-79 \cr \cr
[3] Kelly,J.K. (1997). \emph{A test of neutrality based on interlocus associations.} Genetics 146: 1197-1206 \cr \cr
}

\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- linkage.stats(GENOME.class)
# GENOME.class <- linkage.stats(GENOME.class,list(1:4,5:10),subsites="syn")
# GENOME.class <- linkage.stats(GENOME.class,list(c("seq1","seq5","seq3"),
# c("seq2","seq8")))
# GENOME.class <- linkage.stats(GENOME.class, detail=TRUE)
# show the result:
# get.linkage(GENOME.class)
# GENOME.class@Wall.B --> population specific view
# GENOME.class@region.stats

}


\keyword{methods}
