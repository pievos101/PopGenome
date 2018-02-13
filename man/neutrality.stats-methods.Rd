\name{neutrality.stats-methods}
\docType{methods}
\alias{neutrality.stats,GENOME-method}
\alias{neutrality.stats-methods}

\alias{get.neutrality,GENOME-method}
\alias{get.neutrality-methods}



\title{Neutrality Statistics}

\description{
 This generic function calculates some neutrality statistics.
}


\usage{

\S4method{neutrality.stats}{GENOME}(object,new.populations=FALSE,new.outgroup=FALSE,
subsites=FALSE,detail=FALSE, FAST=FALSE, do.R2=FALSE)
\S4method{get.neutrality}{GENOME}(object,theta=FALSE,stats=TRUE)

}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{new.populations}{list of populations. default:\code{FALSE}}
  \item{new.outgroup}{vector of outgroup sequences. default:\code{FALSE}}
  \item{subsites}{ 

                  \code{"transitions"}: SNPs that are transitions. \cr
                  \code{"transversions"}: SNPs that are transversions. \cr
                  \code{"syn"}: synonymous sites. \cr
                  \code{"nonsyn"}:  non-synonymous sites. \cr
                  \code{"exon"}: SNPs in exon regions. \cr
                  \code{"intron"}: SNPs in intron regions. \cr
                  \code{"coding"}: SNPs in coding regions (CDS). \cr
                  \code{"utr"}:  SNPs in UTR regions. \cr
                  \code{"gene"}: SNPs in genes. \cr
          default:\code{FALSE}}

  
  \item{detail}{default:\code{FALSE}, \code{TRUE} for some detailed statistics. Note: slows down calculations!}  
  \item{FAST}{Fast computation. only works if there is no outgroup defined.}
  \item{do.R2}{Ramos-Onsins' & Rozas' R2}
  \item{stats}{show the results of each statistic. default:\code{TRUE}}
  \item{theta}{show the theta values. default:\code{FALSE}}
}

\value{


  The return value is a modified object of class \code{"GENOME"} \cr \cr
  --------------------------------------------------------- \cr
  The following slots will be modified in the \code{"GENOME"} object \cr
  --------------------------------------------------------- \cr

  \tabular{rlll}{
            \tab Slot                	\tab Reference       \tab Description\cr
          1.\tab \code{n.segregating.sites} 	\tab             \tab Total number of segregating sites\cr
          2.\tab \code{Tajima.D}            	\tab [1]         \tab Tajima's' D statistic 1989\cr
          3.\tab \code{Fu.Li.F}             	\tab [3]         \tab Fu & Li's' F* statistic 1993\cr
          4.\tab \code{Fu.Li.D}             	\tab [3]         \tab Fu & Li's D* statistic 1993\cr
          5.\tab \code{Fay.Wu.H}            	\tab [6]         \tab Fay & Wu's H statistic 2000\cr
          6.\tab \code{Zeng.E}              	\tab [7]         \tab Zeng's E statistic 2006\cr
          7.\tab \code{Strobeck.S}          	\tab [5]         \tab Strobeck's S statistic 1987 (if detail==\code{TRUE})\cr
          8.\tab \code{Fu.F_S}              	\tab [4]         \tab Fu's F$_{S}$ statistic 1997 (if detail==\code{TRUE})   \cr
          9.\tab \code{Rozas.R_2}          	\tab [2]         \tab Ramos-Onsins' & Rozas' $R_{2}$ statistic 2002\cr
          10.\tab \code{theta_Tajima}       	\tab [1]        \tab \cr
          11.\tab \code{theta_Watterson}       \tab          \tab \cr
          12.\tab \code{theta_Fu.Li}           \tab [3]         \tab \cr
          13.\tab \code{theta_Achaz.Watterson}      \tab          \tab \cr
          14.\tab \code{theta_Achaz.Tajima}    \tab          \tab \cr
          15.\tab \code{theta_Fay.Wu}          \tab [6]         \tab \cr
          16.\tab \code{theta_Zeng}            \tab [7]         \tab \cr
       }


  
}

\references{ 

[1]  Tajima, F.(1989) \emph{Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism.} Genetics, 123(3): 585-595.\cr \cr
[2]  Ramos-Onsins, S.E. and J.Rozas (2002). \emph{Statistical Properties of New Neutrality Tests Against Population Growth.} Mol.Biol.Evol.19(12),2092-2100 \cr \cr
[3]  Fu, Y.X. and W.H.Li (1993). \emph{Statistical Tests of Neutrality of Mutations.} Genetics 133(3),693-709 \cr \cr
[4]  Fu, Y.-X.(1997). \emph{Statistical Tests of Neutrality of mutations against population growth, hitchhiking and background selection.} Genetics 147(2),915-925.\cr \cr
[5]  Strobeck, C. (1987). \emph{Average number of nucleotide differences in a sample from a single subpopulation: a test for population subdivision.} Genetics 117, 149-153 \cr \cr
[6]  Fay, J.C. and C.-I. Wu (2000). \emph{Hitchhiking under positive Darwinian selection.} Genetics 155 (3),1405-1413 \cr \cr
[7]  Zeng, K., Y.-X. Fu, S. Shi, and C.-I. Wu (2006). \emph{Statistical tests for detecting positive selection by utilizing high-frequency variants.} Genetics 174, 1431-1439 \cr \cr

}


\examples{

# GENOME.class <- readData("\home\Alignments")
# GENOME.class
# GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE)
# GENOME.class <- neutrality.stats(GENOME.class,list(1:4,5:10),subsites="syn")
# GENOME.class <- neutrality.stats(GENOME.class,list(c("seq1","seq5","seq3"),
# c("seq2","seq8")))
# GENOME.class <- neutrality.stats(GENOME.class,detail=TRUE)
# show the result:
# get.neutrality(GENOME.class)
# GENOME.class@Tajima.D --> population specific view
# detail = TRUE
# GENOME.class@region.stats

}

\keyword{methods}
