\name{weighted.jackknife-methods}

\docType{methods}

\alias{weighted.jackknife,GENOME-method}
\alias{weighted.jackknife-methods}


\title{Weighted Jackknife}

\description{
 Weighted Jackknife calculations.
}

\usage{

\S4method{weighted.jackknife}{GENOME}(object,
do.D=TRUE,
do.df=TRUE,
per.region=FALSE,
block.size=1)
}

 \arguments{

  \item{object}{An object of class \code{"GENOME"}}
  \item{do.D}{Pattersons D and Martin's f statistic}
  \item{do.df}{d-fraction (distance fraction)}
  \item{per.region}{jacknife within regions}
  \item{block.size}{the block size for jackknife (SNPs)}
}


\value{

Returned value is a modified object of class \code{"GENOME"} \cr  \cr
--------------------------------------------------------- \cr
The following slots will be modified in the \code{"GENOME"}  object \cr
--------------------------------------------------------- \cr

\tabular{rlll}{

            \tab Slot                	    \tab Reference   \tab Description \cr
 1.\tab \code{df.z}	            \tab [3]   \tab z values (jackknife) \cr
 2.\tab \code{df.SE}	            \tab [3]   \tab standard error (jackknife) \cr
 3.\tab \code{df.pval}	            \tab [3]   \tab P values (jackknife) \cr
 4.\tab \code{D.z}	            \tab [x]   \tab z values (jackknife) \cr
 5.\tab \code{D.SE}	            \tab [x]   \tab standard error (jackknife) \cr
 6.\tab \code{D.pval}	            \tab [x]   \tab P values (jackknife) \cr
 
          }

}

\details{
To perform the D and f statistic one needs to define 3 populations via the function \code{set.populations}, \cr
where the third population represent the archaic population. In addition, an outgroup is required and have to be
set via the function \code{set.outgroup}. Here, only SNPs where the outgroup is monomorphic are considered.
f is the fraction of the genome that is admixed [2].
}

\references{
[1] Durand, E. Y., Patterson, N. J., Reich, D., & Slatkin, M. (2011). \emph{Testing for ancient admixture between closely related populations.} \cr Molecular Biology and Evolution, 28(8), 2239--2252. doi:10.1093/molbev/msr048 \cr
[2] Simon H Martin, Kanchon K Dasmahapatra, Nicola J Nadeau, et al. (2013). \emph{Genome-wide evidence for speciation with gene flow in Heliconius butterflies.} \cr Genome Res. doi:10.1101/gr.159426.113 \cr
[3] Bastian Pfeifer and Durrell D. Kapan (2019). \emph{Estimates of introgression as a function of pairwise distances.} \cr BMC Bioinformatics. https://doi.org/10.1186/s12859-019-2747-z

}

\examples{
# GENOME.class <- readData("\home\Alignments")
# GENOME.class <- set.populations(GENOME.class,list(1:3,4:8,9:12))
# GENOME.class <- set.outgroup(GENOME.class,13)
# GENOME.class <- introgression.stats(GENOME.class, do.D=TRUE)
# show the result:
# GENOME.class@D
# GENOME.class@f
# GENOME.class <- introgression.stats(GENOME.class, do.df=TRUE)
# show the result:
# GENOME.class@df
}


\keyword{methods}
