\name{test.params-class}
\docType{class}
\alias{test.params}
\alias{test.params-class}

\title{

Set parameters for coalescent simulations with Hudson's MS and Ewing's MSMS.

}

\description{

The object that contains the set parameter values can be passed to the function MS. 
This class simpliifes the process of passing on all necessary values
to the MS function.
 
}

\arguments{
 
  \item{theta}{mutation parameter theta (4Nmu), where N is the diplod effective population size and mu the mutation rate per locus. It needs to be provided as a vector of length n.regions}  
  \item{seeds}{specify 3 random number seeds. a vector of length 3 with positive integer values is expected}
  \item{fixedSegsites}{usually the number of segregating sites varies in each iteration. Please provide a  single numeric value if the number of segregating sites needs to be fixed.}
  \item{recombination}{provide a vector of format: c(p, nsites), p = cross-over parameter rate, nsites is the number of sites between which recombination occurs}
  \item{geneConv}{in addition to recombination, intra-locus non-cross-over exchange gene conversion can be included in the simulation; the expected format is c(f, gamma),  where f denotes the ratio g/r (r is the probability per generation  of crossing-over between adjacent sites (see Wiuf and Hein 2000), and gamma is the mean conversion tract length.}
  \item{growth}{population size is assumed to be $N(t) = N0 exp^{alpha*t}$. Provide alpha as an integer value. Negative values indicate that population was larger in the past than present, positive values indicate that it was smaller.}
  \item{migration}{specify the migration rate between populations. Please provide a single numeric value.}
  \item{demography}{vector of length 3 or 4 with first value denoted as 'type' \cr
					valid 'types' for vectors of length 3 are as follows: \cr
					- \code{1}  set a growth rate change alpha at a certain time t:\cr 
					  		c(1, t, alpha)\cr
\cr
					- \code{2}	set all sub-populations to size $x * N_{0}$ and growth rate to zero:\cr 
					  		c(2, t, x)\cr
\cr
					- \code{3}      set all elements of the migration matrix to $x/(npop-1)$:\cr 
					 		c(3, t, x)					\cr
\cr
							valid 'types' for vectors of length 4 are as follows:\cr	
					- \code{4}      set growth rate of sub-population i to alpha at time z: \cr
					 		 c(4, t, i, alpha)\cr
\cr
					- \code{5}      set sub-population i size to $x * N_{0}$ at time t and growth rate to zero:\cr 
					  		c(5, t, i, x)\cr
\cr
					- \code{6}      split sub-population i into sub-population i and a new sub-population,\cr 
					  		labeled npop + 1. Each ancestral lineage in sub-population i is randomly \cr 
					  		assigned to sub-population i with probability p and sub-population  \cr
					 		npop + 1 with probability 1 - p. The size of sub-population npop + 1 is\cr  
					  		set to $N_{0}$. Migration rates to and from the new sub-population are assumed\cr  
					  		to be zero and the growth rate of the new sub-population is set to zero: \cr
					 		c(6, t, i, p)\cr
\cr
					- \code{7}      move all lineages in sub-population i to sub-population j at time t.\cr 
					  		Migration rates from sub-population i are set to zero: \cr
					  		c(7, t, i, j)\cr
  }

}


\author{
Bastian Pfeifer
}

\seealso{
\link{MS}
}
\examples{

# params           <- new("test.params")
# params@theta     <- rep(5,n.regions)
# params@migration <- 3

}

\keyword{classes}

