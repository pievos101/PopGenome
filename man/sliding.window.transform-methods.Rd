\name{sliding.window.transform-methods}
\docType{methods}

\alias{sliding.window.transform,GENOME-method}
\alias{sliding.window.transform-methods}

\title{Sliding Window Transformation}


\description{
 This generic function transforms an existing object of class \code{"GENOME"} into another object of class \code{"GENOME"}, 
 in which each region corresponds to one window. This allows to apply the full spectrum of PopGenome methods to 
 sliding window data.   
 
}

\usage{

\S4method{sliding.window.transform}{GENOME}(object, 
	width=7, jump=5, 
	type=1,
	start.pos=FALSE,end.pos=FALSE,
	whole.data=TRUE
	)
}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{width}{window size. default:\code{7}}
  \item{jump}{jump size. default:\code{5}}
  \item{type}{\code{1} scan only biallelic positions (SNPs), \code{2} scan the genome. default:\code{1}}
  \item{start.pos}{start position}
  \item{end.pos}{end position}
  \item{whole.data}{scan the complete data by concatenating the regions in "object". If FALSE, each region is scanned seperately.}
}

\value{
The function creates a transformed object of class \code{"GENOME"}.
}

\note{
If you want to scan regions seperately (whole.data=FALSE), you may not use the big.data option in the readData function.
PopGenome will scan the data from position 1 to the last observed SNP if start or end-positions are not specified. 
}


\examples{
# GENOME.class       <- readData("...\Alignments")
# slide.GENOME.class <- sliding.window.transform(GENOME.class)
# slide.GENOME.class <- sliding.window.transform(GENOME.class,100,100)
# slide.GENOME.class <- neutrality.stats(slide.GENOME.class)
# slide.GENOME.class@region.names
# values             <- get.neutrality(slide.GENOME.class)
# GENOME.class       <- readSNP("Arabidopsis", CHR=1)
# GENOME.slide       <- sliding.window.transform(GENOME.split, 10000, 10000, type=2,
# start.pos=10000000, end.pos=12000000)
# GENOME.slide@region.names
}



\keyword{methods}
