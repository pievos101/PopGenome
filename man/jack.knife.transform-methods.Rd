\name{jack.knife.transform}
\docType{methods}

\alias{jack.knife.transform,GENOME-method}
\alias{jack.knife.transform-methods}

\title{Jacknife Transformation}


\description{
 This generic function transforms an existing object of class \code{"GENOME"} into another object of class \code{"GENOME"}, 
 in which each region corresponds to the (JACKNIFE !) window. Each jacknife window will be excluded from the analyses and the 
 calculation will be applied to the union of all other windows.   
}

\usage{

\S4method{jack.knife.transform}{GENOME}(object, 
	width=7, jump=5, 
	type=1,
	start.pos=FALSE,end.pos=FALSE
	)
}

\arguments{

  \item{object}{an object of class \code{"GENOME"}}
  \item{width}{window size. default:\code{7}}
  \item{jump}{jump size. default:\code{5}}
  \item{type}{\code{1} scan only biallelic positions (SNPs), \code{2} scan the genome. default:\code{1}}
  \item{start.pos}{start position}
  \item{end.pos}{end position}
}

\value{
The function creates a transformed object of class \code{"GENOME"}.
}

\note{
This function currently is only available for SNP data formats.
PopGenome will scan the data from position 1 to the last observed SNP if start or end-positions are not specified. 
This mechanism can also be applied to the splitting.data() function. Just
set split.GENOME.class@jack.knife <- TRUE after splitting the data.
}


\examples{
# GENOME.class      <- readData("...", format="VCF")
# jack.GENOME.class <- jack.knife.transform(GENOME.class,100,100)
# jack.GENOME.class <- neutrality.stats(jack.GENOME.class)
# jack.GENOME.class@Tajima.D
}



\keyword{methods}
