\name{PopGplot}

\docType{methods}
\alias{PopGplot}


\title{Smoothed line-plot for multiple populations}


\description{
This function plots values with smoothed lines using spline interpolation.
}

\usage{

PopGplot(values,colors=FALSE,span=0.1,ylab="",xlab="",
	ylim=c(min(values,na.rm=TRUE),max(values,na.rm=TRUE)))

}

\arguments{
\item{values}{the statistical values (matrix); columns=populations}
\item{colors}{the colors for each population (character vector)}
\item{span}{the degree of smoothing} 
\item{ylab}{a title for the y axis} 
\item{xlab}{a title for the x axis} 
\item{ylim}{ranges for the y axis}
}  

\examples{

# GENOME.class <- readSNP("Arabidopsis",CHR=1)
# GENOME.class.slide <- sliding.window.transform(test,1000,1000)
# GENOME.class.slide <- diversity.stats(GENOME.class.slide)
# values <- GENOME.class.slide@nuc.divesity.within
# PopGplot(values)

}



\keyword{methods}
