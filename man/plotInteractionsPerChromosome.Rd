\name{plotInteractionsPerChromosome}
\alias{plotInteractionsPerChromosome}
\alias{plotInteractionsPerChromosome,r3Cseq-method}

\title{Plot interaction regions per each chromosome of interest}

\description{
  Plot the distribution of interaction regions per each chromosome
}

\usage{
plotInteractionsPerChromosome(obj, chromosomeName)
}


\arguments{

	\item{obj}{
		obj is r3Cseq or r3CseqInBatch object. 
	}

	\item{chromosomeName}{
		Character. The input chromosome name (e.g. "chr1")
	}
}


\value{
	Plots of interaction regions per chromosome.
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotInteractionsNearViewpoint}}, \code{\link{plotOverviewInteractions}}, 
	\code{\link{plotDomainogramNearViewpoint}}
}


\examples{

# See the vignette

}
%\keyword{plot interaction regions}
