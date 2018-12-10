\name{plotOverviewInteractions}
\alias{plotOverviewInteractions}
\alias{plotOverviewInteractions,r3Cseq-method}


\title{Plot overview of identified interaction regions for genome-wide}

\description{
  Plot the distribution of identified interaction regions across genome
}

\usage{
plotOverviewInteractions(obj, cutoff.qvalue=0.05)
}


\arguments{

	\item{obj}{
		obj is r3Cseq or r3CseqInBatch object
	}

	\item{cutoff.qvalue}{
		Numeric. The cutoff q-value (default=0.05)
	}
}

\value{
	Plots of identified 3C-Seq interaction regions genome-wide
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotInteractionsNearViewpoint}}, \code{\link{plotInteractionsPerChromosome}}, 
	\code{\link{plotDomainogramNearViewpoint}}
}


\examples{
 # See the vignette
}
%\keyword{plot interaction regions}
