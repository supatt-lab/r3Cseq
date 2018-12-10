\name{plotInteractionsNearViewpoint}
\alias{plotInteractionsNearViewpoint}
\alias{plotInteractionsNearViewpoint,r3Cseq-method}

\title{Plot identified interaction regions near the viewpoint}

\description{
  Plot identified interaction regions near the viewpoint 
}

\usage{
plotInteractionsNearViewpoint(obj,distance=5e5,log2fc_cutoff=1,yLim=0)
}


\arguments{

	\item{obj}{
		obj is r3Cseq or r3CseqInBatch object
	}
	\item{distance}{
		Numeric. The distance relative to the viewpoint (default=5e5)
	}
	\item{log2fc_cutoff}{
		Numeric. The log2 cutoff ratio between the experiment and control (default=1)
	}
	\item{yLim}{
		Numeric. The limited height of y-axis (default=0)
	}
}


\value{
	Plots of identified interaction regions close to the viewpoint
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{plotOverviewInteractions}}, \code{\link{plotInteractionsPerChromosome}}, 
	\code{\link{plotDomainogramNearViewpoint}}
}


\examples{

	# See the vignette
}
%\keyword{plot interaction regions}
