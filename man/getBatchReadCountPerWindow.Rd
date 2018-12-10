\name{getBatchReadCountPerWindow}
\alias{getBatchReadCountPerWindow}
\alias{getBatchReadCountPerWindow,r3CseqInBatch-method}
\title{count reads per window size for replicates analysis}

\description{
  Counts the number of reads from 3C-Seq data per each window size for replicates analysis
}

\usage{
getBatchReadCountPerWindow(object,windowSize=5e3,nFragmentExcludedReadsNearViewpoint=2,mode=c("non-overlapping", "overlapping"))
}


\arguments{

	\item{object}{
		r3CseqInBatch object
	}
	\item{windowSize}{
		Numeric. non-overlapping window size for counting reads (default=5e3)
	}
	\item{nFragmentExcludedReadsNearViewpoint}{
		Numeric. The number of excluded fragments around the viewpoint, reads found in these fragments
		will be removed from the analysis (default=2)
	}
	\item{mode}{
		character. The window-based modes analysis (default="non-overlapping")
	}
}

\value{
   The RangedData represents the number of reads per each window size
}

\seealso{
  \link{getReadCountPerRestrictionFragment},
  \link{getBatchReadCountPerRestrictionFragment},
  \link{getReadCountPerWindow},
}

\author{
	S. Thongjuea
}

\examples{
  #See the vignette	
}
%\keyword{read counts manipulation}
