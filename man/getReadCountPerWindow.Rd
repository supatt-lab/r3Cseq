\name{getReadCountPerWindow}
\alias{getReadCountPerWindow}
\alias{getReadCountPerWindow,r3Cseq-method}
\title{count reads per window size}

\description{
  Counts the number of reads from 3C-Seq data per each window size
}

\usage{
getReadCountPerWindow(object,windowSize=5e3,nFragmentExcludedReadsNearViewpoint=2,mode=c("non-overlapping", "overlapping"))
}


\arguments{

	\item{object}{
		r3Cseq object
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
}

\author{
	S. Thongjuea
}

\examples{
  #See the vignette	
}
%\keyword{read counts manipulation}
