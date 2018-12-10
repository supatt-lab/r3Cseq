\name{getReadCountPerRestrictionFragment}
\alias{getReadCountPerRestrictionFragment}
\alias{getReadCountPerRestrictionFragment,r3Cseq-method}
\title{count reads per resitrcition fragment}

\description{
  Counts the number of reads from 3C-Seq data per each restriction fragment
}

\usage{
getReadCountPerRestrictionFragment(object,getReadsMethod = c("wholeReads", "adjacentFragmentEndsReads"),
			nFragmentExcludedReadsNearViewpoint=2)
}


\arguments{

	\item{object}{
		r3Cseq object
	}
	\item{getReadsMethod}{
		character. To count all reads found in the particular restriction fragment uses wholeReads option.  
		To count reads found around the edge of restriction fragment both 5'utr and 3'utr uses adjacentFragmentEndsReads
		option (default=wholeReads)
	}
	\item{nFragmentExcludedReadsNearViewpoint}{
		Numeric. The number of excluded fragments around the viewpoint, reads found in these fragments
		will be removed from the analysis (default=2)
	}
}

\value{
   The RangedData represents the number of reads per each restriction fragment
}

\seealso{
  \link{getReadCountPerWindow},
  \link{getBatchReadCountPerRestrictionFragment}
}

\author{
	S. Thongjuea
}

\examples{
  #See the vignette	
}
%\keyword{read counts manipulation}
