\name{getRawReads}
\alias{getRawReads}
\alias{getRawReads,r3Cseq-method}
\title{Get aligned reads from the BAM file}

\description{
  Reading in the input BAM file and then store it in the GRanged object
}

\usage{
getRawReads(object)
}


\arguments{

	\item{object}{
		r3Cseq object
	}
}

\value{
   The GRangedData represents the aligned reads from the BAM file
}

\seealso{
  \link{getBatchRawReads},
}

\author{
	S. Thongjuea
}

\examples{
  #See the vignette	
}
%\keyword{raw read manipulation}
