\name{getBatchRawReads}
\alias{getBatchRawReads}
\alias{getBatchRawReads,r3CseqInBatch-method}
\title{Get aligned reads from the replicates BAM files }

\description{
  Reading in the input BAM files from the 3C-Seq replicates analysis and 
  then save files as the local GRanged object .rData files
}

\usage{
getBatchRawReads(object)
}


\arguments{

	\item{object}{
		r3CseqInBatch object
	}
}

\value{
   The GRangedData represents the aligned reads from the BAM file
}

\seealso{
  \link{getRawReads},
}

\author{
	S. Thongjuea
}

\examples{
  #See the vignette	
}
%\keyword{raw read manipulation}
