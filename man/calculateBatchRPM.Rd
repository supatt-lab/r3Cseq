\name{calculateBatchRPM}
\alias{calculateBatchRPM}
\alias{calculateBatchRPM,r3CseqInBatch-method}
\title{calculate read per million (RPM) for replicates analysis}

\description{
  Normalize 3C-Seq data by transforming raw reads to read per million per each region for replication analysis
}

\usage{
calculateBatchRPM(object,normalized_method=c("powerlawFittedRPM","normalRPM"))
}


\arguments{

	\item{object}{
		r3CseqInBatch object
	}
	\item{normalized_method}{
		character. method of normalization (default=powerlawFittedRPM)
	}
}

\seealso{
  \link{calculateRPM},
  \link{expRPM}
  \link{contrRPM}
}

\author{
	S. Thongjuea
}

\examples{
#See the vignette
}
%\keyword{normalization}
