\name{calculateRPM}
\alias{calculateRPM}
\alias{calculateRPM,r3Cseq-method}
\title{calculate read per million (RPM)}

\description{
  Normalize 3C-Seq data by transforming raw reads to read per million per each region 
}

\usage{
calculateRPM(object,normalized_method=c("powerlawFittedRPM","normalRPM"))
}


\arguments{

	\item{object}{
		r3Cseq object
	}
	\item{normalized_method}{
		character. method of normalization (default=powerlawFittedRPM)
	}
}

\seealso{

  \link{contrRPM},
  \link{expRPM},
  \link{calculateBatchRPM}
}

\author{
	S. Thongjuea
}

\examples{
#See the vignette
}
%\keyword{normalization}
