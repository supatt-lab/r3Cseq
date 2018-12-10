\name{expRawData}
\docType{methods}
\alias{expRawData}
\alias{expRawData<-}
\alias{expRawData,r3Cseq-method}
\alias{expRawData<-,r3Cseq-method}

\title{Accessors for the 'expRawData' slot of a r3Cseq object.}

\description{
  The 'expRawData' slot of hold the raw aligned reads data in the GRanges object.
}

\usage{
\S4method{expRawData}{r3Cseq}(object)
\S4method{expRawData}{r3Cseq}(object) <- value
}


\arguments{

	\item{object}{
		r3Cseq object
	}
	\item{value}{
	     a GRanges object of aligned reads
    }
}

\seealso{
  \link{expRawData}
}

\author{
	S. Thongjuea
}

\examples{

#See the vignette	
}
%\keyword{raw reads}
