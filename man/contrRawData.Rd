\name{contrRawData}
\docType{methods}
\alias{contrRawData}
\alias{contrRawData<-}
\alias{contrRawData,r3Cseq-method}
\alias{contrRawData<-,r3Cseq-method}

\title{Accessors for the 'contrRawData' slot of a r3Cseq object.}

\description{
  The 'contrRawData' slot of hold the raw aligned reads data in the GRanges object.
}

\usage{
\S4method{contrRawData}{r3Cseq}(object)
\S4method{contrRawData}{r3Cseq}(object) <- value
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
