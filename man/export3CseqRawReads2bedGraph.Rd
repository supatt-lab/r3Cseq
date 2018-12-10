\name{export3CseqRawReads2bedGraph}
\alias{export3CseqRawReads2bedGraph}
\alias{export3CseqRawReads2bedGraph,r3Cseq-method}

\title{export the interaction signal from the raw reads to the 'bedGraph' format}

\description{
  export interaction regions signal to the bedGraph format, which
  suitable for uploading to the UCSC genome browser
}

\usage{
export3CseqRawReads2bedGraph(object)
}


\arguments{

	\item{object}{
		r3Cseq object
	}
}

\value{
  The text file in 'bedGraph' format
}

\seealso{

  \link{exportInteractions2text},
  \link{export3Cseq2bedGraph},
}

\author{
	S. Thongjuea
}

\examples{
#See the vignette	
}
%\keyword{export identified interaction regions}
