\name{generate3CseqReport}
\alias{generate3CseqReport}
\alias{generate3CseqReport,r3Cseq-method}
\alias{generate3CseqReport,r3CseqInBatch-method}
\title{generate reports for analysis results from r3Cseq}

\description{
  generate reports for analysis results from r3Cseq, the report contains 
  all plots in one pdf file and a text separated out put file.
}

\usage{
generate3CseqReport(obj)
}


\arguments{

	\item{obj}{
		r3Cseq or r3CseqInBatch object
	}
}

\value{
  The text file in the tab separated format and the pdf file of all plots
}

\seealso{

  \link{exportInteractions2text}
  \link{plotOverviewInteractions}, 
  \link{plotInteractionsPerChromosome}, 
  \link{plotInteractionsNearViewpoint}
}

\author{
	S. Thongjuea
}

\examples{
 	
 #See the vignette
}
%\keyword{creat the report for r3Cseq analysis}
