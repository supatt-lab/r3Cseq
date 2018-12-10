\name{getExpInteractionsInRefseq}
\alias{getExpInteractionsInRefseq}
\alias{getExpInteractionsInRefseq,r3Cseq-method}

\title{identified significant interaction regions for RefSeq genes}

\description{
  Get a list of genes that contain strong interaction signals in the experiment
}

\usage{
getExpInteractionsInRefseq(obj,cutoff.qvalue=0.05,expanded_upstream=50e3,expanded_downstream=10e3)
}


\arguments{

	\item{obj}{
		obj is r3Cseq or r3CseqInBatch object
	}

	\item{cutoff.qvalue}{
		Numeric. The cutoff q-value (default=0.05)
	}
	\item{expanded_upstream}{
		Numeric. The expanded distance from the upstream of a gene start (default=50e3)
	}
	\item{expanded_downstream}{
		Numeric. The expanded distance from the downstream of a gene end (default =10e3)
	}
	
}

\value{
	List of identified genes, which contain strong interaction signals
}

\author{
	S. Thongjuea
}

\seealso{
	\code{\link{getContrInteractionsInRefseq}}
}


\examples{
 # See the vignette
}
%\keyword{Interactions in RefGene}
