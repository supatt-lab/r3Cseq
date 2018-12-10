\name{getViewpoint}
\alias{getViewpoint}
\alias{getViewpoint,r3Cseq-method}

\title{get the viewpoint of 3C-seq data}

\description{
  The viewpoint is the bait of 3C method, which can be a promoter region of an interested gene, 
  an enhancer, and a transcrition factor binding region. 
   
}

\usage{
getViewpoint(obj)
}


\arguments{

	\item{obj}{
		r3Cseq or r3CseqInBatch object
	}
}

\value{
   The viewpoint shows in the IRanges
}

\author{
	S. Thongjuea
}

\examples{
 #See the vignette
}
%\keyword{viewpoint}
