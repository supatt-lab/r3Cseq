\name{r3CseqCommon-class}
\docType{class}

% Classes:
\alias{r3CseqCommon}
\alias{r3CseqCommon-class}

\title{r3CseqCommon objects}

\description{
  The r3CseqCommon class is a general container for storing and manipulating a set of input parameters, RangeData of interactions regions from r3Cseq analysis.
  It is a root class for r3Cseq and r3CseqInBatch classes. 
}

\section{Slots}{
  \describe{

    \item{\code{organismName}}{Object of class \code{"character"} the version of 
      particular assembly genome from UCSC (e.g. mm9, hg18, hg19) . The package supports
      three genome assemblies consisting of mouse (mm9), and human (hg18, hg19).}

    \item{\code{restrictionEnzyme}}{Object of class \code{"character"} this is the primary
      restriction enzyme name using in 3C-Seq experiment}

    \item{\code{viewpoint_chromosome}}{Object of class \code{"character"} 
    	chromosome name of where is the viewpoint located eg. chr10, chrX etc.
    }

    \item{\code{viewpoint_primer_forward}}{Object of class \code{"character"}
      the forward primer DNA sequences for the viewpoint amplification}
      
    \item{\code{viewpoint_primer_reverse}}{Object of class \code{"character"}
      the reverse primer DNA sequences for the viewpoint amplification}
    
    \item{\code{expReadCount}}{Object of class \code{"RangedData"}
      the read count in experiment}  
    
    \item{\code{contrReadCount}}{Object of class \code{"RangedData"}
      the read count in control}  
    
    \item{\code{expRPM}}{Object of class \code{"RangedData"}
      the normalized read read per million in experiment}  
   
    \item{\code{contrRPM}}{Object of class \code{"RangedData"}
      the normalized read read per million in control}  
   
    \item{\code{expInteractionRegions}}{Object of class \code{"RangedData"}
      the identified interaction regions in experiment}  
   
    \item{\code{contrInteractionRegions}}{Object of class \code{"RangedData"}
      the identified interaction regions in control}  
      
 	\item{\code{isControlInvolved}}{Object of class \code{"logical"}
      the logical to ask whether the control is involved in the analysis or not}		
  }
}


\author{S. Thongjuea}

\seealso{
  \link{r3Cseq},
  \link{r3CseqInBatch}
}
\examples{
# See the vignette
}

\keyword{classes}
