\name{r3CseqInBatch-class}
\docType{class}

% Classes:
\alias{r3CseqInBatch}
\alias{r3CseqInBatch-class}

\title{r3CseqInBatch objects}

\description{
  The r3CseqInBatch class is the extended class from r3CseqCommon class. It is a general container for storing
  and manipulating a set of input parameters, RangeData of interactions regions from r3Cseq analysis for replicates data sets.
}
\section{Extends}{
Class \code{r3CseqCommon}, directly.
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
  
    
    \item{\code{bamFilesDirectory}}{Object of class \code{"character"}
      the path name of directory that contains BAM files}
      
    \item{\code{BamExpFiles}}{Object of class \code{"vector"}
      the file names of BAM files in the experiment}
        
    \item{\code{BamContrFiles}}{Object of class \code{"vector"}
      the file names of BAM files in the control}
      
    \item{\code{expBatchLabel}}{Object of class \code{"vector"}
      the labeled experiment names}
    
    \item{\code{contrBatchLabel}}{Object of class \code{"vector"}
      the labeled control names}  
      
    \item{\code{readCountTable}}{Object of class \code{"RangedData"}
      the read count table}  
    
    \item{\code{RPMsTable}}{Object of class \code{"RangedData"}
      the normalized read per million table}  
      
    \item{\code{expBatchLibrarySize}}{Object of class \code{"vector"}
      the library size of each experiment} 
    
    \item{\code{contrBatchLibrarySize}}{Object of class \code{"vector"}
      the library size of each control} 
    
    \item{\code{expBatchReadLength}}{Object of class \code{"vector"}
      the read length of experiments}  
      
    \item{\code{contrBatchReadLength}}{Object of class \code{"vector"}
      the read length of controls} 
  }
}

\author{S. Thongjuea}

\seealso{
  \link{r3CseqCommon},
  \link{r3CseqInBatch}
}
\examples{
# See the vignette
}

\keyword{classes}
