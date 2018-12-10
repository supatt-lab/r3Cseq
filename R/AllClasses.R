setClass(
		Class="r3CseqCommon",
		representation(
				organismName="character",
				restrictionEnzyme="character",
				viewpoint_chromosome="character",
				viewpoint_primer_forward="character",
				viewpoint_primer_reverse="character",
				expReadCount="GRanges",
				contrReadCount="GRanges",
				expRPM="GRanges",
				contrRPM="GRanges",
				expInteractionRegions="GRanges",
				contrInteractionRegions="GRanges",
				isControlInvolved="logical"
				
		),
		prototype(
				organismName=character(1),
				restrictionEnzyme=character(1),
				viewpoint_chromosome=character(1),
				viewpoint_primer_forward=character(1),
				viewpoint_primer_reverse=character(1),
				expReadCount=GRanges(),
				contrReadCount=GRanges(),
				expRPM=GRanges(),
				contrRPM=GRanges(),
				expInteractionRegions=GRanges(),
				contrInteractionRegions=GRanges(),
				isControlInvolved=FALSE
		),
		validity = function(object) {
			if(object@organismName==character(1))
				return( "The organism name is empty" )
			if(!object@organismName %in% c("mm9","mm10","hg18","hg19","rn5"))
				return( "This version of r3Cseq supports only for mm9, mm10, hg18, hg19, and rn5 assembly." )
			##############################################
			if(length(object@restrictionEnzyme)==0)
				return( "Please enter the primary_restrictionEnzyme parameter" )
			if(object@restrictionEnzyme==character(1))
				return( "The restrictionEnzyme is empty" )
			#############################################
			if(length(object@viewpoint_chromosome)==0)
				return( "Please enter the viewpoint_chromosome parameter" )
			if(object@viewpoint_chromosome==character(1))
				return( "The chromosome name of the viewpoint is empty. Please enter the chromosome name such as 'chr10'")
			if(grep("^chr[0-9X-Y]",object@viewpoint_chromosome)==0)
				return( "The chromosome name of the viewpoint is incorrect format. Please enter the chromosome name for example: 'chr10'")
			#############################################
			if(length(object@viewpoint_primer_forward)==0)
				return( "Please enter the viewpoint_primer_forward parameter" )
			if(object@viewpoint_primer_forward==character(1))
				return( "The forward primer is empty. Please enter the forward primer sequences")
			##############################################
			if(length(object@viewpoint_primer_reverse)==0)
				return( "Please enter the viewpoint_primer_reverse parameter" )
			if(object@viewpoint_primer_reverse==character(1))
				return( "The reverse primer is empty. Please enter the reverse primer sequences")
		}
)

setClass(
		Class="r3Cseq",contain="r3CseqCommon",
		representation(
				alignedReadsBamExpFile="character",
				alignedReadsBamContrFile="character",
				expLabel ="character",
				contrLabel="character",
				expLibrarySize ="integer",
				contrLibrarySize ="integer",
				expReadLength="integer",
				contrReadLength="integer",
				expRawData="GRanges",
				contrRawData="GRanges"
		),
		prototype(
				
				alignedReadsBamExpFile=character(1),
				alignedReadsBamContrFile=character(1),
				expLabel =character(1),
				contrLabel=character(1),
				expRawData=GRanges(),
				contrRawData=GRanges()
		),
		validity = function(object) {
			if(length(object@expLabel)==0)
				return( "Please enter the experiment name" )
			if(object@expLabel==character(1))
				return( "Please enter the experiment name" )
			##############################################
		}
)

setClass("repbaseEnzyme",
		representation(
				enzymeRestriction="data.frame"
		),
		prototype(
				enzymeRestriction=data.frame()
		)
)

setClass(
		Class="r3CseqInBatch",contain="r3CseqCommon",
		representation(
				bamFilesDirectory="character",
				BamExpFiles="vector",
				BamContrFiles="vector",
				expBatchLabel ="vector",
				contrBatchLabel="vector",
				readCountTable="GRanges",
				RPMsTable="GRanges",
				expBatchLibrarySize="vector",
				contrBatchLibrarySize="vector",
				expBatchReadLength="vector",
				contrBatchReadLength="vector"
		),
		prototype(
				bamFilesDirectory=character(1),
				BamExpFiles=vector(),
				BamContrFiles=vector(),
				expBatchLabel=vector(),
				contrBatchLabel=vector(),
			  readCountTable=GRanges(),
				RPMsTable=GRanges(),
				expBatchLibrarySize=vector(),
				contrBatchLibrarySize=vector(),
				expBatchReadLength=vector(),
				contrBatchReadLength=vector()
		),
		validity = function(object) {
			if(length(object@bamFilesDirectory)==0)
				return( "Please enter the bamFilesDirectory parameter" )
			if(object@bamFilesDirectory==character(1))
				return( "The folder name of the BAM files is empty. Please enter the PATH to the BAM files" )
			##############################################
			if(length(object@BamExpFiles)==0)
				return( "Please enter the list of BAM file names in the experiments" )
			if(length(object@BamContrFiles)==0)
				return( "Please enter the list of BAM file names in the controls" )
			##############################################
			if(length(object@expBatchLabel)==0)
				return( "Please enter the list of label for the experiments" )
			if(length(object@contrBatchLabel)==0)
				return( "Please enter the list of label for the controls" )
		}
)
