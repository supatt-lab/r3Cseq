# Author: Supat Thongjuea
# Contact : supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com
####################################
#########  AllGenerics.R
#########
############Common methods##########
setGeneric(
		name="organismName",
		def=function(object){
			standardGeneric("organismName")
		})
setMethod("organismName",
		signature(object = "r3CseqCommon"),
		function (object){
			object@organismName
		})
########
setGeneric(
		name="restrictionEnzyme",
		def=function(object){
			standardGeneric("restrictionEnzyme")
		})
setMethod("restrictionEnzyme",
		signature(object = "r3CseqCommon"),
		function (object){
			object@restrictionEnzyme
		}
)
########
setGeneric(
		name="viewpoint_chromosome",
		def=function(object){
			standardGeneric("viewpoint_chromosome")
		})
setMethod("viewpoint_chromosome",
		signature(object = "r3CseqCommon"),
		function (object){
			object@viewpoint_chromosome
		}
)
########
setGeneric(
		name="viewpoint_primer_forward",
		def=function(object){
			standardGeneric("viewpoint_primer_forward")
		})
setMethod("viewpoint_primer_forward",
		signature(object = "r3CseqCommon"),
		function (object){
			object@viewpoint_primer_forward
		}
)
########
setGeneric(
		name="viewpoint_primer_reverse",
		def=function(object){
			standardGeneric("viewpoint_primer_reverse")
		})
setMethod("viewpoint_primer_reverse",
		signature(object = "r3CseqCommon"),
		function (object){
			object@viewpoint_primer_reverse
		}
)
#########
setGeneric(
		name="expReadCount",
		def=function(object){
			standardGeneric("expReadCount")
		}
)
setGeneric(
		name="expReadCount<-",
		def=function(object,value){
			standardGeneric("expReadCount<-")
		}
)
setMethod("expReadCount",
		signature(object = "r3CseqCommon"),
		function (object){
			object@expReadCount
		}
)
setReplaceMethod(
		f="expReadCount",
		signature="r3CseqCommon",
		definition=function(object,value){
			initialize(object,expReadCount=value)
		})
#########
setGeneric(
		name="contrReadCount",
		def=function(object){
			standardGeneric("contrReadCount")
		}
)
setGeneric(
		name="contrReadCount<-",
		def=function(object,value){
			standardGeneric("contrReadCount<-")
		}
)
setMethod("contrReadCount",
		signature(object = "r3CseqCommon"),
		function (object){
			object@contrReadCount
		}
)
setReplaceMethod(
		f="contrReadCount",
		signature="r3CseqCommon",
		definition=function(object,value){
			initialize(object,contrReadCount=value)
		})
###########
setGeneric(
		name="expRPM",
		def=function(object){
			standardGeneric("expRPM")
		}
)
setGeneric(
		name="expRPM<-",
		def=function(object,value){
			standardGeneric("expRPM<-")
		}
)
setMethod("expRPM",
		signature(object = "r3CseqCommon"),
		function (object){
			object@expRPM
		}
)
setReplaceMethod(
		f="expRPM",
		signature="r3CseqCommon",
		definition=function(object,value){
			initialize(object,expRPM=value)
		})
###########
setGeneric(
		name="contrRPM",
		def=function(object){
			standardGeneric("contrRPM")
		}
)
setGeneric(
		name="contrRPM<-",
		def=function(object,value){
			standardGeneric("contrRPM<-")
		}
)
setMethod("contrRPM",
		signature(object = "r3CseqCommon"),
		function (object){
			object@contrRPM
		}
)
setReplaceMethod(
		f="contrRPM",
		signature="r3CseqCommon",
		definition=function(object,value){
			initialize(object,contrRPM=value)
		})
############
setGeneric(
		name="expInteractionRegions",
		def=function(object){
			standardGeneric("expInteractionRegions")
		}
)
setGeneric(
		name="expInteractionRegions<-",
		def=function(object,value){
			standardGeneric("expInteractionRegions<-")
		}
)
setMethod("expInteractionRegions",
		signature(object = "r3CseqCommon"),
		function (object){
			object@expInteractionRegions
		}
)
setReplaceMethod(
		f="expInteractionRegions",
		signature="r3CseqCommon",
		definition=function(object,value){
			initialize(object,expInteractionRegions=value)
		})
########
setGeneric(
		name="contrInteractionRegions",
		def=function(object){
			standardGeneric("contrInteractionRegions")
		}
)
setGeneric(
		name="contrInteractionRegions<-",
		def=function(object,value){
			standardGeneric("contrInteractionRegions<-")
		}
)
setMethod("contrInteractionRegions",
		signature(object = "r3CseqCommon"),
		function (object){
			object@contrInteractionRegions
		}
)
setReplaceMethod(
		f="contrInteractionRegions",
		signature="r3CseqCommon",
		definition=function(object,value){
			initialize(object,contrInteractionRegions=value)
		})
##########
setGeneric(
		name="isControlInvolved",
		def=function(object){
			standardGeneric("isControlInvolved")
		})
setMethod("isControlInvolved",
		signature(object = "r3CseqCommon"),
		function (object){
			object@isControlInvolved
		}
)
############Method for r3seq class###########
setGeneric(
		name="alignedReadsBamExpFile",
		def=function(object){
			standardGeneric("alignedReadsBamExpFile")
		})
setMethod("alignedReadsBamExpFile",
		signature(object = "r3Cseq"),
		function (object){
			object@alignedReadsBamExpFile
		}
)
#######
setGeneric(
		name="alignedReadsBamContrFile",
		def=function(object){
			standardGeneric("alignedReadsBamContrFile")
		})
setMethod("alignedReadsBamContrFile",
		signature(object = "r3Cseq"),
		function (object){
			object@alignedReadsBamContrFile
		}
)

########
setGeneric(
		name="expLabel",
		def=function(object){
			standardGeneric("expLabel")
		})
setMethod("expLabel",
		signature(object = "r3Cseq"),
		function (object){
			object@expLabel
		}
)
setGeneric(
		name="contrLabel",
		def=function(object){
			standardGeneric("contrLabel")
		})
setMethod("contrLabel",
		signature(object = "r3Cseq"),
		function (object){
			object@contrLabel
		}
)
########
setGeneric(
		name="expLibrarySize",
		def=function(object){
			standardGeneric("expLibrarySize")
		}
)
setGeneric(
		name="expLibrarySize<-",
		def=function(object,value){
			standardGeneric("expLibrarySize<-")
		}
)
setMethod("expLibrarySize",
		signature(object="r3Cseq"),
		function(object){
			object@expLibrarySize
		}
)
setReplaceMethod(
		f="expLibrarySize",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expLibrarySize=value)
		})
#########
setGeneric(
		name="contrLibrarySize",
		def=function(object){
			standardGeneric("contrLibrarySize")
		}
)
setGeneric(
		name="contrLibrarySize<-",
		def=function(object,value){
			standardGeneric("contrLibrarySize<-")
		}
)
setMethod("contrLibrarySize",
		signature(object="r3Cseq"),
		function(object){
			object@contrLibrarySize
		}
)
setReplaceMethod(
		f="contrLibrarySize",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrLibrarySize=value)
		})

##########
setGeneric(
		name="expReadLength",
		def=function(object){
			standardGeneric("expReadLength")
		}
)
setGeneric(
		name="expReadLength<-",
		def=function(object,value){
			standardGeneric("expReadLength<-")
		}
)
setMethod("expReadLength",
		signature(object = "r3Cseq"),
		function (object){
			object@expReadLength
		}
)
setReplaceMethod(
		f="expReadLength",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expReadLength=value)
		})
##########
setGeneric(
		name="contrReadLength",
		def=function(object){
			standardGeneric("contrReadLength")
		}
)
setGeneric(
		name="contrReadLength<-",
		def=function(object,value){
			standardGeneric("contrReadLength<-")
		}
)
setMethod("contrReadLength",
		signature(object = "r3Cseq"),
		function (object){
			object@contrReadLength
		}
)
setReplaceMethod(
		f="contrReadLength",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrReadLength=value)
		})
##########
setGeneric(
		name="expRawData",
		def=function(object){
			standardGeneric("expRawData")
		}
)
setGeneric(
		name="expRawData<-",
		def=function(object,value){
			standardGeneric("expRawData<-")
		}
)
setMethod("expRawData",
		signature(object = "r3Cseq"),
		function (object){
			object@expRawData
		}
)
setReplaceMethod(
		f="expRawData",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expRawData=value)
		})
##########
setGeneric(
		name="contrRawData",
		def=function(object){
			standardGeneric("contrRawData")
		}
)
setGeneric(
		name="contrRawData<-",
		def=function(object,value){
			standardGeneric("contrRawData<-")
		}
)
setMethod("contrRawData",
		signature(object = "r3Cseq"),
		function (object){
			object@contrRawData
		}
)
setReplaceMethod(
		f="contrRawData",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrRawData=value)
		})

#######Method for r3CseqInBatch#########
setGeneric(
		name="bamFilesDirectory",
		def=function(object){
			standardGeneric("bamFilesDirectory")
		})
setMethod("bamFilesDirectory",
		signature(object = "r3CseqInBatch"),
		function (object){
			object@bamFilesDirectory
		}
)
########
setGeneric(
		name="BamExpFiles",
		def=function(object){
			standardGeneric("BamExpFiles")
		})
setMethod("BamExpFiles",
		signature(object = "r3CseqInBatch"),
		function (object){
			object@BamExpFiles
		}
)
########
setGeneric(
		name="BamContrFiles",
		def=function(object){
			standardGeneric("BamContrFiles")
		})
setMethod("BamContrFiles",
		signature(object = "r3CseqInBatch"),
		function (object){
			object@BamContrFiles
		}
)
########
setGeneric(
		name="expBatchLabel",
		def=function(object){
			standardGeneric("expBatchLabel")
		})
setMethod("expBatchLabel",
		signature(object = "r3CseqInBatch"),
		function (object){
			object@expBatchLabel
		}
)
########
setGeneric(
		name="contrBatchLabel",
		def=function(object){
			standardGeneric("contrBatchLabel")
		})
setMethod("contrBatchLabel",
		signature(object = "r3CseqInBatch"),
		function (object){
			object@contrBatchLabel
		}
)
###########
setGeneric(
		name="readCountTable",
		def=function(object){
			standardGeneric("readCountTable")
		}
)
setGeneric(
		name="readCountTable<-",
		def=function(object,value){
			standardGeneric("readCountTable<-")
		})
setMethod("readCountTable",
		signature(object = "r3CseqInBatch"),
		function (object){
			object@readCountTable
		}
)
setReplaceMethod(
		f="readCountTable",
		signature="r3CseqInBatch",
		definition=function(object,value){
			initialize(object,readCountTable=value)
		}
)

###########
setGeneric(
		name="RPMsTable",
		def=function(object){
			standardGeneric("RPMsTable")
		}
)
setGeneric(
		name="RPMsTable<-",
		def=function(object,value){
			standardGeneric("RPMsTable<-")
		})
setMethod("RPMsTable",
		signature(object = "r3CseqInBatch"),
		function (object){
			object@RPMsTable
		}
)
setReplaceMethod(
		f="RPMsTable",
		signature="r3CseqInBatch",
		definition=function(object,value){
			initialize(object,RPMsTable=value)
		}
)
########
setGeneric(
		name="expBatchLibrarySize",
		def=function(object){
			standardGeneric("expBatchLibrarySize")
		}
)
setGeneric(
		name="expBatchLibrarySize<-",
		def=function(object,value){
			standardGeneric("expBatchLibrarySize<-")
		}
)
setMethod("expBatchLibrarySize",
		signature(object="r3CseqInBatch"),
		function(object){
			object@expBatchLibrarySize
		}
)
setReplaceMethod(
		f="expBatchLibrarySize",
		signature="r3CseqInBatch",
		definition=function(object,value){
			initialize(object,expBatchLibrarySize=value)
		})
########
setGeneric(
		name="contrBatchLibrarySize",
		def=function(object){
			standardGeneric("contrBatchLibrarySize")
		}
)
setGeneric(
		name="contrBatchLibrarySize<-",
		def=function(object,value){
			standardGeneric("contrBatchLibrarySize<-")
		}
)
setMethod("contrBatchLibrarySize",
		signature(object="r3CseqInBatch"),
		function(object){
			object@contrBatchLibrarySize
		}
)
setReplaceMethod(
		f="contrBatchLibrarySize",
		signature="r3CseqInBatch",
		definition=function(object,value){
			initialize(object,contrBatchLibrarySize=value)
		})

########
setGeneric(
		name="expBatchReadLength",
		def=function(object){
			standardGeneric("expBatchReadLength")
		}
)
setGeneric(
		name="expBatchReadLength<-",
		def=function(object,value){
			standardGeneric("expBatchReadLength<-")
		}
)
setMethod("expBatchReadLength",
		signature(object="r3CseqInBatch"),
		function(object){
			object@expBatchReadLength
		}
)
setReplaceMethod(
		f="expBatchReadLength",
		signature="r3CseqInBatch",
		definition=function(object,value){
			initialize(object,expBatchReadLength=value)
		})
########
setGeneric(
		name="contrBatchReadLength",
		def=function(object){
			standardGeneric("contrBatchReadLength")
		}
)
setGeneric(
		name="contrBatchReadLength<-",
		def=function(object,value){
			standardGeneric("contrBatchReadLength<-")
		}
)
setMethod("contrBatchReadLength",
		signature(object="r3CseqInBatch"),
		function(object){
			object@contrBatchReadLength
		}
)
setReplaceMethod(
		f="contrBatchReadLength",
		signature="r3CseqInBatch",
		definition=function(object,value){
			initialize(object,contrBatchReadLength=value)
		})

###Removed functions######
#########
setGeneric(
		name="expCoverage",
		def=function(object){
			standardGeneric("expCoverage")
		}
)
setGeneric(
		name="expCoverage<-",
		def=function(object,value){
			standardGeneric("expCoverage<-")
		})

setMethod("expCoverage",
		signature(object = "r3Cseq"),
		function (object){
			stop("This function has been removed.")
		}
)
setReplaceMethod(
		f="expCoverage",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,expCoverage=value)
		})

###########
setGeneric(
		name="contrCoverage",
		def=function(object){
			standardGeneric("contrCoverage")
		}
)
setGeneric(
		name="contrCoverage<-",
		def=function(object,value){
			standardGeneric("contrCoverage<-")
		})
setMethod("contrCoverage",
		signature(object = "r3Cseq"),
		function (object){
			stop("This function has been removed.")
		}
)
setReplaceMethod(
		f="contrCoverage",
		signature="r3Cseq",
		definition=function(object,value){
			initialize(object,contrCoverage=value)
		})

