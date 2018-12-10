# TODO: The following functions are implemented for exporting the output of 3C-seq results.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com
###############################################################################
setGeneric(
		name="export3CseqRawReads2bedGraph",
		def=function(object){
			standardGeneric("export3CseqRawReads2bedGraph")
		}

)
setMethod("export3CseqRawReads2bedGraph",
		signature(object = "r3Cseq"),
		function (object){
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			
			if(isControlInvolved(object)==FALSE){
				expRawReads.GRanges<-expRawData(object)
				expLabeled<-expLabel(object)	
				orgName  <- organismName(object)
				
				if(length(expRawReads.GRanges)>0){	
					print("making coverage vector......")
					exp.read.cov      <- RleList()
					#***Fix the error in this version by changing from RangeData to GRanges
					exp.read.cov   <- coverage(expRawReads.GRanges)
					
					export.iranges<-as(exp.read.cov,"GRanges")
					export.ucsc<-as(export.iranges,"UCSCData")
					export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=expLabeled,description="Raw reads")		
					file_name<-paste(expLabeled,".bedGraph",".gz",sep="")
					export(export.ucsc,file_name)
					print(paste("File",file_name,"' is created."))
				}else{
					stop("No raw reads found in the object, please run getRawReads function!!!.")
				}
			}
			if(isControlInvolved(object)==TRUE){
				expRawReads.GRanges<-expRawData(object)
				contrRawReads.GRanges<-contrRawData(object)
				expLabeled<-expLabel(object)
				contrLabeled<-contrLabel(object)
				
				orgName  <- organismName(object)
				
				if(length(expRawReads.GRanges)>0){	
					print("making coverage vector......")
					exp.read.cov      <- RleList()
					contr.read.cov      <- RleList()
					
					exp.read.cov   <- coverage(expRawReads.GRanges)	
					contr.read.cov   <- coverage(contrRawReads.GRanges)
					
					#***Fix the error in this version by changing from RangeData to GRanges
					export.iranges<-as(exp.read.cov,"GRanges")
					export.ucsc<-as(export.iranges,"UCSCData")
					export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=expLabeled,description="Raw reads")		
					file_name<-paste(expLabeled,".bedGraph",".gz",sep="")
					export(export.ucsc,file_name)
					print(paste("File",file_name,"' is created."))
					
					#***Fix the error in this version by changing from RangeData to GRanges
					export.iranges<-as(contr.read.cov,"GRanges")
					export.ucsc<-as(export.iranges,"UCSCData")
					export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=contrLabeled,description="Raw reads")		
					file_name<-paste(contrLabeled,".bedGraph",".gz",sep="")
					export(export.ucsc,file_name)
					print(paste("File",file_name,"' is created."))
					
				}else{
					stop("No raw reads found in the object, please run getRawReads function!!!.")
				}
			}
		}
)

setGeneric(
		name="export3Cseq2bedGraph",
		def=function(object,datatype = c("rpm", "read_count")){
			standardGeneric("export3Cseq2bedGraph")
		}
)
setMethod("export3Cseq2bedGraph",
		signature(object = "r3Cseq"),
		function (object,datatype){
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			datatype <- match.arg(datatype)
			if(datatype==""){
				datatype<-"rpm"
			}
			if(isControlInvolved(object)==FALSE){
				expRPMs   <-expRPM(object)
				expLabeled<-expLabel(object)
				
				if(datatype=="rpm"){	
					if(length(expRPMs)>0){
						##Fix export.f by adding "as.character"
						export.f <-data.frame(seqnames=as.character(seqnames(expRPMs)),start=start(expRPMs),end=end(expRPMs),score=expRPMs$RPMs)
						export.GRanges<-GRanges(seqnames=export.f$seqnames,IRanges(start=export.f$start,end=export.f$end),score=export.f$score)
						#export.iranges<-export.iranges[export.iranges$score>0,]<--fix the error in this version
					 	export.GRanges<-export.GRanges[export.GRanges$score>0,]
					 	export.ucsc<-as(export.GRanges,"UCSCData")
					 	export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=paste(expLabeled,"_RPMs",sep=""),description="read per million")		
					 	file_name<-paste(expLabeled,".RPMs.bedGraph.gz",sep="")
					 	export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
					}else{
						stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}
				}else if(datatype=="read_count"){
					if(length(expRPMs)>0){
						export.f <-data.frame(seqnames=as.character(seqnames(expRPMs)),start=start(expRPMs),end=end(expRPMs),score=expRPMs$nReads)
						export.GRanges<-GRanges(seqnames=export.f$seqnames,IRanges(start=export.f$start,end=export.f$end),score=export.f$score)
						#export.iranges<-export.iranges[export.iranges$score>0,]
						export.GRanges<-export.GRanges[export.GRanges$score>0,]
						export.ucsc<-as(export.GRanges,"UCSCData")
						export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=paste(expLabeled,"_ReadCounts",sep=""),description="read count")	
						file_name<-paste(expLabeled,".ReadCounts.bedGraph.gz",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
					}else{
						stop("No read found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}
					
				}else{
					stop("Choose the input datatype parameter: 'read_count' or 'rpm' (reads per million)")
				}
			}
			if(isControlInvolved(object)==TRUE){
				expRPMs   <-expRPM(object)
				contrRPMs <-contrRPM(object)
				expLabeled<-expLabel(object)
				controlLabeled<-contrLabel(object)
				
				if(datatype=="rpm"){

					if(length(expRPMs)>0){
						export.GRanges<-GRanges(seqnames =as.character(seqnames(expRPMs)),IRanges(start=start(expRPMs),end=end(expRPMs)),score=expRPMs$RPMs)
						#export.iranges<-export.iranges[export.iranges$score>0,]
						export.GRanges<-export.GRanges[export.GRanges$score>0,]
						export.ucsc<-as(export.GRanges,"UCSCData")
						export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=paste(expLabeled,"_RPMs",sep=""),description="reads per million")	
						file_name<-paste(expLabeled,".RPMs.bedGraph.gz",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
						

						export.GRanges<-GRanges(seqnames=as.character(seqnames(contrRPMs)),IRanges(start=start(contrRPMs),end=end(contrRPMs)),score=contrRPMs$RPMs)
						#export.iranges<-export.iranges[export.iranges$score>0,]
						export.GRanges<-export.GRanges[export.GRanges$score>0,]
						export.ucsc<-as(export.GRanges,"UCSCData")
						export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=paste(controlLabeled,"_RPMs",sep=""),description="reads per million")	
						file_name<-paste(controlLabeled,".RPMs.bedGraph.gz",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))	
						
					}else{
						stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}
				}else if(datatype=="read_count"){
					if(length(expRPMs)>0){
					  export.GRanges<-GRanges(seqnames=as.character(seqnames(expRPMs)),IRanges(start=start(expRPMs),end=end(expRPMs)),score=expRPMs$nReads)
						#export.iranges<-export.iranges[export.iranges$score>0,]
					  export.GRanges<-export.GRanges[export.GRanges$score>0,]
						export.ucsc<-as(export.GRanges,"UCSCData")
						export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=paste(expLabeled,"_ReadCounts",sep=""),description="read count")	
						file_name<-paste(expLabeled,".ReadCounts.bedGraph.gz",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))
					
						export.GRanges<-GRanges(seqnames=as.character(seqnames(contrRPMs)),IRanges(start=start(contrRPMs),end=end(contrRPMs)),score=contrRPMs$nReads)
						#export.iranges<-export.iranges[export.iranges$score>0,]
						export.GRanges<-export.GRanges[export.GRanges$score>0,]
						export.ucsc<-as(export.GRanges,"UCSCData")
						export.ucsc@trackLine <- new("GraphTrackLine",type="bedGraph",name=paste(controlLabeled,"_ReadCounts",sep=""),description="read count")	
						file_name<-paste(controlLabeled,".ReadCounts.bedGraph.gz",sep="")
						export(export.ucsc,file_name,"bedGraph")
						print(paste("File",file_name,"' is created."))	
					}else{
						stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
					}					
				}else{
					stop("Choose the input datatype parameter: 'read_count' or 'rpm' (reads per million)")
				}
				
			}
	}
)

setGeneric(
		name="exportInteractions2text",
		def=function(object){
			standardGeneric("exportInteractions2text")
		}
)

setMethod("exportInteractions2text",
		signature(object = "r3Cseq"),
		function (object){
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(isControlInvolved(object)==FALSE){
				expInteractions<-expInteractionRegions(object)
				expLabeled<-expLabel(object)
				
					if(length(expInteractions)>0){
						
						file_name<-paste(expLabeled,".interaction.txt",sep="")
						export.data<-data.frame(chromosome=seqnames(expInteractions),
								start=start(expInteractions),end=end(expInteractions),
								nReads=expInteractions$nReads,RPMs=expInteractions$RPMs,
								p.value=expInteractions$p.value,q.value=expInteractions$q.value)
						export.data<-export.data[order(export.data[,7]),]
						export.data$p.value<-format(export.data$p.value,scientific=T)
						export.data$q.value<-format(export.data$q.value,scientific=T)
						write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)					
						print(paste("File",file_name,"' is created."))
					}else{
						stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
					}
			}
			if(isControlInvolved(object)==TRUE){
				expInteractions  <-expInteractionRegions(object)
				contrInteractions <- contrInteractionRegions(object)
				expLabeled<-expLabel(object)
				controlLabeled<-contrLabel(object)
				
				if(length(expInteractions)>0){
					
					file_name<-paste(expLabeled,".interaction.txt",sep="")
					export.data<-data.frame(chromosome=seqnames(expInteractions),
							start=start(expInteractions),end=end(expInteractions),
							nReads=expInteractions$nReads,RPMs=expInteractions$RPMs,
							p.value=expInteractions$p.value,q.value=expInteractions$q.value)
					export.data<-export.data[order(export.data[,7]),]
					export.data$p.value<-format(export.data$p.value,scientific=T)
					export.data$q.value<-format(export.data$q.value,scientific=T)
					write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)					
					print(paste("File",file_name,"' is created."))
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}
				
				if(length(contrInteractions)>0){
					
					file_name<-paste(controlLabeled,".interaction.txt",sep="")
					export.data<-data.frame(chromosome=seqnames(contrInteractions),
							start=start(contrInteractions),end=end(contrInteractions),
							nReads=contrInteractions$nReads,RPMs=contrInteractions$RPMs,
							p.value=contrInteractions$p.value,q.value=contrInteractions$q.value)
					export.data<-export.data[order(export.data[,7]),]
					export.data$p.value<-format(export.data$p.value,scientific=T)
					export.data$q.value<-format(export.data$q.value,scientific=T)
					write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)	
					print(paste("File",file_name,"' is created."))
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}
				
			}
	}
)

setGeneric(
		name="exportBatchInteractions2text",
		def=function(object){
			standardGeneric("exportBatchInteractions2text")
		}
)

setMethod("exportBatchInteractions2text",
		signature(object = "r3CseqInBatch"),
		function (object){
			if(!is(object,"r3CseqInBatch")){
				stop("Need the r3CseqInBatch object")
			}
			if(isControlInvolved(object)==TRUE){
				expInteractions  <-expInteractionRegions(object)
				contrInteractions <- contrInteractionRegions(object)
				
				if(length(expInteractions)>0){
					
					file_name<-"Batch.experiment.interaction.txt"
					export.data<-data.frame(chromosome=seqnames(expInteractions),
							start=start(expInteractions),end=end(expInteractions),
							nReads=expInteractions$nReads,RPMs=expInteractions$RPMs,
							p.value=expInteractions$p.value,q.value=expInteractions$q.value)
					export.data<-export.data[order(export.data[,7]),]
					export.data$p.value<-format(export.data$p.value,scientific=T)
					export.data$q.value<-format(export.data$q.value,scientific=T)
					write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)					
					print(paste("File",file_name,"' is created."))
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}
				
				if(length(contrInteractions)>0){
					
					file_name<-"Batch.control.interaction.txt"
					export.data<-data.frame(chromosome=seqnames(contrInteractions),
							start=start(contrInteractions),end=end(contrInteractions),
							nReads=contrInteractions$nReads,RPMs=contrInteractions$RPMs,
							p.value=contrInteractions$p.value,q.value=contrInteractions$q.value)
					export.data<-export.data[order(export.data[,7]),]
					export.data$p.value<-format(export.data$p.value,scientific=T)
					export.data$q.value<-format(export.data$q.value,scientific=T)
					write.table(export.data,file=file_name,append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)	
					print(paste("File",file_name,"' is created."))
				}else{
					stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
				}
				
			}
		}
)




