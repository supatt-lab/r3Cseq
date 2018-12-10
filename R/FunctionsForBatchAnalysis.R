# TODO: These following functions were implemented to get the interaction regions for replicates data.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com
#####
#####The functions below are completely migrated to BioC 3.9 on R 3.6
#####
########################################
setGeneric(
		name="getBatchRawReads",
		def=function(object){
			standardGeneric("getBatchRawReads")
		}
)
setMethod("getBatchRawReads",
		signature(object = "r3CseqInBatch"),
		function (object){
			if(!is(object,"r3CseqInBatch")){
				stop("Need to initialize the r3CseqInBatch object")
			}
			######check files########################
			bamDir<-bamFilesDirectory(object)
			##########get labels#######################
			expFiles<-BamExpFiles(object)
			contrFiles<-BamContrFiles(object)
			#########checking directory and files####	
			if (file.exists(bamDir)){
				setwd(file.path(bamDir))
			}else{
				stop(paste("Couldn't find path","-->",bamDir))
			}
			objName <- deparse(substitute(object))
			exp.bam.files    <-BamExpFiles(object)
			contr.bam.files  <-BamContrFiles(object)
			
			for(i.file in c(exp.bam.files)){
				if(file.exists(i.file) ==FALSE){
					stop(paste("Couldn't find file","-->",i.file))
				}		
			}
			for(i.file in c(contr.bam.files)){
				if(file.exists(i.file) ==FALSE){
					stop(paste("Couldn't find file","-->",i.file))
				}		
			}
			exp.lib.size<-c()
			exp.read.length<-c()
			for(i in 1:length(exp.bam.files)){
				
					input.bam<-exp.bam.files[i]		    
					print(paste("start reading in...",input.bam, "....file in the experiment"))
					what    <- c("rname", "strand", "pos", "qwidth","seq")
					param   <- ScanBamParam(what=what,flag = scanBamFlag(isUnmappedQuery = FALSE))
					exp.bam <- scanBam(input.bam, param=param)
					exp.length <-width(exp.bam[[1]]$seq[1])
					exp.GRanges<-GRanges(seqnames=as.vector(exp.bam[[1]]$rname),IRanges(start=exp.bam[[1]]$pos,width=exp.length),strand=exp.bam[[1]]$strand)
					fileName<-gsub(".bam","",expFiles[i])
					lib.size <-length(exp.GRanges)
					exp.lib.size<-append(exp.lib.size,lib.size)
					exp.read.length<-append(exp.read.length,exp.length)
					save(exp.GRanges,file=paste(fileName,".rData",sep=""))
			}
			contr.lib.size<-c()
			contr.read.length<-c()
			###########Make count table for each bam file in the controls########
			for(i in 1:length(contr.bam.files)){
					input.bam<-contr.bam.files[i]		    
					print(paste("start reading in...",input.bam, "....file in the control"))
					what    <- c("rname", "strand", "pos", "qwidth","seq")
					param   <- ScanBamParam(what=what,flag = scanBamFlag(isUnmappedQuery = FALSE))
					contr.bam <- scanBam(input.bam, param=param)
					contr.length <-width(contr.bam[[1]]$seq[1])
					contr.GRanges<-GRanges(seqnames=as.vector(contr.bam[[1]]$rname),IRanges(start=contr.bam[[1]]$pos,width=contr.length),strand=contr.bam[[1]]$strand)	
					
					fileName<-gsub(".bam","",contrFiles[i])
					lib.size <-length(contr.GRanges)
					contr.lib.size<-append(contr.lib.size,lib.size)
					contr.read.length<-append(contr.read.length,contr.length)
					save(contr.GRanges,file=paste(fileName,".rData",sep=""))
			}
				expBatchLibrarySize(object)<-exp.lib.size
				contrBatchLibrarySize(object)<-contr.lib.size
				
				expBatchReadLength(object)<-exp.read.length
				contrBatchReadLength(object)<-contr.read.length
				
				assign(objName,object,envir=parent.frame())
				invisible(1)
				print("The raw read files are created!!!.")
		}
)
setGeneric(
		name="getBatchReadCountPerRestrictionFragment",
		def=function(object,getReadsMethod = c("wholeReads", "adjacentFragmentEndsReads"),
				nFragmentExcludedReadsNearViewpoint=2){
			standardGeneric("getBatchReadCountPerRestrictionFragment")
		}
)
setMethod("getBatchReadCountPerRestrictionFragment",
		signature(object = "r3CseqInBatch"),
		function (object,getReadsMethod,nFragmentExcludedReadsNearViewpoint){
			if(!is(object,"r3CseqInBatch")){
				stop("Need to initialize the r3CseqInBatch object")
			}
			selected_methods <- match.arg(getReadsMethod)
			if(selected_methods==""){
				selected_methods<-"wholeReads"
			}
			######check the number of fragment to exclude reads around the viewpoint###
			nExcludedFragments <- nFragmentExcludedReadsNearViewpoint
			if(nExcludedFragments <0 | nExcludedFragments>5){
				stop("The number of fragment that uses for removing reads around the viewpoint must be 0-5 restriction fragments.")
			}
			######check files########################
			objName <- deparse(substitute(object))
			bamDir<-bamFilesDirectory(object)
			#########checking directory and files####	
			if (file.exists(bamDir)){
				setwd(file.path(bamDir))
			}else{
				stop(paste("Couldn't find path","-->",bamDir))
			}
			
			exp.bam.files    <-BamExpFiles(object)
			contr.bam.files  <-BamContrFiles(object)
					
			for(i.file in c(exp.bam.files)){
				i.file<-gsub(".bam","",i.file)
				i.file<-paste(i.file,".rData",sep="")
				if(file.exists(i.file) ==FALSE){
					stop(paste("Couldn't find file","-->",i.file, ". Please run getBatchRawReads function"))
				}		
			}
			for(i.file in c(contr.bam.files)){
				i.file<-gsub(".bam","",i.file)
				i.file<-paste(i.file,".rData",sep="")
				if(file.exists(i.file) ==FALSE){
					stop(paste("Couldn't find file","-->",i.file, ". Please run getBatchRawReads function"))
				}		
			}
			####Get organism#######################
			orgName  <- organismName(object)
			####Initialized Restriction Enzyme Object##
			enzymeDb	<-new("repbaseEnzyme")
			resEnzyme 	<-restrictionEnzyme(object)
			####Build GRanges for restriction site#####
			print("Fragmenting genome..........")
			wholeFragments<-getWholeGenomeRestrictionFragments(enzymeDb,resEnzyme,orgName)
			##########get labels#######################
			expFiles<-BamExpFiles(object)
			contrFiles<-BamContrFiles(object)
			###########################################
			exp.read.length<-expBatchReadLength(object)
			contr.read.length<-contrBatchReadLength(object)
			
			if(selected_methods=="wholeReads"){
				###########Making GRange Object############
				wholeFragments.GRanges<-GRanges(seqnames = as.character(wholeFragments$chromosome),IRanges(start=wholeFragments$start,end=wholeFragments$end))
				###########Make count table for each bam file in the experiments########
				readCountTable<-data.frame(chromosome=seqnames(wholeFragments.GRanges),start=start(wholeFragments.GRanges),end=end(wholeFragments.GRanges))
				for(i in 1:length(exp.bam.files)){   
					print(paste("start counting reads in...",exp.bam.files[i], "....file in the experiment"))
					fileName<-gsub(".bam","",expFiles[i])
					fileName<-paste(fileName,".rData",sep="")
					load(fileName)
					exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
					readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=exp.GRanges)
					reads<-data.frame(nReads=readsPerFragment)
					Name<-gsub(".bam","",expFiles[i])
					colnames(reads)[1]<-Name
					readCountTable<-cbind(readCountTable,reads)
				}
				###########Make count table for each bam file in the controls########
				for(i in 1:length(contr.bam.files)){	    
					print(paste("start counting reads in...",contr.bam.files[i], "....file in the control"))
					fileName<-gsub(".bam","",contrFiles[i])
					fileName<-paste(fileName,".rData",sep="")
					load(fileName)
					contr.GRanges<-excludeReadsNearViewpoint(object,contr.GRanges,nExcludedFragments)
					readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=contr.GRanges)
					reads<-data.frame(nReads=readsPerFragment)
					Name<-gsub(".bam","",contrFiles[i])
					colnames(reads)[1]<-Name
					readCountTable<-cbind(readCountTable,reads)
				}
				#########Filter out fragment with no read###########################
					n<-ncol(readCountTable)
					readCountTable.f<-readCountTable[rowSums(readCountTable[,4:n]) >0,]
					readCountTable.GRanges<-GRanges(seqnames=readCountTable.f$chromosome,
					                                IRanges(start=readCountTable.f$start,end=readCountTable.f$end))
					values(readCountTable.GRanges)<-readCountTable.f[,4:n]
					readCountTable(object)<-readCountTable.GRanges
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("The countTable is created!!!.")
			}
			if(selected_methods=="adjacentFragmentEndsReads"){
				###########Making GRange Object############
				wholeFragments.GRanges<-GRanges(seqnames=as.character(wholeFragments$chromosome),
				                                IRanges(start=wholeFragments$start,end=wholeFragments$end))
				###########Make count table for each bam file in the experiments########
				readCountTable<-data.frame(chromosome=seqnames(wholeFragments.GRanges),
				                           start=start(wholeFragments.GRanges),end=end(wholeFragments.GRanges))
				for(i in 1:length(exp.bam.files)){    
					print(paste("start counting reads in...",exp.bam.files[i], "....file in the experiment"))
					fileName<-gsub(".bam","",expFiles[i])
					fileName<-paste(fileName,".rData",sep="")
					load(fileName)
					exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
					########Make RangedData for the 5' and 3' of RE#######
					wholeFragments_5_prime.GRanges<-GRanges(seqnames=seqnames(wholeFragments.GRanges),
					                                        IRanges(start=start(wholeFragments.GRanges),
					                                        end=start(wholeFragments.GRanges)+exp.read.length[i]))
					
					wholeFragments_3_prime.GRanges<-GRanges(seqnames=seqnames(wholeFragments.GRanges),
					                                        IRanges(start=end(wholeFragments.GRanges)-exp.read.length[i],
					                                        end=end(wholeFragments.GRanges)))
					
					readsPerFragment_5 <- countOverlaps(wholeFragments_5_prime.GRanges,subject=exp.GRanges)
					readsPerFragment_3 <- countOverlaps(wholeFragments_3_prime.GRanges,subject=exp.GRanges)
					combined_reads<-readsPerFragment_5+readsPerFragment_3
					reads<-data.frame(nReads=combined_reads)
					fileName<-gsub(".bam","",expFiles[i])
					colnames(reads)[1]<-fileName
					readCountTable<-cbind(readCountTable,reads)
				}
				for(i in 1:length(contr.bam.files)){
					print(paste("start counting reads in...",contr.bam.files[i], "....file in the control"))
					fileName<-gsub(".bam","",contrFiles[i])
					fileName<-paste(fileName,".rData",sep="")
					load(fileName)	
					contr.GRanges<-excludeReadsNearViewpoint(object,contr.GRanges,nExcludedFragments)
					########Make RangedData for the 5' and 3' of RE#######
					wholeFragments_5_prime.GRanges<-GRanges(seqnames=seqnames(wholeFragments.GRanges),
					                                        IRanges(start=start(wholeFragments.GRanges),
					                                                end=start(wholeFragments.GRanges)+contr.read.length[i]))
					wholeFragments_3_prime.GRanges<-GRanges(seqnames=seqnames(wholeFragments.GRanges),
					                                        IRanges(start=end(wholeFragments.GRanges)-contr.read.length[i],
					                                                end=end(wholeFragments.GRanges)))
					
					readsPerFragment_5 <- countOverlaps(wholeFragments_5_prime.GRanges,subject=contr.GRanges)
					readsPerFragment_3 <- countOverlaps(wholeFragments_3_prime.GRanges,subject=contr.GRanges)
					combined_reads<-readsPerFragment_5+readsPerFragment_3
					reads<-data.frame(nReads=combined_reads)
					fileName<-gsub(".bam","",contrFiles[i])
					colnames(reads)[1]<-fileName
					readCountTable<-cbind(readCountTable,reads)
				}
				#########Filter out fragment with no read###########################
				n<-ncol(readCountTable)
				readCountTable.f<-readCountTable[rowSums(readCountTable[,4:n]) >0,]
				readCountTable.GRanges<-GRanges(seqnames=readCountTable.f$chromosome,
				                                IRanges(start=readCountTable.f$start,end=readCountTable.f$end))
				values(readCountTable.GRanges)<-readCountTable.f[,4:n]
				readCountTable(object)<-readCountTable.GRanges
				assign(objName,object,envir=parent.frame())
				invisible(1)
				print("The countTable is created!!!.")
			}
				
	  }
)

setGeneric(
		name="getBatchReadCountPerWindow",
		def=function(object,windowSize=5e3,nFragmentExcludedReadsNearViewpoint=2,mode=c("non-overlapping", "overlapping")){
			standardGeneric("getBatchReadCountPerWindow")
		}
)

setMethod("getBatchReadCountPerWindow",
		signature(object = "r3CseqInBatch"),
		function (object,windowSize,nFragmentExcludedReadsNearViewpoint,mode){
			objName <- deparse(substitute(object))
			if(!is(object,"r3CseqInBatch")){
				stop("Need to initialize the r3CseqInBatch object")
			}
			######check the number of fragment to exclude reads around the viewpoint###
			nExcludedFragments <- nFragmentExcludedReadsNearViewpoint
			if(nExcludedFragments <0 | nExcludedFragments>5){
				stop("The number of fragment that uses for removing reads around the viewpoint must be 0-5 restriction fragments.")
			}
			#####Select the read count methods#########
			mode <- match.arg(mode)
			if(mode==""){
				mode<-"non-overlapping"
			}
			######check files########################
			objName <- deparse(substitute(object))
			bamDir<-bamFilesDirectory(object)
			#########checking directory and files####	
			if (file.exists(bamDir)){
				setwd(file.path(bamDir))
			}else{
				stop(paste("Couldn't find path","-->",bamDir))
			}
			
			exp.bam.files    <-BamExpFiles(object)
			contr.bam.files  <-BamContrFiles(object)
			
			for(i.file in c(exp.bam.files)){
				i.file<-gsub(".bam","",i.file)
				i.file<-paste(i.file,".rData",sep="")
				if(file.exists(i.file) ==FALSE){
					stop(paste("Couldn't find file","-->",i.file, ". Please run getBatchRawReads function"))
				}		
			}
			for(i.file in c(contr.bam.files)){
				i.file<-gsub(".bam","",i.file)
				i.file<-paste(i.file,".rData",sep="")
				if(file.exists(i.file) ==FALSE){
					stop(paste("Couldn't find file","-->",i.file, ". Please run getBatchRawReads function"))
				}		
			}		
			##########get labels#######################
			expFiles<-BamExpFiles(object)
			contrFiles<-BamContrFiles(object)
			###########################################
			wholeFragments.GRanges<-getFragmentsPerWindow(object,windowSize,mode)
			###########Make count table for each bam file in the experiments########
				readCountTable<-data.frame(chromosome=seqnames(wholeFragments.GRanges),
				                           start=start(wholeFragments.GRanges),
				                           end=end(wholeFragments.GRanges))
				for(i in 1:length(exp.bam.files)){   
					print(paste("start counting reads in...",exp.bam.files[i], "....file in the experiment"))
					fileName<-gsub(".bam","",expFiles[i])
					fileName<-paste(fileName,".rData",sep="")
					load(fileName)
					exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
					readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=exp.GRanges)
					reads<-data.frame(nReads=readsPerFragment)
					Name<-gsub(".bam","",expFiles[i])
					colnames(reads)[1]<-Name
					readCountTable<-cbind(readCountTable,reads)
				}
				###########Make count table for each bam file in the controls########
				for(i in 1:length(contr.bam.files)){	    
					print(paste("start counting reads in...",contr.bam.files[i], "....file in the control"))
					fileName<-gsub(".bam","",contrFiles[i])
					fileName<-paste(fileName,".rData",sep="")
					load(fileName)
					contr.GRanges<-excludeReadsNearViewpoint(object,contr.GRanges,nExcludedFragments)
					readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=contr.GRanges)
					reads<-data.frame(nReads=readsPerFragment)
					Name<-gsub(".bam","",contrFiles[i])
					colnames(reads)[1]<-Name
					readCountTable<-cbind(readCountTable,reads)
				}
				#########Filter out fragment with no read###########################
				n<-ncol(readCountTable)
				readCountTable.f<-readCountTable[rowSums(readCountTable[,4:n]) >0,]
				readCountTable.GRanges<-GRanges(seqnames=readCountTable.f$chromosome,
				                                IRanges(start=readCountTable.f$start,
				                                        end=readCountTable.f$end))
				values(readCountTable.GRanges)<-readCountTable.f[,4:n]
        readCountTable(object)<-readCountTable.GRanges
				assign(objName,object,envir=parent.frame())
				invisible(1)
				print("The countTable is created!!!.")		
		}
)

setGeneric(
		name="calculateBatchRPM",
		def=function(object,normalized_method=c("powerlawFittedRPM","normalRPM")){
			standardGeneric("calculateBatchRPM")
		}
)
setMethod("calculateBatchRPM",
		signature(object = "r3CseqInBatch"),
		function (object,normalized_method){
			
			if(!is(object,"r3CseqInBatch")){
				stop("Need the r3CseqInBatch object")
			}
			#####Select the method#########
			selected_methods <- match.arg(normalized_method)
			if(selected_methods==""){
				selected_methods<-"powerlawFittedRPM"
			}
			##########get labels#######################
			expFiles<-BamExpFiles(object)
			contrFiles<-BamContrFiles(object)
			##########Get library size#################
			expLibs<-expBatchLibrarySize(object)
			contrLibs<-contrBatchLibrarySize(object)
			###########################################
			objName <- deparse(substitute(object))
			count.data<-readCountTable(object)
			if(length(count.data)==0){
				stop("No reads found in the countTable object, please run getBatchReadCountPerRestrictionFragment or getBatchReadCountPerWindow function!!")
			}
			if(selected_methods=="powerlawFittedRPM"){

				RPMsTable<-data.frame(chromosome=seqnames(count.data),start=start(count.data),end=end(count.data))
				for(i in 1:length(expFiles)){
				  
					print(paste("Calculating RPMs reads in...",expFiles[i], "....file in the experiment"))
					fileName<-gsub(".bam","",expFiles[i])
					count.data.f<-mcols(count.data)
					count.i<-count.data.f[,colnames(count.data.f)==fileName]
					values<-count.i
					coef.i<-getPowerLawFittedCoeficient(values)
					RPMs<-data.frame(RPMs=powerLawFittedRPM(values,coef.i)) 
					colnames(RPMs)[1]<-fileName
					RPMsTable<-cbind(RPMsTable,RPMs)
				}	
				for(i in 1:length(contrFiles)){
				  
					print(paste("Calculating RPMs reads in...",contrFiles[i], "....file in the control"))
					fileName<-gsub(".bam","",contrFiles[i])
					count.data.f<-mcols(count.data)
					count.i<-count.data.f[,colnames(count.data.f)==fileName]
					values<-count.i
					coef.i<-getPowerLawFittedCoeficient(values)
					RPMs<-data.frame(RPMs=powerLawFittedRPM(values,coef.i)) 
					colnames(RPMs)[1]<-fileName
					RPMsTable<-cbind(RPMsTable,RPMs)
				}
				####################################
				n<-ncol(RPMsTable)
				RPMsTable.GRanges<-GRanges(seqnames=RPMsTable$chromosome,
				                           IRanges(start=RPMsTable$start,end=RPMsTable$end))
				
				values(RPMsTable.GRanges)<-RPMsTable[,4:n]
				RPMsTable(object)<-RPMsTable.GRanges
				assign(objName,object,envir=parent.frame())
				invisible(1)
				print("The RPMs table is created!!!.")		
				
			}
			if(selected_methods=="normalRPM"){
				RPMsTable<-data.frame(chromosome=seqnames(count.data),start=start(count.data),end=end(count.data))
				for(i in 1:length(expFiles)){
					print(paste("Calculating RPMs reads in...",expFiles[i], "....file in the experiment"))
					fileName<-gsub(".bam","",expFiles[i])
					count.data.f<-mcols(count.data)
					count.i<-count.data.f[,colnames(count.data.f)==fileName]
					values<-count.i
					RPMs<-data.frame(RPMs=normalcalRPM(values,expLibs[i])) 
					colnames(RPMs)[1]<-fileName
					RPMsTable<-cbind(RPMsTable,RPMs)
				}	
				for(i in 1:length(contrFiles)){
					print(paste("Calculating RPMs reads in...",contrFiles[i], "....file in the control"))
					fileName<-gsub(".bam","",contrFiles[i])
					count.data.f<-mcols(count.data)
					count.i<-count.data.f[,colnames(count.data.f)==fileName]
					values<-count.i
					RPMs<-data.frame(RPMs=normalcalRPM(values,contrLibs[i]))
					colnames(RPMs)[1]<-fileName
					RPMsTable<-cbind(RPMsTable,RPMs)
				}
				####################################
				n<-ncol(RPMsTable)
				RPMsTable.GRanges<-GRanges(seqnames=RPMsTable$chromosome,IRanges(start=RPMsTable$start,end=RPMsTable$end))
				values(RPMsTable.GRanges)<-RPMsTable[,4:n]
				RPMsTable(object)<-RPMsTable.GRanges
				assign(objName,object,envir=parent.frame())
				invisible(1)
				print("The RPMs table is created!!!.")		
			}
		}
)

setGeneric(
		name="getBatchInteractions",
		def=function(object,method=c("union","intersection"),smoothing.parameter=0.1,fdr=0.05){
			standardGeneric("getBatchInteractions")
		}
)
setMethod("getBatchInteractions",
		signature(object = "r3CseqInBatch"),
		function (object,method,smoothing.parameter,fdr){
			if(!is(object,"r3CseqInBatch")){
				stop("Need the r3CseqInBatch object")
			}
			#####Select the method#########
			selected_methods <- match.arg(method)
			if(selected_methods==""){
				selected_methods<-"union"
			}
			if(smoothing.parameter <0.06 | smoothing.parameter>0.4){
				stop("The smoothing parameter must be 0.06-0.4 (default=0.1).")
			}
			##########################################
			objName <- deparse(substitute(object))
			##########get labels#######################
			expFiles<-BamExpFiles(object)
			contrFiles<-BamContrFiles(object)
			##########get read count##################
			count.data<-readCountTable(object)
			if(length(count.data)==0){
				stop("No reads found in the countTable object, please run getBatchReadCountPerRestrictionFragment or getBatchReadCountPerWindow function!!")
			}
			##########get RPM##################
			rpm.data<-RPMsTable(object)
			if(length(rpm.data)==0){
				stop("No RPMs found, please run calculateBatchRPM function!!")
			}
			#############get interaction in experiments###########
			expInteraction<-data.frame(chromosome=seqnames(rpm.data),start=start(rpm.data),end=end(rpm.data))
			for(i in 1:length(expFiles)){
			
				print(paste("Analyzing interaction regions ...",expFiles[i], ".... in the experiment"))
				fileName<-gsub(".bam","",expFiles[i])
				count.data.f<-mcols(count.data)
				count.i<-count.data.f[,colnames(count.data.f)==fileName]
				
				rpm.data.f<-mcols(rpm.data)
				rpm.i<-rpm.data.f[,colnames(rpm.data.f)==fileName]
				
				count.i<-GRanges(seqnames = seqnames(count.data),IRanges(start=start(count.data),
				                                                               end=end(count.data)),nReads=count.i,RPMs=rpm.i)
				
				count.i<-count.i[count.i$RPMs>1,]
				expInt.i<-assign3CseqSigContact(object,count.i,smoothing.parameter,fdr)	
				
				expInt.i.frame<-data.frame(chromosome=seqnames(expInt.i),
				                           start=start(expInt.i),
				                           end=end(expInt.i),
				                           nReads=expInt.i$nReads,
				                           RPMs=expInt.i$RPMs,
				                           p.value=expInt.i$p.value,
				                           q.value=expInt.i$q.value)
				
				expInteraction<-merge(expInteraction,expInt.i.frame,by.x=c("chromosome","start","end"),by.y=c("chromosome","start","end"),all=T)
				expInt.i.frame<-data.frame(chromosome=seqnames(expInt.i),
				                           start=start(expInt.i),
				                           end=end(expInt.i),
				                           nReads=expInt.i$nReads,
				                           RPMs=expInt.i$RPMs,
				                           p.value=expInt.i$p.value,
				                           q.value=expInt.i$q.value)
				write.table(expInt.i.frame,file=paste(fileName,".interaction.txt",sep=""),append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)					
				
			}
			n.col<-ncol(expInteraction)
			expInteraction.filtered<-expInteraction[rowSums(is.na(expInteraction[,4:n.col])) !=n.col-3,]
			#############get interaction in controls###############
			contrInteraction<-data.frame(chromosome=seqnames(rpm.data),start=start(rpm.data),end=end(rpm.data))
			
			for(i in 1:length(contrFiles)){
				print(paste("Analyzing interaction regions ...",contrFiles[i], ".... in the control"))
				fileName<-gsub(".bam","",contrFiles[i])
				count.data.f<-mcols(count.data)
				
				count.i<-count.data.f[,colnames(count.data.f)==fileName]
				
				rpm.data.f<-mcols(rpm.data)
				rpm.i<-rpm.data.f[,colnames(rpm.data.f)==fileName]
				

				count.i<-GRanges(seqnames = seqnames(count.data),IRanges(start=start(count.data),
				                                                         end=end(count.data)),nReads=count.i,RPMs=rpm.i)
				count.i<-count.i[count.i$RPMs>1,]
				contrInt.i<-assign3CseqSigContact(object,count.i,smoothing.parameter,fdr)	
				
				contrInt.i.frame<-data.frame(chromosome=seqnames(contrInt.i),
				                             start=start(contrInt.i),
				                             end=end(contrInt.i),
				                             nReads=contrInt.i$nReads,
				                             RPMs=contrInt.i$RPMs,
				                             p.value=contrInt.i$p.value,
				                             q.value=contrInt.i$q.value)
				
				contrInteraction<-merge(contrInteraction,contrInt.i.frame,by.x=c("chromosome","start","end"),by.y=c("chromosome","start","end"),all=T)
				contrInt.i.frame<-data.frame(chromosome=seqnames(contrInt.i),
				                             start=start(contrInt.i),
				                             end=end(contrInt.i),
				                             nReads=contrInt.i$nReads,
				                             RPMs=contrInt.i$RPMs,
				                             p.value=contrInt.i$p.value,
				                             q.value=contrInt.i$q.value)
				write.table(contrInt.i.frame,file=paste(fileName,".interaction.txt",sep=""),append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)					
				
			}
			m.col<-ncol(contrInteraction)
			contrInteraction.filtered<-contrInteraction[rowSums(is.na(contrInteraction[,4:m.col])) !=m.col-3,]
			
			if(selected_methods=="intersection"){
				##############################################
				exp.intersec<-na.omit(expInteraction.filtered)
				exp.readCounts<-exp.intersec[,c(seq(from=4,to=n.col,by=4))]
				exp.RPMs<-exp.intersec[,c(seq(from=5,to=n.col,by=4))]
				exp.pvalue<-exp.intersec[,c(seq(from=6,to=n.col,by=4))]
				exp.nReads.mean<-round(rowMeans(exp.readCounts))
				exp.RPMs.mean<-rowMeans(exp.RPMs)	
				exp.pvalue[exp.pvalue == 0] <- 1e-20
				exp.fisher.pvalue<-apply(exp.pvalue,1,fishersMethod)
				exp.fisher.qvalue<-qvalue(exp.fisher.pvalue, fdr.level=fdr, pi0.method="bootstrap")$qvalues
				exp.intersec.interaction<-GRanges(seqnames=exp.intersec$chromosome,
				                                  IRanges(start=exp.intersec$start,end=exp.intersec$end),
				                                  nReads=exp.nReads.mean,RPMs=exp.RPMs.mean,p.value=exp.fisher.pvalue,q.value=exp.fisher.qvalue)
				##############################################
				contr.intersec<-na.omit(contrInteraction.filtered)
				contr.readCounts<-contr.intersec[,c(seq(from=4,to=n.col,by=4))]
				contr.RPMs<-contr.intersec[,c(seq(from=5,to=n.col,by=4))]
				contr.pvalue<-contr.intersec[,c(seq(from=6,to=n.col,by=4))]
				contr.nReads.mean<-round(rowMeans(contr.readCounts))
				contr.RPMs.mean<-rowMeans(contr.RPMs)
				contr.pvalue[contr.pvalue == 0] <- 1e-20
				contr.fisher.pvalue<-apply(contr.pvalue,1,fishersMethod)
				contr.fisher.qvalue<-qvalue(contr.fisher.pvalue, fdr.level=fdr, pi0.method="bootstrap")$qvalues
				contr.intersec.interaction<-GRanges(seqnames=contr.intersec$chromosome,
				                                    IRanges(start=contr.intersec$start,end=contr.intersec$end),
				                                    nReads=contr.nReads.mean,RPMs=contr.RPMs.mean,
				                                    p.value=contr.fisher.pvalue,q.value=contr.fisher.qvalue)
				##############################################
				expInteractionRegions(object)<-exp.intersec.interaction	
				contrInteractionRegions(object)<-contr.intersec.interaction	
				assign(objName,object,envir=parent.frame())
				invisible(1)
				print("Interactions from the intersection method are created!!!.")		
			}
			if(selected_methods=="union"){
				##############################################
				exp.union<-expInteraction.filtered
				exp.readCounts<-exp.union[,c(seq(from=4,to=n.col,by=4))]
				exp.readCounts[is.na(exp.readCounts)==TRUE]<-0
				exp.RPMs<-exp.union[,c(seq(from=5,to=n.col,by=4))]
				exp.RPMs[is.na(exp.RPMs)==TRUE]<-0
				exp.pvalue<-exp.union[,c(seq(from=6,to=n.col,by=4))]
				exp.pvalue[is.na(exp.pvalue)==TRUE]<-1
				exp.nReads.mean<-round(rowMeans(exp.readCounts))
				exp.RPMs.mean<-rowMeans(exp.RPMs)	
				
				exp.pvalue[exp.pvalue == 0] <- 1e-20
				exp.fisher.pvalue<-apply(exp.pvalue,1,fishersMethod)
				exp.fisher.qvalue<-qvalue(exp.fisher.pvalue, fdr.level=fdr, pi0.method="bootstrap")$qvalues
				
				exp.union.interaction<-GRanges(seqnames=exp.union$chromosome,
				                               IRanges(start=exp.union$start,end=exp.union$end),
				                               nReads=exp.nReads.mean,RPMs=exp.RPMs.mean,
				                               p.value=exp.fisher.pvalue,q.value=exp.fisher.qvalue)
				##############################################
				contr.union<-contrInteraction.filtered
				contr.readCounts<-contr.union[,c(seq(from=4,to=n.col,by=4))]
				contr.readCounts[is.na(contr.readCounts)==TRUE]<-0
				contr.RPMs<-contr.union[,c(seq(from=5,to=n.col,by=4))]
				contr.RPMs[is.na(contr.RPMs)==TRUE]<-0
				contr.pvalue<-contr.union[,c(seq(from=6,to=n.col,by=4))]
				contr.pvalue[is.na(contr.pvalue)==TRUE]<-1
				contr.nReads.mean<-round(rowMeans(contr.readCounts))
				contr.RPMs.mean<-rowMeans(contr.RPMs)
				
				contr.pvalue[contr.pvalue == 0] <- 1e-20
				contr.fisher.pvalue<-apply(contr.pvalue,1,fishersMethod)
				contr.fisher.qvalue<-qvalue(contr.fisher.pvalue, fdr.level=fdr, pi0.method="bootstrap")$qvalues
				contr.union.interaction<-GRanges(seqnames=contr.union$chromosome,
				                                 IRanges(start=contr.union$start,end=contr.union$end),
				                                 nReads=contr.nReads.mean,RPMs=contr.RPMs.mean,
				                                 p.value=contr.fisher.pvalue,q.value=contr.fisher.qvalue)
				##############################################
				expInteractionRegions(object)<-exp.union.interaction	
				contrInteractionRegions(object)<-contr.union.interaction	
				assign(objName,object,envir=parent.frame())
				invisible(1)
				print("Interactions from the union method are created!!!.")		
				
			}
		}
)

