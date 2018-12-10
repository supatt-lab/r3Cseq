# TODO: These following functions were implemented to get interaction regions.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com
#####
#####The functions below are completely migrated to BioC 3.9 on R 3.6
#####
setGeneric(
		name="getCoverage",
		def=function(object){
			standardGeneric("getCoverage")
		}
)
setMethod("getCoverage",
		signature(object = "r3Cseq"),
		function (object){
		stop("This method has been removed!!.")
		}
)
#####
setGeneric(
		name="getRawReads",
		def=function(object){
			standardGeneric("getRawReads")
		}
)
setMethod("getRawReads",
		signature(object = "r3Cseq"),
		function (object){
			if(!is(object,"r3Cseq")){
				stop("Need to initialize the r3Cseq object")
			}
			objName <- deparse(substitute(object))
				
			if(isControlInvolved(object)==TRUE){
					
					exp.bam.file    <-alignedReadsBamExpFile(object)
					contr.bam.file  <-alignedReadsBamContrFile(object)
					
					if(file.exists(exp.bam.file) ==FALSE){
						stop(paste("Couldn't find file","-->",exp.bam.file))
					}			
					if(file.exists(contr.bam.file) ==FALSE){
						stop(paste("Couldn't find file","-->",contr.bam.file))
					}	
					
					expLabeled<-expLabel(object)
					controlLabeled<-contrLabel(object)
					
					if(nchar(expLabeled)==0 & nchar(controlLabeled)==0){
						stop("The package requires input of expLabel and contrLabel.")
					}
					
					print("start reading in ......")
										
					what    <- c("rname", "strand", "pos", "qwidth","seq")
					param   <- ScanBamParam(what=what,flag = scanBamFlag(isUnmappedQuery = FALSE))
					exp.bam <- scanBam(exp.bam.file, param=param)
					contr.bam <- scanBam(contr.bam.file, param=param)
					
					exp.read.length <-width(exp.bam[[1]]$seq[1])
					exp.GRanges<-GRanges(seqnames=as.vector(exp.bam[[1]]$rname),IRanges(start=exp.bam[[1]]$pos,width=exp.read.length),strand=exp.bam[[1]]$strand)	

					contr.read.length <-width(contr.bam[[1]]$seq[1])
					contr.GRanges<-GRanges(seqnames=as.vector(contr.bam[[1]]$rname),IRanges(start=contr.bam[[1]]$pos,width=contr.read.length),strand=contr.bam[[1]]$strand)	
		
					expLibrarySize(object)   <-length(exp.GRanges)
					contrLibrarySize(object) <-length(contr.GRanges)
				
					#######Remove chr.._random chromosomes############################################
					#exp.GRanges<-exp.GRanges[seqnames(exp.GRanges) %in% paste('chr',c(seq(1,50),'X','Y'),sep='')]
					#contr.GRanges<-contr.GRanges[seqnames(contr.GRanges) %in% paste('chr',c(seq(1,50),'X','Y'),sep='')]
					#######assign raw read to the object########
					expReadLength(object)<-exp.read.length
					contrReadLength(object)<-contr.read.length
					expRawData(object)<-exp.GRanges
					contrRawData(object)<-contr.GRanges
					
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("Raw reads processing is done.")
				}
				if(isControlInvolved(object)==FALSE){
					exp.bam.file    <-alignedReadsBamExpFile(object)
					if(file.exists(exp.bam.file) ==FALSE){
						stop(paste("Couldn't find file","-->",exp.bam.file))
					}			
					expLabeled<-expLabel(object)	
					if(nchar(expLabeled)==0){
						stop("The package requires input of expLabel.")
					}
					print("start reading in ......")
					
					what    <- c("rname", "strand", "pos", "qwidth","seq")
					param   <- ScanBamParam(what=what,flag = scanBamFlag(isUnmappedQuery = FALSE))
					exp.bam <- scanBam(exp.bam.file, param=param)
					exp.read.length <-width(exp.bam[[1]]$seq[1])
					exp.GRanges<-GRanges(seqnames=as.vector(exp.bam[[1]]$rname),IRanges(start=exp.bam[[1]]$pos,width=exp.read.length),strand=exp.bam[[1]]$strand)	
					
					#######Remove chr.._random chromosomes############################################
					#exp.GRanges<-exp.GRanges[seqnames(exp.GRanges) %in% paste('chr',c(seq(1,50),'X','Y'),sep='')]
					##########assign raw read to the object######
					expLibrarySize(object)   <-length(exp.GRanges)	
					expReadLength(object)<-exp.read.length
					expRawData(object)<-exp.GRanges
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("Raw reads processing is done.")
					
				}
			
		}
)
########
setGeneric(
		name="getReadCountPerRestrictionFragment",
		def=function(object,getReadsMethod = c("wholeReads", "adjacentFragmentEndsReads"),
			nFragmentExcludedReadsNearViewpoint=2){
			standardGeneric("getReadCountPerRestrictionFragment")
		}
)
setMethod("getReadCountPerRestrictionFragment",
		signature(object = "r3Cseq"),
		function (object,getReadsMethod,nFragmentExcludedReadsNearViewpoint){
			objName <- deparse(substitute(object))
			if(!is(object,"r3Cseq")){
				stop("Need to initialize the r3Cseq object")
			}
			#####Select the read count methods#########
			selected_methods <- match.arg(getReadsMethod)
			if(selected_methods==""){
				selected_methods<-"wholeReads"
			}
			######check the number of fragment to exclude reads around the viewpoint###
			nExcludedFragments <- nFragmentExcludedReadsNearViewpoint
			if(nExcludedFragments <0 | nExcludedFragments>5){
				stop("The number of fragment that uses for removing reads around the viewpoint must be 0-5 restriction fragments.")
			}
			if(length(expRawData(object))==0){
				stop("No reads found, please run getRawReads function to process BAM file.")
			}
			if(isControlInvolved(object)==FALSE){
				if(selected_methods=="wholeReads"){
					####Get organism#######################
					orgName  <- organismName(object)
					####Initialized Restriction Enzyme Object##
					enzymeDb	<-new("repbaseEnzyme")
					resEnzyme 	<-restrictionEnzyme(object)	
					####Build GRanges for restriction site#####
					print(paste("Fragmenting genome by....",resEnzyme,"....",sep=""))
					wholeFragments<-getWholeGenomeRestrictionFragments(enzymeDb,resEnzyme,orgName)
					wholeFragments.GRanges<-GRanges(seqnames=as.character(wholeFragments$chromosome),IRanges(start=wholeFragments$start,end=wholeFragments$end))
					##########Count all reads located inside the fragments######
					exp.GRanges<-expRawData(object)
					exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
					readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=exp.GRanges)
					reads<-data.frame(nReads=readsPerFragment)
					####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
					values(wholeFragments.GRanges)=reads
					countTable<-wholeFragments.GRanges
					countTable.filtered<-countTable[countTable$nReads >0,]
					expReadCount(object)<-countTable.filtered
					print("Count processing is done.")
				}
				if(selected_methods=="adjacentFragmentEndsReads"){
					####Get organism#######################
					orgName  <- organismName(object)
					exp.read.length<-expReadLength(object)
					####Initialized Restriction Enzyme Object##
					enzymeDb	<-new("repbaseEnzyme")
					resEnzyme 	<-restrictionEnzyme(object)	
					####Build GRanges for restriction site#####
					print(paste("Fragmenting genome by....",resEnzyme,"....",sep=""))
					wholeFragments<-getWholeGenomeRestrictionFragments(enzymeDb,resEnzyme,orgName)
					wholeFragments.GRanges<-GRanges(seqnames=as.character(wholeFragments$chromosome),IRanges(start=wholeFragments$start,end=wholeFragments$end))
					
					########Make GRanges for the 5' and 3' of RE#######
					exp.GRanges<-expRawData(object)
					exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
					wholeFragments_5_prime.GRanges<-GRanges(seqnames =seqnames(wholeFragments.GRanges),
					                                        IRanges(start=start(wholeFragments.GRanges),
					                                        end=start(wholeFragments.GRanges)+exp.read.length))
					
					wholeFragments_3_prime.GRanges<-GRanges(seqnames =seqnames(wholeFragments.GRanges),
					                                        IRanges(start=end(wholeFragments.GRanges)-exp.read.length,
					                                        end=end(wholeFragments.GRanges)))
					
					readsPerFragment_5 <- countOverlaps(wholeFragments_5_prime.GRanges,subject=exp.GRanges)
					readsPerFragment_3 <- countOverlaps(wholeFragments_3_prime.GRanges,subject=exp.GRanges)
					combined_reads<-readsPerFragment_5+readsPerFragment_3
					reads<-data.frame(nReads=combined_reads)
					
					####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
					values(wholeFragments.GRanges)=reads
					countTable<-wholeFragments.GRanges
					countTable.filtered<-countTable[countTable$nReads >0,]
					expReadCount(object)<-countTable.filtered
					print("Count processing is done.")
				}
				
				assign(objName,object,envir=parent.frame())
				invisible(1)
			}
			if(isControlInvolved(object)==TRUE){
				if(selected_methods=="wholeReads"){
					####Get organism#######################
					orgName  <- organismName(object)
					####Initialized Restriction Enzyme Object##
					enzymeDb	<-new("repbaseEnzyme")
					resEnzyme 	<-restrictionEnzyme(object)	
					####Build GRanges for restriction site#####
					print(paste("Fragmenting genome by....",resEnzyme,"....",sep=""))
					wholeFragments<-getWholeGenomeRestrictionFragments(enzymeDb,resEnzyme,orgName)
					wholeFragments.GRanges<-GRanges(seqnames =as.character(wholeFragments$chromosome),IRanges(start=wholeFragments$start,end=wholeFragments$end))
					##########Count all reads in the experiment, which are located inside the fragments######
					exp.GRanges<-expRawData(object)
					exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
					readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=exp.GRanges)
					exp.reads<-data.frame(nReads=readsPerFragment)
					
					####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
					values(wholeFragments.GRanges)=exp.reads
					exp.countTable<-wholeFragments.GRanges
					exp.countTable.filtered<-exp.countTable[exp.countTable$nReads >0,]
					expReadCount(object)<-exp.countTable.filtered
					print("Count processing in the experiment is done.")
				  ##########Count all reads in the control, which are located inside the fragments######
					contr.GRanges<-contrRawData(object)
					contr.GRanges<-excludeReadsNearViewpoint(object,contr.GRanges,nExcludedFragments)
					readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=contr.GRanges)
					contr.reads<-data.frame(nReads=readsPerFragment)
					####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
					values(wholeFragments.GRanges)=contr.reads
					contr.countTable<-wholeFragments.GRanges
					contr.countTable.filtered<-contr.countTable[contr.countTable$nReads >0,]
					contrReadCount(object)<-contr.countTable.filtered
					print("Count processing in the control is done.")
					
				}
				if(selected_methods=="adjacentFragmentEndsReads"){
					####Get organism#######################
					orgName  <- organismName(object)
					exp.read.length<-expReadLength(object)
					####Initialized Restriction Enzyme Object##
					enzymeDb	<-new("repbaseEnzyme")
					resEnzyme 	<-restrictionEnzyme(object)	
					####Build GRanges for restriction site#####
					print(paste("Fragmenting genome by....",resEnzyme,"....",sep=""))
					wholeFragments<-getWholeGenomeRestrictionFragments(enzymeDb,resEnzyme,orgName)
					wholeFragments.GRanges<-GRanges(seqnames =as.character(wholeFragments$chromosome),IRanges(start=wholeFragments$start,end=wholeFragments$end))
					
					########Make RangedData for the 5' and 3' of RE#######
					exp.GRanges<-expRawData(object)
					exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
					wholeFragments_5_prime.GRanges<-GRanges(seqnames =seqnames(wholeFragments.GRanges),
					                                        IRanges(start=start(wholeFragments.GRanges),
					                                        end=start(wholeFragments.GRanges)+exp.read.length))
					
					wholeFragments_3_prime.GRanges<-GRanges(seqnames =seqnames(wholeFragments.GRanges),
					                                        IRanges(start=end(wholeFragments.GRanges)-exp.read.length,
					                                        end=end(wholeFragments.GRanges)))
					
					exp.readsPerFragment_5 <- countOverlaps(wholeFragments_5_prime.GRanges,subject=exp.GRanges)
					exp.readsPerFragment_3 <- countOverlaps(wholeFragments_3_prime.GRanges,subject=exp.GRanges)
					exp.combined_reads<-exp.readsPerFragment_5+exp.readsPerFragment_3
					exp.reads<-data.frame(nReads=exp.combined_reads)
					
					####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
					values(wholeFragments.GRanges)=exp.reads
					exp.countTable<-wholeFragments.GRanges
					exp.countTable.filtered<-exp.countTable[exp.countTable$nReads >0,]
					expReadCount(object)<-exp.countTable.filtered
					print("Count processing in the experiment is done.")
					
					contr.GRanges<-contrRawData(object)
					contr.GRanges<-excludeReadsNearViewpoint(object,contr.GRanges,nExcludedFragments)
					contr.readsPerFragment_5 <- countOverlaps(wholeFragments_5_prime.GRanges,subject=contr.GRanges)
					contr.readsPerFragment_3 <- countOverlaps(wholeFragments_3_prime.GRanges,subject=contr.GRanges)
					contr.combined_reads<-contr.readsPerFragment_5+contr.readsPerFragment_3
					contr.reads<-data.frame(nReads=contr.combined_reads)
					####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
					values(wholeFragments.GRanges)=contr.reads
					contr.countTable<-wholeFragments.GRanges
					contr.countTable.filtered<-contr.countTable[contr.countTable$nReads >0,]
					contrReadCount(object)<-contr.countTable.filtered
					print("Count processing in the control is done.")
				}
				assign(objName,object,envir=parent.frame())
				invisible(1)
			}
		}
)

##########
setGeneric(
		name="getReadCountPerWindow",
		def=function(object,windowSize=5e3,nFragmentExcludedReadsNearViewpoint=2,mode=c("non-overlapping")){
			standardGeneric("getReadCountPerWindow")
		}
)

setMethod("getReadCountPerWindow",
		signature(object = "r3Cseq"),
		function (object,windowSize,nFragmentExcludedReadsNearViewpoint,mode){
			objName <- deparse(substitute(object))
			if(!is(object,"r3Cseq")){
				stop("Need to initialize the r3Cseq object")
			}
			######check the number of fragment to exclude reads around the viewpoint###
			nExcludedFragments <- nFragmentExcludedReadsNearViewpoint
			if(nExcludedFragments <0 | nExcludedFragments>5){
				stop("The number of fragment that uses for removing reads around the viewpoint must be 0-5 restriction fragments.")
			}
			if(length(expRawData(object))==0){
				stop("No reads found, please run getRawReads function to process BAM file.")
			}
			#####Select the read count methods#########
			mode <- match.arg(mode)
			if(mode==""){
				mode<-"non-overlapping"
			}
			
			#####Select the read count methods#########			
			if(isControlInvolved(object)==FALSE){
				
				wholeFragments.GRanges<-getFragmentsPerWindow(object,windowSize,mode)
				##########Count all reads located inside the fragments######
				exp.GRanges<-expRawData(object)
				exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
				readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=exp.GRanges)
				reads<-data.frame(nReads=readsPerFragment)
				
				####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
				values(wholeFragments.GRanges)<-reads
				countTable<-wholeFragments.GRanges
				countTable.filtered<-countTable[countTable$nReads >0,]
				expReadCount(object)<-countTable.filtered
				expRPM(object)<-GRanges()
				print("Count processing is done.")
				assign(objName,object,envir=parent.frame())
				invisible(1)
			}
			if(isControlInvolved(object)==TRUE){
				wholeFragments.GRanges<-getFragmentsPerWindow(object,windowSize,mode)
				##########Count all reads located inside the fragments######
				exp.GRanges<-expRawData(object)
				exp.GRanges<-excludeReadsNearViewpoint(object,exp.GRanges,nExcludedFragments)
				exp.readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=exp.GRanges)
				exp.reads<-data.frame(nReads=exp.readsPerFragment)
				
				####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
				values(wholeFragments.GRanges)<-exp.reads
				exp.countTable<-wholeFragments.GRanges
				
				exp.countTable.filtered<-exp.countTable[exp.countTable$nReads >0,]
				expReadCount(object)<-exp.countTable.filtered
				expRPM(object)<-GRanges()
				
				contr.GRanges<-contrRawData(object)
				contr.GRanges<-excludeReadsNearViewpoint(object,contr.GRanges,nExcludedFragments)
				contr.readsPerFragment <- countOverlaps(wholeFragments.GRanges,subject=contr.GRanges)
				contr.reads<-data.frame(nReads=contr.readsPerFragment)
				
				####migrate the countTable to GRanges because BioC3.9 will not support the RagnedData###
				values(wholeFragments.GRanges)<-contr.reads
				contr.countTable<-wholeFragments.GRanges
				contr.countTable.filtered<-contr.countTable[contr.countTable$nReads >0,]
				contrReadCount(object)<-contr.countTable.filtered
				contrRPM(object)<-GRanges()
				
				print("Count processing is done.")
				assign(objName,object,envir=parent.frame())
				invisible(1)
			}
		}
)

##############
setGeneric(
		name="calculateRPM",
		def=function(object,normalized_method=c("powerlawFittedRPM","normalRPM")){
			standardGeneric("calculateRPM")
		}
)
setMethod("calculateRPM",
		signature(object = "r3Cseq"),
		function (object,normalized_method){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			#####Select the method#########
			selected_methods <- match.arg(normalized_method)
			if(selected_methods==""){
				selected_methods<-"powerlawFittedRPM"
			}
			if(selected_methods=="powerlawFittedRPM"){
				if(isControlInvolved(object)==FALSE){
					objName 	 <- deparse(substitute(object))
					expReadCounts   <-expReadCount(object)
					if(length(expReadCounts)>0){
						coef<-getPowerLawFittedCoeficient(expReadCounts$nReads)
						expReadCounts$RPMs<-powerLawFittedRPM(expReadCounts$nReads,coef) 
						expRPM(object) <-expReadCounts
						assign(objName,object,envir=parent.frame())
						invisible(1)
						print(paste("Normal RPM calculation is done."))		
					}else{
						stop("No reads count per regions found in the r3Cseq object.")
					}
				}
				if(isControlInvolved(object)==TRUE){
					expReadCounts   <-expReadCount(object)
					contrReadCounts <-contrReadCount(object)
					if(length(expReadCounts)>0){
						
						objName 	 <- deparse(substitute(object))
						exp.coef	 <- getPowerLawFittedCoeficient(expReadCounts$nReads)
						contr.coef	 <- getPowerLawFittedCoeficient(contrReadCounts$nReads)
						
						expReadCounts$RPMs <- powerLawFittedRPM(expReadCounts$nReads,exp.coef)
						contrReadCounts$RPMs <-powerLawFittedRPM(contrReadCounts$nReads,contr.coef)
						
						expRPM(object) <-expReadCounts
						contrRPM(object) <-contrReadCounts
						
						assign(objName,object,envir=parent.frame())
						invisible(1)
						print(paste("Normal RPM calculation is done."))	
					}else{
						stop("No reads count per regions found in the r3Cseq object.")
					}
				}
			}
			if(selected_methods=="normalRPM"){
				if(isControlInvolved(object)==FALSE){
					objName 	 <- deparse(substitute(object))
					lib.size<-expLibrarySize(object)
					expReadCounts   <-expReadCount(object)
					if(length(expReadCounts)>0){
						expReadCounts$RPMs <- normalcalRPM(expReadCounts$nReads,lib.size)
						expRPM(object) <-expReadCounts
						assign(objName,object,envir=parent.frame())
						invisible(1)
						print(paste("Normal RPM calculation is done."))		
					}else{
						stop("No reads count per regions found in the r3Cseq object.")
					}
				}
				if(isControlInvolved(object)==TRUE){
					expReadCounts   <-expReadCount(object)
					contrReadCounts <-contrReadCount(object)
						if(length(expReadCounts)>0){
							
							objName 	 <- deparse(substitute(object))	
							expLibSize 	 <- expLibrarySize(object)
							contrLibSize <- contrLibrarySize(object)
							expReadCounts<-expReadCount(object)
							contrReadCounts<-contrReadCount(object)
					
							expRPMs   	 <- normalcalRPM(expReadCounts$nReads,expLibSize)
							contrRPMs    <- normalcalRPM(contrReadCounts$nReads,contrLibSize)
					
							expReadCounts$RPMs <- expRPMs
							contrReadCounts$RPMs <-contrRPMs
					
							expRPM(object) <-expReadCounts
							contrRPM(object) <-contrReadCounts
					
							assign(objName,object,envir=parent.frame())
							invisible(1)
							print(paste("Normal RPM calculation is done."))	
						}else{
							stop("No reads count per regions found in the r3Cseq object.")
						}
				}
			}
		}
)

##################
setGeneric(
		name="getInteractions",
		def=function(object,smoothing.parameter=0.1,fdr=0.05){
			standardGeneric("getInteractions")
		}
)
setMethod("getInteractions",
		signature(object = "r3Cseq"),
		function (object,smoothing.parameter,fdr){
			
			if(!is(object,"r3Cseq")){
				stop("Need the r3Cseq object")
			}
			if(smoothing.parameter <0.06 | smoothing.parameter>0.4){
				stop("The smoothing parameter must be 0.06-0.4 (default=0.1).")
			}
			if(isControlInvolved(object)==FALSE){
				expRPM.RangedData <-expRPM(object)
				if(length(expRPM.RangedData)>0){
					objName <- deparse(substitute(object))
					expInteraction<-assign3CseqSigContact(object,expRPM.RangedData,smoothing.parameter,fdr)	
					expInteractionRegions(object) <-expInteraction
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("Calculation is done. Use function 'expInteractionRegions' or 'contrInteractionRegions' to get the result.")
				}else{
					stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")	
				}
			}
			if(isControlInvolved(object)==TRUE){
				expRPM.RangedData   <-expRPM(object)
				contrRPM.RangedData <-contrRPM(object)
				
				if(length(expRPM.RangedData)>0){
					objName <- deparse(substitute(object))
					expInteraction<-assign3CseqSigContact(object,expRPM.RangedData,smoothing.parameter,fdr)	
					expInteractionRegions(object) <-expInteraction
					##calculate in the control##
					contrInteraction<-assign3CseqSigContact(object,contrRPM.RangedData,smoothing.parameter,fdr)	
					contrInteractionRegions(object) <-contrInteraction
					assign(objName,object,envir=parent.frame())
					invisible(1)
					print("Calculation is done. Use function 'expInteractionRegions' or 'contrInteractionRegions' to get the result.")
				}else{
					stop("No RPM found in the r3Cseq object, you have to run the function 'getReadCountPerRestrictionFragment' following by the function 'calculateRPM'.")
				}
			}
		}
)

