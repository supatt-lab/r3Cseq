# TODO: These following functions were implemented to facilitate r3Cseq and r3CseqWithReplocates classes.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com
#####
#####The functions below are completely migrated to BioC 3.9 on R 3.6
#####

getViewpoint<-function (obj){
	stopifnot( is( obj, "r3Cseq" ) |is( obj, "r3CseqInBatch" ) )
	#######get organism name ############
	#######Fixed genome to genome.assembly
	genome.assembly<-organismName(obj)
	#######get forward primer############
	primer_f<-viewpoint_primer_forward(obj)
	primer_r<-viewpoint_primer_reverse(obj)
	#######get viewpoint chromosome######
	viewpoint_chr<-viewpoint_chromosome(obj)
	#####################################
	if('BSgenome.Hsapiens.UCSC.hg19.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Hsapiens.UCSC.hg19.masked,unload=TRUE)
	}
	if('BSgenome.Hsapiens.UCSC.hg18.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Hsapiens.UCSC.hg18.masked,unload=TRUE)
	}
	if('BSgenome.Mmusculus.UCSC.mm9.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Mmusculus.UCSC.mm9.masked,unload=TRUE)
	}
	if('BSgenome.Mmusculus.UCSC.mm10.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Mmusculus.UCSC.mm10.masked,unload=TRUE)
	}
	if('BSgenome.Rnorvegicus.UCSC.rn5.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Rnorvegicus.UCSC.rn5.masked,unload=TRUE)
	}
	
	if(genome.assembly=="hg18"){
		library(BSgenome.Hsapiens.UCSC.hg18.masked)
		genome <- BSgenome.Hsapiens.UCSC.hg18.masked
		hits_f_f<-matchPattern(DNAString(primer_f),genome[[viewpoint_chr]],fixed=FALSE)
		hits_f_r<-matchPattern(reverseComplement(DNAString(primer_f)),genome[[viewpoint_chr]],fixed=FALSE)
		hits_f_positions<-c(start(hits_f_f),end(hits_f_f),start(hits_f_r),end(hits_f_r))
		
		if(length(hits_f_positions)==0){
			stop(paste("Could not find the forward primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input forward primer"))
		}
		
		hits_r_f<-matchPattern(DNAString(primer_r),genome[[viewpoint_chr]],fixed=FALSE)
		hits_r_r<-matchPattern(reverseComplement(DNAString(primer_r)),genome[[viewpoint_chr]],fixed=FALSE)
		
		hits_r_positions<-c(start(hits_r_f),end(hits_r_f),start(hits_r_r),end(hits_r_r))
		
		if(length(hits_r_positions)==0){
			stop(paste("Could not find the reverse primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input reverse primer"))
		}
		
		hit_positions<-c(hits_f_positions,hits_r_positions)
		
		if(length(hit_positions)==4){
			start_p<-min(hit_positions)
			end_p<-max(hit_positions)
			viewpoint<-GRanges(seqnames =as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
			return(viewpoint)
		}
		if(length(hit_positions) > 4){
			
			if(length(hits_r_positions) > 2 & length(hits_f_positions)==2){
				rank.id<-order(abs(hits_r_positions-hits_f_positions))
				hit_positions<-c(hits_f_positions,hits_r_positions[rank.id==1],hits_r_positions[rank.id==2])
				start_p<-min(hit_positions)
				end_p<-max(hit_positions)
				viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
				return(viewpoint)
				
			}
			if(length(hits_f_positions) > 2 & length(hits_r_positions)==2){
				rank.id<-order(abs(hits_f_positions-hits_r_positions))
				hit_positions<-c(hits_f_positions[rank.id==1],hits_f_positions[rank.id==2],hits_r_positions)
				start_p<-min(hit_positions)
				end_p<-max(hit_positions)
				viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
				return(viewpoint)
			}
			if(length(hits_f_positions) > 2 & length(hits_r_positions) > 2){
				stop(paste("The position of your input primers are not unique. We detected multiple positions of your input primers!!!"))
			}
		}
	}else if(genome.assembly=="hg19"){
		library(BSgenome.Hsapiens.UCSC.hg19.masked)
		genome <- BSgenome.Hsapiens.UCSC.hg19.masked
		hits_f_f<-matchPattern(DNAString(primer_f),genome[[viewpoint_chr]],fixed=FALSE)
		hits_f_r<-matchPattern(reverseComplement(DNAString(primer_f)),genome[[viewpoint_chr]],fixed=FALSE)
		hits_f_positions<-c(start(hits_f_f),end(hits_f_f),start(hits_f_r),end(hits_f_r))
		
		if(length(hits_f_positions)==0){
			stop(paste("Could not find the forward primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input forward primer"))
		}
		
		hits_r_f<-matchPattern(DNAString(primer_r),genome[[viewpoint_chr]],fixed=FALSE)
		hits_r_r<-matchPattern(reverseComplement(DNAString(primer_r)),genome[[viewpoint_chr]],fixed=FALSE)
		
		hits_r_positions<-c(start(hits_r_f),end(hits_r_f),start(hits_r_r),end(hits_r_r))
		
		if(length(hits_r_positions)==0){
			stop(paste("Could not find the reverse primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input reverse primer"))
		}
		
		hit_positions<-c(hits_f_positions,hits_r_positions)
		
		if(length(hit_positions)==4){
			start_p<-min(hit_positions)
			end_p<-max(hit_positions)
			viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
			return(viewpoint)
		}
		if(length(hit_positions) > 4){
			
			if(length(hits_r_positions) > 2 & length(hits_f_positions)==2){
				rank.id<-order(abs(hits_r_positions-hits_f_positions))
				hit_positions<-c(hits_f_positions,hits_r_positions[rank.id==1],hits_r_positions[rank.id==2])
				start_p<-min(hit_positions)
				end_p<-max(hit_positions)
				viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
				return(viewpoint)
				
			}
			if(length(hits_f_positions) > 2 & length(hits_r_positions)==2){
				rank.id<-order(abs(hits_f_positions-hits_r_positions))
				hit_positions<-c(hits_f_positions[rank.id==1],hits_f_positions[rank.id==2],hits_r_positions)
				start_p<-min(hit_positions)
				end_p<-max(hit_positions)
				viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
				return(viewpoint)
			}
			if(length(hits_f_positions) > 2 & length(hits_r_positions) > 2){
				stop(paste("The position of your input primers are not unique. We detected multiple positions of your input primers!!!"))
			}
		}
	
	}else if(genome.assembly =="mm9"){
		library(BSgenome.Mmusculus.UCSC.mm9.masked)
		genome <- BSgenome.Mmusculus.UCSC.mm9.masked
		hits_f_f<-matchPattern(DNAString(primer_f),genome[[viewpoint_chr]],fixed=FALSE)
		hits_f_r<-matchPattern(reverseComplement(DNAString(primer_f)),genome[[viewpoint_chr]],fixed=FALSE)
		hits_f_positions<-c(start(hits_f_f),end(hits_f_f),start(hits_f_r),end(hits_f_r))
		
		if(length(hits_f_positions)==0){
			stop(paste("Could not find the forward primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input forward primer"))
		}
		
		hits_r_f<-matchPattern(DNAString(primer_r),genome[[viewpoint_chr]],fixed=FALSE)
		hits_r_r<-matchPattern(reverseComplement(DNAString(primer_r)),genome[[viewpoint_chr]],fixed=FALSE)
		
		hits_r_positions<-c(start(hits_r_f),end(hits_r_f),start(hits_r_r),end(hits_r_r))
		
		if(length(hits_r_positions)==0){
			stop(paste("Could not find the reverse primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input reverse primer"))
		}
		
		hit_positions<-c(hits_f_positions,hits_r_positions)
		
		if(length(hit_positions)==4){
			start_p<-min(hit_positions)
			end_p<-max(hit_positions)
			viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
			return(viewpoint)
		}
		if(length(hit_positions) > 4){
			
			if(length(hits_r_positions) > 2 & length(hits_f_positions)==2){
				rank.id<-order(abs(hits_r_positions-hits_f_positions))
				hit_positions<-c(hits_f_positions,hits_r_positions[rank.id==1],hits_r_positions[rank.id==2])
				start_p<-min(hit_positions)
				end_p<-max(hit_positions)
				viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
				return(viewpoint)
				
			}
			if(length(hits_f_positions) > 2 & length(hits_r_positions)==2){
				rank.id<-order(abs(hits_f_positions-hits_r_positions))
				hit_positions<-c(hits_f_positions[rank.id==1],hits_f_positions[rank.id==2],hits_r_positions)
				start_p<-min(hit_positions)
				end_p<-max(hit_positions)
				viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
				return(viewpoint)
			}
			if(length(hits_f_positions) > 2 & length(hits_r_positions) > 2){
				stop(paste("The position of your input primers are not unique. We detected multiple positions of your input primers!!!"))
			}
		}
	
	}else if(genome.assembly =="mm10"){
			library(BSgenome.Mmusculus.UCSC.mm10.masked)
			genome <- BSgenome.Mmusculus.UCSC.mm10.masked
			hits_f_f<-matchPattern(DNAString(primer_f),genome[[viewpoint_chr]],fixed=FALSE)
			hits_f_r<-matchPattern(reverseComplement(DNAString(primer_f)),genome[[viewpoint_chr]],fixed=FALSE)
			hits_f_positions<-c(start(hits_f_f),end(hits_f_f),start(hits_f_r),end(hits_f_r))
			
			if(length(hits_f_positions)==0){
				stop(paste("Could not find the forward primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input forward primer"))
			}
			
			hits_r_f<-matchPattern(DNAString(primer_r),genome[[viewpoint_chr]],fixed=FALSE)
			hits_r_r<-matchPattern(reverseComplement(DNAString(primer_r)),genome[[viewpoint_chr]],fixed=FALSE)
			
			hits_r_positions<-c(start(hits_r_f),end(hits_r_f),start(hits_r_r),end(hits_r_r))
			
			if(length(hits_r_positions)==0){
				stop(paste("Could not find the reverse primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input reverse primer"))
			}
			
			hit_positions<-c(hits_f_positions,hits_r_positions)
			
			if(length(hit_positions)==4){
				start_p<-min(hit_positions)
				end_p<-max(hit_positions)
				viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
				return(viewpoint)
			}
			if(length(hit_positions) > 4){
				
				if(length(hits_r_positions) > 2 & length(hits_f_positions)==2){
					rank.id<-order(abs(hits_r_positions-hits_f_positions))
					hit_positions<-c(hits_f_positions,hits_r_positions[rank.id==1],hits_r_positions[rank.id==2])
					start_p<-min(hit_positions)
					end_p<-max(hit_positions)
					viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
					return(viewpoint)
					
				}
				if(length(hits_f_positions) > 2 & length(hits_r_positions)==2){
					rank.id<-order(abs(hits_f_positions-hits_r_positions))
					hit_positions<-c(hits_f_positions[rank.id==1],hits_f_positions[rank.id==2],hits_r_positions)
					start_p<-min(hit_positions)
					end_p<-max(hit_positions)
					viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
					return(viewpoint)
				}
				if(length(hits_f_positions) > 2 & length(hits_r_positions) > 2){
					stop(paste("The position of your input primers are not unique. We detected multiple positions of your input primers!!!"))
				}
			}
	}else if(genome.assembly =="rn5"){
				library(BSgenome.Rnorvegicus.UCSC.rn5.masked)
				genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
				hits_f_f<-matchPattern(DNAString(primer_f),genome[[viewpoint_chr]],fixed=FALSE)
				hits_f_r<-matchPattern(reverseComplement(DNAString(primer_f)),genome[[viewpoint_chr]],fixed=FALSE)
				hits_f_positions<-c(start(hits_f_f),end(hits_f_f),start(hits_f_r),end(hits_f_r))
				
				if(length(hits_f_positions)==0){
					stop(paste("Could not find the forward primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input forward primer"))
				}
				
				hits_r_f<-matchPattern(DNAString(primer_r),genome[[viewpoint_chr]],fixed=FALSE)
				hits_r_r<-matchPattern(reverseComplement(DNAString(primer_r)),genome[[viewpoint_chr]],fixed=FALSE)
				
				hits_r_positions<-c(start(hits_r_f),end(hits_r_f),start(hits_r_r),end(hits_r_r))
				
				if(length(hits_r_positions)==0){
					stop(paste("Could not find the reverse primer position in",genome.assembly," :",viewpoint_chr, ". Please check the input reverse primer"))
				}
				
				hit_positions<-c(hits_f_positions,hits_r_positions)
				
				if(length(hit_positions)==4){
					start_p<-min(hit_positions)
					end_p<-max(hit_positions)
					viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
					return(viewpoint)
				}
				if(length(hit_positions) > 4){
					
					if(length(hits_r_positions) > 2 & length(hits_f_positions)==2){
						rank.id<-order(abs(hits_r_positions-hits_f_positions))
						hit_positions<-c(hits_f_positions,hits_r_positions[rank.id==1],hits_r_positions[rank.id==2])
						start_p<-min(hit_positions)
						end_p<-max(hit_positions)
						viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
						return(viewpoint)
						
					}
					if(length(hits_f_positions) > 2 & length(hits_r_positions)==2){
						rank.id<-order(abs(hits_f_positions-hits_r_positions))
						hit_positions<-c(hits_f_positions[rank.id==1],hits_f_positions[rank.id==2],hits_r_positions)
						start_p<-min(hit_positions)
						end_p<-max(hit_positions)
						viewpoint<-GRanges(seqnames=as.character(viewpoint_chr),IRanges(start=start_p,end=end_p))
						return(viewpoint)
					}
					if(length(hits_f_positions) > 2 & length(hits_r_positions) > 2){
						stop(paste("The position of your input primers are not unique. We detected multiple positions of your input primers!!!"))
					}
				}	
	}else{
		stop("Require the selected genome: hg18, hg19, mm9, mm10 or rn5.")
	}
}

excludeReadsNearViewpoint<-function (object,rawReadsGRanges,span_fragment=1){
	if(span_fragment <0 | span_fragment >5){
		stop("Span fragment number from the viewpoint must be in the range of 0-5 fragment.")
	}
	if(length(rawReadsGRanges)==0){
		stop("No reads can be found!!!, please run 'getRawReads' function in order to get raw reads.")
	}
	#####get viewpoint##############
	viewpoint<-getViewpoint(object)
	######get restriction fragment##
	orgName  <- organismName(object)
	####Initialized Restriction Enzyme Object##
	enzymeDb	<-new("repbaseEnzyme")
	resEnzyme 	<-restrictionEnzyme(object)	
	####Build GRanges for restriction site#####
	fragments<-getRestrictionFragments(enzymeDb,resEnzyme,orgName,as.character(seqnames(viewpoint)))
	
	fragments$distance2viewpoint<-abs(fragments$end-start(viewpoint))
	v.index<-which(fragments$distance2viewpoint==min(fragments$distance2viewpoint))
	v.index<-v.index[1]
	f.index<-c((v.index-span_fragment):(v.index+span_fragment))
	selected.f<-fragments[f.index,]
	expected.regions<-GRanges(seqnames =viewpoint_chromosome(object),IRanges(start=min(selected.f$start),end=max(selected.f$end)))
	excluded.index <- subjectHits(findOverlaps(expected.regions,subject=rawReadsGRanges))
	new.GRanges<-rawReadsGRanges[-(excluded.index)]
	return(new.GRanges)
}

getFragmentsPerWindow<-function (obj,windowSize=5e3,mode){
	stopifnot( is( obj, "r3Cseq" ) |is( obj, "r3CseqInBatch" ) )
	if(windowSize <2000 | windowSize >100000){
		stop("Please select window size between 2Kb - 100Kb")
	}
	if(mode==""){
		mode<-"non-overlapping"
	}
	#######get organism name ############
	genome<-organismName(obj)
	#####################################
	viewpoint<-getViewpoint(obj)
	viewpoint.chr<-viewpoint_chromosome(obj)
	#####################################
	if('BSgenome.Hsapiens.UCSC.hg19.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Hsapiens.UCSC.hg19.masked,unload=TRUE)
	}
	if('BSgenome.Hsapiens.UCSC.hg18.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Hsapiens.UCSC.hg18.masked,unload=TRUE)
	}
	if('BSgenome.Mmusculus.UCSC.mm9.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Mmusculus.UCSC.mm9.masked,unload=TRUE)
	}
	if('BSgenome.Mmusculus.UCSC.mm10.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Mmusculus.UCSC.mm10.masked,unload=TRUE)
	}
	if('BSgenome.Rnorvegicus.UCSC.rn5.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Rnorvegicus.UCSC.rn5.masked,unload=TRUE)
	}
	
	if(genome=="hg18"){
		library(BSgenome.Hsapiens.UCSC.hg18.masked)
		genome <- BSgenome.Hsapiens.UCSC.hg18.masked
		fragment<-data.frame()
		for (chr in paste('chr',c(seq(1,22),'X','Y'),sep='')){
			chr.size<-seqlengths(genome)[chr]
			#if(chr==viewpoint.chr){
			#	if(mode=="non-overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#	if(mode=="overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size/2
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		s.starts<-s.starts[-(length(s.starts))]
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		
			#		
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#}else{	
				if(mode=="non-overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
				if(mode=="overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size/2
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=window.size-step.size+1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
			#}
		}
		fragment.GRanges<-GRanges(seqnames =fragment$chromosome,IRanges(start=fragment$start,end=fragment$end))
		return(fragment.GRanges)
	}else if(genome=="hg19"){
		library(BSgenome.Hsapiens.UCSC.hg19.masked)
		genome <- BSgenome.Hsapiens.UCSC.hg19.masked
		fragment<-data.frame()
		for (chr in paste('chr',c(seq(1,22),'X','Y'),sep='')){
			chr.size<-seqlengths(genome)[chr]
				
			#if(chr==viewpoint.chr){
			#	if(mode=="non-overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#	if(mode=="overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size/2
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		s.starts<-s.starts[-(length(s.starts))]
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		
			#		
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#}else{
				if(mode=="non-overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
				if(mode=="overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size/2
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=window.size-step.size+1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
			#}
		}
		fragment.GRanges<-GRanges(seqnames =fragment$chromosome,IRanges(start=fragment$start,end=fragment$end))
		return(fragment.GRanges)
	}else if(genome =="mm9"){
		library(BSgenome.Mmusculus.UCSC.mm9.masked)
		genome <- BSgenome.Mmusculus.UCSC.mm9.masked
		fragment<-data.frame()
		for (chr in paste('chr',c(seq(1,19),'X','Y'),sep='')){
			chr.size<-seqlengths(genome)[chr]
			#if(chr==viewpoint.chr){
			#	if(mode=="non-overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#	if(mode=="overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size/2
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		s.starts<-s.starts[-(length(s.starts))]
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#}else{
				if(mode=="non-overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
				if(mode=="overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size/2
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=window.size-step.size+1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
			#}
		}
		fragment.GRanges<-GRanges(seqnames = fragment$chromosome,IRanges(start=fragment$start,end=fragment$end))
		return(fragment.GRanges)
		
	}else if (genome =="mm10"){
		library(BSgenome.Mmusculus.UCSC.mm10.masked)
		genome <- BSgenome.Mmusculus.UCSC.mm10.masked
		fragment<-data.frame()
		for (chr in paste('chr',c(seq(1,19),'X','Y'),sep='')){
			chr.size<-seqlengths(genome)[chr]
			#if(chr==viewpoint.chr){
			#	if(mode=="non-overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#	if(mode=="overlapping"){
			#		##Now getting the fragment#####		
			#		window.size=windowSize
			#		step.size = window.size/2
			#		start.p<-start(viewpoint)
			#		end.p<-end(viewpoint)
			#		
			#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
			#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
			#		f.starts<-f.starts[-(2)]
			#		
			#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
			#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
			#		s.starts<-s.starts[-(length(s.starts))]
			#		
			#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
			#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
			#		fragment.chr<-rbind(f.fragment,s.fragment)
			#		fragment<-rbind(fragment,fragment.chr)
			#	}
			#}else{
				if(mode=="non-overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
				if(mode=="overlapping"){
					##Now getting the fragment#####		
					window.size=windowSize
					step.size = window.size/2
					ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
					starts <- seq(from=window.size-step.size+1,to=chr.size,by=step.size)
					fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
					fragment<-rbind(fragment,fragment.chr)
				}
			#}
		}
		fragment.GRanges<-GRanges(seqnames =fragment$chromosome,IRanges(start=fragment$start,end=fragment$end))
		return(fragment.GRanges) 
	}else if (genome =="rn5"){
			library(BSgenome.Rnorvegicus.UCSC.rn5.masked)
			genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
			fragment<-data.frame()
			for (chr in paste('chr',c(seq(1,20),'X'),sep='')){
				chr.size<-seqlengths(genome)[chr]
				#if(chr==viewpoint.chr){
				#	if(mode=="non-overlapping"){
				#		##Now getting the fragment#####		
				#		window.size=windowSize
				#		step.size = window.size
				#		start.p<-start(viewpoint)
				#		end.p<-end(viewpoint)
				#		
				#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
				#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
				#		f.starts<-f.starts[-(2)]
				#		
				#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
				#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
				#		
				#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
				#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
				#		fragment.chr<-rbind(f.fragment,s.fragment)
				#		fragment<-rbind(fragment,fragment.chr)
				#	}
				#	if(mode=="overlapping"){
				#		##Now getting the fragment#####		
				#		window.size=windowSize
				#		step.size = window.size/2
				#		start.p<-start(viewpoint)
				#		end.p<-end(viewpoint)
				#		
				#		f.ends<-(rev(seq(-start.p, -window.size,by=step.size)))*-1
				#		f.starts <- c(1,(rev(seq(-start.p+window.size,1,by=step.size)))*-1)
				#		f.starts<-f.starts[-(2)]
				#		
				#		s.ends <- c(seq(from=end.p+window.size,to=chr.size,by=step.size),as.integer(chr.size))
				#		s.starts <- seq(from=end.p,to=chr.size,by=step.size)
				#		s.starts<-s.starts[-(length(s.starts))]
				#		
				#		f.fragment<-data.frame(chromosome=chr,start=f.starts,end=f.ends)
				#		s.fragment<-data.frame(chromosome=chr,start=s.starts,end=s.ends)
				#		fragment.chr<-rbind(f.fragment,s.fragment)
				#		fragment<-rbind(fragment,fragment.chr)
				#	}
				#}else{
					if(mode=="non-overlapping"){
						##Now getting the fragment#####		
						window.size=windowSize
						step.size = window.size
						ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
						starts <- seq(from=1,to=chr.size,by=step.size)
						fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
						fragment<-rbind(fragment,fragment.chr)
					}
					if(mode=="overlapping"){
						##Now getting the fragment#####		
						window.size=windowSize
						step.size = window.size/2
						ends <- c(seq(from=window.size,to=chr.size,by=step.size),as.integer(chr.size))
						starts <- seq(from=window.size-step.size+1,to=chr.size,by=step.size)
						fragment.chr<-data.frame(chromosome=chr,start=starts,end=ends)
						fragment<-rbind(fragment,fragment.chr)
					}
				#}
			}
			fragment.GRanges<-GRanges(seqnames =fragment$chromosome,IRanges(start=fragment$start,end=fragment$end))
			return(fragment.GRanges) 	
	}else{
		stop("Require the selected genome: hg18, hg19, mm9, mm10 or rn5.")
	}
}

normalcalRPM<-function (x, lib.size){
	seqDepth <-10^6
	if(lib.size >=10^6){
		seqDepth <- 10^6
	}else if(lib.size >=10^5 && lib.size< 10^6){
		seqDepth <- 10^5
	}else if(lib.size >=10^4 && lib.size< 10^5 ){
		seqDepth <- 10^4
	}else if(lib.size >=10^3 && lib.size < 10^4){
		seqDepth <-10^3
	}else{
		stop("This is too low sequencing depth!!!!!. We have to stop the analysis.")
	}
	RPM <- x / (lib.size/seqDepth)
	return(RPM)
}

getPowerLawFittedCoeficient <- function(values) {
	values<-values[values>=50]
	v <- data.table(num = 1, nr_reads = as.integer(values))
	v <- v[, sum(num), by = nr_reads]
	setkey(v, nr_reads)
	v$V1 <- rev(cumsum(rev(v$V1)))
	setnames(v, c('nr_reads', 'reverse_cumulative'))
	ranges<-quantile(v$nr_reads, c(.20,0.9998))
	x<-as.vector(round(ranges))
	v <- v[nr_reads >= min(x) & nr_reads <=max(x)]
	lin.m <- lm(log(reverse_cumulative) ~ log(nr_reads), data = v)
	a = coefficients(lin.m)[2]
	b = coefficients(lin.m)[1]
	return(c(a, b))
}

powerLawFittedRPM <- function(values,lin.reg.coef,alpha = 1.35, T = 10^6) {
	a = lin.reg.coef[1]        
	b = lin.reg.coef[2]
	lambda = (T/(zeta(alpha) * exp(b)))^(1/alpha)
	beta = -1 * a/alpha
	values.norm = lambda * (values)^beta
	return(values.norm)
}

assign3CseqSigContact<-function(obj,readcount.ranges,smoothing.parameter,fdr){
	stopifnot( is( obj, "r3Cseq" ) |is( obj, "r3CseqInBatch" ) )
	if(length(readcount.ranges)>0){
		viewpoint<-getViewpoint(obj)
		viewpoint.chr<-viewpoint_chromosome(obj)
		cis.ranges<-readcount.ranges[seqnames(readcount.ranges)==viewpoint.chr,]
		cis<-data.frame(chromosome=seqnames(cis.ranges),
		                start=start(cis.ranges),end=end(cis.ranges),
		                relative.position=start(cis.ranges)-start(viewpoint),
		                observed=cis.ranges$nReads,rpm=cis.ranges$RPMs)
		#######Smooth the data############
		e.smooth<-vsmooth.spline(cis$relative.position,cis$observed,spar=smoothing.parameter)
		e.fitted<-fitted(e.smooth)
		cis$expected<-as.numeric(e.fitted)
		######Calculate z-score and assign p-value in cis
		cis$z<-(cis$observed-cis$expected)/sd(cis$observed-cis$expected)
		cis$p.value<-2*pnorm(-abs(cis$z))
		cis$p.value[cis$z < 0]<-1
		cis$q.value<-qvalue(cis$p.value, fdr.level=fdr, pi0.method="bootstrap")$qvalues
		cis.interaction<-data.frame(chromosome=cis$chromosome,start=cis$start,end=cis$end,nReads=cis$observed,RPMs=cis$rpm,z=cis$z,p.value=cis$p.value,q.value=cis$q.value)
		########Calculate interaction in trans#####
		ref1<-subset(cis,relative.position< -100000)
		ref2<-subset(cis,relative.position> 100000)
		ref_cis<-rbind(ref1,ref2)
		ref_cis<-ref_cis[,c(1:3,5:6)]
		
		ref_cis<-ref_cis[sample(1:nrow(ref_cis),0.5*nrow(ref_cis)),]
		ref.trans.ranges<-readcount.ranges[seqnames(readcount.ranges)!=viewpoint.chr,]
		ref.trans<-data.frame(chromosome=seqnames(ref.trans.ranges),start=start(ref.trans.ranges),end=end(ref.trans.ranges),observed=ref.trans.ranges$nReads,rpm=ref.trans.ranges$RPMs)
		
		trans<-rbind(ref_cis,ref.trans)
		trans$z<-(trans$observed-mean(trans$observed))/sd(trans$observed)
		trans$p.value<-2*pnorm(-abs(trans$z))
		trans$p.value[trans$z< 0]<-1
		trans$q.value<-qvalue(trans$p.value, fdr.level=fdr, pi0.method="bootstrap")$qvalues
		trans<-subset(trans,chromosome !=viewpoint.chr)
		trans.interaction<-data.frame(chromosome=trans$chromosome,start=trans$start,end=trans$end,nReads=trans$observed,RPMs=trans$rpm,z=trans$z,p.value=trans$p.value,q.value=trans$q.value)
		#######combined bothcis and trans results########
		all.int<-rbind(cis.interaction,trans.interaction)
		all.GRangedData<-GRanges(seqnames =all.int$chromosome,
		                         IRanges(start=all.int$start,end=all.int$end),
		                         nReads=all.int$nReads,
		                         RPMs=all.int$RPMs,
		                         z=all.int$z,
		                         p.value=all.int$p.value,
		                         q.value=all.int$q.value)
		###Return RangedData object#######
		return(all.GRangedData)
		
	}else{
		stop("No reads found, please reads the r3Cseq manual to create the read counts!!!.")
	}
}

get3CseqRefGene<-function(obj){
	stopifnot( is( obj, "r3Cseq" ) |is(obj, "r3CseqInBatch" ) )
	####Get organism#######################
	genome  <- organismName(obj)	
	if(genome=="hg18"){
		refGene.file<-system.file("data/hg18refGene.rdata", package="r3Cseq")
		if(file.exists(refGene.file)==TRUE){
			load(file=refGene.file)
			refGene<-hg18refGene
			return(refGene)
		}else{
			stop("Couldn't find hg18refGene.rdata")
		}
	}else if(genome=="hg19"){
		refGene.file<-system.file("data/hg19refGene.rdata", package="r3Cseq")
		if(file.exists(refGene.file)==TRUE){
			load(file=refGene.file)
			refGene<-hg19refGene
			return(refGene)
		}else{
			stop("Couldn't find hg19refGene.rdata")
		}
	}else if(genome =="mm9"){
		refGene.file<-system.file("data/mm9refGene.rdata", package="r3Cseq")
		if(file.exists(refGene.file)==TRUE){
			load(file=refGene.file)
			refGene<-mm9refGene
			return(refGene)
		}else{
			stop("Couldn't find mm9refGene.rdata")
		}
	}else if(genome =="mm10"){
		refGene.file<-system.file("data/mm10refGene.rdata", package="r3Cseq")
		if(file.exists(refGene.file)==TRUE){
			load(file=refGene.file)
			refGene<-mm10refGene
			return(refGene)
		}else{
			stop("Couldn't find mm10refGene.rdata")
		}
	}else if(genome =="rn5"){
		refGene.file<-system.file("data/rn5refGene.rdata", package="r3Cseq")
		if(file.exists(refGene.file)==TRUE){
			load(file=refGene.file)
			refGene<-rn5refGene
			return(refGene)
		}else{
			stop("Couldn't find rn5refGene.rdata")
		}
	}else{
		stop("Require the selected genome: hg18, hg19, mm9, mm10 or rn5.")
	}
}

makeInteractionMatrixNearCisPerWindow<-function(obj,smoothing.parameter,rawReads.Ranged,max.window=25e3,viewpoint,distanceFromViewpoint=5e5){
	stopifnot( is( obj, "r3Cseq" ) |is( obj, "r3CseqInBatch" ) )
	
	if('BSgenome.Hsapiens.UCSC.hg19.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Hsapiens.UCSC.hg19.masked,unload=TRUE)
	}
	if('BSgenome.Hsapiens.UCSC.hg18.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Hsapiens.UCSC.hg18.masked,unload=TRUE)
	}
	if('BSgenome.Mmusculus.UCSC.mm9.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Mmusculus.UCSC.mm9.masked,unload=TRUE)
	}
	if('BSgenome.Mmusculus.UCSC.mm10.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Mmusculus.UCSC.mm10.masked,unload=TRUE)
	}
	if('BSgenome.Rnorvegicus.UCSC.rn5.masked' %in% loadedNamespaces()==TRUE){
		detach(package:BSgenome.Rnorvegicus.UCSC.rn5.masked,unload=TRUE)
	}
	###get raw reads
	if(length(rawReads.Ranged)==0){
		stop("No reads found, please run getRawReads function!!!.")
	}
	if(length(viewpoint)==0){
		stop("No viewpoint, please run getVeiwpoint function to get the viewpoint !!!.")
	}
	viewpoint<-getViewpoint(obj)
	viewpoint.chr <-as.character(seqnames(viewpoint))
	rawReads.chr<-rawReads.Ranged[as.character(seqnames(rawReads.Ranged))==viewpoint.chr]
	rawReads.filted<-excludeReadsNearViewpoint(obj,rawReads.chr,span_fragment=2)
	#######get organism name ############
	genome<-organismName(obj)
	#####################################
	chr.size<-0
	if(genome=="hg19"){
		library(BSgenome.Hsapiens.UCSC.hg19.masked)
		genome <- BSgenome.Hsapiens.UCSC.hg19.masked
		chr.size<-seqlengths(genome)[viewpoint.chr]
	}else if(genome=="hg18"){
		library(BSgenome.Hsapiens.UCSC.hg18.masked)
		genome <- BSgenome.Hsapiens.UCSC.hg18.masked
		chr.size<-seqlengths(genome)[viewpoint.chr]
	}else if(genome =="mm9"){
		library(BSgenome.Mmusculus.UCSC.mm9.masked)
		genome <- BSgenome.Mmusculus.UCSC.mm9.masked
		chr.size<-seqlengths(genome)[viewpoint.chr]
		
	}else if(genome =="mm10"){
		library(BSgenome.Mmusculus.UCSC.mm10.masked)
		genome <- BSgenome.Mmusculus.UCSC.mm10.masked
		chr.size<-seqlengths(genome)[viewpoint.chr]
		
	}else if(genome =="rn5"){
		library(BSgenome.Rnorvegicus.UCSC.rn5.masked)
		genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
		chr.size<-seqlengths(genome)[viewpoint.chr]
	}else{
		stop("Require the selected genome: hg18, hg19, mm9, mm10 or rn5")
	}
	r.start <-start(viewpoint)-distanceFromViewpoint
	r.end 	<-end(viewpoint)+distanceFromViewpoint
	r.size<-r.end-r.start+1
	ref.position<-c(seq(from=1,to=r.size,by=1e3),r.size)
	ref.v<-data.frame(position=ref.position)
	print("Analyzing interaction regions for each window....")
	for(window.size in rev(c(seq(from=2e3,to=max.window,by=1000)))){
		fragment.GRanges<-GRanges()
		step.size = window.size
		ends <- c(seq(from=min(start(rawReads.filted))+window.size,to=chr.size,by=step.size),as.integer(chr.size))
		starts <- seq(from=min(start(rawReads.filted)),to=chr.size,by=step.size)
		fragment<-data.frame(chromosome=viewpoint.chr,start=starts,end=ends-1)
		fragment.GRanges<-GRanges(seqnames =fragment$chromosome,IRanges(start=fragment$start,end=fragment$end))
		
		readsPerFragment <- countOverlaps(fragment.GRanges,subject=rawReads.filted)
		reads<-data.frame(nReads=readsPerFragment)
		#####migrating the code to support the BioC 3.9#####
		values(fragment.GRanges)<-reads
		countTable<-fragment.GRanges
		countTable<-countTable[countTable$nReads >0,]
		
		cis<-data.frame(chromosome=as.character(seqnames(countTable)),
		                start=start(countTable),
		                end=end(countTable),
		                relative.position=start(countTable)-start(viewpoint),observed=countTable$nReads)
		#######Smooth the data############
		e.smooth<-vsmooth.spline(cis$relative.position,cis$observed,spar=smoothing.parameter)
		e.fitted<-fitted(e.smooth)
		cis$expected<-as.numeric(e.fitted)
		
		######Calculate z-score and assign p-value in cis
		cis$z<-(cis$observed-cis$expected)/sd(cis$observed-cis$expected)
		cis$p.value<-2*pnorm(-abs(cis$z))
		cis$p.value[cis$z < 0]<-1
		cis$q.value<-qvalue(cis$p.value, fdr.level=0.05, pi0.method="bootstrap")$qvalues
		
		cis.in<-GRanges(seqnames = as.character(cis$chromosome),
		                IRanges(start=cis$start,end=cis$end),
		                nReads=cis$observed,z=cis$z,p.value=cis$p.value,q.value=cis$q.value)
		
		cis.in<-cis.in[start(cis.in) >=r.start & end(cis.in) <= r.end,]
		cis.in$rel.start<-start(cis.in)-r.start
		cis.in$rel.end<-end(cis.in)-r.start
		
		##Creat vector of q-value##
		my_vector<-rep(1,r.size)
		my.IRanges<-IRanges(start=cis.in$rel.start,end=cis.in$rel.end)
		#fix the IRanges problem here##########
		#cis.p  <-as.integer(my.IRanges)	
		cis.p  <-as.integer(IPos(my.IRanges))
		#######################################	
		my_vector[cis.p]<-c(rep(cis.in$q.value,each=window.size))
		ref.value<-data.frame(q.value=my_vector[ref.position])
		colnames(ref.value)<-paste("window.",window.size,sep="")
		ref.v<-cbind(ref.v,ref.value)
	}
	return(ref.v)
	print("The interaction matrix is created !!!.")
}

fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)

