# TODO: The following functions are implemented for get the list of genes containing high interactions.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com
###############################################################################
getExpInteractionsInRefseq<-function(obj,cutoff.qvalue=0.05,expanded_upstream=50e3,expanded_downstream=10e3){
	stopifnot( is(obj, "r3Cseq") | is(obj,"r3CseqInBatch"))	
	########Get organism#############
	orgName<-organismName(obj)
	expInteractions <-expInteractionRegions(obj)
	expInteractions <-expInteractions[expInteractions$q.value<=cutoff.qvalue,]
	if(length(expInteractions)>0){
		expInteraction.GRanges<-GRanges(seqnames=seqnames(expInteractions),
		                                IRanges(start=start(expInteractions),
		                                        end=end(expInteractions)),
				nReads=expInteractions$nReads,RPMs=expInteractions$RPMs,q.value=expInteractions$q.value)
		expInteraction.frame<-data.frame(chromosome=seqnames(expInteractions),
		                                 start=start(expInteractions),end=end(expInteractions),
				nReads=expInteractions$nReads,RPMs=expInteractions$RPMs)
		#########Get genes###############
		genes<-get3CseqRefGene(obj)
		#################################
		genes$r_start<-ifelse(genes$strand==1,genes$start-expanded_upstream,genes$start-expanded_downstream)
		genes$r_end<-ifelse(genes$strand==1,genes$end+expanded_downstream,genes$end+expanded_upstream)
		genes.GRanges<-GRanges(seqnames=genes$chromosome,IRanges(start=genes$r_start,end=genes$r_end),name=genes$name)
		genes.frame<-data.frame(chromosome_g=genes$chromosome,start=genes$r_start,end=genes$r_end,name=genes$name)
		##########perform overlapping####
		o<- findOverlaps(expInteraction.GRanges,subject=genes.GRanges)
		
		my.interactions<-expInteraction.frame[queryHits(o), c("nReads", "RPMs")]
		my.genes<-genes.frame[subjectHits(o), c("chromosome_g", "name")]
		ol <- cbind(my.interactions,my.genes) 
		annotated<- sqldf("select chromosome_g,name,sum(nReads),sum(RPMs) from ol group by name")
		annotated<- annotated[order(annotated[,4],decreasing=T),]
		colnames(annotated)<-c("chromosome","gene_name","total_nReads","total_RPMs")
		return(annotated)
	}else{
		stop("No interactions found in the control!!!")
	}
}

getContrInteractionsInRefseq<-function(obj,cutoff.qvalue=0.05,expanded_upstream=50e3,expanded_downstream=10e3){
	stopifnot( is(obj, "r3Cseq") | is(obj,"r3CseqInBatch"))	
	########Get organism#############
	orgName<-organismName(obj)
	contrInteractions <-contrInteractionRegions(obj)
	contrInteractions <-contrInteractions[contrInteractions$q.value<=cutoff.qvalue,]
	if(length(contrInteractions)>0){
		contrInteraction.GRanges<-GRanges(seqnames=seqnames(contrInteractions),
		                                  IRanges(start=start(contrInteractions),end=end(contrInteractions)),
				nReads=contrInteractions$nReads,RPMs=contrInteractions$RPMs,q.value=contrInteractions$q.value)
		contrInteraction.frame<-data.frame(chromosome=seqnames(contrInteractions),
		                                   start=start(contrInteractions),end=end(contrInteractions),
				nReads=contrInteractions$nReads,RPMs=contrInteractions$RPMs)
		#########Get genes###############
		genes<-get3CseqRefGene(obj)
		#################################
		genes$r_start<-ifelse(genes$strand==1,genes$start-expanded_upstream,genes$start-expanded_downstream)
		genes$r_end<-ifelse(genes$strand==1,genes$end+expanded_downstream,genes$end+expanded_upstream)
		genes.GRanges<-GRanges(seqnames=genes$chromosome,IRanges(start=genes$r_start,end=genes$r_end),name=genes$name)
		genes.frame<-data.frame(chromosome_g=genes$chromosome,start=genes$r_start,end=genes$r_end,name=genes$name)
		##########perform overlapping####
		o<- findOverlaps(contrInteraction.GRanges,subject=genes.GRanges)
		
		my.interactions<-contrInteraction.frame[queryHits(o),]
		my.genes<-genes.frame[subjectHits(o),]
		
		ol <-cbind(my.interactions,my.genes) 
		annotated<- sqldf("select chromosome_g,name,sum(nReads),sum(RPMs) from ol group by name")
		annotated<- annotated[order(annotated[,4],decreasing=T),]
		colnames(annotated)<-c("chromosome","gene_name","total_nReads","total_RPMs")
		return(annotated)
	}else{
		stop("No interactions found in the control!!!")
	}
}
