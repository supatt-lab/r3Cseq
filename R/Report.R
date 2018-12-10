# TODO: The function below can be used to generate the report.
# Author: Supat Thongjuea
# Contact : supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com
###############################################################################
#####
#####The functions below are completely migrated to BioC 3.9 on R 3.6
#####
generate3CseqReport<-function (obj){
	
	stopifnot( is( obj, "r3Cseq" ) |is( obj, "r3CseqInBatch" ) )
	
	if(is(obj, "r3Cseq" )==TRUE){
		if(isControlInvolved(obj)==FALSE){
			expInteractions<-expInteractionRegions(obj)
			expLabeled<-expLabel(obj)
			viewpoint<-getViewpoint(obj)
			
			if(length(expInteractions)>0){
				
				file_name<-paste(expLabeled,".pdf",sep="")
				pdf(file=file_name, width=12, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)
				plotOverviewInteractions(obj)
				plotInteractionsNearViewpoint(obj)
				plotInteractionsPerChromosome(obj,as.character(seqnames(viewpoint)))
				plotDomainogramNearViewpoint(obj)
				dev.off()
				exportInteractions2text(obj)
				export3Cseq2bedGraph(obj)
				print("Three files are generated : a pdf file of plots, a text file of interaction regions, and a bedGraph file.")
			}else{
				stop("No interaction regions found in the r3Cseq obj, you have to run the r3Cseq pipeline.")	
			}
		}
		if(isControlInvolved(obj)==TRUE){
			expInteractions  <-expInteractionRegions(obj)
			contrInteractions <- contrInteractionRegions(obj)
			expLabeled<-expLabel(obj)
			controlLabeled<-contrLabel(obj)
			viewpoint<-getViewpoint(obj)
			if(length(expInteractions)>0){
				file_name<-paste(expLabeled,"_",controlLabeled,".pdf",sep="")
				pdf(file=file_name, width=12, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)
				plotOverviewInteractions(obj)
				plotInteractionsNearViewpoint(obj)
				plotInteractionsPerChromosome(obj,as.character(seqnames(viewpoint)))
				plotDomainogramNearViewpoint(obj,view="both")
				dev.off()
				exportInteractions2text(obj)
				export3Cseq2bedGraph(obj)
				print("Three files are generated : a pdf file of plots, a text file of interaction regions, and a bedGraph file.")
			}else{
				stop("No interaction regions found in the r3Cseq object, you have to run the r3Cseq pipeline.")	
			}	
		}
	}
	if(is(obj, "r3CseqInBatch" )==TRUE){
		if(isControlInvolved(obj)==TRUE){
			expInteractions  <-expInteractionRegions(obj)
			contrInteractions <- contrInteractionRegions(obj)
			viewpoint<-getViewpoint(obj)
			if(length(expInteractions)>0){
				file_name<-"myBatch-Analysis-Result.pdf"
				pdf(file=file_name, width=12, height=8, onefile=T, bg="transparent",family = "Helvetica",fonts = NULL)
				plotOverviewInteractions(obj)
				plotInteractionsNearViewpoint(obj)
				plotInteractionsPerChromosome(obj,as.character(seqnames(viewpoint)))
				dev.off()
				exportBatchInteractions2text(obj)
				print("Files are generated : a pdf file of plots and a text file of interaction regions from replicates.")
			}else{
				stop("No interaction regions found in the r3CseqInBatch object, you have to run the r3Cseq pipeline.")	
			}
		}
	}
}