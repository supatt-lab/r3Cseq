# These functions are for getting Repbase enzyme information.
###############################################################################
setMethod("initialize", "repbaseEnzyme",
		function(.Object)
		{
			enzyme.file<-system.file("data", "enzymeDb.rda", package="r3Cseq")
			if(file.exists(enzyme.file)==TRUE){
				data<-load(file=enzyme.file)
				.Object@enzymeRestriction = enzyme.db
				.Object
			}else{
				stop("Couldn't find enzymeDb.rda")
			}
		} 
)

setGeneric(
		name="enzymeRestriction",
		def=function(object){
			standardGeneric("enzymeRestriction")
		})

setMethod("enzymeRestriction",
		signature(object = "repbaseEnzyme"),
		function (object){
			object@enzymeRestriction
		}
)
setGeneric(
		name="getEnzymeRestrictionSequences",
		def=function(object,enzyme.name){
			standardGeneric("getEnzymeRestrictionSequences")
		})

setMethod("getEnzymeRestrictionSequences",
		   signature(object = "repbaseEnzyme",enzyme.name="character"),
	 	function(object,enzyme.name){
			if(!is(object,"repbaseEnzyme")){
			   stop("Need the repbaseEnzyme object")
			}
			enzyme.hit<-subset(enzymeRestriction(object),enzyme==enzyme.name)
			if(nrow(enzyme.hit) > 0){
				sequences <- as.character(enzyme.hit$restriction.site.sequences)
			sequences
			}else{
				stop("Not found enzyme's name in the database")
			}
	}
)
setGeneric(
		name="getEnzymeRestrictionPositionInSelectedGenome",
		def=function(object,enzyme.name,genome,chromosome){
			standardGeneric("getEnzymeRestrictionPositionInSelectedGenome")
		})
setMethod("getEnzymeRestrictionPositionInSelectedGenome",
		signature = "repbaseEnzyme",
		definition= function(object,enzyme.name,genome,chromosome){
	
		if(!is(object,"repbaseEnzyme")){
				stop("Need the repbaseEnzyme object")
		}

		if(enzyme.name==character(1)){
			stop("Require the restriction enzyme name for example : 'HindIII'")
		}
		if(genome==character(1)){
			stop("Require the UCSC genome name for example: 'mm9'")
		}
		if(chromosome ==character(1)){
			stop("Require the chromosome name for example : 'chr1'")
		}
		
		if(!chromosome %in% paste('chr',c(seq(1,50),'X','Y'),sep='')){
			stop("Require the correct format chromosome name : 'chr1','chrX','chrY'")
		}
		sequences<-getEnzymeRestrictionSequences(object,enzyme.name)	
		
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
			hits<-matchPattern(sequences,genome[[chromosome]],fixed=FALSE)
			hits.frame<-data.frame(chromosome=chromosome,start=start(hits),end=end(hits))
			return(hits.frame)
		}else if(genome=="hg19"){
			library(BSgenome.Hsapiens.UCSC.hg19.masked)
			genome <- BSgenome.Hsapiens.UCSC.hg19.masked
			hits<-matchPattern(sequences,genome[[chromosome]],fixed=FALSE)
			hits.frame<-data.frame(chromosome=chromosome,start=start(hits),end=end(hits))
			return(hits.frame)
		}else if(genome =="mm9"){
			library(BSgenome.Mmusculus.UCSC.mm9.masked)
			genome <- BSgenome.Mmusculus.UCSC.mm9.masked
			hits<-matchPattern(sequences,genome[[chromosome]],fixed=FALSE)
			hits.frame<-data.frame(chromosome=chromosome,start=start(hits),end=end(hits))
			return(hits.frame)
			
		}else if(genome =="mm10"){
			library(BSgenome.Mmusculus.UCSC.mm10.masked)
			genome <- BSgenome.Mmusculus.UCSC.mm10.masked
			hits<-matchPattern(sequences,genome[[chromosome]],fixed=FALSE)
			hits.frame<-data.frame(chromosome=chromosome,start=start(hits),end=end(hits))
			return(hits.frame)	
		}else if(genome =="rn5"){
			library(BSgenome.Rnorvegicus.UCSC.rn5.masked)
			genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
			hits<-matchPattern(sequences,genome[[chromosome]],fixed=FALSE)
			hits.frame<-data.frame(chromosome=chromosome,start=start(hits),end=end(hits))
			return(hits.frame)	
		}else{
			stop("Require the selected genome: hg18, hg19, mm9, mm10 or rn5.")
		}
	}
)
setGeneric(
		name="getRestrictionFragments",
		def=function(object,enzyme.name,genome,chromosome){
			standardGeneric("getRestrictionFragments")
		})
setMethod("getRestrictionFragments",
		  signature = "repbaseEnzyme",
		  definition= function(object,enzyme.name,genome,chromosome){
			restriction.sites <-getEnzymeRestrictionPositionInSelectedGenome(object,enzyme.name,genome,chromosome)
			restriction.sites <-restriction.sites[order(restriction.sites[,2]),]
			position.x <- restriction.sites[1:nrow(restriction.sites)-1,]
			position.y <- restriction.sites[2:nrow(restriction.sites),]
			fragments  <- data.frame(chromosome=chromosome,start=position.x$end,end=position.y$start)
			good.fragments <- subset(fragments,end-start>=1)
			return(good.fragments)
		  }
)
setGeneric(
		name="getWholeGenomeRestrictionFragments",
		def=function(object,enzyme.name,genome){
			standardGeneric("getWholeGenomeRestrictionFragments")
		})
setMethod("getWholeGenomeRestrictionFragments",
		signature = "repbaseEnzyme",
		definition= function(object,enzyme.name,genome){
			fragments<-data.frame()
			if(genome=="hg18"){
				for (chr in paste('chr',c(seq(1,22),'X','Y'),sep='')){
					chr.fragments<-getRestrictionFragments(object,enzyme.name,"hg18",chr)
					fragments<-rbind(fragments,chr.fragments)
				}
				return(fragments)
			}else if(genome=="hg19"){
				for (chr in paste('chr',c(seq(1,22),'X','Y'),sep='')){
					chr.fragments<-getRestrictionFragments(object,enzyme.name,"hg19",chr)
					fragments<-rbind(fragments,chr.fragments)
				}
				return(fragments)
			}else if(genome =="mm9"){
				for (chr in paste('chr',c(seq(1,19),'X','Y'),sep='')){
					chr.fragments<-getRestrictionFragments(object,enzyme.name,"mm9",chr)
					fragments<-rbind(fragments,chr.fragments)
				}
				return(fragments)
			}else if(genome =="mm10"){
				for (chr in paste('chr',c(seq(1,19),'X','Y'),sep='')){
					chr.fragments<-getRestrictionFragments(object,enzyme.name,"mm10",chr)
					fragments<-rbind(fragments,chr.fragments)
				}
				return(fragments)
			}else if(genome =="rn5"){
				for (chr in paste('chr',c(seq(1,20),'X'),sep='')){
					chr.fragments<-getRestrictionFragments(object,enzyme.name,"rn5",chr)
					fragments<-rbind(fragments,chr.fragments)
				}
				return(fragments)
			}else{
				stop("Require the selected genome: hg18, hg19, mm9, mm10 or rn5.")
			}
		}
)
