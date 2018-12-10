# TODO: These following functions are implemented for visualizing 3C-seq data.
# Author: Supat Thongjuea
# Contact:supat.thongjuea@imm.ox.ac.uk or supat.thongjuea@gmail.com 
#####################################
#####The functions below are completely migrated to BioC 3.9 on R 3.6
#####################################
plotOverviewInteractions<-function (obj,cutoff.qvalue=0.05){
			
			stopifnot( is(obj, "r3Cseq") | is(obj,"r3CseqInBatch"))
			
			if(isControlInvolved(obj)==FALSE){
				orgName<-organismName(obj)
				chr.data=c()
				if(orgName=="hg18"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){
						genome <- BSgenome.Hsapiens.UCSC.hg18.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName=="hg19"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){
						genome <- BSgenome.Hsapiens.UCSC.hg19.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName =="mm9"){
					for (chr in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){
						genome <- BSgenome.Mmusculus.UCSC.mm9.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
					
				}else if(orgName =="mm10"){
					for (chr in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){
						genome <- BSgenome.Mmusculus.UCSC.mm10.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName =="rn5"){
					for (chr in c(paste('chr',seq(1,20),sep=''),'chrX')){
						genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else{
					stop("Your input organism name is not in the list ('mm9','mm10','hg18','hg19', and 'rn5')")
				}
				######check interactions############
				expInteractions <-expInteractionRegions(obj)
				
				if(length(expInteractions) ==0){
					stop("There are no interaction regions found in r3Cseq object. Use 'getInteractions' function to get interaction regions")
				}
				
				exp.filted <-expInteractions[expInteractions$q.value <=cutoff.qvalue,]
				
				if(length(exp.filted) ==0){
					stop("There are no interaction regions pass your qvalue cutoff.")
				}
				
				chr.size.max<-max(chr.data$size)
				max.scale <- floor(chr.size.max/10^6)
				box.x.size=chr.size.max+5e6
				plot(c(1,box.x.size), c(1,100), type= "n", ylab="",yaxt='n',
						xaxt='n',xlab="Chromosomal position (Mbp)",
						main=paste("3C-seq distribution of interaction regions (q-value <=",cutoff.qvalue,")"))
				
				axis(1, at=(seq(0, max.scale*10^6, by=10*10^6)),labels=c(seq(0,max.scale,by=10)),cex.axis=0.8)
				
				#expLabeled<-expLabel(obj)
				polygon(c(chr.size.max-10e6-10e5,chr.size.max-10e6,chr.size.max-10e6+10e5),
						c(52,50,52), col="red")
				text(chr.size.max-2e6,51,"  viewpoint",cex = .8)
				
				y=0;
				
				for(i in 1:nrow(chr.data)){
					y=y+4
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="grey",lwd=16)
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="black")
					text.x<-chr.data$size[i]+5e6
					text(text.x, y+1, chr.data$name[i],cex = .8)
				}
				
				start.exp <-seq(from=3,to=100,by=4)
				end.exp   <-seq(from=6,to=100,by=4)
				
				exp.coor<-data.frame(i=start.exp[1:24],j=end.exp[1:24])
				
				exp.filted$rank<-rank(-exp.filted$RPMs)
				exp.max.rank<-max(exp.filted$rank)
				
				exppalette<-rev(brewer.pal(7,"Reds"))
				exp.filted$col<-""
				exp.filted$col[exp.filted$rank <=round(0.01*exp.max.rank)]<-exppalette[1]
				exp.filted$col[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)]<-exppalette[2]
				exp.filted$col[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)]<-exppalette[3]
				exp.filted$col[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)]<-exppalette[4]
				exp.filted$col[exp.filted$rank > round(0.3*exp.max.rank) & round(0.4*exp.max.rank)]<-exppalette[5]
				exp.filted$col[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)]<-exppalette[6]
				exp.filted$col[exp.filted$rank > round(0.5*exp.max.rank)]<-exppalette[7]
				
				if(round(0.01*exp.max.rank)>0){
					sv1<-min(exp.filted$RPMs[exp.filted$rank <=round(0.01*exp.max.rank)])
					sv2<-min(exp.filted$RPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$RPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$RPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$RPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$RPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
				
					name1<-paste(">",round(sv1)," RPMs")
					name2<-paste(round(sv2),"-",round(sv1)," RPMs")
					name3<-paste(round(sv3),"-",round(sv2)," RPMs")
					name4<-paste(round(sv4),"-",round(sv3)," RPMs")
					name5<-paste(round(sv5),"-",round(sv4)," RPMs")
					name6<-paste(round(sv6),"-",round(sv5)," RPMs")
					name7<-paste("<",round(sv6)," RPMs")
				}else{
					sv1<-max(exp.filted$RPMs)
					sv2<-min(exp.filted$RPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$RPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$RPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$RPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$RPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
					
					name1<-paste(">",round(sv1)," RPMs")
					name2<-paste(round(sv2),"-",round(sv1)," RPMs")
					name3<-paste(round(sv3),"-",round(sv2)," RPMs")
					name4<-paste(round(sv4),"-",round(sv3)," RPMs")
					name5<-paste(round(sv5),"-",round(sv4)," RPMs")
					name6<-paste(round(sv6),"-",round(sv5)," RPMs")
					name7<-paste("<",round(sv6)," RPMs")
				}
				exp.names<-c(name1,name2,name3,name4,name5,name6,name7)
				
				
				viewpoint <-getViewpoint(obj)
				i=0
				for (chri in 1:nrow(chr.data)){
					i=i+1
					if(as.character(chr.data$name[chri]) %in% as.character(seqnames(exp.filted))){
						exp.chr<-exp.filted[seqnames(exp.filted)==as.character(chr.data$name[chri]),]
						if(length(exp.chr) >0){
							rect(start(exp.chr),exp.coor$i[i], end(exp.chr), exp.coor$j[i], border=exp.chr$col)		
						}
						if(as.character(seqnames(viewpoint))==as.character(chr.data$name[chri])){
							polygon(c(start(viewpoint)-10e5,start(viewpoint),start(viewpoint)+10e5),
								c(exp.coor$j[i]+1.5,exp.coor$j[i],exp.coor$j[i]+1.5), col="red")
						}
					}
				}
				legend("topright",legend = exp.names, fill=exppalette, cex=0.55,title="experiment",bty="n")
			}
			if(isControlInvolved(obj)==TRUE){
				
				######check interactions######
				expInteractions <-expInteractionRegions(obj)
				contrInteractions <-contrInteractionRegions(obj)
				
				if(length(expInteractions) ==0){
					stop("There are no interaction regions found in r3Cseq object. Use 'getInteractions' function to get interaction regions")
				}
				
				exp.filted <-expInteractions[expInteractions$q.value <=cutoff.qvalue,]
				contr.filted <-contrInteractions[contrInteractions$q.value <=cutoff.qvalue,]
				if(length(exp.filted) ==0){
					stop("There are no interaction regions pass your input parameters.")
				}
				
				#######draw chromosome########
				
				orgName<-organismName(obj)
				chr.data=c()
				if(orgName=="hg18"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){
						genome <- BSgenome.Hsapiens.UCSC.hg18.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName=="hg19"){
					for (chr in c(paste('chr',seq(1,22),sep=''),'chrX','chrY')){
						genome <- BSgenome.Hsapiens.UCSC.hg19.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName =="mm9"){
					for (chr in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){
						genome <- BSgenome.Mmusculus.UCSC.mm9.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName =="mm10"){
					for (chr in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){
						genome <- BSgenome.Mmusculus.UCSC.mm10.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else if(orgName =="rn5"){
					for (chr in c(paste('chr',seq(1,20),sep=''),'chrX')){
						genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
						chr.size<-seqlengths(genome)[chr]
						chr.info<-data.frame(name=chr,size=chr.size)
						chr.data<-rbind(chr.data,chr.info)
					}
				}else{
					stop("Your input organism name is not in the list ('mm9','mm10','hg18','hg19', and 'rn5')")
				}
	
				chr.size.max<-max(chr.data$size)
				max.scale <- floor(chr.size.max/10^6)
				box.x.size=chr.size.max+5e6
				plot(c(1,box.x.size), c(1,100), type= "n", ylab="",yaxt='n',
						xaxt='n',xlab="Chromosomal position (Mbp)",
						main=paste("3C-seq distribution of interaction regions (q-value <=",cutoff.qvalue,")"))
				
				axis(1, at=(seq(0, max.scale*10^6, by=10*10^6)),labels=c(seq(0,max.scale,by=10)),cex.axis=0.8)
				
				#expLabeled<-expLabel(obj)
				#controlLabeled<-contrLabel(obj)
				polygon(c(chr.size.max-10e6-10e5,chr.size.max-10e6,chr.size.max-10e6+10e5),
						c(52,50,52), col="red")
				text(chr.size.max-2e6,51,"  viewpoint",cex = .8)
				
				y=0;
				
				for(i in 1:nrow(chr.data)){
					y=y+4
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="grey",lwd=16)
					rect(1, y+0.5, chr.data$size[i], y+0.5, border="black")
					text.x<-chr.data$size[i]+5e6
					text(text.x, y+1, chr.data$name[i],cex = .8)
				}
				
				start.contr <-seq(from=3,to=100,by=4)
				end.contr   <-seq(from=4.5,to=100,by=4)
				
				start.exp <-seq(from=4.5,to=100,by=4)
				end.exp   <-seq(from=6,to=100,by=4)
				
				exp.coor<-data.frame(i=start.exp[1:24],j=end.exp[1:24])
				contr.coor<-data.frame(i=start.contr[1:24],j=end.contr[1:24])
				
				exp.filted$rank<-rank(-exp.filted$RPMs)
				contr.filted$rank<-rank(-contr.filted$RPMs)
				
				exp.max.rank<-max(exp.filted$rank)
				contr.max.rank<-max(contr.filted$rank)
				
				exppalette<-rev(brewer.pal(7,"Greens"))
				exp.filted$col<-""
				exp.filted$col[exp.filted$rank <=round(0.01*exp.max.rank)]<-exppalette[1]
				exp.filted$col[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)]<-exppalette[2]
				exp.filted$col[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)]<-exppalette[3]
				exp.filted$col[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)]<-exppalette[4]
				exp.filted$col[exp.filted$rank > round(0.3*exp.max.rank) & round(0.4*exp.max.rank)]<-exppalette[5]
				exp.filted$col[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)]<-exppalette[6]
				exp.filted$col[exp.filted$rank > round(0.5*exp.max.rank)]<-exppalette[7]
				
				if(round(0.01*exp.max.rank)>0){
					sv1<-min(exp.filted$RPMs[exp.filted$rank <=round(0.01*exp.max.rank)])
					sv2<-min(exp.filted$RPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$RPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$RPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$RPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$RPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
	
					name1<-paste(">",round(sv1)," RPMs")
					name2<-paste(round(sv2),"-",round(sv1)," RPMs")
					name3<-paste(round(sv3),"-",round(sv2)," RPMs")
					name4<-paste(round(sv4),"-",round(sv3)," RPMs")
					name5<-paste(round(sv5),"-",round(sv4)," RPMs")
					name6<-paste(round(sv6),"-",round(sv5)," RPMs")
					name7<-paste("<",round(sv6)," RPMs")
				}else{
					sv1<-max(exp.filted$RPMs)
					sv2<-min(exp.filted$RPMs[exp.filted$rank > round(0.01*exp.max.rank) & exp.filted$rank <=round(0.1*exp.max.rank)])
					sv3<-min(exp.filted$RPMs[exp.filted$rank > round(0.1*exp.max.rank) & exp.filted$rank <=round(0.2*exp.max.rank)])
					sv4<-min(exp.filted$RPMs[exp.filted$rank > round(0.2*exp.max.rank) & exp.filted$rank <=round(0.3*exp.max.rank)])
					sv5<-min(exp.filted$RPMs[exp.filted$rank > round(0.3*exp.max.rank) & exp.filted$rank <=round(0.4*exp.max.rank)])
					sv6<-min(exp.filted$RPMs[exp.filted$rank > round(0.4*exp.max.rank) & exp.filted$rank <=round(0.5*exp.max.rank)])
					
					name1<-paste(">",round(sv1)," RPMs")
					name2<-paste(round(sv2),"-",round(sv1)," RPMs")
					name3<-paste(round(sv3),"-",round(sv2)," RPMs")
					name4<-paste(round(sv4),"-",round(sv3)," RPMs")
					name5<-paste(round(sv5),"-",round(sv4)," RPMs")
					name6<-paste(round(sv6),"-",round(sv5)," RPMs")
					name7<-paste("<",round(sv6)," RPMs")
				}
				exp.names<-c(name1,name2,name3,name4,name5,name6,name7)
				
			
				contrpalette<-rev(brewer.pal(7,"Reds"))
				contr.filted$col<-""
				contr.filted$col[contr.filted$rank <=round(0.01*contr.max.rank)]<-contrpalette[1]
				contr.filted$col[contr.filted$rank > round(0.01*contr.max.rank) & contr.filted$rank <=round(0.1*contr.max.rank)]<-contrpalette[2]
				contr.filted$col[contr.filted$rank > round(0.1*contr.max.rank) & contr.filted$rank <=round(0.2*contr.max.rank)]<-contrpalette[3]
				contr.filted$col[contr.filted$rank > round(0.2*contr.max.rank) & contr.filted$rank <=round(0.3*contr.max.rank)]<-contrpalette[4]
				contr.filted$col[contr.filted$rank > round(0.3*contr.max.rank) & round(0.4*contr.max.rank)]<-contrpalette[5]
				contr.filted$col[contr.filted$rank > round(0.4*contr.max.rank) & contr.filted$rank <=round(0.5*contr.max.rank)]<-contrpalette[6]
				contr.filted$col[contr.filted$rank > round(0.5*contr.max.rank)]<-contrpalette[7]
				
				if(round(0.01*contr.max.rank)>0){
					csv1<-min(contr.filted$RPMs[contr.filted$rank <=round(0.01*contr.max.rank)])
					csv2<-min(contr.filted$RPMs[contr.filted$rank > round(0.01*contr.max.rank) & contr.filted$rank <=round(0.1*contr.max.rank)])
					csv3<-min(contr.filted$RPMs[contr.filted$rank > round(0.1*contr.max.rank) & contr.filted$rank <=round(0.2*contr.max.rank)])
					csv4<-min(contr.filted$RPMs[contr.filted$rank > round(0.2*contr.max.rank) & contr.filted$rank <=round(0.3*contr.max.rank)])
					csv5<-min(contr.filted$RPMs[contr.filted$rank > round(0.3*contr.max.rank) & contr.filted$rank <=round(0.4*contr.max.rank)])
					csv6<-min(contr.filted$RPMs[contr.filted$rank > round(0.4*contr.max.rank) & contr.filted$rank <=round(0.5*contr.max.rank)])
					
					tname1<-paste(">",round(csv1)," RPMs")
					tname2<-paste(round(csv2),"-",round(csv1)," RPMs")
					tname3<-paste(round(csv3),"-",round(csv2)," RPMs")
					tname4<-paste(round(csv4),"-",round(csv3)," RPMs")
					tname5<-paste(round(csv5),"-",round(csv4)," RPMs")
					tname6<-paste(round(csv6),"-",round(csv5)," RPMs")
					tname7<-paste("<",round(csv6)," RPMs")
				}else{
					csv1<-max(contr.filted$RPMs)
					csv2<-min(contr.filted$RPMs[contr.filted$rank > round(0.01*contr.max.rank) & contr.filted$rank <=round(0.1*contr.max.rank)])
					csv3<-min(contr.filted$RPMs[contr.filted$rank > round(0.1*contr.max.rank) & contr.filted$rank <=round(0.2*contr.max.rank)])
					csv4<-min(contr.filted$RPMs[contr.filted$rank > round(0.2*contr.max.rank) & contr.filted$rank <=round(0.3*contr.max.rank)])
					csv5<-min(contr.filted$RPMs[contr.filted$rank > round(0.3*contr.max.rank) & contr.filted$rank <=round(0.4*contr.max.rank)])
					csv6<-min(contr.filted$RPMs[contr.filted$rank > round(0.4*contr.max.rank) & contr.filted$rank <=round(0.5*contr.max.rank)])
					
					tname1<-paste(">",round(csv1)," RPMs")
					tname2<-paste(round(csv2),"-",round(csv1)," RPMs")
					tname3<-paste(round(csv3),"-",round(csv2)," RPMs")
					tname4<-paste(round(csv4),"-",round(csv3)," RPMs")
					tname5<-paste(round(csv5),"-",round(csv4)," RPMs")
					tname6<-paste(round(csv6),"-",round(csv5)," RPMs")
					tname7<-paste("<",round(csv6)," RPMs")
				}
				
				contr.names<-c(tname1,tname2,tname3,tname4,tname5,tname6,tname7)
				
				viewpoint <-getViewpoint(obj)
				
				i=0
				for (chri in 1:nrow(chr.data)){
					i=i+1
						if(as.character(chr.data$name[chri]) %in% as.character(seqnames(exp.filted))){  
							exp.chr<-exp.filted[as.character(seqnames(exp.filted))==as.character(chr.data$name[chri]),] 
							contr.chr<-contr.filted[as.character(seqnames(contr.filted))==as.character(chr.data$name[chri]),]
							if(length(exp.chr) >0){
								rect(start(exp.chr),exp.coor$i[i], end(exp.chr), exp.coor$j[i], border=exp.chr$col,lwd=1.5)		
							}
							if(length(contr.chr) >0){
								rect(start(contr.chr), contr.coor$i[i], end(contr.chr), contr.coor$j[i], border=contr.chr$col,lwd=1.5)	
							}
							if(as.character(seqnames(viewpoint))==as.character(chr.data$name[chri])){
								polygon(c(start(viewpoint)-10e5,start(viewpoint),start(viewpoint)+10e5),
								c(exp.coor$j[i]+1.5,exp.coor$j[i],exp.coor$j[i]+1.5), col="red")
							}
						}
				}
				legend(chr.size.max-10e6-10e5,100,legend = exp.names, fill=exppalette, cex=0.55,title="experiment",bty="n")
				legend(chr.size.max-10e6-10e5,80,legend = contr.names, fill=contrpalette, cex=0.55,title="control",bty="n")
			}
}

plotInteractionsNearViewpoint<-function(obj,distance=5e5,log2fc_cutoff=1,yLim=0){
			
			stopifnot( is(obj, "r3Cseq") | is(obj,"r3CseqInBatch"))
			if(distance < 50000 | distance > 500000){
				print("You distance is too high or too low!!!. Please input the distance between 50Kb - 500Kb")
			}
			########Get viewpoint############
			viewpoint <-getViewpoint(obj)
			viewpoint.chr<-as.character(seqnames(viewpoint))
			########Get organism#############
			orgName<-organismName(obj)
			chr.size<-0
			if(orgName=="hg18"){
					genome <- BSgenome.Hsapiens.UCSC.hg18.masked
					chr.size<-seqlengths(genome)[viewpoint.chr]
			}else if(orgName=="hg19"){
					genome <- BSgenome.Hsapiens.UCSC.hg19.masked
					chr.size<-seqlengths(genome)[viewpoint.chr]
			}else if(orgName =="mm9"){
					genome <- BSgenome.Mmusculus.UCSC.mm9.masked
					chr.size<-seqlengths(genome)[viewpoint.chr]	
			}else if(orgName =="mm10"){
					genome <- BSgenome.Mmusculus.UCSC.mm10.masked
					chr.size<-seqlengths(genome)[viewpoint.chr]
			}else if(orgName =="rn5"){
					genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
					chr.size<-seqlengths(genome)[viewpoint.chr]	
			}
			
			if(isControlInvolved(obj)==FALSE){
				expInteractions   <-expInteractionRegions(obj)
				#expLabeled <-expLabel(obj)
				######look around the viewpoint###
				r.start <-start(viewpoint)-distance
				r.end 	<-end(viewpoint)+distance
	
				r.start <-ifelse(r.start <0,1,r.start)
				r.end   <-ifelse(r.end <=chr.size,r.end,chr.size)
				c.vec   <- c(rep(0,r.end-r.start))
				####Get refGenes######
				genes<-get3CseqRefGene(obj)
				g.chr<-subset(genes,chromosome==viewpoint.chr)
				g.r<-subset(g.chr,start>=r.start & end <= r.end)
				
				par(fig=c(0,1,0.7,1))
				par(mar=c(0,5,1,2))
				
				if(nrow(g.r)>0){
					g.r<-g.r[order(g.r$start),]
					g.r$rel.start <-g.r$start-r.start
					g.r$rel.end <-g.r$end-r.start
					g.r$size<-g.r$end-g.r$start+1
					
					plot(c(1,2*distance), c(1,60), type= "n", xlab="", ylab="",xaxt='n',yaxt='n')
					abline(v=start(viewpoint)-r.start, col="red",lty=3,lwd=2)
					y1.start=30
					for(i in 1:nrow(g.r)){
								gx<-g.r[i,]
								gx.start<-gx$rel.start
								gx.end <-gx$rel.end
								gx.size<-gx$size
								gx.strand<-gx$strand
								
								if(c.vec[gx$rel.start]==0){
									y1.start = y1.start
									c.vec[gx$rel.start:gx$rel.end]=1
								}else{
									y1.start = y1.start+8
									if(y1.start >50){
										y1.start=10
									}
								}
								text(gx$rel.start,y1.start+2,gx$name,cex =0.8)
								
								if(gx.size >=100){
									s.q <-ifelse(gx.size <=5000,500,5000)
									s.q <-ifelse(s.q > gx.size, gx.size-10,s.q)
									
									xx<-seq(gx.start,gx.start+gx.size-s.q, by=s.q)
									yy<-seq(gx.start+s.q,gx.start+gx.size, by=s.q)
									
									if(gx.strand==-1){
										for(i in 1:length(xx)){
											lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-1))
											lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-4))
										}
									}else{
										for(i in 1:length(xx)){
											lines(c(xx[i],yy[i]),c(y1.start-1,y1.start-2))
											lines(c(xx[i],yy[i]),c(y1.start-4,y1.start-3))
										}
									}
									rect(gx.start,y1.start-3,gx.start+gx.size,(y1.start-2), col="red")
								}	
					}
				}
				legend("topleft",legend = "Refseq Genes",fill="red", cex=0.55,bty="n")
				####Draw Restriction map####
				enzymeDb	<-new("repbaseEnzyme")
				resEnzyme 	<-restrictionEnzyme(obj)	
				chr.fragment <- getRestrictionFragments(enzymeDb,resEnzyme,orgName,viewpoint.chr)
				
				RE.r <-subset(chr.fragment,start >=r.start & end <=r.end)
				
				RE.r$rel.start <-RE.r$start-r.start
				RE.r$rel.end <-RE.r$end-r.start
				
				par(fig=c(0,1,0.60,0.70),new=T)
				par(mar=c(0,5,0.1,2))
				
				plot(c(1,2*distance), c(1,60), type= "n", xlab="", ylab="",xaxt='n',yaxt='n')
				
				for(i in 1:nrow(RE.r)){
					if((i%%2)==0){
						i.start<-RE.r$rel.start[i]
						i.end  <-RE.r$rel.end[i]
						rect(i.start, 15, i.end, 25, col="blue")
					}else{
						i.start<-RE.r$rel.start[i]
						i.end  <-RE.r$rel.end[i]
						rect(i.start, 25, i.end,35, col="blue")
					}
				}
				legend("topleft",legend = "Restriction Fragments",fill="blue", cex=0.55,bty="n")
				abline(v=start(viewpoint)-r.start, col="red",lty=3,lwd=2)
				#######Draw Interaction regions##########
				exp.data<-expInteractions[start(expInteractions) >=r.start & end(expInteractions) <= r.end & as.character(seqnames(expInteractions))==viewpoint.chr,]
				exp.data$p_start<-start(exp.data)-start(viewpoint)
				
				exppalette<-rev(brewer.pal(9,"Reds"))
				
				exp.data$col[exp.data$q.value <=0.00001]<-exppalette[1]
				exp.data$col[exp.data$q.value > 0.00001 & exp.data$q.value <=0.0001]<-exppalette[2]
				exp.data$col[exp.data$q.value > 0.0001 & exp.data$q.value <=0.001]<-exppalette[3]
				exp.data$col[exp.data$q.value > 0.001 & exp.data$q.value <=0.01]<-exppalette[4]
				exp.data$col[exp.data$q.value > 0.01 & exp.data$q.value <=0.05]<-exppalette[5]
				exp.data$col[exp.data$q.value > 0.05 & exp.data$q.value <=0.1]<-exppalette[6]
				exp.data$col[exp.data$q.value > 0.1 & exp.data$q.value <=0.2]<-exppalette[7]
				exp.data$col[exp.data$q.value > 0.2 & exp.data$q.value <=0.5]<-exppalette[8]
				exp.data$col[exp.data$q.value > 0.5]<-exppalette[9]
				
				t.names<-c("q-value<=0.00001",
						"0.00001 >q-value<= 0.0001",
						"0.0001 >q-value<= 0.001",
						"0.001 >q-value<= 0.01",
						"0.01 >q-value<= 0.05",
						"0.05 >q-value<= 0.1",
						"0.1 >q-value<= 0.2",
						"0.2 >q-value<= 0.5",
						"q-value >0.5"
				)	
				par(fig=c(0,1,0.12,0.60),new=T)
				par(mar=c(4,5,0.1,2))
				y.exp.max <-ifelse(yLim >0,yLim,max(exp.data$nReads))
				
				plot(c(-1*(distance),distance),c(1,y.exp.max), type= "n", ylab="Reads",xlab="Distance (bp) relative to the viewpoint")
				points(exp.data$p_start,exp.data$nReads,col=exp.data$col,pch=19)	
				lines(exp.data$p_start,exp.data$nReads,col='black',lty=1)
				abline(v=0,lty=3,col="red",lwd=2)
				
				legend("topleft",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
				

			}
			if(isControlInvolved(obj)==TRUE){
				expInteractions   <-expInteractionRegions(obj)
				contrInteractions <-contrInteractionRegions(obj)	
				########Get labels################
				#expLabeled <-expLabel(obj)
				#contrLabeled <-contrLabel(obj)
				######look around the viewpoint###
				r.start <-start(viewpoint)-distance
				r.end 	<-end(viewpoint)+distance
				
				r.start <-ifelse(r.start <0,1,r.start)
				r.end   <-ifelse(r.end <=chr.size,r.end,chr.size)
				c.vec   <- c(rep(0,r.end-r.start))
				####Get refGenes######
				genes<-get3CseqRefGene(obj)
				g.chr<-subset(genes,chromosome==viewpoint.chr)
				g.r<-subset(g.chr,start>=r.start & end <= r.end)
				
				par(fig=c(0,1,0.8,1))
				par(mar=c(0,5,1,2))
				
				if(nrow(g.r)>0){
					g.r<-g.r[order(g.r$start),]
					g.r$rel.start <-g.r$start-r.start
					g.r$rel.end <-g.r$end-r.start
					g.r$size<-g.r$end-g.r$start+1
					
					plot(c(1,2*distance), c(1,60), type= "n", xlab="", ylab="",xaxt='n',yaxt='n')
					abline(v=start(viewpoint)-r.start, col="red",lty=3,lwd=2)
					y1.start=30
					for(i in 1:nrow(g.r)){
						gx<-g.r[i,]
						gx.start<-gx$rel.start
						gx.end <-gx$rel.end
						gx.size<-gx$size
						gx.strand<-gx$strand
						
						if(c.vec[gx$rel.start]==0){
							y1.start = y1.start
							c.vec[gx$rel.start:gx$rel.end]=1
						}else{
							y1.start = y1.start+8
							if(y1.start >50){
								y1.start=10
							}
						}
						text(gx$rel.start,y1.start+2,gx$name,cex =0.8)
						
						if(gx.size >=100){
							s.q <-ifelse(gx.size <=5000,500,5000)
							s.q <-ifelse(s.q > gx.size, gx.size-10,s.q)
							
							xx<-seq(gx.start,gx.start+gx.size-s.q, by=s.q)
							yy<-seq(gx.start+s.q,gx.start+gx.size, by=s.q)
							
							if(gx.strand==-1){
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-1))
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-4))
								}
							}else{
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-1,y1.start-2))
									lines(c(xx[i],yy[i]),c(y1.start-4,y1.start-3))
								}
							}
							rect(gx.start,y1.start-3,gx.start+gx.size,(y1.start-2), col="red")
						}	
					}
				}
				legend("topleft",legend = "Refseq Genes",fill="red", cex=0.55,bty="n")
				####Draw Restriction map####
				enzymeDb	<-new("repbaseEnzyme")
				resEnzyme 	<-restrictionEnzyme(obj)	
				chr.fragment <- getRestrictionFragments(enzymeDb,resEnzyme,orgName,viewpoint.chr)
				
				RE.r <-subset(chr.fragment,start >=r.start & end <=r.end)
				
				RE.r$rel.start <-RE.r$start-r.start
				RE.r$rel.end <-RE.r$end-r.start
				
				par(fig=c(0,1,0.75,0.80),new=T)
				par(mar=c(0,5,0.1,2))
				
				plot(c(1,2*distance), c(1,60), type= "n", xlab="", ylab="",xaxt='n',yaxt='n')
				for(i in 1:nrow(RE.r)){
					if((i%%2)==0){
						i.start<-RE.r$rel.start[i]
						i.end  <-RE.r$rel.end[i]
						rect(i.start, 15, i.end, 25, col="blue")
					}else{
						i.start<-RE.r$rel.start[i]
						i.end  <-RE.r$rel.end[i]
						rect(i.start, 25, i.end,35, col="blue")
					}
				}
				legend("topleft",legend = "Restriction Fragments",fill="blue", cex=0.55,bty="n")
				abline(v=start(viewpoint)-r.start, col="red",lty=3,lwd=2)
				#######Draw Interaction regions##########
				exp.data<-expInteractions[start(expInteractions) >=r.start & end(expInteractions) <= r.end & as.character(seqnames(expInteractions))==viewpoint.chr,]
				exp.data$p_start<-start(exp.data)-start(viewpoint)
				
				exppalette<-rev(brewer.pal(9,"Reds"))
				
				exp.data$col[exp.data$q.value <=0.00001]<-exppalette[1]
				exp.data$col[exp.data$q.value > 0.00001 & exp.data$q.value <=0.0001]<-exppalette[2]
				exp.data$col[exp.data$q.value > 0.0001 & exp.data$q.value <=0.001]<-exppalette[3]
				exp.data$col[exp.data$q.value > 0.001 & exp.data$q.value <=0.01]<-exppalette[4]
				exp.data$col[exp.data$q.value > 0.01 & exp.data$q.value <=0.05]<-exppalette[5]
				exp.data$col[exp.data$q.value > 0.05 & exp.data$q.value <=0.1]<-exppalette[6]
				exp.data$col[exp.data$q.value > 0.1 & exp.data$q.value <=0.2]<-exppalette[7]
				exp.data$col[exp.data$q.value > 0.2 & exp.data$q.value <=0.5]<-exppalette[8]
				exp.data$col[exp.data$q.value > 0.5]<-exppalette[9]
				
				
				contr.data<-contrInteractions[start(contrInteractions) >=r.start & end(contrInteractions) <= r.end & as.character(seqnames(contrInteractions))==viewpoint.chr,]
				contr.data$p_start<-start(contr.data)-start(viewpoint)
				
				contrpalette<-rev(brewer.pal(9,"Blues"))
				
				contr.data$col[contr.data$q.value <=0.00001]<-contrpalette[1]
				contr.data$col[contr.data$q.value > 0.00001 & contr.data$q.value <=0.0001]<-contrpalette[2]
				contr.data$col[contr.data$q.value > 0.0001 & contr.data$q.value <=0.001]<-contrpalette[3]
				contr.data$col[contr.data$q.value > 0.001 & contr.data$q.value <=0.01]<-contrpalette[4]
				contr.data$col[contr.data$q.value > 0.01 & contr.data$q.value <=0.05]<-contrpalette[5]
				contr.data$col[contr.data$q.value > 0.05 & contr.data$q.value <=0.1]<-contrpalette[6]
				contr.data$col[contr.data$q.value > 0.1 & contr.data$q.value <=0.2]<-contrpalette[7]
				contr.data$col[contr.data$q.value > 0.2 & contr.data$q.value <=0.5]<-contrpalette[8]
				contr.data$col[contr.data$q.value > 0.5]<-contrpalette[9]
							
				t.names<-c("q-value<=0.00001",
						"0.00001 >q-value<= 0.0001",
						"0.0001 >q-value<= 0.001",
						"0.001 >q-value<= 0.01",
						"0.01 >q-value<= 0.05",
						"0.05 >q-value<= 0.1",
						"0.1 >q-value<= 0.2",
						"0.2 >q-value<= 0.5",
						"q-value >0.5"
				)	
				par(fig=c(0,1,0.50,0.75),new=T)
				par(mar=c(0,5,0.1,2))
				y.exp.max <-ifelse(yLim >0,yLim,max(exp.data$RPMs))
		
				plot(c(-1*(distance),distance),c(1,y.exp.max), type= "n", ylab="RPMs",xaxt='n')
				points(exp.data$p_start,exp.data$RPMs,col=exp.data$col,pch=19)
				lines(exp.data$p_start,exp.data$RPMs,col='black',lty=1)
				
				legend("topleft",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
				legend("topright",legend = "experiment", cex=1,bty="n")
				abline(v=0,lty=3,col="red",lwd=2)
				
				par(fig=c(0,1,0.25,0.50),new=T)
				par(mar=c(0,5,0.1,2))
				y.exp.max <-ifelse(yLim >0,yLim,max(exp.data$RPMs))
				plot(c(-1*(distance),distance),c(1,y.exp.max), type= "n", ylab="RPMs",xaxt='n')				
				points(contr.data$p_start,contr.data$RPMs,col=contr.data$col,pch=19)	
				lines(contr.data$p_start,contr.data$RPMs,col='black',lty=1)
				legend("topleft",legend = t.names, fill=contrpalette, cex=0.55,title="q-value",bty="n")
				legend("topright",legend = "control", cex=1,bty="n")
				abline(v=0,lty=3,col="red",lwd=2)
				
				#####Draw log2 of subtraction#####
				exp.rpm<-data.frame(start=start(exp.data),end=end(exp.data),exp_RPMs=exp.data$RPMs)
				exp.rpm<-subset(exp.rpm,exp_RPMs>=1)
				contr.rpm<-data.frame(start=start(contr.data),end=end(contr.data),contr_RPMs=contr.data$RPMs)
				contr.rpm<-subset(contr.rpm,contr_RPMs>=1)
				combined.rpm<-merge(exp.rpm,contr.rpm,by.x=c("start","end"),by.y=c("start","end"),all=T)
				combined.rpm[is.na(combined.rpm)==TRUE]<-0
				combined.rpm$log2fold<-log2(combined.rpm$exp_RPMs+1)-log2(combined.rpm$contr_RPMs+1)
				combined.rpm<-subset(combined.rpm,abs(log2fold) >=log2fc_cutoff)
				combined.rpm$p_start<-combined.rpm$start-start(viewpoint)
				combined.rpm$p_end<-combined.rpm$end-start(viewpoint)
				
				par(fig=c(0,1,0,0.25),new=T)
				par(mar=c(4,5,0.1,2))
				
				plot(c(-1*(distance),distance),c(min(combined.rpm$log2fold),max(combined.rpm$log2fold)), type= "n", ylab="log2(Exp/Cont)",xlab="Distance (bp) from the viewpoint")
				
				for(i in 1:nrow(combined.rpm)){
					start<-combined.rpm$p_start[i]
					end  <-combined.rpm$p_end[i]
					log2.fold <-combined.rpm$log2fold[i]
					if(log2.fold >0){
						rect(start,0,end,log2.fold, col="red")
					}else{
						rect(start,0,end,log2.fold, col="green")
					}
				}
				abline(v=0,lty=3,col="red",lwd=2)
				abline(h=0,lty=1,col="black",lwd=1)
				abline(h=log2fc_cutoff,lty=3,col="black",lwd=1)
				abline(h=-log2fc_cutoff,lty=3,col="black",lwd=1)
				legend("bottomleft",legend =paste("log2fc >=",log2fc_cutoff,"& log2fc <=",-log2fc_cutoff), cex=0.7,bty="n")
			}
}

plotInteractionsPerChromosome<-function(obj,chromosomeName){
	
			stopifnot( is(obj, "r3Cseq") | is(obj,"r3CseqInBatch"))
			if(chromosomeName ==character(1)){
				stop("Require the chromosome name for example : 'chr1'")
			}
			if(!chromosomeName %in% paste('chr',c(seq(1,50),'X','Y','M'),sep='')){
				stop("Require the correct format chromosome name : 'chr1','chrX','chrY'")
			}
			
			if(isControlInvolved(obj)==FALSE){
				orgName<-organismName(obj)
				chr.size=c()
				
				if(orgName=="hg18"){
					genome <- BSgenome.Hsapiens.UCSC.hg18.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName=="hg19"){
					genome <- BSgenome.Hsapiens.UCSC.hg19.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName =="mm9"){
					genome <- BSgenome.Mmusculus.UCSC.mm9.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName =="mm10"){
					genome <- BSgenome.Mmusculus.UCSC.mm10.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName =="rn5"){
					genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else{
					stop("Your input organism name is not in the list ('mm9','mm10','hg18','hg19', and 'rn5')")
				}
				
				#expLabeled<-expLabel(obj)
				
				expInteractions   <-expInteractionRegions(obj)
				
				if(length(expInteractions) >0){
					viewpoint <-getViewpoint(obj)
					
					if(as.character(seqnames(viewpoint))==chromosomeName){
						chr.exp   <-expInteractions[as.character(seqnames(expInteractions))==chromosomeName,]
						if(length(chr.exp ) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}
						chr.exp.data  <-data.frame(start=start(chr.exp),nReads=chr.exp$nReads,q.value=chr.exp$q.value)
						e.smooth<-smooth.spline(chr.exp.data$nReads,spar=0.5)
						e.fitted<-fitted(e.smooth)
						chr.exp.data$expected<-e.fitted
						
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$q.value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$q.value > 0.00001 & chr.exp.data$q.value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$q.value > 0.0001 & chr.exp.data$q.value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$q.value > 0.001 & chr.exp.data$q.value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$q.value > 0.01 & chr.exp.data$q.value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$q.value > 0.05 & chr.exp.data$q.value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$q.value > 0.1 & chr.exp.data$q.value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$q.value > 0.2 & chr.exp.data$q.value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$q.value > 0.5]<-exppalette[9]
						
						t.names<-c("q-value<=0.00001",
								"0.00001 >q-value<= 0.0001",
								"0.0001 >q-value<= 0.001",
								"0.001 >q-value<= 0.01",
								"0.01 >q-value<= 0.05",
								"0.05 >q-value<= 0.1",
								"0.1 >q-value<= 0.2",
								"0.2 >q-value<= 0.5",
								"q-value >0.5"
						)							
						plot(chr.exp.data$start,chr.exp.data$nReads,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (bp)",
								ylab="Reads",
								main = paste(chromosomeName,":","experiment")
						)
						abline(v=start(viewpoint), col="black",lty=3)
						lines(chr.exp.data$start,chr.exp.data$expected,col="blue",lwd=2)
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
						legend("topleft",legend = "expected reads",lty=c(1),col="blue", cex=0.55,bty="n")
						
					}else{
						chr.exp <-expInteractions[as.character(seqnames(expInteractions))==chromosomeName,]
						if(length(chr.exp) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}	
						chr.exp.data  <-data.frame(start=start(chr.exp),nReads=chr.exp$nReads,q.value=chr.exp$q.value)
					
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$q.value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$q.value > 0.00001 & chr.exp.data$q.value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$q.value > 0.0001 & chr.exp.data$q.value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$q.value > 0.001 & chr.exp.data$q.value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$q.value > 0.01 & chr.exp.data$q.value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$q.value > 0.05 & chr.exp.data$q.value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$q.value > 0.1 & chr.exp.data$q.value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$q.value > 0.2 & chr.exp.data$q.value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$q.value > 0.5]<-exppalette[9]
						
						t.names<-c("q-value<=0.00001",
								"0.00001 >q-value<= 0.0001",
								"0.0001 >q-value<= 0.001",
								"0.001 >q-value<= 0.01",
								"0.01 >q-value<= 0.05",
								"0.05 >q-value<= 0.1",
								"0.1 >q-value<= 0.2",
								"0.2 >q-value<= 0.5",
								"q-value >0.5"
						)							
						plot(chr.exp.data$start,chr.exp.data$nReads,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (bp)",
								ylab="Reads",
								main = paste(chromosomeName,":","experiment")
						)
						
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
					
					}
					
				}else{
					stop("No interactions were found in your r3Cseq obj!!!")
				}
			}			
			if(isControlInvolved(obj)==TRUE){
				orgName<-organismName(obj)
				chr.size=c()
				
				if(orgName=="hg18"){
					genome <- BSgenome.Hsapiens.UCSC.hg18.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName=="hg19"){
					genome <- BSgenome.Hsapiens.UCSC.hg19.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName =="mm9"){
					genome <- BSgenome.Mmusculus.UCSC.mm9.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName =="mm10"){
					genome <- BSgenome.Mmusculus.UCSC.mm10.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else if(orgName =="rn5"){
					genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
					chr.size<-seqlengths(genome)[chromosomeName]
				}else{
					stop("Your input organism name is not in the list ('mm9','mm10','hg18','hg19', and 'rn5')")
				}
				
				#expLabeled<-expLabel(obj)
				#controlLabeled<-contrLabel(obj)
				
				expInteractions   <-expInteractionRegions(obj)
				contrInteractions <-contrInteractionRegions(obj)
				
				if(length(expInteractions) >0){
					
					viewpoint <-getViewpoint(obj)
					
					if(as.character(seqnames(viewpoint))==chromosomeName){
						chr.exp   <-expInteractions[as.character(seqnames(expInteractions))==chromosomeName,]
						if(length(chr.exp ) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}
						chr.contr <-contrInteractions[as.character(seqnames(contrInteractions))==chromosomeName,]
						
						chr.exp.data  <-data.frame(start=start(chr.exp),RPMs=chr.exp$RPMs,q.value=chr.exp$q.value)
						chr.contr.data <-data.frame(start=start(chr.contr),RPMs=chr.contr$RPMs,q.value=chr.contr$q.value)
						
						e.smooth<-smooth.spline(chr.exp.data$RPMs,spar=0.5)
						e.fitted<-fitted(e.smooth)
						chr.exp.data$expected<-e.fitted
						
						e.smooth<-smooth.spline(chr.contr.data$RPMs,spar=0.5)
						e.fitted<-fitted(e.smooth)
						chr.contr.data$expected<-e.fitted
						par(mfrow=c(2,1))
						
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$q.value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$q.value > 0.00001 & chr.exp.data$q.value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$q.value > 0.0001 & chr.exp.data$q.value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$q.value > 0.001 & chr.exp.data$q.value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$q.value > 0.01 & chr.exp.data$q.value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$q.value > 0.05 & chr.exp.data$q.value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$q.value > 0.1 & chr.exp.data$q.value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$q.value > 0.2 & chr.exp.data$q.value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$q.value > 0.5]<-exppalette[9]
						
						t.names<-c("q-value<=0.00001",
								"0.00001 >q-value<= 0.0001",
								"0.0001 >q-value<= 0.001",
								"0.001 >q-value<= 0.01",
								"0.01 >q-value<= 0.05",
								"0.05 >q-value<= 0.1",
								"0.1 >q-value<= 0.2",
								"0.2 >q-value<= 0.5",
								"q-value >0.5"
								)							
						y.max<-max(chr.exp.data$RPMs)
						plot(chr.exp.data$start,chr.exp.data$RPMs,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (bp)",
								ylab="Reads/Million",
								ylim=c(0,y.max),
								main = paste(chromosomeName,":","experiment")
								)
						abline(v=start(viewpoint), col="black",lty=3)
						lines(chr.exp.data$start,chr.exp.data$expected,col="blue",lwd=2)
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
						legend("topleft",legend = "expected reads",lty=c(1),col="blue", cex=0.55,bty="n")
						
						#####plot control###
						contrpalette<-rev(brewer.pal(9,"Blues"))
						
						chr.contr.data$col[chr.contr.data$q.value <=0.00001]<-contrpalette[1]
						chr.contr.data$col[chr.contr.data$q.value > 0.00001 & chr.contr.data$q.value <=0.0001]<-contrpalette[2]
						chr.contr.data$col[chr.contr.data$q.value > 0.0001 & chr.contr.data$q.value <=0.001]<-contrpalette[3]
						chr.contr.data$col[chr.contr.data$q.value > 0.001 & chr.contr.data$q.value <=0.01]<-contrpalette[4]
						chr.contr.data$col[chr.contr.data$q.value > 0.01 & chr.contr.data$q.value <=0.05]<-contrpalette[5]
						chr.contr.data$col[chr.contr.data$q.value > 0.05 & chr.contr.data$q.value <=0.1]<-contrpalette[6]
						chr.contr.data$col[chr.contr.data$q.value > 0.1 & chr.contr.data$q.value <=0.2]<-contrpalette[7]
						chr.contr.data$col[chr.contr.data$q.value > 0.2 & chr.contr.data$q.value <=0.5]<-contrpalette[8]
						chr.contr.data$col[chr.contr.data$q.value > 0.5]<-contrpalette[9]
						
						t.names<-c("q-value<=0.00001",
								"0.00001 >q-value<= 0.0001",
								"0.0001 >q-value<= 0.001",
								"0.001 >q-value<= 0.01",
								"0.01 >q-value<= 0.05",
								"0.05 >q-value<= 0.1",
								"0.1 >q-value<= 0.2",
								"0.2 >q-value<= 0.5",
								"q-value >0.5"
						)							
						plot(chr.contr.data$start,chr.contr.data$RPMs,pch=19,
								col=chr.contr.data$col,
								xlab="Chromosomal position (bp)",
								ylab="Reads/Million",
								ylim=c(0,y.max),
								main = paste(chromosomeName,":","control")
						)
						abline(v=start(viewpoint), col="black",lty=3)
						lines(chr.contr.data$start,chr.contr.data$expected,col="blue",lwd=2)
						legend("topright",legend = t.names, fill=contrpalette, cex=0.55,title="q-value",bty="n")
						legend("topleft",legend = "expected reads",lty=c(1),col="blue", cex=0.55,bty="n")
						
					}else{
						chr.exp   <-expInteractions[as.character(seqnames(expInteractions))==chromosomeName,]
						chr.contr <-contrInteractions[as.character(seqnames(contrInteractions))==chromosomeName,]
						
						if(length(chr.exp) ==0){
							stop("There is no interaction regions found in your selected chromosome!!.")
						}
						
						chr.exp.data  <-data.frame(start=start(chr.exp),RPMs=chr.exp$RPMs,q.value=chr.exp$q.value)
						chr.contr.data <-data.frame(start=start(chr.contr),RPMs=chr.contr$RPMs,q.value=chr.contr$q.value)
						
						par(mfrow=c(2,1))
						
						exppalette<-rev(brewer.pal(9,"Reds"))
						
						chr.exp.data$col[chr.exp.data$q.value <=0.00001]<-exppalette[1]
						chr.exp.data$col[chr.exp.data$q.value > 0.00001 & chr.exp.data$q.value <=0.0001]<-exppalette[2]
						chr.exp.data$col[chr.exp.data$q.value > 0.0001 & chr.exp.data$q.value <=0.001]<-exppalette[3]
						chr.exp.data$col[chr.exp.data$q.value > 0.001 & chr.exp.data$q.value <=0.01]<-exppalette[4]
						chr.exp.data$col[chr.exp.data$q.value > 0.01 & chr.exp.data$q.value <=0.05]<-exppalette[5]
						chr.exp.data$col[chr.exp.data$q.value > 0.05 & chr.exp.data$q.value <=0.1]<-exppalette[6]
						chr.exp.data$col[chr.exp.data$q.value > 0.1 & chr.exp.data$q.value <=0.2]<-exppalette[7]
						chr.exp.data$col[chr.exp.data$q.value > 0.2 & chr.exp.data$q.value <=0.5]<-exppalette[8]
						chr.exp.data$col[chr.exp.data$q.value > 0.5]<-exppalette[9]
						
						t.names<-c("q-value<=0.00001",
								"0.00001 >q-value<= 0.0001",
								"0.0001 >q-value<= 0.001",
								"0.001 >q-value<= 0.01",
								"0.01 >q-value<= 0.05",
								"0.05 >q-value<= 0.1",
								"0.1 >q-value<= 0.2",
								"0.2 >q-value<= 0.5",
								"q-value >0.5"
						)							
						plot(chr.exp.data$start,chr.exp.data$RPMs,pch=19,
								col=chr.exp.data$col,
								xlab="Chromosomal position (bp)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":","experiment")
						)
						
						legend("topright",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
						
						#####plot control###
						contrpalette<-rev(brewer.pal(9,"Blues"))
						
						chr.contr.data$col[chr.contr.data$q.value <=0.00001]<-contrpalette[1]
						chr.contr.data$col[chr.contr.data$q.value > 0.00001 & chr.contr.data$q.value <=0.0001]<-contrpalette[2]
						chr.contr.data$col[chr.contr.data$q.value > 0.0001 & chr.contr.data$q.value <=0.001]<-contrpalette[3]
						chr.contr.data$col[chr.contr.data$q.value > 0.001 & chr.contr.data$q.value <=0.01]<-contrpalette[4]
						chr.contr.data$col[chr.contr.data$q.value > 0.01 & chr.contr.data$q.value <=0.05]<-contrpalette[5]
						chr.contr.data$col[chr.contr.data$q.value > 0.05 & chr.contr.data$q.value <=0.1]<-contrpalette[6]
						chr.contr.data$col[chr.contr.data$q.value > 0.1 & chr.contr.data$q.value <=0.2]<-contrpalette[7]
						chr.contr.data$col[chr.contr.data$q.value > 0.2 & chr.contr.data$q.value <=0.5]<-contrpalette[8]
						chr.contr.data$col[chr.contr.data$q.value > 0.5]<-contrpalette[9]
						
						t.names<-c("q-value<=0.00001",
								"0.00001 >q-value<= 0.0001",
								"0.0001 >q-value<= 0.001",
								"0.001 >q-value<= 0.01",
								"0.01 >q-value<= 0.05",
								"0.05 >q-value<= 0.1",
								"0.1 >q-value<= 0.2",
								"0.2 >q-value<= 0.5",
								"q-value >0.5"
						)							
						plot(chr.contr.data$start,chr.contr.data$RPMs,pch=19,
								col=chr.contr.data$col,
								xlab="Chromosomal position (bp)",
								ylab="Reads/Million",
								main = paste(chromosomeName,":","control")
						)
						
						legend("topright",legend = t.names, fill=contrpalette, cex=0.55,title="q-value",bty="n")
					}
						
				}else{
					stop("No interactions were found in your r3Cseq object!!!")
				}
				
			}
}

setGeneric(
		name="plotDomainogramNearViewpoint",
		def=function(object,smoothing.parameter=0.1,distance=5e5,maximum_window=25e3, view=c("experiment","control","both")){
			standardGeneric("plotDomainogramNearViewpoint")
		}
)

setMethod("plotDomainogramNearViewpoint",
		signature(object = "r3Cseq"),
		
		function (object,smoothing.parameter,distance,maximum_window,view){
			
			if(!is(object,"r3Cseq")){
				stop("Need to initialize the r3Cseq object")
				
			}
			
			if(distance < 100000 | distance > 500000){
				stop("Please select distance between 100Kb - 500Kb")
			}
			if(maximum_window > 50e3){
				stop("The maximum allow for windowing is 50Kb.")
			}
			########Get input################################
			selected_view <- match.arg(view)
			if(selected_view==""){
				selected_view<-"experiment"
			}
			########Get viewpoint############
			viewpoint <-getViewpoint(object)
			viewpoint.chr<-as.character(seqnames(viewpoint))
			########Get organism#############
			orgName<-organismName(object)
			chr.size<-0
			if(orgName=="hg18"){
				genome <- BSgenome.Hsapiens.UCSC.hg18.masked
				chr.size<-seqlengths(genome)[viewpoint.chr]
			
			}else if(orgName=="hg19"){
				genome <- BSgenome.Hsapiens.UCSC.hg19.masked
				chr.size<-seqlengths(genome)[viewpoint.chr]
			
			}else if(orgName =="mm9"){
				genome <- BSgenome.Mmusculus.UCSC.mm9.masked
				chr.size<-seqlengths(genome)[viewpoint.chr]
				
			}else if(orgName =="mm10"){
					genome <- BSgenome.Mmusculus.UCSC.mm10.masked
					chr.size<-seqlengths(genome)[viewpoint.chr]
					
			}else if(orgName =="rn5"){
					genome <- BSgenome.Rnorvegicus.UCSC.rn5.masked
					chr.size<-seqlengths(genome)[viewpoint.chr]
			}
			if(selected_view=="experiment"){
				expInteractions   <-expInteractionRegions(object)
				######look around the viewpoint###
				r.start <-start(viewpoint)-distance
				r.end 	<-end(viewpoint)+distance
				
				r.start <-ifelse(r.start <0,1,r.start)
				r.end   <-ifelse(r.end <=chr.size,r.end,chr.size)
				c.vec   <- c(rep(0,r.end-r.start))
				####Get refGenes######
				genes<-get3CseqRefGene(object)
				g.chr<-subset(genes,chromosome==viewpoint.chr)
				g.r<-subset(g.chr,start>=r.start & end <= r.end)
				
				par(fig=c(0,1,0.75,1))
				par(mar=c(0,5,1,2))
				
				if(nrow(g.r)>0){
					g.r<-g.r[order(g.r$start),]
					g.r$rel.start <-g.r$start-r.start
					g.r$rel.end <-g.r$end-r.start
					g.r$size<-g.r$end-g.r$start+1
					
					plot(c(1,2*distance), c(1,60), type= "n", xlab="", ylab="",xaxt='n',yaxt='n')
					abline(v=start(viewpoint)-r.start, col="red",lty=3,lwd=2)
					y1.start=30
					for(i in 1:nrow(g.r)){
						gx<-g.r[i,]
						gx.start<-gx$rel.start
						gx.end <-gx$rel.end
						gx.size<-gx$size
						gx.strand<-gx$strand
						
						if(c.vec[gx$rel.start]==0){
							y1.start = y1.start
							c.vec[gx$rel.start:gx$rel.end]=1
						}else{
							y1.start = y1.start+8
							if(y1.start >50){
								y1.start=10
							}
						}
						text(gx$rel.start,y1.start+2,gx$name,cex =0.8)
						
						if(gx.size >=100){
							s.q <-ifelse(gx.size <=5000,500,5000)
							s.q <-ifelse(s.q > gx.size, gx.size-10,s.q)
							
							xx<-seq(gx.start,gx.start+gx.size-s.q, by=s.q)
							yy<-seq(gx.start+s.q,gx.start+gx.size, by=s.q)
							
							if(gx.strand==-1){
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-1))
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-4))
								}
							}else{
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-1,y1.start-2))
									lines(c(xx[i],yy[i]),c(y1.start-4,y1.start-3))
								}
							}
							rect(gx.start,y1.start-3,gx.start+gx.size,(y1.start-2), col="red")
						}	
					}
				}
				legend("topleft",legend = "Refseq Genes",fill="red", cex=0.55,bty="n")
				#######Draw Interaction regions##########
				exp.data<-expInteractions[start(expInteractions) >=r.start & end(expInteractions) <= r.end & as.character(seqnames(expInteractions))==viewpoint.chr,]
				exp.data$p_start<-start(exp.data)-start(viewpoint)
				
				exppalette<-rev(brewer.pal(9,"Reds"))
				
				exp.data$col[exp.data$q.value <=0.00001]<-exppalette[1]
				exp.data$col[exp.data$q.value > 0.00001 & exp.data$q.value <=0.0001]<-exppalette[2]
				exp.data$col[exp.data$q.value > 0.0001 & exp.data$q.value <=0.001]<-exppalette[3]
				exp.data$col[exp.data$q.value > 0.001 & exp.data$q.value <=0.01]<-exppalette[4]
				exp.data$col[exp.data$q.value > 0.01 & exp.data$q.value <=0.05]<-exppalette[5]
				exp.data$col[exp.data$q.value > 0.05 & exp.data$q.value <=0.1]<-exppalette[6]
				exp.data$col[exp.data$q.value > 0.1 & exp.data$q.value <=0.2]<-exppalette[7]
				exp.data$col[exp.data$q.value > 0.2 & exp.data$q.value <=0.5]<-exppalette[8]
				exp.data$col[exp.data$q.value > 0.5]<-exppalette[9]
				
				t.names<-c("q-value<=0.00001",
						"0.00001 >q-value<= 0.0001",
						"0.0001 >q-value<= 0.001",
						"0.001 >q-value<= 0.01",
						"0.01 >q-value<= 0.05",
						"0.05 >q-value<= 0.1",
						"0.1 >q-value<= 0.2",
						"0.2 >q-value<= 0.5",
						"q-value >0.5"
				)	
				par(fig=c(0,1,0.35,0.75),new=T)
				par(mar=c(4,5,0.1,2))
				y.exp.max <-max(exp.data$nReads)
				plot(c(-1*(distance),distance),c(1,y.exp.max), type= "n", ylab="Reads",xaxt="n",xlab="")
				points(exp.data$p_start,exp.data$nReads,col=exp.data$col,pch=19)	
				lines(exp.data$p_start,exp.data$nReads,col='black',lty=1)
				abline(v=0,lty=3,col="red",lwd=2)
				
				legend("topleft",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
				
				###########Draw Domainogram################
				expRawReads<-expRawData(object)
				expInteractionMatrix<-makeInteractionMatrixNearCisPerWindow(object,smoothing.parameter,expRawReads,max.window=maximum_window,viewpoint,distanceFromViewpoint=distance)
				col.n<-ncol(expInteractionMatrix)
				spec = colorRampPalette(c(rev(brewer.pal(5,"Reds")),"gray"))
				shades = spec(300)
				
				par(fig=c(0,1,0,0.44),new=T)
				par(mar=c(4,5,0.1,2))
				
				s.names<-c("very high interaction",
						"high interaction",
						"moderate interaction",
						"moderate-low interaction",
						"low interaction",
						"no interaction"
				)	
				image(as.matrix(log2(1e-10+expInteractionMatrix[,2:col.n])),xaxt="n",yaxt="n",ylab="Window (Kb)",col=shades,xlab="Distance (bp) relative to the viewpoint")
				abline(v=0.5,lty=3,col="red",lwd=2)
				legend("topleft",legend = s.names, fill=c(rev(brewer.pal(5,"Reds")),"gray"), cex=0.55,title="",bty="n")
				axis(1, at=c(seq(0,1,by=0.1)),labels=c(seq(-distance,distance,distance/5)),cex.axis=0.8,las=2)
				axis(2, at=c(0, 1),labels=c(paste(maximum_window/1000,"Kb"),"2 Kb"),cex.axis=0.8,las=2)
			}
			if(selected_view=="control"){
				contrInteractions <-contrInteractionRegions(object)	
				######look around the viewpoint###
				r.start <-start(viewpoint)-distance
				r.end 	<-end(viewpoint)+distance
				
				r.start <-ifelse(r.start <0,1,r.start)
				r.end   <-ifelse(r.end <=chr.size,r.end,chr.size)
				c.vec   <- c(rep(0,r.end-r.start))
				####Get refGenes######
				genes<-get3CseqRefGene(object)
				g.chr<-subset(genes,chromosome==viewpoint.chr)
				g.r<-subset(g.chr,start>=r.start & end <= r.end)
				
				par(fig=c(0,1,0.75,1))
				par(mar=c(0,5,1,2))
				
				if(nrow(g.r)>0){
					g.r<-g.r[order(g.r$start),]
					g.r$rel.start <-g.r$start-r.start
					g.r$rel.end <-g.r$end-r.start
					g.r$size<-g.r$end-g.r$start+1
					
					plot(c(1,2*distance), c(1,60), type= "n", xlab="", ylab="",xaxt='n',yaxt='n')
					abline(v=start(viewpoint)-r.start, col="red",lty=3,lwd=2)
					y1.start=30
					for(i in 1:nrow(g.r)){
						gx<-g.r[i,]
						gx.start<-gx$rel.start
						gx.end <-gx$rel.end
						gx.size<-gx$size
						gx.strand<-gx$strand
						
						if(c.vec[gx$rel.start]==0){
							y1.start = y1.start
							c.vec[gx$rel.start:gx$rel.end]=1
						}else{
							y1.start = y1.start+8
							if(y1.start >50){
								y1.start=10
							}
						}
						text(gx$rel.start,y1.start+2,gx$name,cex =0.8)
						
						if(gx.size >=100){
							s.q <-ifelse(gx.size <=5000,500,5000)
							s.q <-ifelse(s.q > gx.size, gx.size-10,s.q)
							
							xx<-seq(gx.start,gx.start+gx.size-s.q, by=s.q)
							yy<-seq(gx.start+s.q,gx.start+gx.size, by=s.q)
							
							if(gx.strand==-1){
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-1))
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-4))
								}
							}else{
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-1,y1.start-2))
									lines(c(xx[i],yy[i]),c(y1.start-4,y1.start-3))
								}
							}
							rect(gx.start,y1.start-3,gx.start+gx.size,(y1.start-2), col="red")
						}	
					}
				}
				legend("topleft",legend = "Refseq Genes",fill="red", cex=0.55,bty="n")
				#######Draw Interaction regions##########
				exp.data<-expInteractions[start(expInteractions) >=r.start & end(expInteractions) <= r.end & as.character(seqnames(expInteractions))==viewpoint.chr,]
				exp.data$p_start<-start(exp.data)-start(viewpoint)
				
				exppalette<-rev(brewer.pal(9,"Blues"))
				
				exp.data$col[exp.data$q.value <=0.00001]<-exppalette[1]
				exp.data$col[exp.data$q.value > 0.00001 & exp.data$q.value <=0.0001]<-exppalette[2]
				exp.data$col[exp.data$q.value > 0.0001 & exp.data$q.value <=0.001]<-exppalette[3]
				exp.data$col[exp.data$q.value > 0.001 & exp.data$q.value <=0.01]<-exppalette[4]
				exp.data$col[exp.data$q.value > 0.01 & exp.data$q.value <=0.05]<-exppalette[5]
				exp.data$col[exp.data$q.value > 0.05 & exp.data$q.value <=0.1]<-exppalette[6]
				exp.data$col[exp.data$q.value > 0.1 & exp.data$q.value <=0.2]<-exppalette[7]
				exp.data$col[exp.data$q.value > 0.2 & exp.data$q.value <=0.5]<-exppalette[8]
				exp.data$col[exp.data$q.value > 0.5]<-exppalette[9]
				
				t.names<-c("q-value<=0.00001",
						"0.00001 >q-value<= 0.0001",
						"0.0001 >q-value<= 0.001",
						"0.001 >q-value<= 0.01",
						"0.01 >q-value<= 0.05",
						"0.05 >q-value<= 0.1",
						"0.1 >q-value<= 0.2",
						"0.2 >q-value<= 0.5",
						"q-value >0.5"
				)	
				par(fig=c(0,1,0.35,0.75),new=T)
				par(mar=c(4,5,0.1,2))
				y.exp.max <-max(exp.data$nReads)
				plot(c(-1*(distance),distance),c(1,y.exp.max), type= "n", ylab="Reads",xaxt="n",xlab="")
				points(exp.data$p_start,exp.data$nReads,col=exp.data$col,pch=19)	
				lines(exp.data$p_start,exp.data$nReads,col='black',lty=1)
				abline(v=0,lty=3,col="red",lwd=2)
				
				legend("topleft",legend = t.names, fill=exppalette, cex=0.55,title="q-value",bty="n")
				
				###########Draw Domainogram################
				contrRawReads<-contrRawData(object)
				contrInteractionMatrix<-makeInteractionMatrixNearCisPerWindow(object,smoothing.parameter,contrRawReads,max.window=maximum_window,viewpoint,distanceFromViewpoint=distance)
				col.n<-ncol(contrInteractionMatrix)
				
				spec = colorRampPalette(c(rev(brewer.pal(5,"Blues")),"gray"))
				shades = spec(300)
				
				par(fig=c(0,1,0,0.44),new=T)
				par(mar=c(4,5,0.1,2))
				
				s.names<-c("very high interaction",
						"high interaction",
						"moderate interaction",
						"moderate-low interaction",
						"low interaction",
						"no interaction"
				)	
				image(as.matrix(log2(1e-10+contrInteractionMatrix[,2:col.n])),xaxt="n",yaxt="n",ylab="Window (Kb)",col=shades,xlab="Distance (bp) relative to the viewpoint")
				abline(v=0.5,lty=3,col="red",lwd=2)
				legend("topleft",legend = s.names, fill=c(rev(brewer.pal(5,"Blues")),"gray"), cex=0.55,title="",bty="n")
				axis(1, at=c(seq(0,1,by=0.1)),labels=c(seq(-distance,distance,distance/5)),cex.axis=0.8,las=2)
				axis(2, at=c(0, 1),labels=c(paste(maximum_window/1000,"Kb"),"2 Kb"),cex.axis=0.8,las=2)
			}
			if(selected_view=="both"){
				######look around the viewpoint###
				r.start <-start(viewpoint)-distance
				r.end 	<-end(viewpoint)+distance
				
				r.start <-ifelse(r.start <0,1,r.start)
				r.end   <-ifelse(r.end <=chr.size,r.end,chr.size)
				c.vec   <- c(rep(0,r.end-r.start))
				####Get refGenes######
				genes<-get3CseqRefGene(object)
				g.chr<-subset(genes,chromosome==viewpoint.chr)
				g.r<-subset(g.chr,start>=r.start & end <= r.end)
				
				par(fig=c(0,1,0.78,1))
				par(mar=c(0,5,1,2))
				
				if(nrow(g.r)>0){
					g.r<-g.r[order(g.r$start),]
					g.r$rel.start <-g.r$start-r.start
					g.r$rel.end <-g.r$end-r.start
					g.r$size<-g.r$end-g.r$start+1
					
					plot(c(1,2*distance), c(1,60), type= "n", xlab="", ylab="",xaxt='n',yaxt='n')
					abline(v=start(viewpoint)-r.start, col="red",lty=3,lwd=2)
					y1.start=30
					for(i in 1:nrow(g.r)){
						gx<-g.r[i,]
						gx.start<-gx$rel.start
						gx.end <-gx$rel.end
						gx.size<-gx$size
						gx.strand<-gx$strand
						
						if(c.vec[gx$rel.start]==0){
							y1.start = y1.start
							c.vec[gx$rel.start:gx$rel.end]=1
						}else{
							y1.start = y1.start+8
							if(y1.start >50){
								y1.start=10
							}
						}
						text(gx$rel.start,y1.start+2,gx$name,cex =0.8)
						
						if(gx.size >=100){
							s.q <-ifelse(gx.size <=5000,500,5000)
							s.q <-ifelse(s.q > gx.size, gx.size-10,s.q)
							
							xx<-seq(gx.start,gx.start+gx.size-s.q, by=s.q)
							yy<-seq(gx.start+s.q,gx.start+gx.size, by=s.q)
							
							if(gx.strand==-1){
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-1))
									lines(c(xx[i],yy[i]),c(y1.start-2,y1.start-4))
								}
							}else{
								for(i in 1:length(xx)){
									lines(c(xx[i],yy[i]),c(y1.start-1,y1.start-2))
									lines(c(xx[i],yy[i]),c(y1.start-4,y1.start-3))
								}
							}
							rect(gx.start,y1.start-3,gx.start+gx.size,(y1.start-2), col="red")
						}	
					}
				}
				legend("topleft",legend = "Refseq Genes",fill="red", cex=0.55,bty="n")
				###########Draw Domainogram in the experiment################
				expRawReads<-expRawData(object)
				expInteractionMatrix<-makeInteractionMatrixNearCisPerWindow(object,smoothing.parameter,expRawReads,max.window=maximum_window,viewpoint,distanceFromViewpoint=distance)
				col.n<-ncol(expInteractionMatrix)
				spec = colorRampPalette(c(rev(brewer.pal(5,"Reds")),"gray"))
				shades = spec(300)
				
				par(fig=c(0,1,0.45,0.78),new=T)
				par(mar=c(0,5,0.1,2))
				
				s.names<-c("very high interaction",
						"high interaction",
						"moderate interaction",
						"moderate-low interaction",
						"low interaction",
						"no interaction"
				)	
				image(as.matrix(log2(1e-10+expInteractionMatrix[,2:col.n])),xaxt="n",yaxt="n",ylab="Window (Kb)",col=shades,xlab="")
				abline(v=0.5,lty=3,col="red",lwd=2)
				legend("topleft",legend = s.names, fill=c(rev(brewer.pal(5,"Reds")),"gray"), cex=0.55,title="experiment",bty="n")
				axis(2, at=c(0, 1),labels=c(paste(maximum_window/1000,"Kb"),"2 Kb"),cex.axis=0.8,las=2)
				
				
				###########Draw Domainogram in the control################
				contrRawReads<-contrRawData(object)
				contrInteractionMatrix<-makeInteractionMatrixNearCisPerWindow(object,smoothing.parameter,contrRawReads,max.window=maximum_window,viewpoint,distanceFromViewpoint=distance)
				col.n<-ncol(contrInteractionMatrix)
				
				spec = colorRampPalette(c(rev(brewer.pal(5,"Blues")),"gray"))
				shades = spec(300)
				
				par(fig=c(0,1,0,0.44),new=T)
				par(mar=c(4,5,0.1,2))
				
				image(as.matrix(log2(1e-10+contrInteractionMatrix[,2:col.n])),xaxt="n",yaxt="n",ylab="Window (Kb)",col=shades,xlab="Distance (bp) relative to the viewpoint")
				abline(v=0.5,lty=3,col="red",lwd=2)
				legend("topleft",legend = s.names, fill=c(rev(brewer.pal(5,"Blues")),"gray"), cex=0.55,title="control",bty="n")
				axis(1, at=c(seq(0,1,by=0.1)),labels=c(seq(-distance,distance,distance/5)),cex.axis=0.8,las=2)
				axis(2, at=c(0, 1),labels=c(paste(maximum_window/1000,"Kb"),"2 Kb"),cex.axis=0.8,las=2)	
			}	
		}
)

###########
setGeneric(
		name="plot3Cecdf",
		def=function(object){
			standardGeneric("plot3Cecdf")
		}
)
setMethod("plot3Cecdf",
		signature(object = "r3Cseq"),
		function(object){
			stop( "The function 'plot3Cecdf' has been removed." )
			
		}
)
