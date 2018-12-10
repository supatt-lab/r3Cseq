### R code from vignette source 'r3Cseq.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 60)
olocale=Sys.setlocale(locale="C")


###################################################
### code chunk number 2: r3Cseq.Rnw:208-209
###################################################
library(r3Cseq)


###################################################
### code chunk number 3: r3Cseq.Rnw:212-214
###################################################
data(Myb_prom_FL)
data(Myb_prom_FB)


###################################################
### code chunk number 4: r3Cseq.Rnw:233-237
###################################################
my3Cseq.obj<-new("r3Cseq",organismName='mm9',isControlInvolved=TRUE,
viewpoint_chromosome='chr10',viewpoint_primer_forward='TCTTTGTTTGATGGCATCTGTT',
viewpoint_primer_reverse='AAAGGGGAGGAGAAGGAGGT',expLabel="Myb_prom_FL",
contrLabel="MYb_prom_FB",restrictionEnzyme='HindIII')


###################################################
### code chunk number 5: r3Cseq.Rnw:241-243
###################################################
expRawData(my3Cseq.obj)<-exp.GRanges
contrRawData(my3Cseq.obj)<-contr.GRanges	


###################################################
### code chunk number 6: r3Cseq.Rnw:246-247
###################################################
my3Cseq.obj	


###################################################
### code chunk number 7: r3Cseq.Rnw:252-253
###################################################
getReadCountPerRestrictionFragment(my3Cseq.obj)


###################################################
### code chunk number 8: r3Cseq.Rnw:261-262
###################################################
calculateRPM(my3Cseq.obj)	


###################################################
### code chunk number 9: r3Cseq.Rnw:266-267
###################################################
getInteractions(my3Cseq.obj,fdr=0.05)	


###################################################
### code chunk number 10: r3Cseq.Rnw:272-274
###################################################
fetal.liver.interactions<-expInteractionRegions(my3Cseq.obj)
fetal.liver.interactions


###################################################
### code chunk number 11: r3Cseq.Rnw:277-279
###################################################
fetal.brain.interactions<-contrInteractionRegions(my3Cseq.obj)
fetal.brain.interactions


###################################################
### code chunk number 12: r3Cseq.Rnw:284-286
###################################################
viewpoint<-getViewpoint(my3Cseq.obj)
viewpoint


###################################################
### code chunk number 13: plotOverviewInteractions
###################################################
plotOverviewInteractions(my3Cseq.obj)	


###################################################
### code chunk number 14: plotInteractionsNearViewpoint
###################################################
plotInteractionsNearViewpoint(my3Cseq.obj)	


###################################################
### code chunk number 15: plotInteractionsPerChromosome
###################################################
plotInteractionsPerChromosome(my3Cseq.obj,"chr10")	


###################################################
### code chunk number 16: r3Cseq.Rnw:335-336
###################################################
#plotDomainogramNearViewpoint(my3Cseq.obj)


###################################################
### code chunk number 17: r3Cseq.Rnw:341-343
###################################################
detected_genes<-getExpInteractionsInRefseq(my3Cseq.obj)
head(detected_genes)


###################################################
### code chunk number 18: r3Cseq.Rnw:349-350
###################################################
#export3Cseq2bedGraph(my3Cseq.obj)	


###################################################
### code chunk number 19: r3Cseq.Rnw:358-359
###################################################
#generate3CseqReport(my3Cseq.obj)	


###################################################
### code chunk number 20: r3Cseq.Rnw:376-377
###################################################
sessionInfo()


