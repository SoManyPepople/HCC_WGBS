#All figures
setwd("/data/huangp/Methy/Mixed/New")
options(stringsAsFactors = F)
library(fastSave)
################### part 1 Age-Dependent Global hypomethylation in HCC ################
# 1.1 The depth distribution of all CpGs 
if(1<0) {fileNames <- list.files("/data/huangp/Methy/methylation")
library(foreach)
library(doParallel)
#registerDoParallel(cores=16)
system.time(Meth.list <- foreach(i = 1:length(fileNames)) %do%
              {
                file = paste("/data/huangp/Methy/methylation/", fileNames[i], sep = "")
                methy <- read.table(file, sep="\t", stringsAsFactors = F, header=T)
                methy <- methy[, c(1,2,8,7)]
                names(methy) <- c("chr", "pos", "N", "X") 
                methy
                #methy$pos <- as.integer(methy$pos)
                #methy$N <- as.numeric(methy$N)
                #methy$X <- as.integer(methy$X)
                #methy[!is.na(methy$pos),]
              })

library(DSS)
load.lbzip2("Meth.list.unQC.RDataFS", n.cores=64)
system.time(BSobj <- makeBSseqData(Meth.list, names(Meth.list)))

library(DSS)
load.lbzip2("../BSobj.unsm.unQC.RDataFS", n.cores=24)
coverage <- as.matrix(getCoverage(BSobj, type="Cov"))
ranges <- as.data.frame(granges(BSobj))
rownames(coverage) <- paste(ranges$seqnames, ranges$start, sep="_")
save.lbzip2(coverage, file="Coverage_of_all_CpGs_in_each_sample.RDataFS", n.cores=32)
## calculate the total CpGs, mean depth of each WGBS samples
NoCpGs <- colSums(coverage != 0)
MeanDepth <- colSums(coverage)/NoCpGs
WGBS_info <- data.frame(SampleID = colnames(coverage),NoCpGs = NoCpGs, MeanDepth = MeanDepth)
write.csv(WGBS_info, file="Supplementary_tables_WGBS_info_summary.csv", row.names = F)

Meth.all <- as.data.frame(getMeth(BSobj, type="raw"))
rownames(Meth.all) <- paste(ranges$seqnames, ranges$start, sep="_")
save.lbzip2(Meth.all, file="NewRaw.Methylation.unmasked.RDataFS", n.cores=24)
#set methylation level to NA for coverage <5 in each sample
library(foreach)
library(doParallel)
registerDoParallel(cl=24, cores = 4)
Meth.raw <- foreach(i = 1:nrow(Meth.all), .combine = 'rbind')%dopar%
  {
    methy <- Meth.all[i,]
    cov <- coverage[i,]
    methy[cov<5] <- NA
    methy
  }
stopImplicitCluster()
rownames(Meth.raw) <- rownames(Meth.all)
save.lbzip2(Meth.raw, file="NewRaw.Methylation.masked.RDataFS", n.cores=24)
}
#load coverage of all CpGs in all samples
load.lbzip2("Coverage_of_all_CpGs_in_each_sample.RDataFS", n.cores=24)
NA_count <- function(depth){
  return(sum(depth<5))
}
#NA_count_tumor <- apply(coverage[,grep(pattern = "T", x = colnames(coverage), fixed=T)], 1, NA_count)
#NA_count_adjacent <- apply(coverage[,grep(pattern = "N", x = colnames(coverage), fixed=T)], 1, NA_count)
system.time(NA_count_all <- apply(coverage, 1, NA_count))#length(NA_count_all) = 28,978,826
# distribution of depth of all CpGs with missing rate < 100%
coverage_100 <- coverage[NA_count_all < 66,] #28,355,551
mean_depth_100 <- rowMeans(coverage_100)
#ECDF_mean_depth_100 <- ecdf(mean_depth_100)
#1-ECDF_mean_depth_100(5)#0.8812
#coverage_50 <- coverage[NA_count_all < 66*0.5,] #25,193,103
#mean_depth_50 <- rowMeans(coverage_50)
#coverage_20 <- coverage[NA_count_all < 66*0.2,] #21,927,292
#mean_depth_20 <- rowMeans(coverage_20)
#coverage_0 <- coverage[NA_count_all == 0,] #10,949,657
#mean_depth_0 <- rowMeans(coverage_0)
# remove CpGs with depth > 50 in plotting (N = 42,188)
df <- data.frame(MeanDepth = mean_depth_100[mean_depth_100<=50],
                 MissingRateCutoff = rep("100%", sum(mean_depth_100<=50)))


library(ggplot2)
library(ggpubr)
library(RColorBrewer)
p_1.1<- ggplot(df, aes(x = MeanDepth, colour = MissingRateCutoff)) + 
  stat_ecdf(geom = "step", size = 1.2, show.legend = F)+
  scale_color_manual(values="Red") +
  geom_vline(xintercept = 12.76,color = "blue", lty = 2)+
  labs(title="Empirical Cumulative Density",
       y = "F(AverageDepth)", x="Average depth of each CpG site")+
  theme_pubr()+theme(text=element_text(size=18))
p_1.1
dev.off()
save.lbzip2(p_1.1, file="./Figures/Figure_1.1.RDataFS", n.cores=24)

#1.2 The PCA of all samples base on smoothed DNA methylation level of all CpGs 
setwd("/data/huangp/Methy/Mixed/New")
options(stringsAsFactors = F)
library(fastSave)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(fastSave)
library(ggpubr)
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores =64)#
load("../All66HCC.phenotypes.RData")
system.time(AllCpG.pca <- prcomp(Meth.all.qc))
save.lbzip2(AllCpG.pca, file="pca.smooth.AllCpG.RDataFS", n.cores=64)
df <- HCC.phenotypes[match(rownames(Meth.all.qc), HCC.phenotypes$SampleID),c("SampleID","Condition", "age")]
df$age_group <- factor(df$age < 60, levels = c(TRUE, FALSE), 
                       labels = c("age<60", "age>=60"))
df$Condition <- factor(df$Condition, levels = c("T", "N"), labels = c("Tumor", "Non-Tumor"))
df$Group <- paste(df$Condition, "(", df$age_group, ")", sep="")
df$Group <- factor(df$Group, levels = c("Tumor(age>=60)", "Tumor(age<60)", "Non-Tumor(age>=60)", "Non-Tumor(age<60)"))
df$PC1 <- AllCpG.pca$x[,1]#58.1% Variance
df$PC2 <- AllCpG.pca$x[,2]#4.2% Variance

p_1.2 <- ggplot(data=df, aes(x = PC1, y = PC2, color = Condition))+
  geom_point(shape=19, size = 3)+
  scale_color_brewer(palette = "Set1")+
  guides(color = guide_legend(title = NULL))+
  xlab("PC1 (58.1%)")+ ylab("PC2 (4.2%)")+
  theme_pubr()+ theme(text = element_text(size=18))
p_1.2
dev.off()
save.lbzip2(p_1.2, file="./Figures/Figure_1.2.RDataFS", n.cores=24)
#1.3 The overall average methylation between paired tumor and adjacent non-tumor tissues (CpGs in different region version as supplementary)
setwd("/data/huangp/Methy/Mixed/New")#29
options(stringsAsFactors = F)
library(fastSave)
load("../All66HCC.phenotypes.RData")#load phenotypes
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores=64)#load smoothed methylation
#ggplot2
df <- HCC.phenotypes[,c("SampleID", "age", "Condition")]
df$MethyMean <- rowMeans(Meth.all.qc)
df$MethySD <- apply(Meth.all.qc, 1, sd)
df$MethyMedian <- apply(Meth.all.qc, 1, median)
df$MethyUpquarter <- apply(Meth.all.qc, 1, quantile, probs = 0.75)
df$MethyDownquarter <- apply(Meth.all.qc, 1, quantile, probs = 0.25)
df$Condition <- factor(df$Condition, levels = c("T", "N"), labels = c("Tumor", "Non-Tumor"))
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
pd <- position_dodge(0.3) # move them .05 to the left and right

p_1.3 <- ggplot(df, aes(x=age, y=MethyMedian, colour=Condition, group=Condition)) + 
  geom_errorbar(aes(ymin=MethyDownquarter, ymax=MethyUpquarter), width=.1,position=pd, lty = 1, alpha = 0.6) +
  geom_point(position=pd, fill = "white", shape = 21, size = 3.0)+
  geom_smooth(method="lm",se = F, lty = 2, show.legend = F)+
  scale_color_brewer(palette = "Set1")+
  xlab("Age (year)")+
  ylab("Average DNA methylation level")+
  guides(color = guide_legend(title=NULL))+
  theme_pubr()+ theme(text = element_text(size=20))
p_1.3
dev.off()
save.lbzip2(p_1.3, file="./Figures/Figure_1.3.RDataFS", n.cores=24)
#1.4 The physical distribution and enrichment of hypomethylated DML and hypermethylated DMLs
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")#29
library(fastSave)
load.lbzip2("AllCpG.annotation.RDataFS", n.cores=48)#load All 29M CpGs annotation
load.lbzip2("SigDMLs_from_combined_BWS_based_DML_detection_summary.RDataFS", n.cores=48)#load 9.8M SigDMLs
HyperDMLs.annotation <- AllCpG.annotation[rownames(SigDMLs)[SigDMLs$Median.delta > 0],]
HypoDMLs.annotation <- AllCpG.annotation[rownames(SigDMLs)[SigDMLs$Median.delta < 0],]
table(AllCpG.annotation$Promoter)#promoter
table(AllCpG.annotation$Intron)#Intron
table(AllCpG.annotation$Intergenic)#intergenic
table(AllCpG.annotation$Exon)#Exon
table(AllCpG.annotation$utr5)#5'utr
table(AllCpG.annotation$utr3)#3'utr
# Modify the overlap region accroding to Promoter > utr5 > utr3 > Exon > Intron > Intergenic
AllCpG.annotation$utr5[AllCpG.annotation$Promoter] <- FALSE
AllCpG.annotation$utr3[AllCpG.annotation$Promoter] <- FALSE
AllCpG.annotation$utr3[AllCpG.annotation$utr5] <- FALSE
table(AllCpG.annotation$Promoter)#promoter #4240443
table(AllCpG.annotation$Intron)#Intron #13381633
table(AllCpG.annotation$Intergenic)#intergenic #9824444
table(AllCpG.annotation$Exon)#Exon #974016
table(AllCpG.annotation$utr5)#5'utr #68123
table(AllCpG.annotation$utr3)#3'utr #465294 
table(HyperDMLs.annotation$Promoter)#promoter #62211
table(HyperDMLs.annotation$Intron)#Intron #55919
table(HyperDMLs.annotation$Intergenic)#intergenic #16285
table(HyperDMLs.annotation$Exon)#Exon #17023
table(HyperDMLs.annotation$utr5)#5'utr #9863
table(HyperDMLs.annotation$utr3)#3'utr #7061 
table(HypoDMLs.annotation$Promoter)#promoter #534272
table(HypoDMLs.annotation$Intron)#Intron #4217913
table(HypoDMLs.annotation$Intergenic)#intergenic #4655293
table(HypoDMLs.annotation$Exon)#Exon #225486
table(HypoDMLs.annotation$utr5)#5'utr #21583
table(HypoDMLs.annotation$utr3)#3'utr #79369 
df <- data.frame(Group = rep(c("All 28,978,826 CpGs", "157,320 Hyper-DMLs", "9,710,380 Hypo-DMLs"), each =6),
                 Region = rep(c("Promoter region", "Intron", "Intergenic region", "Exon", "5'UTR", "3'UTR"), 3),
                 Count = c(4240443,13381633,9824444,974016,68123,465294, 
                           62211,55919,16285,17023,9863,7061, 
                           534272,4217913,4655293,225486,21583,79369))
df$Group <- factor(df$Group, levels=c("All 28,978,826 CpGs", "157,320 Hyper-DMLs", "9,710,380 Hypo-DMLs"))
df$Ratio <- rep(0, nrow(df))
df$Ratio[1:6] <- df$Count[1:6]/sum(df$Count[1:6])
df$Ratio[7:12] <- df$Count[7:12]/sum(df$Count[7:12])
df$Ratio[13:18] <- df$Count[13:18]/sum(df$Count[13:18])
df$Ratio <- round(df$Ratio*100, digits = 2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
p_1.4 <- ggplot(df, aes(x = Region, y = Ratio, fill=Group))+
         geom_bar(stat = "identity", position = "dodge")+
         scale_fill_manual(values = brewer.pal(9,"Set1")[c(9,1,2)])+
         xlab(NULL)+ylab("Ratio(%)")+
         theme_pubr()+
         theme(axis.text.x = element_text(angle = 30, hjust = 1.0))+
         guides(fill = guide_legend( ncol = 1, byrow = TRUE, title = NULL))+
         theme(legend.position = c(0,1), legend.justification = c(0,1))+
         theme(text = element_text(size=20))

p_1.4
dev.off()
save(p_1.4, file="./Figures/Figure_1.4.RData")


################## part 2 Genes and pathways regulated by tumor asscoiated DMRs #############
#2.0 An example of repressed promoter in HCC using genome track via SeqMonk
#prepare input for SeqMonk
#input file list
#1.gencode.v31.gtf 2.Normal promoter bed files 3. Tumor promoter bed files 
#4.Raw DMR methylation level in adjacent bed file 5. Raw DMR methylation level in tumor bed file 
#6.Adjacent gene expression bed files #7.Tumor gene expression bed files.
setwd("/data/huangp/HCC/Methy/Mixed/New")#51
options(stringsAsFactors = F)
library(fastSave)
#load SigDMR
load.lbzip2("preDMRs.summary.RDataFS", n.cores=24)
DMR.Tumor.summary <- preDMR.summary
load("Granges.chroHMM.LiverAdult.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Adjacent_1.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Adjacent_2.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Tumor_1.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Tumor_2.Hepatology.hg38.RData")
load("Granges.HepG2.cRE.hg38.RData")
load("Granges.Hepatocyte.cRE.hg38.RData")
load("Granges.FemaleAdultLiver.cRE.hg38.RData")
#calculate the deltaScore of PromoterScore.T-PromoterScore.N
DMR.Tumor.summary$PromoterScoreDelta <- DMR.Tumor.summary$PromoterLike.T - DMR.Tumor.summary$PromoterLike.N
#pick out the top repressed promoter and activated promoter
CandidatePromoters <- DMR.Tumor.summary[abs(DMR.Tumor.summary$PromoterScoreDelta) >=2,]
DMR.Tumor.summary$EnhancerScoreDelta <- DMR.Tumor.summary$EnhancerLike.T - DMR.Tumor.summary$EnhancerLike.N
CandidateEnhancers <- DMR.Tumor.summary[abs(DMR.Tumor.summary$EnhancerScoreDelta) >=2,]
write.csv(CandidatePromoters, file="Candidate top activated or repressed Promoter-like DMR for plotting.csv", row.names=F)
write.csv(CandidateEnhancers, file="Candidate top activated or repressed Enhancer-like DMR for plotting.csv", row.names=F)
#build bed file for promoter like cRE in HepG2, H9, Liver, NormalLiver, Adjacent1-2, Tumor1-2
HepG2.promoter <- as.data.frame(Granges.HepG2.cRE)
HepG2.promoter$seqnames <- as.character(HepG2.promoter$seqnames)
HepG2.promoter <- HepG2.promoter[HepG2.promoter$annotRgb == "Promoter",]
#modify the HepG2.promoter to a bed format
HepG2.promoter <- HepG2.promoter[,-4]
HepG2.promoter <- HepG2.promoter[,c(1:3,5,6,4,7:9)]
HepG2.promoter$annotRgb <- rep("255,0,0",nrow(HepG2.promoter))
HepG2.promoter$thickStart <- HepG2.promoter$start
HepG2.promoter$thickEnd <- HepG2.promoter$start
HepG2.promoter$name <- rep("PromoterLike", nrow(HepG2.promoter))
write.table(HepG2.promoter, file="HepG2.promoter.bed", row.names=F, col.names=F, sep="\t", quote=F)
Hepatocyte.promoter <- as.data.frame(Granges.Hepatocyte.cRE)
Hepatocyte.promoter$seqnames <- as.character(Hepatocyte.promoter$seqnames)
Hepatocyte.promoter <- Hepatocyte.promoter[Hepatocyte.promoter$annotRgb == "Promoter",]
#modify the Hepatocyte.promoter to a bed format
Hepatocyte.promoter <- Hepatocyte.promoter[,-4]
Hepatocyte.promoter <- Hepatocyte.promoter[,c(1:3,5,6,4,7:9)]
Hepatocyte.promoter$annotRgb <- rep("255,0,0",nrow(Hepatocyte.promoter))
Hepatocyte.promoter$thickStart <- Hepatocyte.promoter$start
Hepatocyte.promoter$thickEnd <- Hepatocyte.promoter$start
Hepatocyte.promoter$name <- rep("PromoterLike", nrow(Hepatocyte.promoter))
write.table(Hepatocyte.promoter, file="Hepatocyte.promoter.bed", row.names=F, col.names=F, sep="\t", quote=F)
FemaleAdultLiver.promoter <- as.data.frame(Granges.AdultLiver.cRE)
FemaleAdultLiver.promoter$seqnames <- as.character(FemaleAdultLiver.promoter$seqnames)
FemaleAdultLiver.promoter <- FemaleAdultLiver.promoter[FemaleAdultLiver.promoter$annotRgb == "Promoter",]
#modify the FemaleAdultLiver.promoter to a bed format
FemaleAdultLiver.promoter <- FemaleAdultLiver.promoter[,-4]
FemaleAdultLiver.promoter <- FemaleAdultLiver.promoter[,c(1:3,5,6,4,7:9)]
FemaleAdultLiver.promoter$annotRgb <- rep("255,0,0",nrow(FemaleAdultLiver.promoter))
FemaleAdultLiver.promoter$thickStart <- FemaleAdultLiver.promoter$start
FemaleAdultLiver.promoter$thickEnd <- FemaleAdultLiver.promoter$start
FemaleAdultLiver.promoter$name <- rep("PromoterLike",nrow(FemaleAdultLiver.promoter))
write.table(FemaleAdultLiver.promoter, file="FemaleAdultLiver.promoter.bed", row.names=F, col.names=F, sep="\t", quote=F)
LiverAdult.promoter <- as.data.frame(Granges.Normal)
LiverAdult.promoter$seqnames <- as.character(LiverAdult.promoter$seqnames)
LiverAdult.promoter <- LiverAdult.promoter[,-4]
LiverAdult.promoter <- LiverAdult.promoter[LiverAdult.promoter$X %in% c("5_state5", "10_state10"),]
LiverAdult.promoter$names <- rep("", nrow(LiverAdult.promoter))
LiverAdult.promoter$annotRgb <- rep("", nrow(LiverAdult.promoter))
LiverAdult.promoter$annotRgb[LiverAdult.promoter$X == "5_state5"] <- "178,34,34"
LiverAdult.promoter$annotRgb[LiverAdult.promoter$X == "10_state10"] <- "255,0,0"
LiverAdult.promoter$thickStart <- LiverAdult.promoter$thickEnd <- LiverAdult.promoter$start
LiverAdult.promoter$score <- rep(0, nrow(LiverAdult.promoter))
LiverAdult.promoter <- LiverAdult.promoter[,c(1:3,6,10,4,8,9,7)]
LiverAdult.promoter$names[LiverAdult.promoter$annotRgb=="255,0,0"] <- "activePromoter"
LiverAdult.promoter$names[LiverAdult.promoter$annotRgb=="178,34,34"] <- "activeTSS"
write.table(LiverAdult.promoter, file="NormalLiver.promoter.bed", col.names=F, row.names=F, sep = "\t", quote=F)
Adjacent_1.promoter <- as.data.frame(Granges.Adjacent_1)
Adjacent_1.promoter$seqnames <- as.character(Adjacent_1.promoter$seqnames)
Adjacent_1.promoter <- Adjacent_1.promoter[,-4]
Adjacent_1.promoter <- Adjacent_1.promoter[Adjacent_1.promoter$X %in% c("5_state5", "10_state10"),]
Adjacent_1.promoter$names <- rep("", nrow(Adjacent_1.promoter))
Adjacent_1.promoter$annotRgb <- rep("", nrow(Adjacent_1.promoter))
Adjacent_1.promoter$annotRgb[Adjacent_1.promoter$X == "5_state5"] <- "178,34,34"
Adjacent_1.promoter$annotRgb[Adjacent_1.promoter$X == "10_state10"] <- "255,0,0"
Adjacent_1.promoter$thickStart <- Adjacent_1.promoter$thickEnd <- Adjacent_1.promoter$start
Adjacent_1.promoter$score <- rep(0, nrow(Adjacent_1.promoter))
Adjacent_1.promoter <- Adjacent_1.promoter[,c(1:3,6,10,4,8,9,7)]
Adjacent_1.promoter$names[Adjacent_1.promoter$annotRgb=="255,0,0"] <- "activePromoter"
Adjacent_1.promoter$names[Adjacent_1.promoter$annotRgb=="178,34,34"] <- "activeTSS"
write.table(Adjacent_1.promoter, file="Adjacent_1.promoter.bed", col.names=F, row.names=F, sep = "\t", quote=F)
Adjacent_2.promoter <- as.data.frame(Granges.Adjacent_2)
Adjacent_2.promoter$seqnames <- as.character(Adjacent_2.promoter$seqnames)
Adjacent_2.promoter <- Adjacent_2.promoter[,-4]
Adjacent_2.promoter <- Adjacent_2.promoter[Adjacent_2.promoter$X %in% c("5_state5", "10_state10"),]
Adjacent_2.promoter$names <- rep("", nrow(Adjacent_2.promoter))
Adjacent_2.promoter$annotRgb <- rep("", nrow(Adjacent_2.promoter))
Adjacent_2.promoter$annotRgb[Adjacent_2.promoter$X == "5_state5"] <- "178,34,34"
Adjacent_2.promoter$annotRgb[Adjacent_2.promoter$X == "10_state10"] <- "255,0,0"
Adjacent_2.promoter$thickStart <- Adjacent_2.promoter$thickEnd <- Adjacent_2.promoter$start
Adjacent_2.promoter$score <- rep(0, nrow(Adjacent_2.promoter))
Adjacent_2.promoter <- Adjacent_2.promoter[,c(1:3,6,10,4,8,9,7)]
Adjacent_2.promoter$names[Adjacent_2.promoter$annotRgb=="255,0,0"] <- "activePromoter"
Adjacent_2.promoter$names[Adjacent_2.promoter$annotRgb=="178,34,34"] <- "activeTSS"
write.table(Adjacent_2.promoter, file="Adjacent_2.promoter.bed", col.names=F, row.names=F, sep = "\t", quote=F)
Tumor_1.promoter <- as.data.frame(Granges.Tumor_1)
Tumor_1.promoter$seqnames <- as.character(Tumor_1.promoter$seqnames)
Tumor_1.promoter <- Tumor_1.promoter[,-4]
Tumor_1.promoter <- Tumor_1.promoter[Tumor_1.promoter$X %in% c("5_state5", "10_state10"),]
Tumor_1.promoter$names <- rep("", nrow(Tumor_1.promoter))
Tumor_1.promoter$annotRgb <- rep("", nrow(Tumor_1.promoter))
Tumor_1.promoter$annotRgb[Tumor_1.promoter$X == "5_state5"] <- "178,34,34"
Tumor_1.promoter$annotRgb[Tumor_1.promoter$X == "10_state10"] <- "255,0,0"
Tumor_1.promoter$thickStart <- Tumor_1.promoter$thickEnd <- Tumor_1.promoter$start
Tumor_1.promoter$score <- rep(0, nrow(Tumor_1.promoter))
Tumor_1.promoter <- Tumor_1.promoter[,c(1:3,6,10,4,8,9,7)]
Tumor_1.promoter$names[Tumor_1.promoter$annotRgb=="255,0,0"] <- "activePromoter"
Tumor_1.promoter$names[Tumor_1.promoter$annotRgb=="178,34,34"] <- "activeTSS"
write.table(Tumor_1.promoter, file="Tumor_1.promoter.bed", col.names=F, row.names=F, sep = "\t", quote=F)
Tumor_2.promoter <- as.data.frame(Granges.Tumor_2)
Tumor_2.promoter$seqnames <- as.character(Tumor_2.promoter$seqnames)
Tumor_2.promoter <- Tumor_2.promoter[,-4]
Tumor_2.promoter <- Tumor_2.promoter[Tumor_2.promoter$X %in% c("5_state5", "10_state10"),]
Tumor_2.promoter$names <- rep("", nrow(Tumor_2.promoter))
Tumor_2.promoter$annotRgb <- rep("", nrow(Tumor_2.promoter))
Tumor_2.promoter$annotRgb[Tumor_2.promoter$X == "5_state5"] <- "178,34,34"
Tumor_2.promoter$annotRgb[Tumor_2.promoter$X == "10_state10"] <- "255,0,0"
Tumor_2.promoter$thickStart <- Tumor_2.promoter$thickEnd <- Tumor_2.promoter$start
Tumor_2.promoter$score <- rep(0, nrow(Tumor_2.promoter))
Tumor_2.promoter <- Tumor_2.promoter[,c(1:3,6,10,4,8,9,7)]
Tumor_2.promoter$names[Tumor_2.promoter$annotRgb=="255,0,0"] <- "activePromoter"
Tumor_2.promoter$names[Tumor_2.promoter$annotRgb=="178,34,34"] <- "activeTSS"
write.table(Tumor_2.promoter, file="Tumor_2.promoter.bed", col.names=F, row.names=F, sep = "\t", quote=F)

options(stringsAsFactors = F)
library(fastSave)
library(GenomicFeatures)
library(genomation)

setwd("/data/huangp/HCC/Methy/Mixed/New")
load.lbzip2("47.5k_DMR.promoter.methy.matrix.RDataFS", n.cores=24)
PromoterDMR.Tumor <- colMeans(DMR.promoter.methy.matrix[1:33,])
PromoterDMR.Adjacent <- colMeans(DMR.promoter.methy.matrix[34:66,])
#build bedGraph file
#track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20
PromoterDMR.Tumor.bedGraph <- data.frame(chr = sapply(strsplit(names(PromoterDMR.Tumor), split="_"), "[", 1),
                                         start = as.integer(sapply(strsplit(names(PromoterDMR.Tumor), split="_"), "[", 2)),
                                         end = as.integer(sapply(strsplit(names(PromoterDMR.Tumor), split="_"), "[", 3)),
                                         methyMean = PromoterDMR.Tumor 
                                         )
PromoterDMR.Adjacent.bedGraph <- data.frame(chr = sapply(strsplit(names(PromoterDMR.Adjacent), split="_"), "[", 1),
                                         start = as.integer(sapply(strsplit(names(PromoterDMR.Adjacent), split="_"), "[", 2)),
                                         end = as.integer(sapply(strsplit(names(PromoterDMR.Adjacent), split="_"), "[", 3)),
                                         methyMean = PromoterDMR.Adjacent 
                                            )
write.table(PromoterDMR.Tumor.bedGraph, file="PromoterDMR.Tumor.bedGraph", sep="\t", col.names=F, row.names=F, quote=F)
write.table(PromoterDMR.Adjacent.bedGraph, file="PromoterDMR.Adjacent.bedGraph", sep="\t", col.names=F, row.names=F, quote=F)

#load gene expr
load.lbzip2("Gene.TPM.filtered.RDataFS", n.cores=24)
Granges.gtf <- gffToGRanges("gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
Granges.gene <- Granges.gene[Granges.gene$gene_id %in% rownames(gene.TPM),]
range_gene_DF <- as.data.frame(Granges.gene)


#2.1 The distribution of hyper DMR and hypo DMR in genomic region.
setwd("/data/huangp/Methy/Mixed/New")#29
options(stringsAsFactors = F)
library(fastSave)
#load DMR summar
load.lbzip2("preDMRs.summary.RDataFS", n.cores=32)
HyperDMR <- preDMR.summary[preDMR.summary$Median.delta > 0,]
# Modify the overlap region accroding to Promoter > utr5 > utr3 > Exon > Intron > Intergenic
HyperDMR$utr5[HyperDMR$Promoter] <- FALSE
HyperDMR$utr3[HyperDMR$Promoter] <- FALSE
HyperDMR$utr3[HyperDMR$utr5] <- FALSE
attach(HyperDMR)
table(Promoter)#2315
table(Intron)#3024
table(Intergenic)#641
table(Exon)#512
table(utr5)#80
table(utr3)#352
detach(HyperDMR)
HypoDMR <- preDMR.summary[preDMR.summary$Median.delta < 0,]
# Modify the overlap region accroding to Promoter > utr5 > utr3 > Exon > Intron > Intergenic
HypoDMR$utr5[HypoDMR$Promoter] <- FALSE
HypoDMR$utr3[HypoDMR$Promoter] <- FALSE
HypoDMR$utr3[HypoDMR$utr5] <- FALSE
attach(HypoDMR)
table(Promoter)#39,140
table(Intergenic)#282,990
table(Intron)#250,932
table(Exon)#20,521
table(utr5)#1325
table(utr3)#6445
detach(HypoDMR)
df <- data.frame(Group = rep(c("Hyper-DMR", "Hypo-DMR"), each = 6),
                 Region = rep(c("Promoter region", "Intron", "Intergenic region", "Exon", "5'UTR", "3'UTR"),2),
                 Count = c(2315,3024,641,512,80,352, 39140,282990,250932,20521,1325,6445))
df$Ratio <- rep(0, nrow(df))
df$Ratio[1:6] <- df$Count[1:6]/sum(df$Count[1:6])*100
df$Ratio[7:12] <- df$Count[7:12]/sum(df$Count[7:12])*100
df$Ratio <- round(df$Ratio,2)
df$Group <- factor(df$Group, levels = c("Hyper-DMR", "Hypo-DMR"), labels = c("6,924 Hyper-DMRs", "601,353 Hypo-DMRs"))
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
p_2.1 <- ggplot(df, aes(x = Region, y = Ratio, fill=Group))+
         geom_bar(stat = "identity", position = "dodge")+
         scale_fill_manual(values = brewer.pal(9,"Set1")[c(1,2)])+
         xlab(NULL)+ylab("Ratio (%)")+
         theme_pubr()+guides(fill= guide_legend(title = NULL, ncol = 1, byrow = T))+
        theme(legend.position = c(0,1), legend.justification = c(0,1))+
        theme(text = element_text(size=20))+theme(axis.text.x = element_text(angle = 30, hjust = 1.0))
p_2.1
dev.off()
save(p_2.1, file="./Figures/Figure_2.1.RData")
#2.2 The distribution of liver active promoter score.T and score.N of promoter-like DMRs
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")#29
load.lbzip2(file="600k.preDMRs.summary.RDataFS",n.cores=24)
if(1<0){
preDMR.summary$Promoter.Normal <- (preDMR.summary$activeTSS.Normal | preDMR.summary$activePromoter.Normal)
preDMR.summary$Promoter.Adjacent_1 <- (preDMR.summary$activeTSS.Adjacent_1 | preDMR.summary$activePromoter.Adjacent_1)
preDMR.summary$Promoter.Adjacent_2 <- (preDMR.summary$activeTSS.Adjacent_2 | preDMR.summary$activePromoter.Adjacent_2)
preDMR.summary$Promoter.Tumor_1 <- (preDMR.summary$activeTSS.Tumor_1 | preDMR.summary$activePromoter.Tumor_1)
preDMR.summary$Promoter.Tumor_2 <- (preDMR.summary$activeTSS.Tumor_2 | preDMR.summary$activePromoter.Tumor_2)
preDMR.summary$PromoterLike.T <- rowSums(preDMR.summary[,c("Promoter.cRE.HepG2", "Promoter.Tumor_1", "Promoter.Tumor_2")])
preDMR.summary$PromoterLike.N <- rowSums(preDMR.summary[,c("Promoter.cRE.AdultLiver", "Promoter.cRE.Hepatocyte", "Promoter.Normal", "Promoter.Adjacent_1", "Promoter.Adjacent_2")])
preDMR.summary$PromoterLike.score <- preDMR.summary$PromoterLike.T + preDMR.summary$PromoterLike.N
preDMR.summary$PromoterLike <- preDMR.summary$PromoterLike.score > 0     
}
PromoterDMRs <- preDMR.summary[preDMR.summary$Promoter | preDMR.summary$PromoterLike.score > 0, 
                               c("Median.delta", "Promoter", "PromoterLike.score", "PromoterLike.T", "PromoterLike.N")]
DF <- cbind(PromoterDMRs, DM = PromoterDMRs$Median.delta > 0, Active = PromoterDMRs$PromoterLike.score > 0)
DF$DM <- factor(DF$DM, levels = c(TRUE, FALSE), labels = c("Hyper-DMR", "Hypo-DMR"))
DF$Active <- factor(DF$Active, levels = c(TRUE, FALSE), labels = c("Active promoter", "Inactive promoter"))
DF$PromoterLike <- paste(DF$PromoterLike.N, DF$PromoterLike.T, sep="|") 


data <- as.matrix(table(DF$PromoterLike, DF$DM))
data <- data.frame(Score.N = as.integer(sapply(strsplit(rownames(data), split="|", fixed=T), "[", 1)),
                   Score.T = as.integer(sapply(strsplit(rownames(data), split="|", fixed=T), "[", 2)),
                   HyperDMR = data[,1], HypoDMR = data[,2])
data$DMR <- data$HyperDMR+data$HypoDMR
data$Type <- (data$Score.T/3)>(data$Score.N/5)
data$Type[data$Type==T] <- "Activated in Tumor"
data$Type[data$Type==F] <- "Repressed in Tumor"
data$Type[data$Score.N == 0 & data$Score.T == 0] <- "Non-detected"
df <- data[data$Type != "Non-detected",]
data$Type <- factor(data$Type, levels=c("Activated in Tumor", "Repressed in Tumor", "Non-detected"))
df$Type <- factor(df$Type, levels = c("Activated in Tumor", "Repressed in Tumor"))
library(patchwork)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
p_2.2.1 <- ggplot(data, aes(Score.N, Score.T, label = HyperDMR))+
  geom_point(aes(size=HyperDMR), fill = "Firebrick1", color = "Firebrick1")+
  #scale_shape_manual(values = c(24,25,21))+
  scale_size_area(max_size=25,guide=F)+
  geom_text(vjust = -1)+
  xlab("Promoter activity score in five nontumor samples")+
  ylab("Promoter activity score in three tumor samples")+ 
  ggtitle("Promoter-like Hyper-DMRs (N=2,882)")+
  geom_abline(slope = 3/5, intercept = 0, lty = 3, color ="black")+
  ylim(c(-0.2,3.2))+
  xlim(c(-0.2, 5.2))+
  theme_pubr()+theme(text = element_text(size=18))
postscript("Fig 2.2.1.eps", width = 8, height = 8)
p_2.2.1
dev.off()


p_2.2.2 <- ggplot(data, aes(Score.N, Score.T, label = HypoDMR))+
  geom_point(aes(size=HypoDMR), fill = "Firebrick1", color = "Firebrick1")+
  #scale_shape_manual(values = c(24,25,21))+
  scale_size_area(max_size=25,guide=F)+
  geom_text(vjust = -1)+
  xlab("Promoter activity score in five nontumor samples")+
  ylab("Promoter activity score in three tumor samples")+ 
  ggtitle("Promoter-like Hypo-DMRs (N=44,611)")+
  geom_abline(slope = 3/5, intercept = 0, lty = 3, color ="black")+
  ylim(c(-0.2,3.2))+
  xlim(c(-0.2, 5.2))+
  theme_pubr()+theme(text = element_text(size=18))
postscript("Fig 2.2.2.eps", width = 8, height = 8)
p_2.2.2
dev.off()

df3 <- data.frame(group = rep(c("Hyper-DMR", "Hypo-DMR"), each=2),
                  type = rep(c("Activated in Tumor", "Repressed in Tumor"), 2),
                  count = c(sum(df$HyperDMR[df$Type == "Activated in Tumor"]), sum(df$HyperDMR[df$Type == "Repressed in Tumor"]),
                            sum(df$HypoDMR[df$Type == "Activated in Tumor"]), sum(df$HypoDMR[df$Type == "Repressed in Tumor"]))
)
p_2.2.3 <- ggplot(data=df3, aes(x=group, y=count, fill = type))+ geom_bar(stat = "identity", position = "fill", width = 0.8, size=0.25)+
  scale_fill_manual(values = brewer.pal(8,"Paired")[c(2,1)])+ xlab(NULL)+ylab("Ratio")+
  coord_flip()+
  theme_classic()+
  guides(fill = guide_legend(title = NULL, nrow = 1))+
  theme(legend.position = "top")+theme(text = element_text(size=18))
postscript("Fig 2.2.3.eps", width = 8, height = 8/3)
p_2.2.3
dev.off()


#2.3 The distribution of liver active enhancer score.T and score.N of genic enhancer-like DMRs
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")#29
load.lbzip2(file="600k.preDMRs.summary.RDataFS",n.cores=24)
EnhancerDMRs <- preDMR.summary[preDMR.summary$EnhancerLike.score > 0, 
                               c("Median.delta", "Intergenic", "EnhancerLike.score", "EnhancerLike.T", "EnhancerLike.N")]
DF <- cbind(EnhancerDMRs, DM = EnhancerDMRs$Median.delta > 0, Active = EnhancerDMRs$EnhancerLike.score > 0)
DF$DM <- factor(DF$DM, levels = c(TRUE, FALSE), labels = c("Hyper-DMR", "Hypo-DMR"))
DF$Active <- factor(DF$Active, levels = c(TRUE, FALSE), labels = c("Active enhancer", "Inactive enhancer"))
DF$EnhancerLike <- paste(DF$EnhancerLike.N, DF$EnhancerLike.T, sep="|") 

data <- as.matrix(table(DF$EnhancerLike, DF$DM))
data <- data.frame(Score.N = as.integer(sapply(strsplit(rownames(data), split="|", fixed=T), "[", 1)),
                   Score.T = as.integer(sapply(strsplit(rownames(data), split="|", fixed=T), "[", 2)),
                   HyperDMR = data[,1], HypoDMR = data[,2])
data$DMR <- data$HyperDMR+data$HypoDMR
data$Type <- (data$Score.T/3)>(data$Score.N/5)
data$Type[data$Type==T] <- "Activated in Tumor"
data$Type[data$Type==F] <- "Repressed in Tumor"
data$Type <- factor(data$Type, levels = c("Activated in Tumor", "Repressed in Tumor"))
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
p_2.3.1 <- ggplot(data, aes(Score.N, Score.T, label = HyperDMR))+
  geom_point(aes(size=HyperDMR), fill="gold", color="gold")+
  scale_size_area(max_size=25,guide=F)+
  geom_text(vjust = -1.5)+
  xlab("Enhancer activity score in five nontumor samples")+
  ylab("Enhancer activity score in three tumor samples")+ 
  ggtitle("Enhancer-like Hyper-DMRs (N=3,232)")+
  geom_abline(slope = 3/5, intercept = 0, lty = 3, color ="black")+
  ylim(c(-0.2,3.2))+
  xlim(c(-0.2, 5.2))+
  theme_pubr()+
  theme(text = element_text(size=18))
p_2.3.1
dev.off()
postscript("Fig 2.3.1.eps", width = 8, height = 8)
p_2.3.1
dev.off()

p_2.3.2 <- ggplot(data, aes(Score.N, Score.T, label = HypoDMR))+
  geom_point(aes(size=HypoDMR), fill="gold", color="gold")+
  scale_size_area(max_size=25,guide=F)+
  geom_text(vjust = -1.5)+
  xlab("Enhancer activity score in five nontumor samples")+
  ylab("Enhancer activity score in three tumor samples")+ 
  ggtitle("Enhancer-like Hypo-DMRs (N=20,568)")+
  geom_abline(slope = 3/5, intercept = 0, lty = 3, color ="black")+
  ylim(c(-0.2,3.2))+
  xlim(c(-0.2, 5.2))+
  theme_pubr()+
  theme(text = element_text(size=18))
p_2.3.2
dev.off()
postscript("Fig 2.3.2.eps", width = 8, height = 8)
p_2.3.2
dev.off()

df3 <- data.frame(group = rep(c("Hyper-DMR", "Hypo-DMR"), each=2),
                  type = rep(c("Activated in Tumor", "Repressed in Tumor"), 2),
                  count = c(sum(data$HyperDMR[data$Type == "Activated in Tumor"]), sum(data$HyperDMR[data$Type == "Repressed in Tumor"]),
                            sum(data$HypoDMR[data$Type == "Activated in Tumor"]), sum(data$HypoDMR[data$Type == "Repressed in Tumor"]))
)
p_2.3.3 <- ggplot(data=df3, aes(x=group, y=count, fill = type))+ geom_bar(stat = "identity", position = "fill", width = 0.8, size=0.25)+
  scale_fill_manual(values = brewer.pal(8,"Paired")[c(2,1)])+ xlab(NULL)+ylab("Ratio")+
  coord_flip()+guides(fill = guide_legend(title = NULL))+
  theme_classic()+theme(legend.position = "top")+theme(text = element_text(size=18))
postscript("Fig 2.3.3.eps", width = 8, height = 8/3)
p_2.3.3
dev.off()

ggsave("Rplots.pdf", width = 15, height = 5, units = "cm")

p_2.3 <- p_2.3.1 / p_2.3.2/p_2.3.3+plot_layout(guides = "collect", heights = c(5,5,1))
p_2.3
save(p_2.3, file="./Figures/Figure_2.3.RData")
#2.4 The distribution of liver active enhancer score of intergenic enhancer-like DMRs
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")#29
load.lbzip2("IntergenicDMR.eRNAandMethy.correlation.RDataFS", n.cores=24)#load eRNAandMethy correlation
#table(Sig = IntergenicDMR.eRNA.methy.correlation$eRNAandMethy.padj < 0.05, 
#       Classic = IntergenicDMR.eRNA.methy.correlation$eRNAandMethy.rho < 0)
#         Classic
#Sig     FALSE  TRUE
#FALSE  9215 20711
#TRUE   2178  4856
#scatter plot
df <- IntergenicDMR.eRNA.methy.correlation[,c("Median.delta", "eRNAandMethy.rho", "eRNAandMethy.padj")]
df$DM <- factor(df$Median.delta > 0, levels = c(TRUE, FALSE), labels = c("Hyper-DMRs", "Hypo-DMRs"))
df$Regulation[df$eRNAandMethy.padj > 0.05] <- "N.S." 
df$Regulation[df$eRNAandMethy.padj <= 0.05 & df$eRNAandMethy.rho > 0] <- "Positive"
df$Regulation[df$eRNAandMethy.padj <= 0.05 & df$eRNAandMethy.rho < 0] <- "Negative"
df$Regulation <- factor(df$Regulation, levels=c("Positive", "Negative", "N.S."))
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
#barplot 
DF <- data.frame(DM = rep(c("Hyper-DMR","Hypo-DMR"), each =3), Regulation = rep(c("Positive", "Negative", "N.S."), 2),
                 Count = c(52, 23, 234, 2126, 4823, 29692))
DF$Regulation <- factor(DF$Regulation, levels=c("N.S.", "Positive", "Negative"))
DF$DM <- factor(DF$DM, levels = c("Hyper-DMR", "Hypo-DMR"), labels = c("Hyper-DMRs(N=309)", "Hypo-DMRs(N=36,651)"))
p_2.4 <- ggplot(DF, aes(x = DM, y = Count, fill = Regulation))+
  geom_bar(stat = "identity", position = "fill", color = "black", width=0.5, size =0.25)+
  scale_fill_manual(values = brewer.pal(n=9, "Set1")[c(9,1,2)])+
  coord_flip()+xlab(NULL)+ylab("Ratio of intergenic DMRs with active eRNA expression")+
  theme_pubr()
p_2.4
save(p_2.4, file="./Figures/Figure_2.4.RData")
#2.5 The distribution of four types of promoter-like/genic/intergenic enhancer-like DMR associated DEGs
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")
if(1<0)
{
  #load HCC66 DEG result
  promoter_genes <- read.csv("DMR.promoterLike_Gene_pairs.csv", header=T)
  load("/data/huangp/HCC_Enhancer/HCC66_DEG_res.RData")
  promoter_genes <- data.frame(promoter_genes[,1:10], DEG.LFC = HCC66_DEG_res$log2FoldChange[match(promoter_genes$gene_name, HCC66_DEG_res$geneName)],
                               DEG.padj = HCC66_DEG_res$padj[match(promoter_genes$gene_name, HCC66_DEG_res$geneName)],
                               promoter_genes[,11:21])
  library(fastSave)
  load.lbzip2(file="600k.preDMRs.summary.RDataFS",n.cores=24)
  promoter_genes[,c("Promoter", "PromoterLike.T", "PromoterLike.N", "PromoterLike.score")] <- preDMR.summary[match(promoter_genes$DMR.ID, rownames(preDMR.summary)),c("PromoterLike", "PromoterLike.T", "PromoterLike.N", "PromoterLike.score")]
  promoter_genes <- promoter_genes[promoter_genes$DEG.padj < 0.05 & !is.na(promoter_genes$DEG.padj),]#filter out non-significant DEGs
  load.lbzip2("47.5k_DMR.promoter.methy.matrix.RDataFS", n.cores=32)
  DMR.promoter.delta.methy <- colMeans(DMR.promoter.methy.matrix[grep(rownames(DMR.promoter.methy.matrix), pattern="T"), ]) - colMeans(DMR.promoter.methy.matrix[grep(rownames(DMR.promoter.methy.matrix), pattern="N"), ])
  promoter_genes$Median.delta <- DMR.promoter.delta.methy[match(promoter_genes$DMR.ID, names(DMR.promoter.delta.methy))]
  colnames(promoter_genes)[2] <- "Delta.Methy"
  promoter_genes <- promoter_genes[abs(promoter_genes$Delta.Methy) >= 0.15,]
  promoter_genes <- promoter_genes[abs(promoter_genes$DEG.LFC) >= 0.5,]
  write.csv(promoter_genes, file="FilteredDMR.promoterLike_Gene_pairs.csv", row.names=F, quote=F)
}
promoter_genes <- read.csv("FilteredDMR.promoterLike_Gene_pairs.csv", header=T)#3046 pairs
table(Hyper = promoter_genes$Median.delta >0, Classic = promoter_genes$rho < 0)
#         Classic
#Hyper   FALSE TRUE
#FALSE   659 2138
#TRUE    145  104


if(1<0)
{
  #load HCC66 DEG result
  load("/data/huangp/HCC_Enhancer/HCC66_DEG_res.RData")
  genic_enhancer_genes <- data.frame(genic_enhancer_genes[,1:7], DEG.LFC = HCC66_DEG_res$log2FoldChange[match(genic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                     DEG.padj = HCC66_DEG_res$padj[match(genic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                     genic_enhancer_genes[,8:19])
  genic_enhancer_genes <- genic_enhancer_genes[genic_enhancer_genes$DEG.padj <= 0.05 & !is.na(genic_enhancer_genes$DEG.padj),]
  load.lbzip2("23.8k_DMR.enhancer.methy.matrix.RDataFS", n.cores=32)
  DMR.enhancer.delta.methy <- colMeans(DMR.enhancer.methy.matrix[grep(rownames(DMR.enhancer.methy.matrix), pattern="T"), ]) - colMeans(DMR.enhancer.methy.matrix[grep(rownames(DMR.enhancer.methy.matrix), pattern="N"), ])
  genic_enhancer_genes$Median.delta <- DMR.enhancer.delta.methy[match(genic_enhancer_genes$DMR.ID, names(DMR.enhancer.delta.methy))]
  colnames(genic_enhancer_genes)[1] <- "Delta.Methy"
  genic_enhancer_genes <- genic_enhancer_genes[abs(genic_enhancer_genes$Delta.Methy)>=0.15,]
  genic_enhancer_genes <- genic_enhancer_genes[abs(genic_enhancer_genes$DEG.LFC) >= 0.5,]
  write.csv(genic_enhancer_genes, file="FilteredGenicDMR.enhancerLike_Gene_pairs.csv", row.names=F, quote=F)
  
}
genic_enhancer_genes <- read.csv("FilteredGenicDMR.enhancerLike_Gene_pairs.csv", header=T)
if(1<0)
{
  #load HCC66 DEG result
  load("/data/huangp/HCC_Enhancer/HCC66_DEG_res.RData")
  intergenic_enhancer_genes <- data.frame(intergenic_enhancer_genes[,1:13], DEG.LFC = HCC66_DEG_res$log2FoldChange[match(intergenic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                     DEG.padj = HCC66_DEG_res$padj[match(intergenic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                     intergenic_enhancer_genes[,14:27])
  intergenic_enhancer_genes <- intergenic_enhancer_genes[intergenic_enhancer_genes$DEG.padj <= 0.05 & !is.na(intergenic_enhancer_genes$DEG.padj),]
  load.lbzip2("IntergenicDMR.eRNA.methy.matrix.RDataFS", n.cores=32)
  IntergenicDMR.delta.methy <- colMeans(IntergenicDMR.eRNA.methy.matrix[grep(rownames(IntergenicDMR.eRNA.methy.matrix), pattern="T"), ]) - colMeans(IntergenicDMR.eRNA.methy.matrix[grep(rownames(IntergenicDMR.eRNA.methy.matrix), pattern="N"), ])
  intergenic_enhancer_genes$DM.Median <- IntergenicDMR.delta.methy[match(intergenic_enhancer_genes$DMR.ID, names(IntergenicDMR.delta.methy))]
  colnames(intergenic_enhancer_genes)[1] <- "Delta.Methy"
  intergenic_enhancer_genes <- intergenic_enhancer_genes[abs(intergenic_enhancer_genes$Delta.Methy)>=0.15,]
  intergenic_enhancer_genes <- intergenic_enhancer_genes[abs(intergenic_enhancer_genes$DEG.LFC) >= 0.5,]
  write.csv(intergenic_enhancer_genes, file="FilteredConsistent_eRNARegulatedDMRs_Gene_pairs.csv", row.names=F, quote=F)
}
intergenic_enhancer_genes <- read.csv("FilteredConsistent_eRNARegulatedDMRs_Gene_pairs.csv", header=T)

if(1<0)
{
  hc_promoter_genes <- read.csv("HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", header = T)
  load("/data/huangp/HCC_Enhancer/HCC66_DEG_res.RData")
  hc_promoter_genes <- data.frame(promoter_genes[,1:10], DEG.LFC = HCC66_DEG_res$log2FoldChange[match(promoter_genes$gene_name, HCC66_DEG_res$geneName)],
                               DEG.padj = HCC66_DEG_res$padj[match(promoter_genes$gene_name, HCC66_DEG_res$geneName)],
                               promoter_genes[,11:25])
  library(fastSave)
  load.lbzip2(file="600k.preDMRs.summary.RDataFS",n.cores=24)
  hc_promoter_genes[,c("Promoter", "PromoterLike.T", "PromoterLike.N", "PromoterLike.score")] <- preDMR.summary[match(hc_promoter_genes$DMR.ID, rownames(preDMR.summary)),c("PromoterLike", "PromoterLike.T", "PromoterLike.N", "PromoterLike.score")]
  hc_promoter_genes <- promoter_genes[hc_promoter_genes$DEG.padj < 0.05 & !is.na(hc_promoter_genes$DEG.padj),]#filter out non-significant DEGs
  hc_promoter_genes <- hc_promoter_genes[hc_promoter_genes$DMR.ID %in% promoter_genes$DMR.ID,]
  hc_promoter_genes$Median.delta <- promoter_genes$Delta.Methy[match(hc_promoter_genes$DMR.ID, promoter_genes$DMR.ID)]
  colnames(hc_promoter_genes)[2] <- "Delta.Methy"
  hc_promoter_genes <- hc_promoter_genes[abs(hc_promoter_genes$DEG.LFC)>=0.5,]
  write.csv(hc_promoter_genes, file="HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", row.names=F, quote=F)
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC < 0]))#hyper_down 12
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC < 0 & hc_promoter_genes$LIHC.validation == "TypeI"]))#4
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC < 0 & hc_promoter_genes$LIHC.validation == "TypeII"]))#1
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC < 0 & hc_promoter_genes$LIHC.validation == "Validation"]))#7
  
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC > 0]))#hypo_up 108
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC > 0 & hc_promoter_genes$LIHC.validation == "TypeI"]))#61
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC > 0 & hc_promoter_genes$LIHC.validation == "TypeII"]))#15
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC > 0 & hc_promoter_genes$LIHC.validation == "Validation"]))#32
  
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC > 0]))#hyper_up 21
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC > 0 & hc_promoter_genes$LIHC.validation == "TypeI"]))#1
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC > 0 & hc_promoter_genes$LIHC.validation == "TypeII"]))#10
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy > 0 & hc_promoter_genes$DEG.LFC > 0 & hc_promoter_genes$LIHC.validation == "Validation"]))#10
  
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC < 0]))#hyper_up 30
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC < 0 & hc_promoter_genes$LIHC.validation == "TypeI"]))
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC < 0 & hc_promoter_genes$LIHC.validation == "TypeII"]))
  length(unique(hc_promoter_genes$gene_name[hc_promoter_genes$Delta.Methy < 0 & hc_promoter_genes$DEG.LFC < 0 & hc_promoter_genes$LIHC.validation == "Validation"]))
  
  
  hc_genic_enhancer_genes <- read.csv("HighConfidence_genic_Enhancer_Like_DMR_Gene_pairs_0.8.csv", header=T)
  hc_genic_enhancer_genes <- data.frame(hc_genic_enhancer_genes[,1:2], DEG.LFC = HCC66_DEG_res$log2FoldChange[match(hc_genic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                        DEG.padj = HCC66_DEG_res$padj[match(hc_genic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                        hc_genic_enhancer_genes[,3:8])
  hc_genic_enhancer_genes <- hc_genic_enhancer_genes[hc_genic_enhancer_genes$DMR.ID %in% genic_enhancer_genes$DMR.ID,]
  hc_genic_enhancer_genes <- data.frame(Delta.Methy = genic_enhancer_genes$Delta.Methy[match(hc_genic_enhancer_genes$DMR.ID, genic_enhancer_genes$DMR.ID)],
                                        hc_genic_enhancer_genes)
  hc_genic_enhancer_genes <- hc_genic_enhancer_genes[!is.na(hc_genic_enhancer_genes$DEG.padj) & hc_genic_enhancer_genes$DEG.padj<= 0.05 & abs(hc_genic_enhancer_genes$DEG.LFC) >= 0.5,]
  write.csv(hc_genic_enhancer_genes, file="HighConfidence_genic_Enhancer_Like_DMR_Gene_pairs_0.8.csv", row.names=F, quote=F)
  length(unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC < 0]))#hyper_down 14
  X <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC < 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeI"])#4
  Y <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC < 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeII"])#1
  Z <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC < 0 & hc_genic_enhancer_genes$LIHC.Replication == "Replicated"])#7
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)  

  length(unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC > 0]))#hypo_up 108
  X <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC > 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeI"])
  Y <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC > 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeII"])
  Z <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC > 0 & hc_genic_enhancer_genes$LIHC.Replication == "Replicated"])
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)
  
  length(unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC > 0]))#hyper_up 21
  X <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC > 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeI"])#1
  Y <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC > 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeII"])#10
  Z <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy > 0 & hc_genic_enhancer_genes$DEG.LFC > 0 & hc_genic_enhancer_genes$LIHC.Replication == "Replicated"])#10
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)
  
  length(unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC < 0]))#hyper_up 30
  X <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC < 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeI"])
  Y <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC < 0 & hc_genic_enhancer_genes$LIHC.Replication == "TypeII"])
  Z <- unique(hc_genic_enhancer_genes$gene_name[hc_genic_enhancer_genes$Delta.Methy < 0 & hc_genic_enhancer_genes$DEG.LFC < 0 & hc_genic_enhancer_genes$LIHC.Replication == "Replicated"])
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)
  
  hc_intergenic_enhancer_genes <- read.csv("HighConfidence_Intergenic_Enhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)
  hc_intergenic_enhancer_genes <- data.frame(hc_intergenic_enhancer_genes[,1:13], DEG.LFC = HCC66_DEG_res$log2FoldChange[match(hc_intergenic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                             DEG.padj = HCC66_DEG_res$padj[match(hc_intergenic_enhancer_genes$gene_name, HCC66_DEG_res$geneName)],
                                             hc_intergenic_enhancer_genes[,14:30])
  hc_intergenic_enhancer_genes <- hc_intergenic_enhancer_genes[hc_intergenic_enhancer_genes$DMR.ID %in% intergenic_enhancer_genes$DMR.ID,]
  hc_intergenic_enhancer_genes$DM.Median <- intergenic_enhancer_genes$Delta.Methy[match(hc_intergenic_enhancer_genes$DMR.ID, intergenic_enhancer_genes$DMR.ID)]
  colnames(hc_intergenic_enhancer_genes)[1] <- "Delta.Methy"
  hc_intergenic_enhancer_genes <- hc_intergenic_enhancer_genes[!is.na(hc_intergenic_enhancer_genes$DEG.padj) & hc_intergenic_enhancer_genes$DEG.padj<= 0.05 & abs(hc_intergenic_enhancer_genes$DEG.LFC) >= 0.5,]
  write.csv(hc_intergenic_enhancer_genes, file="HighConfidence_Intergenic_Enhancer_Like_DMR_Gene_pairs_0.7.csv", row.names=F, quote=F)
  
  length(unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0]))#hypo_up 108
  X <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeI"])#4
  Y <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeII"])#1
  Z <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0 & hc_intergenic_enhancer_genes$LIHC.validation == "Validation"])#7
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)  
  
  length(unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0]))#hypo_up 108
  X <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeI"])
  Y <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeII"])
  Z <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0 & hc_intergenic_enhancer_genes$LIHC.validation == "Validation"])
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)
  
  length(unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0]))#hyper_up 21
  X <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeI"])#1
  Y <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeII"])#10
  Z <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy > 0 & hc_intergenic_enhancer_genes$DEG.LFC > 0 & hc_intergenic_enhancer_genes$LIHC.validation == "Validation"])#10
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)
  
  length(unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0]))#hyper_up 30
  X <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeI"])
  Y <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0 & hc_intergenic_enhancer_genes$LIHC.validation == "TypeII"])
  Z <- unique(hc_intergenic_enhancer_genes$gene_name[hc_intergenic_enhancer_genes$Delta.Methy < 0 & hc_intergenic_enhancer_genes$DEG.LFC < 0 & hc_intergenic_enhancer_genes$LIHC.validation == "Validation"])
  X <- X[!(X %in% Y) & !(X %in% Z)]
  Y <- Y[!(Y %in% Z)]
  length(X)
  length(Y)
  length(Z)
}





hyper_down <- unique(promoter_genes$gene_name[promoter_genes$Delta.Methy > 0 & promoter_genes$rho < 0])#75
hypo_up <- unique(promoter_genes$gene_name[promoter_genes$Delta.Methy < 0 & promoter_genes$rho < 0])#851
hyper_up <- unique(promoter_genes$gene_name[promoter_genes$Delta.Methy > 0 & promoter_genes$rho > 0])#77
hypo_down <- unique(promoter_genes$gene_name[promoter_genes$Delta.Methy < 0 & promoter_genes$rho > 0])#320
length(hyper_down)
length(hypo_up)
length(hyper_up)
length(hypo_down)


genic.hyper_down <- unique(genic_enhancer_genes$gene_name[genic_enhancer_genes$Delta.Methy > 0 & genic_enhancer_genes$rho < 0])#147
genic.hypo_up <- unique(genic_enhancer_genes$gene_name[genic_enhancer_genes$Delta.Methy < 0 & genic_enhancer_genes$rho < 0])#991
genic.hyper_up <- unique(genic_enhancer_genes$gene_name[genic_enhancer_genes$Delta.Methy > 0 & genic_enhancer_genes$rho > 0])#291
genic.hypo_down <- unique(genic_enhancer_genes$gene_name[genic_enhancer_genes$Delta.Methy < 0 & genic_enhancer_genes$rho > 0])#322
length(genic.hyper_down)
length(genic.hypo_up)
length(genic.hyper_up)
length(genic.hypo_down)

intergenic.hyper_down <- unique(intergenic_enhancer_genes$gene_name[intergenic_enhancer_genes$Delta.Methy > 0 & intergenic_enhancer_genes$MethyandGene.rho < 0])#10
intergenic.hypo_up <- unique(intergenic_enhancer_genes$gene_name[intergenic_enhancer_genes$Delta.Methy < 0 & intergenic_enhancer_genes$MethyandGene.rho < 0])#434
intergenic.hyper_up <- unique(intergenic_enhancer_genes$gene_name[intergenic_enhancer_genes$Delta.Methy > 0 & intergenic_enhancer_genes$MethyandGene.rho > 0])#19
intergenic.hypo_down <- unique(intergenic_enhancer_genes$gene_name[intergenic_enhancer_genes$Delta.Methy < 0 & intergenic_enhancer_genes$MethyandGene.rho > 0])#99
length(intergenic.hyper_down)
length(intergenic.hypo_up)
length(intergenic.hyper_up)
length(intergenic.hypo_down)

enhancer_hyper_down <- unique(c(genic.hyper_down, intergenic.hyper_down))#154
enhancer_hypo_up <- unique(c(genic.hypo_up, intergenic.hypo_up))#1236
enhancer_hyper_up <- unique(c(genic.hyper_up, intergenic.hyper_up))#308
enhancer_hypo_down <- unique(c(genic.hypo_down, intergenic.hypo_down))#395
length(enhancer_hyper_down)
length(enhancer_hypo_up)
length(enhancer_hyper_up)
length(enhancer_hypo_down)

df <- data.frame(DMR = rep(c("Promoter-like DMRs", "GenicEnhancer-like DMRs", "IntergenicEnhancer-like DMRs"), each =4),
                 Type = rep(c("Hyper & Down", "Hypo & Up", "Hyper & Up", "Hypo & Down"), 3),
                 Count = c(75,851,77,320,  147, 991, 291, 322,  10,434,19,99))
df$DMR <- factor(df$DMR, levels=c("Promoter-like DMRs", "GenicEnhancer-like DMRs", "IntergenicEnhancer-like DMRs"))
df$Type <- factor(df$Type, levels=c("Hyper & Down", "Hypo & Up", "Hyper & Up", "Hypo & Down"))
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
p_2.5 <- ggplot(df, aes(x = DMR, y = Count, fill = Type))+
  geom_bar(stat="identity", position = position_dodge())+scale_fill_brewer(palette = "Paired")+
  coord_flip()+
  xlab(NULL)+ylab("Count of DMR-DEGs")+
  theme_pubr()
pdf("Fig 3a.pdf", width = 11.69, height = 8.268)
p_2.5
dev.off()

#2.6 The pathway enrichment of promoter-like/enhancer-like DMR associated DEGs
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")
hypo_promoter_pathways <- read.csv("PathwayEnrichmentofHypomethylatedPromoterAssociatedActivatedDEGs.csv", header=T)
hypo_promoter_pathways$Name <- paste(hypo_promoter_pathways$Category, hypo_promoter_pathways$Description, sep=":")
hypo_promoter_pathways$logFDR <- -1*hypo_promoter_pathways$Log10.q.
hypo_promoter_pathways$Name <- factor(hypo_promoter_pathways$Name, levels = hypo_promoter_pathways$Name)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
p_2.6.1 <- ggplot(hypo_promoter_pathways, aes(x = Name, y = Count, fill = logFDR))+ 
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_gradient(low = brewer.pal(9,"Reds")[4], high = brewer.pal(9,"Reds")[9])+xlab(NULL)+ylab("Count of genes")+
  ggtitle("Hypomethylated promoter associated activated DEGs")+
  guides(fill=guide_legend(direction = "vertical",title = "-logFDR"))+
  coord_flip()+theme_pubr()
pdf("Fig 3b.pdf", width = 11.69, height = 8.268)
p_2.6.1
dev.off()

hypo_enhancer_pathways <- read.csv("PathwayEnrichmentofHypomethylatedEnhancerAssociatedActivatedDEGs.csv", header=T)
hypo_enhancer_pathways$Name <- paste(hypo_enhancer_pathways$Category, hypo_enhancer_pathways$Description, sep=":")
hypo_enhancer_pathways$logFDR <- -1*hypo_enhancer_pathways$Log10.q.
hypo_enhancer_pathways$Name <- factor(hypo_enhancer_pathways$Name, levels = hypo_enhancer_pathways$Name)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
p_2.6.2 <- ggplot(hypo_enhancer_pathways, aes(x = Name, y = Count, fill = logFDR))+ 
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_gradient(low = brewer.pal(9,"Reds")[4], high = brewer.pal(9,"Reds")[9])+xlab(NULL)+ylab("Count of genes")+
  ggtitle("Hypomethylated enhancer associated activated DEGs")+
  guides(fill=guide_legend(direction = "vertical",title = "-logFDR"))+
  coord_flip()+theme_pubr()
pdf("Fig 3c.pdf", width = 11.69, height = 8.268)
p_2.6.2
dev.off()

hyper_enhancer_pathways <- read.csv("PathwayEnrichmentofHyperMethylatedEnhancerAssociatedRepressedDEGs.csv", header=T)
hyper_enhancer_pathways$Name <- paste(hyper_enhancer_pathways$Category, hyper_enhancer_pathways$Description, sep=":")
hyper_enhancer_pathways$logFDR <- -1*hyper_enhancer_pathways$Log10.q.
hyper_enhancer_pathways$Name <- factor(hyper_enhancer_pathways$Name, levels = hyper_enhancer_pathways$Name)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
p_2.6.3 <- ggplot(hyper_enhancer_pathways, aes(x = Name, y = Count, fill = logFDR))+ 
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_gradient(low = brewer.pal(9,"Blues")[4], high = brewer.pal(9,"Blues")[9])+xlab(NULL)+ylab("Count of genes")+
  ggtitle("hypermethylated enhancer associated repressed DEGs")+
  guides(fill=guide_legend(direction = "vertical",title = "-logFDR"))+
  coord_flip()+theme_pubr()
pdf("Fig 3d.pdf", width = 11.69, height = 8.268)
p_2.6.3
dev.off()

library(patchwork)
p_2.6 <- plot_spacer()+p_2.6.3+p_2.6.2+p_2.6.1+plot_layout(guides = "collect", heights = c(7,5))
p_2.6
save(p_2.6, file="./Figures/Figure_2.6.RData")


## The bar plot of pathway hub genes validated by 5aza-dcR
#pathway1 - cell cycle 
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("CDC20", "CDKN3", "DLGAP5", "KIF2C", "BUB1B"), each =2), 
                 Condition = rep(c("0 ¦ÌM", "100 ¦ÌM"), 5),
                 Mean = c(1,1.29, 1, 3.20, 1, 1.36, 1, 1.49, 1, 2.66),
                 SD = c(0.079, 0.024, 0.12, 0.62, 0.24, 0.17, 0.046, 0.0021, 0.14, 0.68),
                 Pvalue = rep(c(0.00368, 0.003721, 0.09954, 0.0000549, 0.014233), each = 2)
                )


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the Cell Cycle pathway", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubr()+
  theme(text = element_text(color = "black"))
ggsave(filename = "Rplots.pdf", width = 9.5, height =8, units = "cm")
#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("CDC20", "CDKN3", "DLGAP5", "KIF2C", "BUB1B"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 5),
                 Mean = c(1,2.44, 1, 2.43, 1, 1.53, 1, 2.55, 1, 1.84),
                 SD = c(0.17, 0.59, 1.21, 1.08, 0.37, 1.80, 0.15, 0.55, 0.35, 0.25),
                 Pvalue = rep(c(0.01548, 0.20, 0.64, 0.008979, 0.027536), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the Cell Cycle pathway", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()

#pathway2 - DNA repair
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("POLR2B", "POLR2G", "POLR2D", "UBC", "LIG1", "RFC1"), each =2), 
                 Condition = rep(c("0 ¦ÌM", "100 ¦ÌM"), 6),
                 Mean = c(1, 1.54, 1, 2.1, 1, 1.44, 1, 1.69, 1, 1.96, 1, 2.62),
                 SD = c(0.1, 0.28, 0.21, 0.63, 0.04, 0.24, 0.32, 0.12, 0.08, 0.27, 0.49, 0.38),
                 Pvalue = rep(c(0.035, 0.046, 0.034, 0.024, 0.0039, 0.01), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the DNA Repair pathway", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubr()+
  theme(text = element_text(color = "black"))
ggsave(filename = "Rplots.pdf", width = 11.5, height =8, units = "cm")
#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("POLR2B", "POLR2G", "POLR2D", "UBC", "LIG1", "RFC1"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 6),
                 Mean = c(1,1.36, 1, 0.93, 1, 1.5, 1, 0.97, 1, 2.52, 1, 1.08),
                 SD = c(0.13, 0.1, 0.29, 0.1, 0.03, 0.15, 0.17, 0.09, 0.29, 0.3, 0.24, 0.17),
                 Pvalue = rep(c(0.0081, 0.72, 0.0044, 0.81, 0.0032, 0.65), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the DNA Repair pathway", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#pathway3 TP53 regulation
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("POLR2B", "POLR2G", "POLR2D", "UBC", "TAF5"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 5),
                 Mean = c(1, 1.54, 1, 2.1, 1, 1.44, 1, 1.69, 1, 2.34),
                 SD = c(0.1, 0.28, 0.21, 0.63, 0.04, 0.24, 0.32, 0.12, 0.16, 0.60),
                 Pvalue = rep(c(0.035, 0.046, 0.034, 0.024, 0.021), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the Transcriptional Regulation of TP53 pathway", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("POLR2B", "POLR2G", "POLR2D", "UBC", "TAF5"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 5),
                 Mean = c(1,1.36, 1, 0.93, 1, 1.5, 1, 0.97, 1, 0.97),
                 SD = c(0.13, 0.1, 0.29, 0.1, 0.03, 0.15, 0.17, 0.09, 0.41, 0.15),
                 Pvalue = rep(c(0.0081, 0.72, 0.0044, 0.81, 0.92), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the Transcriptional Regulation of TP53 pathway", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#pathway4 HIV infection
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("POLR2B", "POLR2G", "UBC", "NCBP2"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 4),
                 Mean = c(1, 1.54, 1, 2.1, 1, 1.69, 1, 1.81),
                 SD = c(0.1, 0.28, 0.21, 0.63, 0.32, 0.12, 0.16, 0.33),
                 Pvalue = rep(c(0.035, 0.046, 0.024, 0.018), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the HIV Infection pathway", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("POLR2B", "POLR2G", "UBC", "NCBP2"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 4),
                 Mean = c(1,1.36, 1, 0.93, 1, 0.97, 1, 0.95),
                 SD = c(0.13, 0.1, 0.29, 0.1, 0.17, 0.09, 0.34, 0.14),
                 Pvalue = rep(c(0.0081, 0.72, 0.81, 0.83), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the HIV Infection pathway", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#pathway5 Metabolic pathway
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("COMT", "ADH1A"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 2),
                 Mean = c(1, 0.72, 1, 1.49),
                 SD = c(0.19, 0.15, 0.046, 0.0021),
                 Pvalue = rep(c(0.11, 5.5E-5), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the Metabolic pathway", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("COMT", "ADH1A"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 2),
                 Mean = c(1, 2.74, 1, 1.86),
                 SD = c(0.36, 0.31, 0.066, 0.42),
                 Pvalue = rep(c(0.0031, 0.026), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Hub genes of the Metabolic pathway", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()

################# part 3 In silico validation of DMR associated DEGs #####################
#(Table) 3.1 The number of DMR associated DEGs passed in silico validation (TCGA-LIHC replication and 5-aza treated HCC lines validation)
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")
promoterGenes <- read.csv("HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", header=T)
library(fastSave)
load.lbzip2("Delta.gene.logTPM.RDataFS", n.cores = 32)

library(foreach)
LFC <- foreach(i = 1:nrow(promoterGenes), .combine = 'c')%do%{
  gene.id <- range_gene_DF$gene_id[range_gene_DF$gene_name==promoterGenes$gene_name[i]]
  LFC.median <- median(delta.expr[gene.id,])
  LFC.median
}
promoterGenes$LFC.median <- round(LFC,2)
write.csv(promoterGenes, file="HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", row.names=F)

promoterGenes$type <- rep("NA", nrow(promoterGenes))
promoterGenes$type[promoterGenes$Median.delta>0 & promoterGenes$rho < 0] <- "HyperDown"
promoterGenes$type[promoterGenes$Median.delta<0 & promoterGenes$rho < 0] <- "HypoUp"
promoterGenes$type[promoterGenes$Median.delta>0 & promoterGenes$rho > 0] <- "HyperUp"
promoterGenes$type[promoterGenes$Median.delta<0 & promoterGenes$rho > 0] <- "HypoDown"
length(unique(promoterGenes$gene_name[promoterGenes$type=="HyperDown" & promoterGenes$lfc.up > 0]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HypoUp" & promoterGenes$lfc.up > 0]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HyperUp" & promoterGenes$lfc.down > 0]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HypoDown" & promoterGenes$lfc.down > 0]))

length(unique(promoterGenes$gene_name[promoterGenes$type=="HyperDown" & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HypoUp" & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HyperUp" & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HypoDown" & promoterGenes$LIHCValidation]))

length(unique(promoterGenes$gene_name[promoterGenes$type=="HyperDown" & promoterGenes$lfc.up > 0 & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HypoUp" & promoterGenes$lfc.up > 0 & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HyperUp" & promoterGenes$lfc.down > 0 & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$type=="HypoDown" & promoterGenes$lfc.down > 0 & promoterGenes$LIHCValidation]))
genicEnhancerGenes <- read.csv("HighConfidence_GenicEnhancer_Like_DMR_Gene_pairs_0.8.csv", header = T)
library(foreach)
LFC <- foreach(i = 1:nrow(genicEnhancerGenes), .combine = 'c')%do%{
  gene.id <- range_gene_DF$gene_id[range_gene_DF$gene_name==genicEnhancerGenes$gene_name[i]]
  LFC.median <- median(delta.expr[gene.id,])
  LFC.median
}
genicEnhancerGenes$LFC.median <- round(LFC,2)
write.csv(genicEnhancerGenes, file="HighConfidence_GenicEnhancer_Like_DMR_Gene_pairs_0.8.csv", row.names=F)
genicEnhancerGenes$type <- rep("NA", nrow(genicEnhancerGenes))
genicEnhancerGenes$type[genicEnhancerGenes$Median.delta>0 & genicEnhancerGenes$rho < 0] <- "HyperDown"
genicEnhancerGenes$type[genicEnhancerGenes$Median.delta<0 & genicEnhancerGenes$rho < 0] <- "HypoUp"
genicEnhancerGenes$type[genicEnhancerGenes$Median.delta>0 & genicEnhancerGenes$rho > 0] <- "HyperUp"
genicEnhancerGenes$type[genicEnhancerGenes$Median.delta<0 & genicEnhancerGenes$rho > 0] <- "HypoDown"
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HyperDown" & genicEnhancerGenes$lfc.up > 0]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HypoUp" & genicEnhancerGenes$lfc.up > 0]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HyperUp" & genicEnhancerGenes$lfc.down > 0]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HypoDown" & genicEnhancerGenes$lfc.down > 0]))

length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HyperDown" & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HypoUp" & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HyperUp" & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HypoDown" & genicEnhancerGenes$LIHCValidation]))

length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HyperDown" & genicEnhancerGenes$lfc.up > 0 & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HypoUp" & genicEnhancerGenes$lfc.up > 0 & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HyperUp" & genicEnhancerGenes$lfc.down > 0 & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$type=="HypoDown" & genicEnhancerGenes$lfc.down > 0 & genicEnhancerGenes$LIHCValidation]))
intergenicEnhancerGenes <- read.csv("HighConfidence_InterGenicEnhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)
library(foreach)
LFC <- foreach(i = 1:nrow(intergenicEnhancerGenes), .combine = 'c')%do%{
  gene.id <- range_gene_DF$gene_id[range_gene_DF$gene_name==intergenicEnhancerGenes$gene_name[i]]
  LFC.median <- median(delta.expr[gene.id,])
  LFC.median
}
intergenicEnhancerGenes$LFC.median <- round(LFC,2)
write.csv(intergenicEnhancerGenes, file="HighConfidence_InterGenicEnhancer_Like_DMR_Gene_pairs_0.7.csv", row.names=F)

intergenicEnhancerGenes$type <- rep("NA", nrow(intergenicEnhancerGenes))
intergenicEnhancerGenes$type[intergenicEnhancerGenes$DM.Median>0 & intergenicEnhancerGenes$MethyandGene.rho < 0] <- "HyperDown"
intergenicEnhancerGenes$type[intergenicEnhancerGenes$DM.Median<0 & intergenicEnhancerGenes$MethyandGene.rho < 0] <- "HypoUp"
intergenicEnhancerGenes$type[intergenicEnhancerGenes$DM.Median>0 & intergenicEnhancerGenes$MethyandGene.rho > 0] <- "HyperUp"
intergenicEnhancerGenes$type[intergenicEnhancerGenes$DM.Median<0 & intergenicEnhancerGenes$MethyandGene.rho > 0] <- "HypoDown"
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HyperDown" & intergenicEnhancerGenes$lfc.up > 0]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HypoUp" & intergenicEnhancerGenes$lfc.up > 0]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HyperUp" & intergenicEnhancerGenes$lfc.down > 0]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HypoDown" & intergenicEnhancerGenes$lfc.down > 0]))

length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HyperDown" & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HypoUp" & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HyperUp" & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HypoDown" & intergenicEnhancerGenes$LIHCValidation]))

length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HyperDown" & intergenicEnhancerGenes$lfc.up > 0 & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HypoUp" & intergenicEnhancerGenes$lfc.up > 0 & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HyperUp" & intergenicEnhancerGenes$lfc.down > 0 & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$type=="HypoDown" & intergenicEnhancerGenes$lfc.down > 0 & intergenicEnhancerGenes$LIHCValidation]))

#High confidence DMR-genes with overlap ratio >= 80% and gene |LFC| > 1
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta > 0 & promoterGenes$rho < 0 & promoterGenes$LFC.median <= -1]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta < 0 & promoterGenes$rho < 0 & promoterGenes$LFC.median >= 1]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta > 0 & promoterGenes$rho > 0 & promoterGenes$LFC.median >= 1]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta < 0 & promoterGenes$rho > 0 & promoterGenes$LFC.median <= -1]))
#The LIHC validation ratio
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta > 0 & promoterGenes$rho < 0 & promoterGenes$LFC.median <= -1 & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta < 0 & promoterGenes$rho < 0 & promoterGenes$LFC.median >= 1 & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta > 0 & promoterGenes$rho > 0 & promoterGenes$LFC.median >= 1 & promoterGenes$LIHCValidation]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta < 0 & promoterGenes$rho > 0 & promoterGenes$LFC.median <= -1 & promoterGenes$LIHCValidation]))
#The 5-aza-dcR demethylation validation
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta > 0 & promoterGenes$rho < 0 & promoterGenes$LFC.median <= -1 & promoterGenes$lfc.up > 0]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta < 0 & promoterGenes$rho < 0 & promoterGenes$LFC.median >= 1 & promoterGenes$lfc.up > 0]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta > 0 & promoterGenes$rho > 0 & promoterGenes$LFC.median >= 1 & promoterGenes$lfc.down > 0]))
length(unique(promoterGenes$gene_name[promoterGenes$Median.delta < 0 & promoterGenes$rho > 0 & promoterGenes$LFC.median <= -1 & promoterGenes$lfc.down > 0]))

#High confidence DMR-genes with overlap ratio >= 80% and gene |LFC| > 1
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta > 0 & genicEnhancerGenes$rho < 0 & genicEnhancerGenes$LFC.median <= -1]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta < 0 & genicEnhancerGenes$rho < 0 & genicEnhancerGenes$LFC.median >= 1]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta > 0 & genicEnhancerGenes$rho > 0 & genicEnhancerGenes$LFC.median >= 1]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta < 0 & genicEnhancerGenes$rho > 0 & genicEnhancerGenes$LFC.median <= -1]))
#The LIHC validation ratio
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta > 0 & genicEnhancerGenes$rho < 0 & genicEnhancerGenes$LFC.median <= -1 & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta < 0 & genicEnhancerGenes$rho < 0 & genicEnhancerGenes$LFC.median >= 1 & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta > 0 & genicEnhancerGenes$rho > 0 & genicEnhancerGenes$LFC.median >= 1 & genicEnhancerGenes$LIHCValidation]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta < 0 & genicEnhancerGenes$rho > 0 & genicEnhancerGenes$LFC.median <= -1 & genicEnhancerGenes$LIHCValidation]))
#The 5-aza-dcR demethylation validation
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta > 0 & genicEnhancerGenes$rho < 0 & genicEnhancerGenes$LFC.median <= -1 & genicEnhancerGenes$lfc.up > 0]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta < 0 & genicEnhancerGenes$rho < 0 & genicEnhancerGenes$LFC.median >= 1 & genicEnhancerGenes$lfc.up > 0]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta > 0 & genicEnhancerGenes$rho > 0 & genicEnhancerGenes$LFC.median >= 1 & genicEnhancerGenes$lfc.down > 0]))
length(unique(genicEnhancerGenes$gene_name[genicEnhancerGenes$Median.delta < 0 & genicEnhancerGenes$rho > 0 & genicEnhancerGenes$LFC.median <= -1 & genicEnhancerGenes$lfc.down > 0]))

#High confidence DMR-genes with overlap ratio >= 80% and gene |LFC| > 1
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median > 0 & intergenicEnhancerGenes$MethyandGene.rho < 0 & intergenicEnhancerGenes$LFC.median <= -1]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median < 0 & intergenicEnhancerGenes$MethyandGene.rho < 0 & intergenicEnhancerGenes$LFC.median >= 1]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median > 0 & intergenicEnhancerGenes$MethyandGene.rho > 0 & intergenicEnhancerGenes$LFC.median >= 1]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median < 0 & intergenicEnhancerGenes$MethyandGene.rho > 0 & intergenicEnhancerGenes$LFC.median <= -1]))
#The LIHC validation ratio
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median > 0 & intergenicEnhancerGenes$MethyandGene.rho < 0 & intergenicEnhancerGenes$LFC.median <= -1 & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median < 0 & intergenicEnhancerGenes$MethyandGene.rho < 0 & intergenicEnhancerGenes$LFC.median >= 1 & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median > 0 & intergenicEnhancerGenes$MethyandGene.rho > 0 & intergenicEnhancerGenes$LFC.median >= 1 & intergenicEnhancerGenes$LIHCValidation]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median < 0 & intergenicEnhancerGenes$MethyandGene.rho > 0 & intergenicEnhancerGenes$LFC.median <= -1 & intergenicEnhancerGenes$LIHCValidation]))
#The 5-aza-dcR demethylation validation
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median > 0 & intergenicEnhancerGenes$MethyandGene.rho < 0 & intergenicEnhancerGenes$LFC.median <= -1 & intergenicEnhancerGenes$lfc.up > 0]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median < 0 & intergenicEnhancerGenes$MethyandGene.rho < 0 & intergenicEnhancerGenes$LFC.median >= 1 & intergenicEnhancerGenes$lfc.up > 0]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median > 0 & intergenicEnhancerGenes$MethyandGene.rho > 0 & intergenicEnhancerGenes$LFC.median >= 1 & intergenicEnhancerGenes$lfc.down > 0]))
length(unique(intergenicEnhancerGenes$gene_name[intergenicEnhancerGenes$DM.Median < 0 & intergenicEnhancerGenes$MethyandGene.rho > 0 & intergenicEnhancerGenes$LFC.median <= -1 & intergenicEnhancerGenes$lfc.down > 0]))
#(Table) 3.2 The summarized information of top candidate DMR-Gene pairs that passed in silico validation 
#Finished
#(Table) 3.2 Significant independent prognostic biomarkers in DMR-associated DEGs.
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")
promoter.prognostic <- read.csv("ClinicalRelevance.LIHC_validated_promoterCpG_Genes.csv", header=T)
#OS
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$OS.pvalue.methy < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$OS.pvalue.methy < 0.05])
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$OS.pvalue.expr < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$OS.pvalue.expr < 0.05])
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$OS.pvalue.methy < 0.05 & promoter.prognostic$OS.pvalue.expr < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$OS.pvalue.methy < 0.05 & promoter.prognostic$OS.pvalue.expr < 0.05])
#PFS
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$PFS.pvalue.methy < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$PFS.pvalue.methy < 0.05])
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$PFS.pvalue.expr < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$PFS.pvalue.expr < 0.05])
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$PFS.pvalue.methy < 0.05 & promoter.prognostic$PFS.pvalue.expr < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$PFS.pvalue.methy < 0.05 & promoter.prognostic$PFS.pvalue.expr < 0.05])
#stage
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$stage.pvalue.methy < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$stage.pvalue.methy < 0.05])
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$stage.pvalue.expr < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$stage.pvalue.expr < 0.05])
length(unique(promoter.prognostic$Gene_name[promoter.prognostic$stage.pvalue.methy < 0.05 & promoter.prognostic$stage.pvalue.expr < 0.05]))
table(promoter.prognostic$Gene_name[promoter.prognostic$stage.pvalue.methy < 0.05 & promoter.prognostic$stage.pvalue.expr < 0.05])

genicEnhancer.prognostic <- read.csv("ClinicalRelevance.LIHC_validated_GenicEnhancerCpG_Genes.csv", header=T)
#OS
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$OS.pvalue.methy < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$OS.pvalue.methy < 0.05])
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$OS.pvalue.expr < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$OS.pvalue.expr < 0.05])
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$OS.pvalue.methy < 0.05 & genicEnhancer.prognostic$OS.pvalue.expr < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$OS.pvalue.methy < 0.05 & genicEnhancer.prognostic$OS.pvalue.expr < 0.05])
#PFS
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$PFS.pvalue.methy < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$PFS.pvalue.methy < 0.05])
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$PFS.pvalue.expr < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$PFS.pvalue.expr < 0.05])
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$PFS.pvalue.methy < 0.05 & genicEnhancer.prognostic$PFS.pvalue.expr < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$PFS.pvalue.methy < 0.05 & genicEnhancer.prognostic$PFS.pvalue.expr < 0.05])
#stage
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$stage.pvalue.methy < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$stage.pvalue.methy < 0.05])
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$stage.pvalue.expr < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$stage.pvalue.expr < 0.05])
length(unique(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$stage.pvalue.methy < 0.05 & genicEnhancer.prognostic$stage.pvalue.expr < 0.05]))
table(genicEnhancer.prognostic$Gene_name[genicEnhancer.prognostic$stage.pvalue.methy < 0.05 & genicEnhancer.prognostic$stage.pvalue.expr < 0.05])

intergenicEnhancer.prognostic <- read.csv("ClinicalRelevance.LIHC_validated_InterGenicEnhancerCpG_Genes.csv", header=T)
#OS
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$OS.pvalue.methy < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$OS.pvalue.methy < 0.05])
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$OS.pvalue.expr < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$OS.pvalue.expr < 0.05])
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$OS.pvalue.methy < 0.05 & intergenicEnhancer.prognostic$OS.pvalue.expr < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$OS.pvalue.methy < 0.05 & intergenicEnhancer.prognostic$OS.pvalue.expr < 0.05])
#PFS
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$PFS.pvalue.methy < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$PFS.pvalue.methy < 0.05])
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$PFS.pvalue.expr < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$PFS.pvalue.expr < 0.05])
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$PFS.pvalue.methy < 0.05 & intergenicEnhancer.prognostic$PFS.pvalue.expr < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$PFS.pvalue.methy < 0.05 & intergenicEnhancer.prognostic$PFS.pvalue.expr < 0.05])
#stage
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$stage.pvalue.methy < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$stage.pvalue.methy < 0.05])
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$stage.pvalue.expr < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$stage.pvalue.expr < 0.05])
length(unique(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$stage.pvalue.methy < 0.05 & intergenicEnhancer.prognostic$stage.pvalue.expr < 0.05]))
table(intergenicEnhancer.prognostic$Gene_name[intergenicEnhancer.prognostic$stage.pvalue.methy < 0.05 & intergenicEnhancer.prognostic$stage.pvalue.expr < 0.05])

####make venn disagram for the OS, PFS, and tumor stage
#load methylation and expression level ClinicalRelevance of 1335 high confident methylation driven genes
ClinicalRelevance.methy <- read.csv("/data/huangp/Methy/Mixed/New/ClinicalRelevance_highconfident_methylation_driven_genes_at_methylation_level.csv", header=T)
#220genes
ClinicalRelevance.expr <- read.csv("/data/huangp/Methy/Mixed/New/ClinicalRelevance_highconfident_methylation_driven_genes_at_expression_level.csv", header=T)
#941genes
OS.expr <- ClinicalRelevance.expr$Gene[ClinicalRelevance.expr$OS.effect %in% c("Favorable", "Unfavorable")]
#168
PFS.expr <- ClinicalRelevance.expr$Gene[ClinicalRelevance.expr$PFS.effect %in% c("Favorable", "Unfavorable")]
#135
Stage.expr <- ClinicalRelevance.expr$Gene[ClinicalRelevance.expr$Stage.effect %in% c("Inverse", "Positive")]
#217
hc_promoter_genes <- read.csv("HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", header = T)
hc_genic_enhancer_genes <- read.csv("HighConfidence_genic_Enhancer_Like_DMR_Gene_pairs_0.8.csv", header=T)
hc_intergenic_enhancer_genes <- read.csv("HighConfidence_Intergenic_Enhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)
OS.expr <- OS.expr[OS.expr %in% hc_promoter_genes$gene_name | 
                     OS.expr %in% hc_genic_enhancer_genes$gene_name | 
                     OS.expr %in% hc_intergenic_enhancer_genes$gene_name]#108
PFS.expr <- PFS.expr[PFS.expr %in% hc_promoter_genes$gene_name | 
                     PFS.expr %in% hc_genic_enhancer_genes$gene_name | 
                     PFS.expr %in% hc_intergenic_enhancer_genes$gene_name]#82
Stage.expr <- Stage.expr[Stage.expr %in% hc_promoter_genes$gene_name | 
                     Stage.expr %in% hc_genic_enhancer_genes$gene_name | 
                     Stage.expr %in% hc_intergenic_enhancer_genes$gene_name]#140
library(foreach)
hc_promoter_genes_pairs <- foreach(i = 1:nrow(hc_promoter_genes), .combine = 'c') %do%
  {
    DMRID = hc_promoter_genes$DMR.ID[i]
    DMR.chr = strsplit(DMRID, split="_")[[1]][1]
    DMR.str = as.integer(strsplit(DMRID, split="_")[[1]][2])
    DMR.end = as.integer(strsplit(DMRID, split="_")[[1]][3])
    DMR.CpGs = paste(DMR.chr, DMR.str:DMR.end, sep="_")
    pairs <- paste(DMR.CpGs, hc_promoter_genes$gene_name[i], sep="|")
  }
hc_genic_enhancer_genes_pairs <- foreach(i = 1:nrow(hc_genic_enhancer_genes), .combine = 'c') %do%
  {
    DMRID = hc_genic_enhancer_genes$DMR.ID[i]
    DMR.chr = strsplit(DMRID, split="_")[[1]][1]
    DMR.str = as.integer(strsplit(DMRID, split="_")[[1]][2])
    DMR.end = as.integer(strsplit(DMRID, split="_")[[1]][3])
    DMR.CpGs = paste(DMR.chr, DMR.str:DMR.end, sep="_")
    pairs <- paste(DMR.CpGs, hc_genic_enhancer_genes$gene_name[i], sep="|")
  }
hc_intergenic_enhancer_genes_pairs <- foreach(i = 1:nrow(hc_intergenic_enhancer_genes), .combine = 'c') %do%
  {
    DMRID = hc_intergenic_enhancer_genes$DMR.ID[i]
    DMR.chr = strsplit(DMRID, split="_")[[1]][1]
    DMR.str = as.integer(strsplit(DMRID, split="_")[[1]][2])
    DMR.end = as.integer(strsplit(DMRID, split="_")[[1]][3])
    DMR.CpGs = paste(DMR.chr, DMR.str:DMR.end, sep="_")
    pairs <- paste(DMR.CpGs, hc_intergenic_enhancer_genes$gene_name[i], sep="|")
  }
ClinicalRelevance.methy$pairs <- paste(ClinicalRelevance.methy$CpG, ClinicalRelevance.methy$Gene, sep="|")
table(ClinicalRelevance.methy$pairs %in% c(hc_promoter_genes_pairs, hc_genic_enhancer_genes_pairs, hc_intergenic_enhancer_genes_pairs))
ClinicalRelevance.methy <- ClinicalRelevance.methy[ClinicalRelevance.methy$pairs %in% c(hc_promoter_genes_pairs, hc_genic_enhancer_genes_pairs, hc_intergenic_enhancer_genes_pairs),]
OS.methy <- unique(ClinicalRelevance.methy$Gene[ClinicalRelevance.methy$OS.effect %in% c("Favorable", "Unfavorable")])
#27
PFS.methy <- unique(ClinicalRelevance.methy$Gene[ClinicalRelevance.methy$PFS.effect %in% c("Favorable", "Unfavorable")])
#30
Stage.methy <- unique(ClinicalRelevance.methy$Gene[ClinicalRelevance.methy$Stage.effect %in% c("Inverse", "Positive")])
#18

OS.both <- intersect(OS.expr, OS.methy)
#6
PFS.both <- intersect(PFS.expr, PFS.methy)
#8
Stage.both <- intersect(Stage.expr, Stage.methy)
#8
write.table(OS.expr, file="/data/huangp/Methy/Mixed/New/OS.expr.txt", sep="\t", col.names=F, row.names=F)
write.table(PFS.expr, file="/data/huangp/Methy/Mixed/New/PFS.expr.txt", sep="\t", col.names=F, row.names=F)
write.table(Stage.expr, file="/data/huangp/Methy/Mixed/New/Stage.expr.txt", sep="\t", col.names=F, row.names=F)
write.table(OS.methy, file="/data/huangp/Methy/Mixed/New/OS.methy.txt", sep="\t", col.names=F, row.names=F)
write.table(PFS.methy, file="/data/huangp/Methy/Mixed/New/PFS.methy.txt", sep="\t", col.names=F, row.names=F)
write.table(Stage.methy, file="/data/huangp/Methy/Mixed/New/Stage.methy.txt", sep="\t", col.names=F, row.names=F)

options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/LIHC")
#load Overall Survival, Progression-free Survival phenotypes of TCGA
load("LIHC_TCGA.phenotypes.RData")
### load methy
library(fastSave)
load.lbzip2("methy.matrix.LIHC.450k.RDataFS", n.cores=32)
### transform hg19 location to hg38 location
hg19tohg38 <- read.table("LIHC.TCGA.450k.hg38.bed", header=F, sep="\t")
colnames(hg19tohg38) <- c("Chr", "Start", "End", "hg19_pos")
hg19tohg38$hg38_pos <- paste(hg19tohg38$Chr, hg19tohg38$Start, sep="_")
LIHC.methy <- LIHC.methy[rownames(LIHC.methy) %in% hg19tohg38$hg19_pos,]
rownames(LIHC.methy) <- hg19tohg38$hg38_pos[match(rownames(LIHC.methy), hg19tohg38$hg19_pos)]
### load expr 
load("LIHC.rnaseq.expr.RData")
dim(LIHC.rnaseq) #423  20532
LIHC.rnaseq$bcr_patient_barcode <- as.character(LIHC.rnaseq$bcr_patient_barcode)
LIHC.rnaseq$bcr_patient_barcode <- sapply(LIHC.rnaseq$bcr_patient_barcode, substr, 1, 15)
x <- t(LIHC.rnaseq)
y <- x[-1,]
colnames(y) <- x[1,]
z <- as.data.frame(apply(y, 2, as.numeric))
rownames(z) <- rownames(y)
z <- cbind(Gene_name = sapply(strsplit(rownames(z), split="|", fixed=T), "[", 1), 
           Gene_ID = sapply(strsplit(rownames(z), split="|", fixed=T), "[", 2),
           z)
LIHC.expr <- z 
# CDC20 chr1_43360233
# OS
CpGID = "chr1_43360233"
Gene = "CDC20"

#UCK2 chr1_165857033
CpGID = "chr1_165857033"
Gene = "UCK2"

#HEART6
#SLC9A3R2

expr = unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name),-(1:2)])
expr <- log2(expr+0.0001)
methy = unlist(LIHC.methy[CpGID,])
library(survival)
library(survminer)
#cox for methy
df.methy <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
df.methy <- df.methy[df.methy$sampleID %in% names(methy),]#keep the sample with both methylation and phenotypes
methy.scaled <- scale(methy)[,1]# transform the methylation data to Z score
df.methy$methy <- round(methy[match(df.methy$sampleID, names(methy))],2)
df.methy$methy.scaled <- round(methy.scaled[match(df.methy$sampleID, names(methy.scaled))],2)
df.methy <- df.methy[df.methy$condition == "Tumor" & df.methy$stage2 != "Stage IV",]#keep only tumor samples & stage I-III
df.methy$stage2 <- factor(df.methy$stage2, levels = c("Stage I", "Stage II", "Stage III"))
df.methy <- df.methy[!is.na(df.methy$age) & !is.na(df.methy$race) & !is.na(df.methy$gender) & !is.na(df.methy$stage2) & !is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                       !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy) & !is.na(df.methy$methy.scaled),]
# set Z-scale cut-offs for high and low expression
#highExpr <- 1.0
#lowExpr <- -1.0
#df.methy$methy.group <- ifelse(df.methy$methy.scaled >= highExpr, 'High', ifelse(df.methy$methy.scaled <= lowExpr, 'Low', 'Mid'))
#df.methy$methy.group <- ifelse(df.methy$methy.scaled >= median(df.methy$methy.scaled), 'High','Low')
# relevel the factors to have mid as the ref level
#df.methy$methy.group <- factor(df.methy$methy.group, levels = c("High", "Low", "Mid"))
df.methy$methy.group <- ifelse(df.methy$methy >= median(df.methy$methy), "HighMethy", "LowMethy")
df.methy$methy.group <- factor(df.methy$methy.group, levels = c("HighMethy", "LowMethy"))

#cox for expr
df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
df.expr <- df.expr[df.expr$condition == "Tumor" & df.expr$stage2 != "Stage IV",]#keep only tumor samples & stage I-III
df.expr$stage2 <- factor(df.expr$stage2, levels = c("Stage I", "Stage II", "Stage III"))
df.expr <- df.expr[!is.na(df.expr$age) & !is.na(df.expr$race) & !is.na(df.expr$gender) & !is.na(df.expr$stage2) & !is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                     !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr) & !is.na(df.expr$expr.scaled),]
# set Z-scale cut-offs for high and low expression
#highExpr <- 1.0
#lowExpr <- -1.0
#df.expr$expr.group <- ifelse(df.expr$expr.scaled >= highExpr, 'High',
#                             ifelse(df.expr$expr.scaled <= lowExpr, 'Low', 'Mid'))

# relevel the factors to have mid as the ref level
#df.expr$expr.group <- factor(df.expr$expr.group,
#                             levels = c('High', 'Low', 'Mid'))
df.expr$expr.group <- ifelse(df.expr$expr >= median(df.expr$expr), "HighExpr", "LowExpr")
df.expr$expr.group <- factor(df.expr$expr.group, levels = c("HighExpr", "LowExpr"))

df.methy <- df.methy[df.methy$sampleID %in% df.expr$sampleID,]
df <- cbind(df.expr, df.methy[match(df.expr$sampleID, df.methy$sampleID),11:13])
df$group <- paste(df$methy.group, df$expr.group, sep=" & ")
table(df$group)
df.new <-  df[df$group %in% c("HighMethy & LowExpr", "LowMethy & HighExpr"),]
colnames(df.new)[17] <- "UCK2"

library(survminer)
library(ggpubr)
library(ggplot2)
OS <- ggsurvplot(survfit(Surv(OS.time, OS) ~ UCK2,
                              data = df.new),
                      data = df.new,
                      risk.table = F,
                      pval = TRUE,
                      palette = "npg",
                      break.time.by = 500,
                      ggtheme = theme(text=element_text(size=18),
                                      legend.direction = "vertical",
                                      panel.grid.major =element_blank(), 
                                      panel.grid.minor = element_blank(), 
                                      panel.background = element_blank(),
                                      axis.line = element_line(colour = "black")),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE
)+ylab("Overall survival probability")
OS
dev.off()
ggsave(file="Rplots.pdf", height =12 , width =19.2 , units = "cm")

PFS <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ UCK2,
                         data = df.new),
                 data = df.new,
                 risk.table = F,
                 pval = TRUE,
                 palette = "npg",
                 break.time.by = 500,
                 ggtheme = theme(text=element_text(size=18),
                                 legend.direction = "vertical",
                                 panel.grid.major =element_blank(), 
                                 panel.grid.minor = element_blank(), 
                                 panel.background = element_blank(),
                                 axis.line = element_line(colour = "black")),                 risk.table.y.text.col = TRUE,
                 risk.table.y.text = FALSE
)+ylab("Progression-free survival probability")
PFS
dev.off()
ggsave(file="Rplots.pdf", height =12 , width =19.2 , units = "cm")


# boxplot of methy in different tumor stage
stage.methy <- ggplot(data=df.methy, aes(x = stage2, y = methy))+
  geom_jitter(aes(fill = stage2), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_brewer(palette = "Reds", guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("Tumor stage")+
  ylab(paste("Methylation of ",Gene,  "(CpG:", CpGID, ")", sep=""))+
  labs(subtitle = "p = 2.5e-2")+
  theme_pubr()+  
  theme(text = element_text(size = 24))
stage.methy
dev.off()

# boxplot of expr in different tumor stage
stage.expr <- ggplot(data=df.expr, aes(x = stage2, y = expr))+
  geom_jitter(aes(fill = stage2), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_brewer(palette = "Reds", guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("Tumor stage")+
  ylab(paste("mRNA expression of ", Gene, sep=""))+
  labs(subtitle = "p = 4.3e-3")+
  theme_pubr()+
  theme(text = element_text(size = 24))

stage.expr  
dev.off()


## The bar plot of prognostic markers


# Tumor stages associated genes
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("KIF15", "SLC22A15", "CDKN2C","CDC20", "UCK2", "HEATR6", "SLC9A3R2", "KIAA0100", "SLC16A3"), each =2), 
                 Condition = rep(c("0 ¦ÌM", "100 ¦ÌM"),9),
                 Mean = c(1,1.49, 1,2.61, 1,1.57, 
                          1,1.29, 1,1.61, 1,1.76, 1,1.51,
                          1, 1.4, 1, 1.27),
                 SD = c( 0.046,0.0021, 0.47,0.17, 0.18,0.23, 
                         0.08,0.02, 0.13,0.19, 0.092,0.32, 0.14,0.14, 
                         0.13, 0.11, 0.33, 0.33),
                 Pvalue = rep(c(5e-5, 0.0051, 0.028, 
                                0.0037, 0.011, 0.016, 0.011,
                                0.015, 0.38), each = 2)
)
df$Gene <- factor(df$Gene, levels = c("KIF15", "SLC22A15", "CDKN2C","CDC20", "UCK2", "HEATR6", "SLC9A3R2", "KIAA0100", "SLC16A3"))

p <- ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x = "", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubr()+
  theme(text = element_text(color = "black"))
ggsave(filename = "Rplots.pdf", width = 20, height = 8, units = "cm")


#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("HEATR6", "CDC20", "SLC9A3R2", "UCK2", "KIAA0100", "SLC16A3"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 6),
                 Mean = c(1, 2.54, 1, 2.44, 1, 1.25, 1, 1.78, 1, 1.61, 1, 2.23),
                 SD = c(0.2, 0.53, 0.17, 0.59, 0.06, 0.09, 0.22, 0.4, 0.29, 0.21, 0.57, 0.96),
                 Pvalue = rep(c(0.0094, 0.015, 0.015, 0.044, 0.042, 0.13), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Tumor stages associated genes", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
# PFS associated genes
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("CDC20", "SLC9A3R2", "UCK2", "HEATR6", "CDKN2C"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 5),
                 Mean = c(1, 1.29, 1, 1.51, 1, 1.61, 1, 1.76, 1, 1.57),
                 SD = c(0.08, 0.02, 0.14, 0.14, 0.13, 0.19, 0.092, 0.32, 0.18, 0.23),
                 Pvalue = rep(c(0.0037, 0.011, 0.011, 0.016, 0.028), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Progression-free Survival associated genes", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("CDC20", "SLC9A3R2", "UCK2", "HEATR6", "CDKN2C"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 5),
                 Mean = c(1, 2.44, 1, 1.25, 1, 1.78, 1, 1.2, 1, 2.16),
                 SD = c(0.17, 0.59, 0.06, 0.09, 0.22, 0.4, 0.17, 0.19, 0.12, 0.3),
                 Pvalue = rep(c(0.015, 0.015, 0.044, 0.25, 0.0032), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Progression-free Survival associated genes", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
# OS associated genes
#LO2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("CDC20", "SLC9A3R2", "UCK2", "KIF15", "SLC22A15"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 5),
                 Mean = c(1, 1.29, 1, 1.51, 1, 1.61, 1, 1.49, 1, 2.61),
                 SD = c(0.08, 0.02, 0.14, 0.14, 0.13, 0.19, 0.046, 0.0021, 0.47, 0.17),
                 Pvalue = rep(c(0.0037, 0.011, 0.011, 5.5E-5, 0.0051), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Overall Survival associated genes", y = "Relative expression", title = "LO2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()
#HepG2
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
options(stringsAsFactors = F)
df <- data.frame(Gene = rep(c("CDC20", "SLC9A3R2", "UCK2", "KIF15", "SLC22A15"), each =2), 
                 Condition = rep(c("0¦ÌM", "100¦ÌM"), 5),
                 Mean = c(1, 2.44, 1, 1.25, 1, 1.78, 1, 2.63, 1, 2.28),
                 SD = c(0.17, 0.59, 0.06, 0.09, 0.22, 0.4, 0.07, 0.55, 0.37, 0.44),
                 Pvalue = rep(c(0.015, 0.015, 0.044, 0.0071, 0.018), each = 2)
)


ggplot(data = df, aes(Gene, Mean, fill = Condition))+
  geom_bar(stat = "identity", color = "black", position = "dodge", width = 0.6)+
  scale_fill_manual(values = c("white", "grey50"))+
  geom_errorbar(aes(ymin = Mean, ymax = (Mean+SD)), position = position_dodge(0.6), width = 0.2)+
  labs(x= "Overall Survival associated genes", y = "Relative expression", title = "HepG2")+
  guides(fill = guide_legend(title = "5aza-dcR"))+
  theme_pubclean()

# the prioritation of genes for RT-PCR validation
library(fastSave)
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")
load.lbzip2("../Gene.TPM.filtered.RDataFS", n.cores=32)
library(genomation)
library(GenomicFeatures)
Granges.gtf <- gffToGRanges("../gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
Granges.gene <- Granges.gene[match(rownames(gene.TPM), Granges.gene$gene_id),]
range_gene_DF <- as.data.frame(Granges.gene)
load("../All66HCC.phenotypes.RData")
tumor.gene.TPM <- gene.TPM[,match(HCC.phenotypes$SampleID[1:33], colnames(gene.TPM))]
normal.gene.TPM <- gene.TPM[,match(HCC.phenotypes$SampleID[34:66], colnames(gene.TPM))]
delta.expr <- log2(tumor.gene.TPM+0.001) - log2(normal.gene.TPM+0.001)
colnames(delta.expr) <- HCC.phenotypes$patientID[1:33]
save.lbzip2(delta.expr, file="Delta.gene.logTPM.RDataFS", n.cores = 32)

promoterGenes <- read.csv("HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", header=T)
genicEnhancerGenes <- read.csv("HighConfidence_GenicEnhancer_Like_DMR_Gene_pairs_0.8.csv", header = T)
intergenicEnhancerGenes <- read.csv("HighConfidence_InterGenicEnhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)

pcr_genes <- read.csv("GenesforRTPCRValidation.csv", header=T)
library(foreach)
pcr_info <- foreach(i = 1:nrow(pcr_genes), .combine = "rbind")%do%{
  gene = pcr_genes$Gene[i]
  gene.ID <- range_gene_DF$gene_id[range_gene_DF$gene_name == gene]
  Median.LFC <- median(delta.expr[gene.ID,])
  if(gene %in% promoterGenes$gene_name)
  {
    deltaMethy.median <- mean(promoterGenes$Median.delta[promoterGenes$gene_name == gene])
    cor <- mean(promoterGenes$rho[promoterGenes$gene_name == gene])
  }
  else if(gene %in% genicEnhancerGenes$gene_name)
  {
    deltaMethy.median <- mean(genicEnhancerGenes$Median.delta[genicEnhancerGenes$gene_name == gene])
    cor <- mean(genicEnhancerGenes$rho[genicEnhancerGenes$gene_name == gene])
  } else {
    deltaMethy.median <- mean(intergenicEnhancerGenes$DM.Median[intergenicEnhancerGenes$gene_name == gene])
    cor <- mean(intergenicEnhancerGenes$MethyandGene.rho[intergenicEnhancerGenes$gene_name == gene])
  }
  c(Median.LFC, deltaMethy.median, cor)
}
pcr_info <- as.data.frame(pcr_info)
pcr_info <- cbind(pcr_genes, pcr_info)
colnames(pcr_info)[3:5] <- c("LFC.median", "DelatMethy.Median", "Rho")
write.csv(pcr_info, file="GenesforRTPCRValidation.csv", row.names = F)
