############### BWS test based differential methylation loci detection ##################
setwd("/data/huangp/Methy/Mixed/New")
options(stringsAsFactors = F)
library(fastSave)
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores =64)
#load.lbzip2("Delta.Methylation.Matrix.RDataFS", n.cores=64)
load("../All66HCC.phenotypes.RData")
Meth.all.qc <- Meth.all.qc[HCC.phenotypes$SampleID,]

#load.lbzip2("Paired_Multivariable_linear_regression_adjusted_for_Age_summary.RDataFS", n.cores=32)
load.lbzip2("SigDML_from_Paired_Multivariable_linear_regression_adjusted_for_Age_summary.RDataFS", n.cores=32)

# CpGs with median of delta methylation >= 15%
#obatain median of delta methylation of each CpG
library(foreach)
library(doParallel)
registerDoParallel(cl=3, cores = 16)
system.time(median.delta.methy <- foreach(i = 1:ncol(median.delta.methy), .combine = 'c')%dopar%
  {
    median(delta.methy[,i])
  })
stopImplicitCluster()
names(median.delta.methy) <- colnames(delta.methy)

#median.delta.methy <- apply(delta.methy, 2, median)
table(median.delta.methy >= 0.15)
table(median.delta.methy <= -0.15)
Candidate_DML <- delta.methy[,abs(median.delta.methy)>=0.15]
table(colnames(Candidate_DML) %in% rownames(SigDML))

HCC.phenotypes$age_group <- rep(NA, nrow(HCC.phenotypes))
HCC.phenotypes$age_group[HCC.phenotypes$age<55] <- "young"
HCC.phenotypes$age_group[HCC.phenotypes$age>=55 & HCC.phenotypes$age<65] <- "medium"
HCC.phenotypes$age_group[HCC.phenotypes$age>=65] <- "old"

phenotypes <- HCC.phenotypes[,c("SampleID", "patientID", "Condition", "age", "age_group")]
phenotypes$Group <- paste(phenotypes$age_group, phenotypes$Condition, sep="|")

library(BWStest)
CDF_bws <- function(B)
{
  if(B<0) {
    nB <- -1*B
    if(nB < 1.5) P = 1-exp(-0.699-1.255*nB)
    else if(nB<9) P = 1-exp(-0.895-1.153*nB+0.0173*nB*nB)
    else P = 1-exp(-2.895-0.786*nB)
  }
  else if(B<1.5) P = exp(-0.699-1.255*B)
  else if(B<9) P = exp(-0.895-1.153*B+0.0173*B*B)
  else P = exp(-2.895-0.786*B)
  
  return(1-P)
}
obtain_T <- function(methy, onlyTvalue = T)
{
  bws_stat_less_young <- bws_test(x = methy[names(methy) == "young|T"], y = methy[names(methy) == "young|N"],
                                  method = "Neuhauser", alternative = "less")[[1]]
  bws_stat_greater_young <- bws_test(x = methy[names(methy) == "young|T"], y = methy[names(methy) == "young|N"],
                                     method = "Neuhauser", alternative = "greater")[[1]]
  bws_stat_less_medium <- bws_test(x = methy[names(methy) == "medium|T"], y = methy[names(methy) == "medium|N"],
                                   method = "Neuhauser", alternative = "less")[[1]]
  bws_stat_greater_medium <- bws_test(x = methy[names(methy) == "medium|T"], y = methy[names(methy) == "medium|N"],
                                      method = "Neuhauser", alternative = "greater")[[1]]
  bws_stat_less_old <- bws_test(x = methy[names(methy) == "old|T"], y = methy[names(methy) == "old|N"],
                                method = "Neuhauser", alternative = "less")[[1]]
  bws_stat_greater_old <- bws_test(x = methy[names(methy) == "old|T"], y = methy[names(methy) == "old|N"],
                                   method = "Neuhauser", alternative = "greater")[[1]]
  bws_logP_less_young <- log(CDF_bws(bws_stat_less_young), 10)
  bws_logP_greater_young <- log(CDF_bws(bws_stat_greater_young), 10)
  bws_logP_less_medium <- log(CDF_bws(bws_stat_less_medium), 10)
  bws_logP_greater_medium <- log(CDF_bws(bws_stat_greater_medium), 10)
  bws_logP_less_old <- log(CDF_bws(bws_stat_less_old), 10)
  bws_logP_greater_old <- log(CDF_bws(bws_stat_greater_old), 10)
  
  T_less <- -2*sum(c(bws_logP_less_young, bws_logP_less_medium, bws_logP_less_old))
  T_greater <- -2*sum(c(bws_logP_greater_young, bws_logP_greater_medium, bws_logP_greater_old))
  T_statistic <- max(T_less, T_greater) 
  sixlogP <- c(-1*bws_logP_less_young, -1*bws_logP_greater_young, 
                  -1*bws_logP_less_medium, -1*bws_logP_greater_medium,
                  -1*bws_logP_less_old, -1*bws_logP_greater_old)
  if(onlyTvalue == T) return(T_statistic)
  else return(c(sixlogP, T_statistic))
}


#calculate combined BWS pvalue
#randomly generated 10^8 permutation T 
setwd("/data/huangp/HCC/Methy/Mixed/New")#51
options(stringsAsFactors = F)
library(fastSave)
load.lbzip2("./NewSmoothed.Methylation.qc.RDataFS", n.cores =64)
load("./All66HCC.phenotypes.RData")
phenotypes <- HCC.phenotypes[,c("SampleID", "patientID", "Condition", "age", "age_group")]
phenotypes$Group <- paste(phenotypes$age_group, phenotypes$Condition, sep="|")

library(foreach)
library(doParallel)
library(BWStest)

system.time(
Tvalue_ECDF <- foreach(i = 1:300, .combine = 'c')%do%
  {
    DML <- sample(x=colnames(Meth.all.qc), size = 1)
    methy <- Meth.all.qc[,DML]
    names(methy) <- phenotypes$Group
    #Tvalues <- obtain_T(methy, onlyTvalue = F)
    ## permutation to obtain empirical distribution of T for each CpG
    registerDoParallel(cl = 24, cores = 4)
    permutation_T <- foreach(i = 1:100000, .combine = 'c')%dopar%
                  {
                    permutated_group <- sample(x = names(methy), size = length(methy), replace = F)
                    permutated_methy <- methy
                    names(permutated_methy) <- permutated_group
                    obtain_T(permutated_methy, onlyTvalue = T)
                  }
    stopImplicitCluster()
    message(paste("Finished the ", i, "th random CpG", sep=""))
    permutation_T
  }
)
save.lbzip2(Tvalue_ECDF, file="300RandomCpGx10wPermutation_1.RDataFS", n.cores=24)



### comparison of different ECDF of T statistics
# for each sigDML, calculate their pvalue according to those 300 ECDF 
# obtained from 10^5 permutation of those first 300 sigDMLs
setwd("/data/huangp/HCC/Methy/Mixed/New")#51 serve
options(stringsAsFactors = F)
library(fastSave)
load.lbzip2("BWS_pvalue.RDataFS", n.cores=32)
BWS_pvalue$logP <- -1*log10(BWS_pvalue$BWS.Pvalue+0.00001)
cor.test(BWS_pvalue$Tvalue, BWS_pvalue$logP, method="spearman")#pvalue<2.2e-16 rho = 0.71
library(fastSave)
library(doParallel)
system.time(BWSPvalue_matrix <- foreach(i = 1:300, .combine = 'rbind')%do%
  {
    load.lbzip2(paste("./ECDF/ECDF_SigDML_", i, ".RDataFS", sep=""))
    1-ECDF(BWS_pvalue$Tvalue)
  })
str(BWSPvalue_matrix)
colnames(BWSPvalue_matrix) <- BWS_pvalue$pos
rownames(BWSPvalue_matrix) <- paste("ECDF", BWS_pvalue$pos, sep="|")
BWSlogP_matrix <- -1*apply((BWSPvalue_matrix+0.00001), 2, log10)

Means <- colMeans(BWSlogP_matrix)
Medians <- apply(BWSlogP_matrix, 2, median)
cor.test(Means, Medians) # cor = 0.99977871
cor.test(Means, BWS_pvalue$logP) # cor =0.9711889

### 


#### calcluate the pvalue and ECDF of those candidate DMRs in Table 2
#ongoing in 51 server
#result: all CpGs were significant in combined BWS test based DML analyses
#### examine the overlap between sigDML and important DMRs in Table 2
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")

DMRs <- read.csv("Important_DMRs_for_comparison.csv", header=T)
DMRs$DMR.chr <- sapply(strsplit(DMRs$DMR_ID, split="_", fixed=T), "[", 1)
DMRs$DMR.start <- as.integer(sapply(strsplit(DMRs$DMR_ID, split="_", fixed=T), "[", 2))
DMRs$DMR.end <- as.integer(sapply(strsplit(DMRs$DMR_ID, split="_", fixed=T), "[", 3))
# choose the start and stop CpGs of each DMR for comparison
start_CpGs <- cbind(chr = DMRs$DMR.chr, pos = DMRs$DMR.start, ID = paste(DMRs$DMR.chr, DMRs$DMR.start, sep="_"))
end_CpGs <- cbind(chr = DMRs$DMR.chr, pos =DMRs$DMR.end, ID = paste(DMRs$DMR.chr, DMRs$DMR.end, sep="_"))
all_CpGs <- as.data.frame(rbind(start_CpGs, end_CpGs))
all_CpGs$pos <- as.integer(all_CpGs$pos)
all_CpGs <- all_CpGs[order(all_CpGs$chr, all_CpGs$pos, decreasing = F),]
#save(all_CpGs, file="Important_CpGs_for_comparison.RData")
load("Important_CpGs_for_comparison.RData")
load.lbzip2("SigDML_from_Paired_Multivariable_linear_regression_adjusted_for_Age_summary.RDataFS", n.cores=32)
load.lbzip2("Paired_Multivariable_linear_regression_adjusted_for_Age_summary.RDataFS", n.cores=64)
table(all_CpGs$ID %in% rownames(SigDML))#FASLE 52 TRUE 6
table(SigDML$shapiro.p < 0.05)#TRUE 957,497 36.7%
table(LM.summary$shapiro.p < 0.05)#TRUE 11,958,240 40%
all_CpGs_LM.summary <- LM.summary[all_CpGs$ID,]
#There were indeed some CpGs show quite huge difference between tumor and adjacent but become insignificant in new paired-MLR
#Thus make a slope chart (?¶?ͼ) to determine whether these DMLs are false positive
#load methylation all
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores =64)
load("../All66HCC.phenotypes.RData")
slope_chart <- function(DML)
{
  #DML = "chr1_58576756"
  library(ggplot2)
  library(RColorBrewer)
  methy <- Meth.all.qc[,DML]
  df <- data.frame(patientID = HCC.phenotypes$patientID[1:33], TumorMethy = methy[1:33], AdjacentMethy = methy[34:66])
  df$delta <- df$TumorMethy - df$AdjacentMethy
  df$class <- ifelse((df$TumorMethy - df$AdjacentMethy) > 0, "red", "green")
  left_label <- paste(df$patientID, df$AdjacentMethy, sep=",")
  right_label <- paste(df$patientID, df$TumorMethy, sep=",")
  p <- ggplot(df)+
    geom_segment(aes(x=1,xend=2,y=AdjacentMethy, yend =TumorMethy, col=class),size=0.75,show.legend=F)+ #��????
    geom_vline(xintercept=1, linetype="solid", size=1)+#Adjacent?Ĵ?ֱ??
    geom_vline(xintercept=2, linetype="solid", size=1)+
    geom_point(aes(x=1,y=AdjacentMethy),size=3, shape=21, fill = "grey80", color="black")+
    geom_point(aes(x=2,y=TumorMethy),size=3, shape=21, fill = "grey80", color="black")+
    scale_color_manual(labels=c("Up", "Down"), values=c("green"="#A6D854", "red" = "#FC4E07"))
  ggsave(p, file=paste("slope_plot_", DML, ".pdf", sep=""))
}
slope_chart(DML = "chr1_58575901")
slope_chart(DML="chr1_161224143")
slope_chart(DML="chr1_161224295")
slope_chart(DML="chr10_10793485")
slope_chart(DML="chr10_10793810")

#bws_test(x, y, method = c("default", "BWS", "Neuhauser", "B1", "B2", "B3",
#                          "B4", "B5"), alternative = c("two.sided", "greater", "less"))

############################# BWS test using 300x10w for all CpGs
# to see the power and potential false positive of combined BWS test based DML identification
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")#29
load.lbzip2("150RandomCpGx10wPermutation_1.RDataFS", n.cores=24)
Tvalue_1 <- Tvalue_ECDF
load.lbzip2("300RandomCpGx10wPermutation_2.RDataFS", n.cores=24)
Tvalue_2 <- Tvalue_ECDF
load.lbzip2("300RandomCpGx10wPermutation_3.RDataFS", n.cores=24)
Tvalue_3 <- Tvalue_ECDF
load.lbzip2("300RandomCpGx10wPermutation_29_1.RDataFS", n.cores=24)
Tvalue_4 <- Tvalue_ECDF
load.lbzip2("300RandomCpGx10wPermutation_29_2.RDataFS", n.cores=24)
Tvalue_5 <- Tvalue_ECDF
load.lbzip2("300RandomCpGx10wPermutation_29_3.RDataFS", n.cores=24)
Tvalue_6 <- Tvalue_ECDF

Tvalue_total <- c(Tvalue_1,Tvalue_2,Tvalue_3, Tvalue_4, Tvalue_5, Tvalue_6)#150x10w + 5x300x10w = 1.65x10^8
ECDF_T <- ecdf(Tvalue_total)
save.lbzip2(Tvalue_total, file="1.65x10e8_permutation_T_statistic.RDataFS", n.cores=24)
save.lbzip2(ECDF_T, file="Function_1.65x10e8_ECDF_of_combined_BWS_T_statistic.RDataFS", n.cores=24)

# Combined BWS test based DML detection
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")#29

library(foreach)
library(doParallel)
library(BWStest)
##load methy
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores =64)
#load phenotypes
load("../All66HCC.phenotypes.RData")

#load function CDF_bws, obtain_T, ECDF_T
load("Function_CDF_bws.RData")
load("Function_obtain_T.RData")
load.lbzip2("Function_1.65x10e8_ECDF_of_combined_BWS_T_statistic.RDataFS", n.cores=24)

chunkSize = 100000
chunkNum = ceiling(ncol(Meth.all.qc)/chunkSize) -1

system.time(BWS_DML_summary <- foreach(j = 1:chunNum, .combine = 'rbind')%do%
  {
    registerDoParallel(cl=24, cores=8)
    system.time(BWS_DML_summary_chunk <- foreach(i = ((j-1)*chunkSize+1):(j*chunkSize), .combine = 'rbind')%dopar%
                  {
                    # obtain six one-side BWS pvalue and combined T value
                    Methy <- Meth.all.qc[,i]
                    names(Methy) <- HCC.phenotypes$Group
                    Tvalues <- obtain_T(methy = Methy, onlyTvalue=F)
                    # obtain combined BWS Pvalue according to ECDF of T statistics
                    combinedPvalue <- 1-ECDF_T(Tvalues[7])
                    c(Tvalues, combinedPvalue)
                  })
    
    message(paste("Finished the ", j, "th chunk", sep=""))
    BWS_DML_summary_chunk
  })
BWS_DML_summary_chunk <- foreach(i = (chunkNum*chunkSize+1):ncol(Meth.all.qc), .combine = 'rbind')%dopar%
  {
    # obtain six one-side BWS pvalue and combined T value
    Methy <- Meth.all.qc[,i]
    names(Methy) <- HCC.phenotypes$Group
    Tvalues <- obtain_T(methy = Methy, onlyTvalue=F)
    # obtain combined BWS Pvalue according to ECDF of T statistics
    combinedPvalue <- 1-ECDF_T(Tvalues[7])
    c(Tvalues, combinedPvalue)
  }
BWS_DML_summary <- rbind(BWS_DML_summary, BWS_DML_summary_chunk)
stopImplicitCluster()

BWS_DML_summary <- apply(BWS_DML_summary, 2, as.numeric)
colnames(BWS_DML_summary) <- c("-logP.hypo.young", "-logP.hyper.young", 
                               "-logP.hypo.medium", "-logP.hyper.medium", 
                               "-logP.hypo.old", "-logP.hyper.old", 
                               "Tvalue", "combined.BWS.Pvalue")
load.lbzip2("delta.methy.median.RDataFS", n.cores=32)
BWS_DML_summary <- data.frame(BWS_DML_summary, Median.delta = median_deltaMethy)
row.names(BWS_DML_summary) <- colnames(Meth.all.qc)
BWS_DML_summary <- data.frame(chr = sapply(strsplit(rownames(BWS_DML_summary), split="_", fixed = T), "[", 1),
                              pos = as.integer(sapply(strsplit(rownames(BWS_DML_summary), split="_", fixed = T), "[", 2)),
                              BWS_DML_summary)
save.lbzip2(BWS_DML_summary, file="Combined_BWS_based_DML_detection_summary.RDataFS", n.cores=32)
table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 0.05, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 0.01, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 0.001, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 0.0001, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 0.00001, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 1E-6, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 1E-7, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

table(Significance = BWS_DML_summary$combined.BWS.Pvalue < 1E-8, 
      Median = abs(BWS_DML_summary$Median.delta) >= 0.15)

SigDMLs <- BWS_DML_summary[BWS_DML_summary$combined.BWS.Pvalue < 1E-5,]
SigDMLs <- data.frame(chr = sapply(strsplit(rownames(SigDMLs), split="_", fixed = T), "[", 1),
                      pos = as.integer(sapply(strsplit(rownames(SigDMLs), split="_", fixed = T), "[", 2)),
                      SigDMLs)
save.lbzip2(SigDMLs, file="SigDMLs_from_combined_BWS_based_DML_detection_summary.RDataFS", n.cores=32)
table(SigDMLs$Median.delta < 0) #FALSE 157,341 TRUE 9,710,359
#load previous DML_Tumor
load.lbzip2("../Multivariable_linear_regression_withDNAmAge_summary_for_unQC_qc_CpG.RDataFS", n.cores=64)
SigDMLs_LM <- LM.summary[abs(LM.summary$Estimate.Tumor) >= 0.15 & LM.summary$pvalue.Tumor < 0.05,]
table(rownames(SigDMLs_LM) %in% rownames(SigDMLs)) # TRUE 2960052/3404029


## as a comparison, perform age-unadjusted BWS test based DML detection 
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")#29

library(foreach)
library(doParallel)
library(BWStest)
##load methy
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores =64)
#load phenotypes
load("../All66HCC.phenotypes.RData")

chunkSize = 100000
chunkNum = ceiling(ncol(Meth.all.qc)/chunkSize) -1

system.time(BWS_DML_pvalue <- foreach(j = 1:2, .combine = 'c')%do%
              {
                registerDoParallel(cl=48, cores=8)
                system.time(BWS_DML_pvalue_chunk <- foreach(i = ((j-1)*chunkSize+1):(j*chunkSize), .combine = 'c')%dopar%
                              {
                                Methy <- Meth.all.qc[,i]
                                bws_test(x = Methy[1:33], y = Methy[34:66], 
                                              method="default", alternative = "two.sided")[[2]]
                              })
                
                message(paste("Finished the ", j, "th chunk", sep=""))
                BWS_DML_pvalue_chunk
              })
BWS_DML_pvalue_chunk <- foreach(i = (chunkNum*chunkSize+1):ncol(Meth.all.qc), .combine = 'c')%dopar%
  {
    Methy <- Meth.all.qc[,i]
    bws_test(x = Methy[1:33], y = Methy[34:36], 
             method="default", alternative = "two.sided")[[2]]
  }
BWS_DML_pvalue <- c(BWS_DML_pvalue, BWS_DML_pvalue_chunk)
stopImplicitCluster()
names(BWS_DML_pvalue) <- colnames(Meth.all.qc)
save.lbzip2(BWS_DML_pvalue, file="Age_unadjusted_BWS_based_DML_detection_pvals.RDataFS", n.cores=32)

############### DMR detection #####################################################
DMRCalling <- function(InputDML, allCpGs = BWS_DML_summary, Distance=200, DML_count=5, DML_ratio=0.8)
{
  ###load LM.summary_all
  #library(fastSave)
  #load.lbzip2("/data/huangp/HCC/Methy/Mixed/MLR.summary_all.RDataFS", n.cores=64)
  #message("Finished loading All CpG")
  all_DMR <- list()
  #all_pre_DMR <- list()
  for(chr in 1:length(table(InputDML$chr)))
  {
    ##DML_HBV_Both in chr1 as a example
    CHR = names(table(InputDML$chr))[chr]
    DML_example <- InputDML[InputDML$chr==CHR,]
    DML_example$pos <- as.integer(DML_example$pos)
    if(nrow(DML_example)!=0)
    {
      DML_example <- DML_example[order(DML_example$pos, decreasing = F),]
      distance <- c(0, diff(DML_example$pos, lag = 1))
      names(distance) <- c(1:nrow(DML_example))
      ###define pre_DMR
      intersect <- as.integer(names(distance)[distance>Distance])
      pre_DMR <- vector(mode="list", length=length(intersect)+1)
      for(j in 1:length(pre_DMR))
      {
        if(j ==1) {pre_DMR[[j]] <- DML_example[1:(intersect[j]-1),]}
        else  if(j <= length(intersect)) 
          pre_DMR[[j]] <-  DML_example[intersect[j-1]:(intersect[j]-1),]
        else pre_DMR[[j]] <- DML_example[intersect[j-1]:nrow(DML_example), ]
      }
      ###filtering pre_DMRs with DML less than DML_count
      pre_DMR_count <- unlist(lapply(pre_DMR, nrow))
      pre_DMR <- pre_DMR[pre_DMR_count>=DML_count]
      pre_DMR_count <- pre_DMR_count[pre_DMR_count>=DML_count]
      
      if(length(pre_DMR)>0)
      {
        ##get the range of pre_DMRs
        for(j in 1:length(pre_DMR))
        {
          names(pre_DMR)[j] <- paste(CHR, pre_DMR[[j]]$pos[1],  pre_DMR[[j]]$pos[nrow(pre_DMR[[j]])], sep="_")
        }
        ###expand the pre_DMR with non-significant CpGs
        exp_DMR <- vector(mode="list", length=length(pre_DMR))
        intra_chr_CpG <- allCpGs[allCpGs$chr == CHR,]
        intra_chr_CpG <- intra_chr_CpG[order(intra_chr_CpG$pos, decreasing = F),]
        
        for(j in 1:length(exp_DMR))
        {
          str <- as.integer(strsplit(names(pre_DMR)[j], "_")[[1]][2])
          end <- as.integer(strsplit(names(pre_DMR)[j], "_")[[1]][3])
          str_ix <- match(str, intra_chr_CpG$pos)
          end_ix <- match(end, intra_chr_CpG$pos)
          exp_DMR[[j]] <- intra_chr_CpG[str_ix:end_ix,]
        }
        ##check the CpG count of each exP_DMR
        exp_DMR_count <- unlist(lapply(exp_DMR, nrow))
        ##name the exp_DMR with its location
        for(i in 1:length(exp_DMR))
        {
          names(exp_DMR)[i] <- paste(exp_DMR[[i]]$chr[1], exp_DMR[[i]]$pos[1], exp_DMR[[i]]$pos[nrow(exp_DMR[[i]])], sep="_")
        }
        ##filtering the exp_DMR with sig_DML ratio > DML_ratio
        exp_DMR <- exp_DMR[pre_DMR_count/exp_DMR_count>=DML_ratio]
        all_DMR<- c(all_DMR, exp_DMR)
      }
    }
    message(paste("Finished the DMR calling of", chr, "of", length(table(InputDML$chr)), sep=" "))
  }
  return(all_DMR)
}
save(DMRCalling, file="Function_DMRCalling.RData")

DMR_summary <- function(DMR)
{
  chr <- sapply(strsplit(names(DMR), "_", fixed=T), "[", 1)
  str <- as.integer(sapply(strsplit(names(DMR), "_", fixed=T), "[", 2))
  end <- as.integer(sapply(strsplit(names(DMR), "_", fixed=T), "[", 3))
  nCG <- unlist(lapply(DMR, nrow))
  
  logP.hypo.young <- unlist(lapply(sapply(DMR, "[", 3), mean))
  logP.hyper.young <- unlist(lapply(sapply(DMR, "[", 4), mean))
  logP.hypo.medium <- unlist(lapply(sapply(DMR, "[", 5), mean))
  logP.hyper.medium <- unlist(lapply(sapply(DMR, "[", 6), mean))
  logP.hypo.old <- unlist(lapply(sapply(DMR, "[", 7), mean))
  logP.hyper.old <- unlist(lapply(sapply(DMR, "[", 8), mean))
  Tvalue <- unlist(lapply(sapply(DMR, "[", 9), mean))
  Median.delta <- unlist(lapply(sapply(DMR, "[", 11), mean))
  DMR.summary <- data.frame(chr, str, end, width = end-str+1, nCG, 
                            logP.hypo.young, logP.hyper.young, 
                            logP.hypo.medium, logP.hyper.medium,
                            logP.hypo.old, logP.hyper.old, 
                            Tvalue, Median.delta)
  return(DMR.summary)
}

options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New")#29

#load SigDML, BWS_DML_summary
load.lbzip2("SigDMLs_from_combined_BWS_based_DML_detection_summary.RDataFS", n.cores=24)
load.lbzip2("Combined_BWS_based_DML_detection_summary.RDataFS", n.cores=24)

system.time(AllDMRs <- DMRCalling(InputDML = SigDMLs, allCpGs = BWS_DML_summary))
save.lbzip2(AllDMRs, file="All_pre_DMRslist.RDataFS", n.cores=24)
#summarizing DMRs
AllDMR.summary <- DMR_summary(AllDMRs)
load.lbzip2("Function_1.65x10e8_ECDF_of_combined_BWS_T_statistic.RDataFS", n.cores=24)
AllDMR.summary$Pvalue <- 1-ECDF_T(AllDMR.summary$Tvalue)
save.lbzip2(AllDMR.summary, file="All_pre_DMRs.summary.RDataFS", n.cores=24)
# SigDMRs
table(Pvalue = AllDMR.summary$Pvalue < 1E-5, Median = abs(AllDMR.summary$Median.delta)>=0.15)
#           Median
#Pvalue   FALSE   TRUE
#  FALSE      1      1
#  TRUE   62718 545559
table(Pvalue = AllDMR.summary$Pvalue < 1E-6, Median = abs(AllDMR.summary$Median.delta)>=0.20)
#          Median
#Pvalue   FALSE   TRUE
#  FALSE  55756  46605
#  TRUE   71933 433985
table(Pvalue = AllDMR.summary$Pvalue < 1E-7, Median = abs(AllDMR.summary$Median.delta)>=0.25)
#         Median
#Pvalue   FALSE   TRUE
#  FALSE 158532  86549
#  TRUE   61832 301366

SigDMR.summary <- AllDMR.summary[AllDMR.summary$Pvalue < 1E-7 & abs(AllDMR.summary$Median.delta) >= 0.25,]
save.lbzip2(SigDMR.summary, file="Significant_DMRs.summary.RDataFS", n.cores=24)
#obtain smoothed and raw group level and perBase level methylation level for each SigDMR
#smoothed methylation level
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores=24)
getDMRMethy_Smooth <- function(DMR, methy.sm=Meth.all.qc, what=c("perBase", "perRegion"))
{
  chr = DMR[1]
  start = as.integer(DMR[2])
  end = as.integer(DMR[3])
  possible_CpGs <- paste(chr, start:end, sep="_")
  avail_CpGs <- possible_CpGs[possible_CpGs %in% colnames(Meth.all.qc)]
  methy.matrix <- Meth.all.qc[,avail_CpGs]
  if(what == "perBase")
    return(methy.matrix)
  else return(rowMeans(methy.matrix))
}

library(foreach)
library(doParallel)
#get group level MethyMatrix
registerDoParallel(cl = 24, cores = 4)
system.time(SigDMR.methy.matrix <- foreach(i = 1:nrow(SigDMR.summary), .combine = 'cbind')%dopar%
  {
    getDMRMethy_Smooth(DMR=SigDMR.summary[i,1:3], methy.sm = Meth.all.qc, what = "perRegion")
  })
stopImplicitCluster()
colnames(SigDMR.methy.matrix) <- rownames(SigDMR.summary)
rownames(SigDMR.methy.matrix) <- rownames(Meth.all.qc)
save.lbzip2(SigDMR.methy.matrix, file="SigDMR_NewSmoothed_MethyMatrix.RDataFS", n.cores=24)
#get perbase level MethyMatrixList
registerDoParallel(cl =24, cores = 4)
SigDMR.methy.matrix.list <- foreach(i = 1:nrow(SigDMR.summary), .combine = 'list')%dopar%
  {
    getDMRMethy_Smooth(DMR=SigDMR.summary[i,1:3], methy.sm = Meth.all.qc, what = "perBase")
  }
stopImplicitCluster()
names(SigDMR.methy.matrix.list) <- rownames(SigDMR.summary)
save.lbzip2(SigDMR.methy.matrix.list, file="SigDMR_NewSmoothed_MethyMatrixList.RDataFS", n.cores=24)

#raw methylation level 
#latter

############### DMR annotation ###########################################
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/HCC/Methy/Mixed/New")#51
#load SigDMR
load.lbzip2("preDMR.summary.RDataFS", n.cores=24)
DMR.Tumor.summary <- SigDMR.summary
### Granges file buliding
library(genomation)
library(GenomicFeatures)
Granges.DMR.Tumor <- makeGRangesFromDataFrame(DMR.Tumor.summary, start.field = 'str', keep.extra.columns = T)
############################################################## part 1 DMR physical annotation ############################
load("Granges.CpGIsland.hg38.RData")
load("Granges.CpGShelf.hg38.RData")
load("Granges.CpGShore.hg38.RData")
load("Granges.Promoter.hg38.RData")
load("Granges.Exon.hg38.gencode.v31.RData")
load("Granges.gene.hg38.gencode.v31.RData")
load("Granges.3UTR.hg38.gencode.v31.RData")
load("Granges.5UTR.hg38.gencode.v31.RData")
load("Granges.FANTOM5.enhancer.hg38.RData")
### CpG island
DMR_island_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.island))
DMR_island <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_island_hit$queryHits,1:4],
                    as.data.frame(Granges.island)[DMR_island_hit$subjectHits,c(2:3,6)])
colnames(DMR_island) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "island_str", "island_end", "island_pos_hg19")
DMR_island$overlap_str <- apply(DMR_island[,c("DMR_str", "island_str")], 1, max)
DMR_island$overlap_end <- apply(DMR_island[,c("DMR_end", "island_end")], 1, min)
DMR_island$ratio.DMR <- (DMR_island$overlap_end - DMR_island$overlap_str +1)/(DMR_island$DMR_width +1)
DMR_island$ratio.island <- (DMR_island$overlap_end - DMR_island$overlap_str +1)/(DMR_island$island_end - DMR_island$island_str +1)
save(DMR_island, file="preDMR_CpG_island_overlap.RData")
### summarize 
X <- paste(DMR_island$chr, DMR_island$DMR_str, DMR_island$DMR_end, sep="_")
DMR.Tumor.summary$CpGIsland <- rownames(DMR.Tumor.summary) %in% X

### CpG shore
DMR_shore_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.shore))
DMR_shore <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_shore_hit$queryHits,1:4],
                   as.data.frame(Granges.shore)[DMR_shore_hit$subjectHits,c(2:3,6)])
colnames(DMR_shore) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "shore_str", "shore_end", "shore_pos_hg19")
DMR_shore$overlap_str <- apply(DMR_shore[,c("DMR_str", "shore_str")], 1, max)
DMR_shore$overlap_end <- apply(DMR_shore[,c("DMR_end", "shore_end")], 1, min)
DMR_shore$ratio.DMR <- (DMR_shore$overlap_end - DMR_shore$overlap_str +1)/(DMR_shore$DMR_width +1)
DMR_shore$ratio.shore <- (DMR_shore$overlap_end - DMR_shore$overlap_str +1)/(DMR_shore$shore_end - DMR_shore$shore_str +1)
save(DMR_shore, file="preDMR_CpG_shore_overlap.RData")
### summarize 
X <- paste(DMR_shore$chr, DMR_shore$DMR_str, DMR_shore$DMR_end, sep="_")
DMR.Tumor.summary$CpGShore <- rownames(DMR.Tumor.summary) %in% X

### CpG shelf
DMR_shelf_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.shelf))
DMR_shelf <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_shelf_hit$queryHits,1:4],
                   as.data.frame(Granges.shelf)[DMR_shelf_hit$subjectHits,c(2:3,6)])
colnames(DMR_shelf) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "shelf_str", "shelf_end", "shelf_pos_hg19")
DMR_shelf$overlap_str <- apply(DMR_shelf[,c("DMR_str", "shelf_str")], 1, max)
DMR_shelf$overlap_end <- apply(DMR_shelf[,c("DMR_end", "shelf_end")], 1, min)
DMR_shelf$ratio.DMR <- (DMR_shelf$overlap_end - DMR_shelf$overlap_str +1)/(DMR_shelf$DMR_width +1)
DMR_shelf$ratio.shelf <- (DMR_shelf$overlap_end - DMR_shelf$overlap_str +1)/(DMR_shelf$shelf_end - DMR_shelf$shelf_str +1)
save(DMR_shelf, file="preDMR_CpG_shelf_overlap.RData")
### summarize 
X <- paste(DMR_shelf$chr, DMR_shelf$DMR_str, DMR_shelf$DMR_end, sep="_")
DMR.Tumor.summary$CpGShelf <- rownames(DMR.Tumor.summary) %in% X
## part 2  functional annotation
### Intergenic
DMR_gene_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.gene))
DMR_gene <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_gene_hit$queryHits,1:4],
                  as.data.frame(Granges.gene)[DMR_gene_hit$subjectHits,c(2,3,5,10,11,12)])
colnames(DMR_gene) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "gene_str", "gene_end", "gene.strand", "gene_id", "gene_type", "gene_name")

#save(DMR_gene, file="preDMR_gene_overlap.RData")
### summarize
X <- paste(DMR_gene$chr, DMR_gene$DMR_str, DMR_gene$DMR_end, sep="_")
DMR.Tumor.summary$Intergenic <- !(rownames(DMR.Tumor.summary) %in% X)
### colapse function
COMBINE <- function(x, sep)
{
  x <- as.character(x)
  if(length(x)!=0)
  {
    y <- x[1]
    if(length(x)>1)
    {
      for(i in 2:length(x)) y <- paste(y, x[i], sep = sep)
    }
  }
  return(y)
}

Collapse <- function(groupVar,collapseVar)
{
  library(doParallel)
  library(foreach)
  registerDoParallel(8)
  if(length(groupVar)==length(collapseVar))
  {
    group <- unique(groupVar)
    collapsed.var <- vector(mode="character", length=length(group))
    collapsed.var <- foreach(i = 1:length(group), .combine = 'c') %dopar%{
      COMBINE(collapseVar[groupVar %in% group[i]], sep="|")
    }
    stopImplicitCluster()
    return(as.data.frame(cbind(group, collapsed.var)))
  }
}
#DMRID <- paste(DMR_gene$chr, DMR_gene$DMR_str, DMR_gene$DMR_end, sep="_")
#system.time(collapse.gene <- Collapse(groupVar= DMRID, collapseVar = DMR_gene$gene_name))
#colnames(collapse.gene) <- c("DMRID", "Gene")
#save(collapse.gene, file="preDMR_gene_collapsed_overlap.RData")
#DMR.Tumor.summary$Gene <- collapse.gene$Gene[match(rownames(DMR.Tumor.summary), collapse.gene$DMRID)]
### exon
DMR_exon_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.exon))
DMR_exon <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_exon_hit$queryHits,1:4],
                  as.data.frame(Granges.exon)[DMR_exon_hit$subjectHits,2:3])
colnames(DMR_exon) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "exon_str", "exon_end")
DMR_exon$overlap_str <- apply(DMR_exon[,c("DMR_str", "exon_str")], 1, max)
DMR_exon$overlap_end <- apply(DMR_exon[,c("DMR_end", "exon_end")], 1, min)
DMR_exon$ratio.DMR <- (DMR_exon$overlap_end - DMR_exon$overlap_str +1)/(DMR_exon$DMR_width +1)
DMR_exon$ratio.exon <- (DMR_exon$overlap_end - DMR_exon$overlap_str +1)/(DMR_exon$exon_end - DMR_exon$exon_str +1)
save(DMR_exon, file="preDMR_exon_overlap.RData")
#### summarize
X <- paste(DMR_exon$chr, DMR_exon$DMR_str, DMR_exon$DMR_end, sep="_")
DMR.Tumor.summary$Exon <- rownames(DMR.Tumor.summary) %in% X
### promoter
DMR_promoter_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.promoter))
DMR_promoter <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_promoter_hit$queryHits,1:4],
                      as.data.frame(Granges.promoter)[DMR_promoter_hit$subjectHits,c(2,3,6)])
colnames(DMR_promoter) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "promoter_str", "promoter_end", "transcript_name")
DMR_promoter$overlap_str <- apply(DMR_promoter[,c("DMR_str", "promoter_str")], 1, max)
DMR_promoter$overlap_end <- apply(DMR_promoter[,c("DMR_end", "promoter_end")], 1, min)
DMR_promoter$ratio.DMR <- (DMR_promoter$overlap_end - DMR_promoter$overlap_str +1)/(DMR_promoter$DMR_width +1)
DMR_promoter$ratio.promoter <- (DMR_promoter$overlap_end - DMR_promoter$overlap_str +1)/(DMR_promoter$promoter_end - DMR_promoter$promoter_str +1)
save(DMR_promoter, file="preDMR_promoter_overlap.RData")
#### summarize
X <- paste(DMR_promoter$chr, DMR_promoter$DMR_str, DMR_promoter$DMR_end, sep="_")
DMR.Tumor.summary$Promoter <- rownames(DMR.Tumor.summary) %in% X
### 5'UTR
DMR_utr5_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.utr5))
DMR_utr5 <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_utr5_hit$queryHits,1:4],
                  as.data.frame(Granges.utr5)[DMR_utr5_hit$subjectHits,c(2,3,7)])
colnames(DMR_utr5) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "utr5_str", "utr5_end", "transcript_name")
DMR_utr5$overlap_str <- apply(DMR_utr5[,c("DMR_str", "utr5_str")], 1, max)
DMR_utr5$overlap_end <- apply(DMR_utr5[,c("DMR_end", "utr5_end")], 1, min)
DMR_utr5$ratio.DMR <- (DMR_utr5$overlap_end - DMR_utr5$overlap_str +1)/(DMR_utr5$DMR_width +1)
DMR_utr5$ratio.utr5 <- (DMR_utr5$overlap_end - DMR_utr5$overlap_str +1)/(DMR_utr5$utr5_end - DMR_utr5$utr5_str +1)
save(DMR_utr5, file="preDMR_utr5_overlap.RData")
#### summarize
X <- paste(DMR_utr5$chr, DMR_utr5$DMR_str, DMR_utr5$DMR_end, sep="_")
DMR.Tumor.summary$utr5 <- rownames(DMR.Tumor.summary) %in% X
### 3'UTR
DMR_utr3_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.utr3))
DMR_utr3 <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_utr3_hit$queryHits,1:4],
                  as.data.frame(Granges.utr3)[DMR_utr3_hit$subjectHits,c(2,3,7)])
colnames(DMR_utr3) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "utr3_str", "utr3_end", "transcript_name")
DMR_utr3$overlap_str <- apply(DMR_utr3[,c("DMR_str", "utr3_str")], 1, max)
DMR_utr3$overlap_end <- apply(DMR_utr3[,c("DMR_end", "utr3_end")], 1, min)
DMR_utr3$ratio.DMR <- (DMR_utr3$overlap_end - DMR_utr3$overlap_str +1)/(DMR_utr3$DMR_width +1)
DMR_utr3$ratio.utr3 <- (DMR_utr3$overlap_end - DMR_utr3$overlap_str +1)/(DMR_utr3$utr3_end - DMR_utr3$utr3_str +1)
save(DMR_utr3, file="preDMR_utr3_overlap.RData")
#### summarize
X <- paste(DMR_utr3$chr, DMR_utr3$DMR_str, DMR_utr3$DMR_end, sep="_")
DMR.Tumor.summary$utr3 <- rownames(DMR.Tumor.summary) %in% X
### Intron (located among gene but not overlap with any exon, promoter, 5UTR or 3UTR)
DMR.Tumor.summary$Intron <- !(DMR.Tumor.summary$Intergenic | DMR.Tumor.summary$Exon | DMR.Tumor.summary$Promoter | DMR.Tumor.summary$utr5 | DMR.Tumor.summary$utr3)

### modify the overlap according to the priority of physical annotation (promoter = 5utr = 3utr > Exon > Intron > Intergenic)
#### Intergenic & Promoter -> Promoter
DMR.Tumor.summary$Intergenic[DMR.Tumor.summary$Promoter] <- FALSE
#### Exon & Promoter -> Promoter
DMR.Tumor.summary$Exon[DMR.Tumor.summary$Promoter] <-F
#### Exon & utr5 -> utr5
DMR.Tumor.summary$Exon[DMR.Tumor.summary$utr5] <- F
#### Exon & utr3 -> utr3
DMR.Tumor.summary$Exon[DMR.Tumor.summary$utr3] <- F

### modify the overlap according to the priority of pysical annotation (CpG island > CpG shore > CpG shelf)
#### CpG Shelf & CpG Shore -> CpG Shore
#### CpG Shelf & CpG Island -> CpG Island
DMR.Tumor.summary$CpGShelf[DMR.Tumor.summary$CpGShore | DMR.Tumor.summary$CpGIsland] <- FALSE
#### CpG Shore & CpG Island -> CpG Island
DMR.Tumor.summary$CpGShore[DMR.Tumor.summary$CpGIsland] <- FALSE


load("Granges.chroHMM.LiverAdult.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Adjacent_1.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Adjacnt_2.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Tumor_1.Hepatology.hg38.RData")
load("Granges.chroHMM.HCC_Tumor_2.Hepatology.hg38.RData")
load("Granges.chroHMM.LiverAdult.RoadMap.hg38.RData")
load("Granges.chroHMM.HepG2.RoadMap.hg38.RData")
#load("Granges.FANTOM5.enhancer.hg38.RData")
#HepG2.cRE <- read.table(file="/data/huangp/Methy/Mixed/AnnotationFiles/HepG2_hg38.cRE", sep="\t", header=T, stringsAsFactors = F)
#Granges.HepG2.cRE <- makeGRangesFromDataFrame(HepG2.cRE, keep.extra.columns = T,seqnames.field = "chrom", start.field = "chromStart", end.field = "chromEnd")
#save(Granges.HepG2.cRE, file="/data/huangp/Methy/Mixed/AnnotationFiles/Granges.HepG2.cRE.hg38.RData")
load("Granges.HepG2.cRE.hg38.RData")
#Hepatocyte.cRE <- read.table(file="/data/huangp/Methy/Mixed/AnnotationFiles/Hepatocyte_hg38.cRE", sep="\t", header=T, stringsAsFactors = F)
#Granges.Hepatocyte.cRE <- makeGRangesFromDataFrame(Hepatocyte.cRE, keep.extra.columns = T,seqnames.field = "chrom", start.field = "chromStart", end.field = "chromEnd")
#save(Granges.Hepatocyte.cRE, file="/data/huangp/Methy/Mixed/AnnotationFiles/Granges.Hepatocyte.cRE.hg38.RData")
load("Granges.Hepatocyte.cRE.hg38.RData")
#AdultLiver.cRE <- read.table(file="/data/huangp/Methy/Mixed/AnnotationFiles/FemaleAdultLiver_hg38.cRE", sep="\t", header=T, stringsAsFactors = F)
#Granges.AdultLiver.cRE <- makeGRangesFromDataFrame(AdultLiver.cRE, keep.extra.columns = T,seqnames.field = "chrom", start.field = "chromStart", end.field = "chromEnd")
#save(Granges.AdultLiver.cRE, file="/data/huangp/Methy/Mixed/AnnotationFiles/Granges.FemaleAdultLiver.cRE.hg38.RData")
load("Granges.FemaleAdultLiver.cRE.hg38.RData")
### NormalLiver Hepatology
DMR_Normal_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.Normal))
DMR_Normal <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_Normal_hit$queryHits,1:4],
                    as.data.frame(Granges.Normal)[DMR_Normal_hit$subjectHits,c(2,3,6)])
colnames(DMR_Normal) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "chroHMM_str", "chroHMM_end", "state")
DMR_Normal$overlap_str <- apply(DMR_Normal[,c("DMR_str", "chroHMM_str")], 1, max)
DMR_Normal$overlap_end <- apply(DMR_Normal[,c("DMR_end", "chroHMM_end")], 1, min)
DMR_Normal$ratio.DMR <- (DMR_Normal$overlap_end - DMR_Normal$overlap_str +1)/(DMR_Normal$DMR_width +1)
DMR_Normal$ratio.chroHMM <- (DMR_Normal$overlap_end - DMR_Normal$overlap_str +1)/(DMR_Normal$chroHMM_end - DMR_Normal$chroHMM_str +1)
save(DMR_Normal, file="preDMR_chromHMM_Normal_Liver_overlap.RData")
#### collapse
DMR_Normal$DMRID <- paste(DMR_Normal$chr, DMR_Normal$DMR_str, DMR_Normal$DMR_end, sep="_")
DMR_Normal$chroHMMID <- paste(DMR_Normal$chr, DMR_Normal$chroHMM_str, DMR_Normal$chroHMM_end, sep="_")
system.time(DMR_Normal_collapsed <- Collapse(groupVar = DMR_Normal$DMRID, collapseVar = DMR_Normal$state))
colnames(DMR_Normal_collapsed) <- c("DMRID", "state")
DMR.Tumor.summary$chroHMM.Normal <- DMR_Normal_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_Normal_collapsed$DMRID)]

### HCC Adjacent 1 Hepatology
DMR_Adjacent_1_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.Adjacent_1))
DMR_Adjacent_1 <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_Adjacent_1_hit$queryHits,1:4],
                        as.data.frame(Granges.Adjacent_1)[DMR_Adjacent_1_hit$subjectHits,c(2,3,6)])
colnames(DMR_Adjacent_1) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "chroHMM_str", "chroHMM_end", "state")
DMR_Adjacent_1$overlap_str <- apply(DMR_Adjacent_1[,c("DMR_str", "chroHMM_str")], 1, max)
DMR_Adjacent_1$overlap_end <- apply(DMR_Adjacent_1[,c("DMR_end", "chroHMM_end")], 1, min)
DMR_Adjacent_1$ratio.DMR <- (DMR_Adjacent_1$overlap_end - DMR_Adjacent_1$overlap_str +1)/(DMR_Adjacent_1$DMR_width +1)
DMR_Adjacent_1$ratio.chroHMM <- (DMR_Adjacent_1$overlap_end - DMR_Adjacent_1$overlap_str +1)/(DMR_Adjacent_1$chroHMM_end - DMR_Adjacent_1$chroHMM_str +1)
save(DMR_Adjacent_1, file="preDMR_chromHMM_Adjacent_1_Liver_overlap.RData")
#### collapse
DMR_Adjacent_1$DMRID <- paste(DMR_Adjacent_1$chr, DMR_Adjacent_1$DMR_str, DMR_Adjacent_1$DMR_end, sep="_")
DMR_Adjacent_1$chroHMMID <- paste(DMR_Adjacent_1$chr, DMR_Adjacent_1$chroHMM_str, DMR_Adjacent_1$chroHMM_end, sep="_")
system.time(DMR_Adjacent_1_collapsed <- Collapse(groupVar = DMR_Adjacent_1$DMRID, collapseVar = DMR_Adjacent_1$state))
colnames(DMR_Adjacent_1_collapsed) <- c("DMRID", "state")
#save(DMR_Adjacent_1_collapsed, file="DMR_Adjacent_1_collapsed.RData")
DMR.Tumor.summary$chroHMM.Adjacent_1 <- DMR_Adjacent_1_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_Adjacent_1_collapsed$DMRID)]

### HCC Adjacent 2 Hepatology
DMR_Adjacent_2_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.Adjacent_2))
DMR_Adjacent_2 <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_Adjacent_2_hit$queryHits,1:4],
                        as.data.frame(Granges.Adjacent_2)[DMR_Adjacent_2_hit$subjectHits,c(2,3,6)])
colnames(DMR_Adjacent_2) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "chroHMM_str", "chroHMM_end", "state")
DMR_Adjacent_2$overlap_str <- apply(DMR_Adjacent_2[,c("DMR_str", "chroHMM_str")], 1, max)
DMR_Adjacent_2$overlap_end <- apply(DMR_Adjacent_2[,c("DMR_end", "chroHMM_end")], 1, min)
DMR_Adjacent_2$ratio.DMR <- (DMR_Adjacent_2$overlap_end - DMR_Adjacent_2$overlap_str +1)/(DMR_Adjacent_2$DMR_width +1)
DMR_Adjacent_2$ratio.chroHMM <- (DMR_Adjacent_2$overlap_end - DMR_Adjacent_2$overlap_str +1)/(DMR_Adjacent_2$chroHMM_end - DMR_Adjacent_2$chroHMM_str +1)
save(DMR_Adjacent_2, file="preDMR_chromHMM_Adjacent_2_Liver_overlap.RData")
#### collapse
DMR_Adjacent_2$DMRID <- paste(DMR_Adjacent_2$chr, DMR_Adjacent_2$DMR_str, DMR_Adjacent_2$DMR_end, sep="_")
DMR_Adjacent_2$chroHMMID <- paste(DMR_Adjacent_2$chr, DMR_Adjacent_2$chroHMM_str, DMR_Adjacent_2$chroHMM_end, sep="_")
system.time(DMR_Adjacent_2_collapsed <- Collapse(groupVar = DMR_Adjacent_2$DMRID, collapseVar = DMR_Adjacent_2$state))
colnames(DMR_Adjacent_2_collapsed) <- c("DMRID", "state")
#save(DMR_Adjacent_2_collapsed, file="DMR_Adjacent_2_collapsed.RData")
DMR.Tumor.summary$chroHMM.Adjacent_2 <- DMR_Adjacent_2_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_Adjacent_2_collapsed$DMRID)]

### HCC Tumor 1 Hepatology
DMR_Tumor_1_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.Tumor_1))
DMR_Tumor_1 <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_Tumor_1_hit$queryHits,1:4],
                     as.data.frame(Granges.Tumor_1)[DMR_Tumor_1_hit$subjectHits,c(2,3,6)])
colnames(DMR_Tumor_1) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "chroHMM_str", "chroHMM_end", "state")
DMR_Tumor_1$overlap_str <- apply(DMR_Tumor_1[,c("DMR_str", "chroHMM_str")], 1, max)
DMR_Tumor_1$overlap_end <- apply(DMR_Tumor_1[,c("DMR_end", "chroHMM_end")], 1, min)
DMR_Tumor_1$ratio.DMR <- (DMR_Tumor_1$overlap_end - DMR_Tumor_1$overlap_str +1)/(DMR_Tumor_1$DMR_width +1)
DMR_Tumor_1$ratio.chroHMM <- (DMR_Tumor_1$overlap_end - DMR_Tumor_1$overlap_str +1)/(DMR_Tumor_1$chroHMM_end - DMR_Tumor_1$chroHMM_str +1)
save(DMR_Tumor_1, file="preDMR_chromHMM_Tumor_1_Liver_overlap.RData")
#### collapse
DMR_Tumor_1$DMRID <- paste(DMR_Tumor_1$chr, DMR_Tumor_1$DMR_str, DMR_Tumor_1$DMR_end, sep="_")
DMR_Tumor_1$chroHMMID <- paste(DMR_Tumor_1$chr, DMR_Tumor_1$chroHMM_str, DMR_Tumor_1$chroHMM_end, sep="_")
system.time(DMR_Tumor_1_collapsed <- Collapse(groupVar = DMR_Tumor_1$DMRID, collapseVar = DMR_Tumor_1$state))
colnames(DMR_Tumor_1_collapsed) <- c("DMRID", "state")
#save(DMR_Tumor_1_collapsed, file="DMR_Tumor_1_collapsed.RData")
DMR.Tumor.summary$chroHMM.Tumor_1 <- DMR_Tumor_1_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_Tumor_1_collapsed$DMRID)]

### HCC Tumor 2 Hepatology
DMR_Tumor_2_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.Tumor_2))
DMR_Tumor_2 <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_Tumor_2_hit$queryHits,1:4],
                     as.data.frame(Granges.Tumor_2)[DMR_Tumor_2_hit$subjectHits,c(2,3,6)])
colnames(DMR_Tumor_2) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "chroHMM_str", "chroHMM_end", "state")
DMR_Tumor_2$overlap_str <- apply(DMR_Tumor_2[,c("DMR_str", "chroHMM_str")], 1, max)
DMR_Tumor_2$overlap_end <- apply(DMR_Tumor_2[,c("DMR_end", "chroHMM_end")], 1, min)
DMR_Tumor_2$ratio.DMR <- (DMR_Tumor_2$overlap_end - DMR_Tumor_2$overlap_str +1)/(DMR_Tumor_2$DMR_width +1)
DMR_Tumor_2$ratio.chroHMM <- (DMR_Tumor_2$overlap_end - DMR_Tumor_2$overlap_str +1)/(DMR_Tumor_2$chroHMM_end - DMR_Tumor_2$chroHMM_str +1)
save(DMR_Tumor_2, file="preDMR_chromHMM_Tumor_2_Liver_overlap.RData")
#### collapse
DMR_Tumor_2$DMRID <- paste(DMR_Tumor_2$chr, DMR_Tumor_2$DMR_str, DMR_Tumor_2$DMR_end, sep="_")
DMR_Tumor_2$chroHMMID <- paste(DMR_Tumor_2$chr, DMR_Tumor_2$chroHMM_str, DMR_Tumor_2$chroHMM_end, sep="_")
system.time(DMR_Tumor_2_collapsed <- Collapse(groupVar = DMR_Tumor_2$DMRID, collapseVar = DMR_Tumor_2$state))
colnames(DMR_Tumor_2_collapsed) <- c("DMRID", "state")
#save(DMR_Tumor_2_collapsed, file="DMR_Tumor_2_collapsed.RData")
DMR.Tumor.summary$chroHMM.Tumor_2 <- DMR_Tumor_2_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_Tumor_2_collapsed$DMRID)]

### Normal Liver Roadmap
DMR_Liver_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.Liver))
DMR_Liver <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_Liver_hit$queryHits,1:4],
                   as.data.frame(Granges.Liver)[DMR_Liver_hit$subjectHits,c(2,3,6)])
colnames(DMR_Liver) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "chroHMM_str", "chroHMM_end", "state")
DMR_Liver$overlap_str <- apply(DMR_Liver[,c("DMR_str", "chroHMM_str")], 1, max)
DMR_Liver$overlap_end <- apply(DMR_Liver[,c("DMR_end", "chroHMM_end")], 1, min)
DMR_Liver$ratio.DMR <- (DMR_Liver$overlap_end - DMR_Liver$overlap_str +1)/(DMR_Liver$DMR_width +1)
DMR_Liver$ratio.chroHMM <- (DMR_Liver$overlap_end - DMR_Liver$overlap_str +1)/(DMR_Liver$chroHMM_end - DMR_Liver$chroHMM_str +1)
save(DMR_Liver, file="preDMR_chromHMM_Roadmap_Adult_Liver_overlap.RData")
#### collapse
DMR_Liver$DMRID <- paste(DMR_Liver$chr, DMR_Liver$DMR_str, DMR_Liver$DMR_end, sep="_")
DMR_Liver$chroHMMID <- paste(DMR_Liver$chr, DMR_Liver$chroHMM_str, DMR_Liver$chroHMM_end, sep="_")
system.time(DMR_Liver_collapsed <- Collapse(groupVar = DMR_Liver$DMRID, collapseVar = DMR_Liver$state))
colnames(DMR_Liver_collapsed) <- c("DMRID", "state")
#save(DMR_Liver_collapsed, file="DMR_Liver_collapsed.RData")
DMR.Tumor.summary$chroHMM.Liver <- DMR_Liver_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_Liver_collapsed$DMRID)]


### HepG2 Roadmap
DMR_HepG2_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.HepG2))
DMR_HepG2 <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_HepG2_hit$queryHits,1:4],
                   as.data.frame(Granges.HepG2)[DMR_HepG2_hit$subjectHits,c(2,3,6)])
colnames(DMR_HepG2) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "chroHMM_str", "chroHMM_end", "state")
DMR_HepG2$overlap_str <- apply(DMR_HepG2[,c("DMR_str", "chroHMM_str")], 1, max)
DMR_HepG2$overlap_end <- apply(DMR_HepG2[,c("DMR_end", "chroHMM_end")], 1, min)
DMR_HepG2$ratio.DMR <- (DMR_HepG2$overlap_end - DMR_HepG2$overlap_str +1)/(DMR_HepG2$DMR_width +1)
DMR_HepG2$ratio.chroHMM <- (DMR_HepG2$overlap_end - DMR_HepG2$overlap_str +1)/(DMR_HepG2$chroHMM_end - DMR_HepG2$chroHMM_str +1)
save(DMR_HepG2, file="preDMR_chromHMM_Roadmap_HepG2_overlap.RData")
#### collapse
DMR_HepG2$DMRID <- paste(DMR_HepG2$chr, DMR_HepG2$DMR_str, DMR_HepG2$DMR_end, sep="_")
DMR_HepG2$chroHMMID <- paste(DMR_HepG2$chr, DMR_HepG2$chroHMM_str, DMR_HepG2$chroHMM_end, sep="_")
system.time(DMR_HepG2_collapsed <- Collapse(groupVar = DMR_HepG2$DMRID, collapseVar = DMR_HepG2$state))
colnames(DMR_HepG2_collapsed) <- c("DMRID", "state")
#save(DMR_HepG2_collapsed, file="DMR_HepG2_collapsed.RData")
DMR.Tumor.summary$chroHMM.HepG2 <- DMR_HepG2_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_HepG2_collapsed$DMRID)]

### HepG2 cRE ENCODE
DMR_HepG2.cRE_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.HepG2.cRE))
DMR_HepG2.cRE <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_HepG2.cRE_hit$queryHits,1:4],
                       as.data.frame(Granges.HepG2.cRE)[DMR_HepG2.cRE_hit$subjectHits,c(2,3,10)])
colnames(DMR_HepG2.cRE) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "cRE_str", "cRE_end", "state")
DMR_HepG2.cRE$overlap_str <- apply(DMR_HepG2.cRE[,c("DMR_str", "cRE_str")], 1, max)
DMR_HepG2.cRE$overlap_end <- apply(DMR_HepG2.cRE[,c("DMR_end", "cRE_end")], 1, min)
DMR_HepG2.cRE$ratio.DMR <- (DMR_HepG2.cRE$overlap_end - DMR_HepG2.cRE$overlap_str +1)/(DMR_HepG2.cRE$DMR_width +1)
DMR_HepG2.cRE$ratio.chroHMM <- (DMR_HepG2.cRE$overlap_end - DMR_HepG2.cRE$overlap_str +1)/(DMR_HepG2.cRE$cRE_end - DMR_HepG2.cRE$cRE_str +1)
save(DMR_HepG2.cRE, file="preDMR_cRE_HepG2_overlap.RData")
#### collapse
DMR_HepG2.cRE$DMRID <- paste(DMR_HepG2.cRE$chr, DMR_HepG2.cRE$DMR_str, DMR_HepG2.cRE$DMR_end, sep="_")
system.time(DMR_HepG2.cRE_collapsed <- Collapse(groupVar = DMR_HepG2.cRE$DMRID, collapseVar = DMR_HepG2.cRE$state))
colnames(DMR_HepG2.cRE_collapsed) <- c("DMRID", "state")
#save(DMR_HepG2.cRE_collapsed, file="DMR_HepG2.cRE_collapsed.RData")
DMR.Tumor.summary$cRE.HepG2 <- DMR_HepG2.cRE_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_HepG2.cRE_collapsed$DMRID)]
### FemaleAdultLiver cRE ENCODE
DMR_AdultLiver.cRE_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.AdultLiver.cRE))
DMR_AdultLiver.cRE <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_AdultLiver.cRE_hit$queryHits,1:4],
                            as.data.frame(Granges.AdultLiver.cRE)[DMR_AdultLiver.cRE_hit$subjectHits,c(2,3,10)])
colnames(DMR_AdultLiver.cRE) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "cRE_str", "cRE_end", "state")
DMR_AdultLiver.cRE$overlap_str <- apply(DMR_AdultLiver.cRE[,c("DMR_str", "cRE_str")], 1, max)
DMR_AdultLiver.cRE$overlap_end <- apply(DMR_AdultLiver.cRE[,c("DMR_end", "cRE_end")], 1, min)
DMR_AdultLiver.cRE$ratio.DMR <- (DMR_AdultLiver.cRE$overlap_end - DMR_AdultLiver.cRE$overlap_str +1)/(DMR_AdultLiver.cRE$DMR_width +1)
DMR_AdultLiver.cRE$ratio.chroHMM <- (DMR_AdultLiver.cRE$overlap_end - DMR_AdultLiver.cRE$overlap_str +1)/(DMR_AdultLiver.cRE$cRE_end - DMR_AdultLiver.cRE$cRE_str +1)
save(DMR_AdultLiver.cRE, file="preDMR_cRE_AdultLiver_overlap.RData")
#### collapse
DMR_AdultLiver.cRE$DMRID <- paste(DMR_AdultLiver.cRE$chr, DMR_AdultLiver.cRE$DMR_str, DMR_AdultLiver.cRE$DMR_end, sep="_")
system.time(DMR_AdultLiver.cRE_collapsed <- Collapse(groupVar = DMR_AdultLiver.cRE$DMRID, collapseVar = DMR_AdultLiver.cRE$state))
colnames(DMR_AdultLiver.cRE_collapsed) <- c("DMRID", "state")
#save(DMR_AdultLiver.cRE_collapsed, file="DMR_AdultLiver.cRE_collapsed.RData")
DMR.Tumor.summary$cRE.AdultLiver <- DMR_AdultLiver.cRE_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_AdultLiver.cRE_collapsed$DMRID)]
### Hepatocyte cRE ENCODE
DMR_Hepatocyte.cRE_hit <- as.data.frame(findOverlaps(Granges.DMR.Tumor, Granges.Hepatocyte.cRE))
DMR_Hepatocyte.cRE <- cbind(as.data.frame(Granges.DMR.Tumor)[DMR_Hepatocyte.cRE_hit$queryHits,1:4],
                            as.data.frame(Granges.Hepatocyte.cRE)[DMR_Hepatocyte.cRE_hit$subjectHits,c(2,3,10)])
colnames(DMR_Hepatocyte.cRE) <- c("chr", "DMR_str", "DMR_end", "DMR_width", "cRE_str", "cRE_end", "state")
DMR_Hepatocyte.cRE$overlap_str <- apply(DMR_Hepatocyte.cRE[,c("DMR_str", "cRE_str")], 1, max)
DMR_Hepatocyte.cRE$overlap_end <- apply(DMR_Hepatocyte.cRE[,c("DMR_end", "cRE_end")], 1, min)
DMR_Hepatocyte.cRE$ratio.DMR <- (DMR_Hepatocyte.cRE$overlap_end - DMR_Hepatocyte.cRE$overlap_str +1)/(DMR_Hepatocyte.cRE$DMR_width +1)
DMR_Hepatocyte.cRE$ratio.chroHMM <- (DMR_Hepatocyte.cRE$overlap_end - DMR_Hepatocyte.cRE$overlap_str +1)/(DMR_Hepatocyte.cRE$cRE_end - DMR_Hepatocyte.cRE$cRE_str +1)
save(DMR_Hepatocyte.cRE, file="preDMR_cRE_Hepatocyte_overlap.RData")
#### collapse
DMR_Hepatocyte.cRE$DMRID <- paste(DMR_Hepatocyte.cRE$chr, DMR_Hepatocyte.cRE$DMR_str, DMR_Hepatocyte.cRE$DMR_end, sep="_")
system.time(DMR_Hepatocyte.cRE_collapsed <- Collapse(groupVar = DMR_Hepatocyte.cRE$DMRID, collapseVar = DMR_Hepatocyte.cRE$state))
colnames(DMR_Hepatocyte.cRE_collapsed) <- c("DMRID", "state")
#save(DMR_Hepatocyte.cRE_collapsed, file="DMR_Hepatocyte.cRE_collapsed.RData")
DMR.Tumor.summary$cRE.Hepatocyte <- DMR_Hepatocyte.cRE_collapsed$state[match(rownames(DMR.Tumor.summary), DMR_Hepatocyte.cRE_collapsed$DMRID)]

### pull out active TSS, active Enhancer, active Promoter
#### Normal Liver 
DMR.Tumor.summary$activeTSS.Normal <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeTSS.Normal[grep("5_state5", DMR.Tumor.summary$chroHMM.Normal, fixed=T)] <- TRUE
DMR.Tumor.summary$activeEnhancer.Normal <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeEnhancer.Normal[grep("7_state7", DMR.Tumor.summary$chroHMM.Normal, fixed=T)] <- TRUE
DMR.Tumor.summary$activePromoter.Normal <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activePromoter.Normal[grep("10_state10", DMR.Tumor.summary$chroHMM.Normal, fixed=T)] <- TRUE
#### Adjacent_1
DMR.Tumor.summary$activeTSS.Adjacent_1 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeTSS.Adjacent_1[grep("5_state5", DMR.Tumor.summary$chroHMM.Adjacent_1, fixed=T)] <- TRUE
DMR.Tumor.summary$activeEnhancer.Adjacent_1 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeEnhancer.Adjacent_1[grep("7_state7", DMR.Tumor.summary$chroHMM.Adjacent_1, fixed=T)] <- TRUE
DMR.Tumor.summary$activePromoter.Adjacent_1 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activePromoter.Adjacent_1[grep("10_state10", DMR.Tumor.summary$chroHMM.Adjacent_1, fixed=T)] <- TRUE

#### Adjacent_2
DMR.Tumor.summary$activeTSS.Adjacent_2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeTSS.Adjacent_2[grep("5_state5", DMR.Tumor.summary$chroHMM.Adjacent_2, fixed=T)] <- TRUE
DMR.Tumor.summary$activeEnhancer.Adjacent_2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeEnhancer.Adjacent_2[grep("7_state7", DMR.Tumor.summary$chroHMM.Adjacent_2, fixed=T)] <- TRUE
DMR.Tumor.summary$activePromoter.Adjacent_2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activePromoter.Adjacent_2[grep("10_state10", DMR.Tumor.summary$chroHMM.Adjacent_2, fixed=T)] <- TRUE
#### Tumor_1
DMR.Tumor.summary$activeTSS.Tumor_1 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeTSS.Tumor_1[grep("5_state5", DMR.Tumor.summary$chroHMM.Tumor_1, fixed=T)] <- TRUE
DMR.Tumor.summary$activeEnhancer.Tumor_1 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeEnhancer.Tumor_1[grep("7_state7", DMR.Tumor.summary$chroHMM.Tumor_1, fixed=T)] <- TRUE
DMR.Tumor.summary$activePromoter.Tumor_1 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activePromoter.Tumor_1[grep("10_state10", DMR.Tumor.summary$chroHMM.Tumor_1, fixed=T)] <- TRUE
#### Tumor_2
DMR.Tumor.summary$activeTSS.Tumor_2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeTSS.Tumor_2[grep("5_state5", DMR.Tumor.summary$chroHMM.Tumor_2, fixed=T)] <- TRUE
DMR.Tumor.summary$activeEnhancer.Tumor_2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activeEnhancer.Tumor_2[grep("7_state7", DMR.Tumor.summary$chroHMM.Tumor_2, fixed=T)] <- TRUE
DMR.Tumor.summary$activePromoter.Tumor_2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$activePromoter.Tumor_2[grep("10_state10", DMR.Tumor.summary$chroHMM.Tumor_2, fixed=T)] <- TRUE
### HepG2 cRE
DMR.Tumor.summary$Promoter.cRE.HepG2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$Promoter.cRE.HepG2[grep("Promoter", DMR.Tumor.summary$cRE.HepG2, fixed=T)] <-TRUE
DMR.Tumor.summary$Enhancer.cRE.HepG2 <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$Enhancer.cRE.HepG2[grep("Enhancer", DMR.Tumor.summary$cRE.HepG2, fixed=T)] <-TRUE
### Hepatocyte cRE
DMR.Tumor.summary$Promoter.cRE.Hepatocyte <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$Promoter.cRE.Hepatocyte[grep("Promoter", DMR.Tumor.summary$cRE.Hepatocyte, fixed=T)] <-TRUE
DMR.Tumor.summary$Enhancer.cRE.Hepatocyte <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$Enhancer.cRE.Hepatocyte[grep("Enhancer", DMR.Tumor.summary$cRE.Hepatocyte, fixed=T)] <-TRUE
### AdultLiver cRE
DMR.Tumor.summary$Promoter.cRE.AdultLiver <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$Promoter.cRE.AdultLiver[grep("Promoter", DMR.Tumor.summary$cRE.AdultLiver, fixed=T)] <-TRUE
DMR.Tumor.summary$Enhancer.cRE.AdultLiver <- rep(FALSE, nrow(DMR.Tumor.summary))
DMR.Tumor.summary$Enhancer.cRE.AdultLiver[grep("Enhancer", DMR.Tumor.summary$cRE.AdultLiver, fixed=T)] <-TRUE
attach(DMR.Tumor.summary)
table(activeTSS.Normal)
table(activeTSS.Adjacent_1)
table(activeTSS.Adjacent_2)
table(activeTSS.Tumor_1)
table(activeTSS.Tumor_2)
table(Promoter.cRE.HepG2)
table(Promoter.cRE.Hepatocyte)
table(Promoter.cRE.AdultLiver)

# overlaped with promoter region
table(Promoter)
table(activePromoter.Normal)
table(activePromoter.Adjacent_1)
table(activePromoter.Adjacent_2)
table(activePromoter.Tumor_1)
table(activePromoter.Tumor_2)
table(activeTSS.Normal|activeTSS.Adjacent_1|activeTSS.Adjacent_2|activeTSS.Tumor_1|activeTSS.Tumor_2|
        Promoter|activePromoter.Normal|activePromoter.Adjacent_1|activePromoter.Adjacent_2|activePromoter.Tumor_1|activePromoter.Tumor_2|
        Promoter.cRE.HepG2|Promoter.cRE.Hepatocyte|Promoter.cRE.AdultLiver)

table(activeTSS.Normal|activeTSS.Adjacent_1|activeTSS.Adjacent_2|activeTSS.Tumor_1|activeTSS.Tumor_2)
table(activePromoter.Normal|activePromoter.Adjacent_1|activePromoter.Adjacent_2|activePromoter.Tumor_1|activePromoter.Tumor_2|
        Promoter.cRE.HepG2|Promoter.cRE.Hepatocyte|Promoter.cRE.AdultLiver)
table(Promoter)
DMR.Tumor.summary$PromoterLike <- activeTSS.Normal|activeTSS.Adjacent_1|activeTSS.Adjacent_2|activeTSS.Tumor_1|activeTSS.Tumor_2|
  activePromoter.Normal|activePromoter.Adjacent_1|activePromoter.Adjacent_2|activePromoter.Tumor_1|activePromoter.Tumor_2|
  Promoter.cRE.HepG2|Promoter.cRE.Hepatocyte|Promoter.cRE.AdultLiver

table(activeEnhancer.Normal)
table(activeEnhancer.Adjacent_1)
table(activeEnhancer.Adjacent_2)
table(activeEnhancer.Tumor_1)
table(activeEnhancer.Tumor_2)

table(Enhancer.cRE.HepG2)
table(Enhancer.cRE.Hepatocyte)
table(Enhancer.cRE.AdultLiver)
table(activeEnhancer.Normal|activeEnhancer.Adjacent_1|activeEnhancer.Adjacent_2|activeEnhancer.Tumor_1|activeEnhancer.Tumor_2|
        Enhancer.cRE.HepG2|Enhancer.cRE.Hepatocyte|Enhancer.cRE.AdultLiver)

DMR.Tumor.summary$EnhancerLike <- activeEnhancer.Normal|activeEnhancer.Adjacent_1|activeEnhancer.Adjacent_2|activeEnhancer.Tumor_1|activeEnhancer.Tumor_2|
  Enhancer.cRE.HepG2|Enhancer.cRE.Hepatocyte|Enhancer.cRE.AdultLiver
detach(DMR.Tumor.summary)

preDMR.summary <- DMR.Tumor.summary
save.lbzip2(preDMR.summary, file="preDMRs.summary.RDataFS", n.cores=32)
PromoterDMR.Tumor <- DMR.Tumor.summary[DMR.Tumor.summary$PromoterLike|DMR.Tumor.summary$Promoter,]
EnhancerDMR.Tumor <- DMR.Tumor.summary[DMR.Tumor.summary$EnhancerLike,]
save.lbzip2(PromoterDMR.Tumor, file="PromoterLike_preDMRs.RDataFS", n.cores=24)
save.lbzip2(EnhancerDMR.Tumor, file="EnhancerLike_preDMRs.RDataFS", n.cores=24)

#### Overlap ratio filtering (cut off = 0%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 100%) to hightlight promoter-/enhancer-like DMRs with high confidence
#### for long range features with multiple overlaped DMRs, combined their overlap ratio of feature together 
options(sstringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")#29
load("preDMR_promoter_overlap.RData")
load("preDMR_chromHMM_Normal_Liver_overlap.RData")
load("preDMR_chromHMM_Adjacent_1_Liver_overlap.RData")
load("preDMR_chromHMM_Adjacent_2_Liver_overlap.RData")
load("preDMR_chromHMM_Tumor_1_Liver_overlap.RData")
load("preDMR_chromHMM_Tumor_2_Liver_overlap.RData")
load("preDMR_cRE_AdultLiver_overlap.RData")
load("preDMR_cRE_Hepatocyte_overlap.RData")
load("preDMR_cRE_HepG2_overlap.RData")
load.lbzip2("preDMRs.summary.RDataFS", n.cores=24)
PromoterDMR.summary <- preDMR.summary[preDMR.summary$PromoterLike.score>0|preDMR.summary$Promoter,]#47493
summary(PromoterDMR.summary$width)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#9.0   280.0   443.0   615.2   718.0 40440.0
EnhancerDMR.summary <- preDMR.summary[preDMR.summary$EnhancerLike.score > 0,]#23800
summary(EnhancerDMR.summary$width)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.0   253.0   381.0   479.2   566.0 16803.0 
table(rownames(EnhancerDMR.summary) %in% rownames(PromoterDMR.summary))#16538 7262
#keep only promoter or enhancer-like chromHMM state or cRE state
DMR_promoter$FeatureID <- paste(DMR_promoter$chr, DMR_promoter$promoter_str, DMR_promoter$promoter_end, sep="_")
DMR_promoter$DMRID <- paste(DMR_promoter$chr, DMR_promoter$DMR_str, DMR_promoter$DMR_end, sep="_")
DMR_Normal$FeatureID <- paste(DMR_Normal$chr, DMR_Normal$chroHMM_str, DMR_Normal$chroHMM_end, DMR_Normal$state, sep="_")
DMR_Normal$DMRID <- paste(DMR_Normal$chr, DMR_Normal$DMR_str, DMR_Normal$DMR_end, sep="_")
DMR_Normal <- DMR_Normal[DMR_Normal$state %in% c("5_state5", "7_state7", "10_state10"),]#9050
DMR_Adjacent_1$FeatureID <- paste(DMR_Adjacent_1$chr, DMR_Adjacent_1$chroHMM_str, DMR_Adjacent_1$chroHMM_end, DMR_Adjacent_1$state, sep="_")
DMR_Adjacent_1$DMRID <- paste(DMR_Adjacent_1$chr, DMR_Adjacent_1$DMR_str, DMR_Adjacent_1$DMR_end, sep="_")
DMR_Adjacent_1 <- DMR_Adjacent_1[DMR_Adjacent_1$state %in% c("5_state5", "7_state7", "10_state10"),]#16826
DMR_Adjacent_2$FeatureID <- paste(DMR_Adjacent_2$chr, DMR_Adjacent_2$chroHMM_str, DMR_Adjacent_2$chroHMM_end, DMR_Adjacent_2$state, sep="_")
DMR_Adjacent_2$DMRID <- paste(DMR_Adjacent_2$chr, DMR_Adjacent_2$DMR_str, DMR_Adjacent_2$DMR_end, sep="_")
DMR_Adjacent_2 <- DMR_Adjacent_2[DMR_Adjacent_2$state %in% c("5_state5", "7_state7", "10_state10"),]#15928
DMR_Tumor_1$FeatureID <- paste(DMR_Tumor_1$chr, DMR_Tumor_1$chroHMM_str, DMR_Tumor_1$chroHMM_end, DMR_Tumor_1$state, sep="_")
DMR_Tumor_1$DMRID <- paste(DMR_Tumor_1$chr, DMR_Tumor_1$DMR_str, DMR_Tumor_1$DMR_end, sep="_")
DMR_Tumor_1 <- DMR_Tumor_1[DMR_Tumor_1$state %in% c("5_state5", "7_state7", "10_state10"),]#11879
DMR_Tumor_2$FeatureID <- paste(DMR_Tumor_2$chr, DMR_Tumor_2$chroHMM_str, DMR_Tumor_2$chroHMM_end, DMR_Tumor_2$state, sep="_")
DMR_Tumor_2$DMRID <- paste(DMR_Tumor_2$chr, DMR_Tumor_2$DMR_str, DMR_Tumor_2$DMR_end, sep="_")
DMR_Tumor_2 <- DMR_Tumor_2[DMR_Tumor_2$state %in% c("5_state5", "7_state7", "10_state10"),]#14194
DMR_AdultLiver.cRE$FeatureID <- paste(DMR_AdultLiver.cRE$chr, DMR_AdultLiver.cRE$cRE_str, DMR_AdultLiver.cRE$cRE_end, DMR_AdultLiver.cRE$state, sep="_")
DMR_AdultLiver.cRE$DMRID <- paste(DMR_AdultLiver.cRE$chr, DMR_AdultLiver.cRE$DMR_str, DMR_AdultLiver.cRE$DMR_end, sep="_")
DMR_AdultLiver.cRE <- DMR_AdultLiver.cRE[DMR_AdultLiver.cRE$state %in% c("Promoter", "Enhancer"),]#5201
DMR_Hepatocyte.cRE$FeatureID <- paste(DMR_Hepatocyte.cRE$chr, DMR_Hepatocyte.cRE$cRE_str, DMR_Hepatocyte.cRE$cRE_end, DMR_Hepatocyte.cRE$state, sep="_")
DMR_Hepatocyte.cRE$DMRID <- paste(DMR_Hepatocyte.cRE$chr, DMR_Hepatocyte.cRE$DMR_str, DMR_Hepatocyte.cRE$DMR_end, sep="_")
DMR_Hepatocyte.cRE <- DMR_Hepatocyte.cRE[DMR_Hepatocyte.cRE$state %in% c("Promoter", "Enhancer"),]#7691
DMR_HepG2.cRE$FeatureID <- paste(DMR_HepG2.cRE$chr, DMR_HepG2.cRE$cRE_str, DMR_HepG2.cRE$cRE_end, DMR_HepG2.cRE$state, sep="_")
DMR_HepG2.cRE$DMRID <- paste(DMR_HepG2.cRE$chr, DMR_HepG2.cRE$DMR_str, DMR_HepG2.cRE$DMR_end, sep="_")
DMR_HepG2.cRE <- DMR_HepG2.cRE[DMR_HepG2.cRE$state %in% c("Promoter", "Enhancer"),]#5421

summary(DMR_promoter$ratio.promoter)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0004998 0.0949525 0.1644178 0.2069264 0.2668666 1.0000000
summary(DMR_Normal$ratio.chroHMM)
#0.0003845 0.0815987 0.1690610 0.2456590 0.3277870 1.0000000
summary(DMR_Adjacent_1$ratio.chroHMM)
#0.0000188 0.0539820 0.1249519 0.1952870 0.2556203 1.0000000
summary(DMR_Adjacent_2$ratio.chroHMM)
#0.0000187 0.0584695 0.1253670 0.1905766 0.2498216 1.0000000
summary(DMR_Tumor_1$ratio.chroHMM)
#0.0000187 0.0754879 0.1610487 0.2377862 0.3171036 1.0000000
summary(DMR_Tumor_2$ratio.chroHMM)
#0.0000187 0.0620293 0.1369190 0.2095160 0.2770683 1.0000000
summary(DMR_AdultLiver.cRE$ratio.chroHMM)
#0.0005058 0.1572165 0.3756906 0.4391865 0.6877637 1.0000000 
summary(DMR_Hepatocyte.cRE$ratio.chroHMM)
#0.0005058 0.1462224 0.3609756 0.4392148 0.7134110 1.0000000
summary(DMR_HepG2.cRE$ratio.chroHMM)
#0.0005058 0.1586052 0.3814757 0.4508474 0.7270955 1.0000000


summary(DMR_Normal$ratio.chroHMM[DMR_Normal$state=="7_state7"])#enhancer
#0.0003845 0.0699636 0.1449525 0.2148702 0.2802499 1.0000000
summary(DMR_Normal$ratio.chroHMM[DMR_Normal$state!="7_state7"])#promoter
#0.0004543 0.1099695 0.2241499 0.2991020 0.4040740 1.0000000
summary(DMR_Adjacent_1$ratio.chroHMM[DMR_Adjacent_1$state=="7_state7"])#enhancer
#0.0001921 0.0452016 0.1003748 0.1623159 0.2059209 1.0000000
summary(DMR_Adjacent_1$ratio.chroHMM[DMR_Adjacent_1$state!="7_state7"])#promoter
#0.0000188 0.0898933 0.1865630 0.2582639 0.3518915 1.0000000
summary(DMR_Adjacent_2$ratio.chroHMM[DMR_Adjacent_2$state=="7_state7"])#enhancer
#0.0000187 0.0472804 0.1028738 0.1604788 0.2060362 1.0000000
summary(DMR_Adjacent_2$ratio.chroHMM[DMR_Adjacent_2$state!="7_state7"])#promoter
#0.0006246 0.1038312 0.1979010 0.2686218 0.3605123 1.0000000
summary(DMR_Tumor_1$ratio.chroHMM[DMR_Tumor_1$state=="7_state7"])#enhancer
#0.0000187 0.0592080 0.1307244 0.2030918 0.2676660 1.0000000
summary(DMR_Tumor_1$ratio.chroHMM[DMR_Tumor_1$state!="7_state7"])#promoter
#0.0005552 0.1023720 0.2027123 0.2785467 0.3676324 1.0000000
summary(DMR_Tumor_2$ratio.chroHMM[DMR_Tumor_2$state=="7_state7"])#enhancer
#0.0000187 0.0517166 0.1103072 0.1652928 0.2144461 1.0000000
summary(DMR_Tumor_2$ratio.chroHMM[DMR_Tumor_2$state!="7_state7"])#promoter
#0.0007689 0.1184753 0.2436010 0.3179447 0.4393136 1.0000000
summary(DMR_AdultLiver.cRE$ratio.chroHMM[DMR_AdultLiver.cRE$state=="Enhancer"])#enhancer
#0.0007559 0.2042921 0.4378049 0.4761926 0.7295462 1.0000000
summary(DMR_AdultLiver.cRE$ratio.chroHMM[DMR_AdultLiver.cRE$state=="Promoter"])#promoter
#0.0005058 0.0938740 0.2543103 0.3623731 0.5736670 1.0000000
summary(DMR_Hepatocyte.cRE$ratio.chroHMM[DMR_Hepatocyte.cRE$state=="Enhancer"])#enhancer
#0.0009681 0.1916673 0.4191280 0.4801477 0.7629101 1.0000000
summary(DMR_Hepatocyte.cRE$ratio.chroHMM[DMR_Hepatocyte.cRE$state=="Promoter"])#promoter
#0.0005058 0.0865874 0.2194030 0.3502736 0.5627738 1.0000000
summary(DMR_HepG2.cRE$ratio.chroHMM[DMR_HepG2.cRE$state=="Enhancer"])#enhancer
#0.0007559 0.2062589 0.4417040 0.4871748 0.7642974 1.0000000
summary(DMR_HepG2.cRE$ratio.chroHMM[DMR_HepG2.cRE$state=="Promoter"])#promoter
#0.0005058 0.1050426 0.2767553 0.3955792 0.6612999 1.0000000
### cutoff set as 10%, 20%, 30%
### colapse function
COMBINE <- function(x, sep)
{
  x <- as.character(x)
  if(length(x)!=0)
  {
    y <- x[1]
    if(length(x)>1)
    {
      for(i in 2:length(x)) y <- paste(y, x[i], sep = sep)
    }
  }
  return(y)
}

Collapse <- function(groupVar,collapseVar, CL =16, Cores=2)
{
  library(doParallel)
  library(foreach)
  registerDoParallel(cl=CL, cores = Cores)
  if(length(groupVar)==length(collapseVar))
  {
    group <- unique(groupVar)
    collapsed.var <- vector(mode="character", length=length(group))
    collapsed.var <- foreach(i = 1:length(group), .combine = 'c') %dopar%{
      COMBINE(collapseVar[groupVar %in% group[i]], sep="|")
    }
    stopImplicitCluster()
    return(as.data.frame(cbind(group, collapsed.var)))
  }
}
HighConfidentDMR <- function(DMR.summary = preDMR.summary[preDMR.summary$PromoterLike.score > 0 | preDMR.summary$Promoter | preDMR.summary$EnhancerLike.score > 0, 1:14], Feature.ratio = 0.1)
{
  #filtered feature overlap ratio lower than cutoff
  promoterOverlap <- DMR_promoter[DMR_promoter$ratio.promoter >= Feature.ratio,]
  system.time(promoterOverlap_collapsed <- Collapse(groupVar = promoterOverlap$DMRID, collapseVar = promoterOverlap$FeatureID))
  colnames(promoterOverlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$promoter <- promoterOverlap_collapsed$Features[match(rownames(DMR.summary), promoterOverlap_collapsed$DMRID)]
  #table(!is.na(DMR.summary$Promoter)) # 31765
  NormalOverlap <- DMR_Normal[DMR_Normal$ratio.chroHMM >= Feature.ratio,]
  system.time(NormalOverlap_collapsed <- Collapse(groupVar = NormalOverlap$DMRID, collapseVar = NormalOverlap$FeatureID))
  colnames(NormalOverlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$Normal <- NormalOverlap_collapsed$Features[match(rownames(DMR.summary), NormalOverlap_collapsed$DMRID)]
  
  Adjacent_1Overlap <- DMR_Adjacent_1[DMR_Adjacent_1$ratio.chroHMM >= Feature.ratio,]
  system.time(Adjacent_1Overlap_collapsed <- Collapse(groupVar = Adjacent_1Overlap$DMRID, collapseVar = Adjacent_1Overlap$FeatureID))
  colnames(Adjacent_1Overlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$Adjacent_1 <- Adjacent_1Overlap_collapsed$Features[match(rownames(DMR.summary), Adjacent_1Overlap_collapsed$DMRID)]
  
  Adjacent_2Overlap <- DMR_Adjacent_2[DMR_Adjacent_2$ratio.chroHMM >= Feature.ratio,]
  system.time(Adjacent_2Overlap_collapsed <- Collapse(groupVar = Adjacent_2Overlap$DMRID, collapseVar = Adjacent_2Overlap$FeatureID))
  colnames(Adjacent_2Overlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$Adjacent_2 <- Adjacent_2Overlap_collapsed$Features[match(rownames(DMR.summary), Adjacent_2Overlap_collapsed$DMRID)]
  
  Tumor_1Overlap <- DMR_Tumor_1[DMR_Tumor_1$ratio.chroHMM >= Feature.ratio,]
  system.time(Tumor_1Overlap_collapsed <- Collapse(groupVar = Tumor_1Overlap$DMRID, collapseVar = Tumor_1Overlap$FeatureID))
  colnames(Tumor_1Overlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$Tumor_1 <- Tumor_1Overlap_collapsed$Features[match(rownames(DMR.summary), Tumor_1Overlap_collapsed$DMRID)]
  
  Tumor_2Overlap <- DMR_Tumor_2[DMR_Tumor_2$ratio.chroHMM >= Feature.ratio,]
  system.time(Tumor_2Overlap_collapsed <- Collapse(groupVar = Tumor_2Overlap$DMRID, collapseVar = Tumor_2Overlap$FeatureID))
  colnames(Tumor_2Overlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$Tumor_2 <- Tumor_2Overlap_collapsed$Features[match(rownames(DMR.summary), Tumor_2Overlap_collapsed$DMRID)]
  
  AdultLiver.cREOverlap <- DMR_AdultLiver.cRE[DMR_AdultLiver.cRE$ratio.chroHMM >= Feature.ratio,]
  system.time(AdultLiver.cREOverlap_collapsed <- Collapse(groupVar = AdultLiver.cREOverlap$DMRID, collapseVar = AdultLiver.cREOverlap$FeatureID))
  colnames(AdultLiver.cREOverlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$AdultLiver.cRE <- AdultLiver.cREOverlap_collapsed$Features[match(rownames(DMR.summary), AdultLiver.cREOverlap_collapsed$DMRID)]
  
  Hepatocyte.cREOverlap <- DMR_Hepatocyte.cRE[DMR_Hepatocyte.cRE$ratio.chroHMM >= Feature.ratio,]
  system.time(Hepatocyte.cREOverlap_collapsed <- Collapse(groupVar = Hepatocyte.cREOverlap$DMRID, collapseVar = Hepatocyte.cREOverlap$FeatureID))
  colnames(Hepatocyte.cREOverlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$Hepatocyte.cRE <- Hepatocyte.cREOverlap_collapsed$Features[match(rownames(DMR.summary), Hepatocyte.cREOverlap_collapsed$DMRID)]
  
  HepG2.cREOverlap <- DMR_HepG2.cRE[DMR_HepG2.cRE$ratio.chroHMM >= Feature.ratio,]
  system.time(HepG2.cREOverlap_collapsed <- Collapse(groupVar = HepG2.cREOverlap$DMRID, collapseVar = HepG2.cREOverlap$FeatureID))
  colnames(HepG2.cREOverlap_collapsed) <- c("DMRID", "Features")
  DMR.summary$HepG2.cRE <- HepG2.cREOverlap_collapsed$Features[match(rownames(DMR.summary), HepG2.cREOverlap_collapsed$DMRID)]
  #summarize 
  NormalE_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Normal, pattern = "7_state7", fixed = T)
  NormalP_status <- ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Normal, pattern = "5_state5", fixed = T)) |  
                    ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Normal, pattern = "10_state10", fixed = T))
  Adjacent_1E_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Adjacent_1, pattern = "7_state7", fixed = T)
  Adjacent_1P_status <- ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Adjacent_1, pattern = "5_state5", fixed = T)) |  
                    ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Adjacent_1, pattern = "10_state10", fixed = T))
  Adjacent_2E_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Adjacent_2, pattern = "7_state7", fixed = T)
  Adjacent_2P_status <- ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Adjacent_2, pattern = "5_state5", fixed = T)) |  
                    ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Adjacent_2, pattern = "10_state10", fixed = T))
  Tumor_1E_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Tumor_1, pattern = "7_state7", fixed = T)
  Tumor_1P_status <- ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Tumor_1, pattern = "5_state5", fixed = T)) |  
                    ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Tumor_1, pattern = "10_state10", fixed = T))
  Tumor_2E_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Tumor_2, pattern = "7_state7", fixed = T)
  Tumor_2P_status <- ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Tumor_2, pattern = "5_state5", fixed = T)) |  
                    ((1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Tumor_2, pattern = "10_state10", fixed = T))
  AdultLiver.cREE_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$AdultLiver.cRE, pattern = "Enhancer", fixed = T)
  AdultLiver.cREP_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$AdultLiver.cRE, pattern = "Promoter", fixed = T) 
  HepatocyteE_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Hepatocyte, pattern = "Enhancer", fixed = T)
  HepatocyteP_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$Hepatocyte, pattern = "Promoter", fixed = T) 
  HepG2E_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$HepG2, pattern = "Enhancer", fixed = T)
  HepG2P_status <- (1:nrow(DMR.summary)) %in% grep(x = DMR.summary$HepG2, pattern = "Promoter", fixed = T) 
  DMR.summary$PromoterLike.score <- rowSums(cbind(NormalP_status,Adjacent_1P_status, Adjacent_2P_status, AdultLiver.cREP_status, HepatocyteP_status,
                                        Tumor_1P_status, Tumor_2P_status, HepG2P_status))
  DMR.summary$PromoterLike.T <- rowSums(cbind(Tumor_1P_status, Tumor_2P_status, HepG2P_status))
  DMR.summary$PromoterLike.N <- rowSums(cbind(NormalP_status, Adjacent_1P_status, Adjacent_2P_status, AdultLiver.cREP_status, HepatocyteP_status))
  
  DMR.summary$EnhancerLike.score <- rowSums(cbind(NormalE_status, Adjacent_1E_status, Adjacent_2E_status, AdultLiver.cREE_status, HepatocyteE_status,
                                        Tumor_1E_status, Tumor_2E_status, HepG2E_status))
  DMR.summary$EnhancerLike.T <- rowSums(cbind(Tumor_1E_status, Tumor_2E_status, HepG2E_status))
  DMR.summary$EnhancerLike.N <- rowSums(cbind(NormalE_status, Adjacent_1E_status, Adjacent_2E_status, AdultLiver.cREE_status, HepatocyteE_status))
  
  return(DMR.summary)
}
system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.1))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#18294
table(PromoterLike.score > 0)#9029
table(!is.na(promoter))#31765
table(!is.na(promoter) | PromoterLike.score > 0)#37618
table(!is.na(promoter) & PromoterLike.score > 0)#3176
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#1102 258
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#3071 4598
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#1955 494
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#8988 6857
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.2))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#13294
table(PromoterLike.score > 0)#7206
table(!is.na(promoter))#17836
table(!is.na(promoter) | PromoterLike.score > 0)#23557
table(!is.na(promoter) & PromoterLike.score > 0)#1485
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#901 193
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#2493 3619
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#1353 339
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#6899 4703
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.3))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#9716
table(PromoterLike.score > 0)#5679
table(!is.na(promoter))#9176
table(!is.na(promoter) | PromoterLike.score > 0)#14160
table(!is.na(promoter) & PromoterLike.score > 0)#695
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#765 130
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#1984 2801
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#970 209
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#5243 3294
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.4))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#7331
table(PromoterLike.score > 0)#4407
table(!is.na(promoter))#4958
table(!is.na(promoter) | PromoterLike.score > 0)#9021
table(!is.na(promoter) & PromoterLike.score > 0)#344
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#619 106
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#1563 2119
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#690 142
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#4098 2401
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.5))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#5692
table(PromoterLike.score > 0)#3457
table(!is.na(promoter))#2806
table(!is.na(promoter) | PromoterLike.score > 0)#6082
table(!is.na(promoter) & PromoterLike.score > 0)#181
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#510 89
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#1204 1654
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#548 92
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#3249 1803
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.6))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#4381
table(PromoterLike.score > 0)#2706
table(!is.na(promoter))#1708
table(!is.na(promoter) | PromoterLike.score > 0)#4317
table(!is.na(promoter) & PromoterLike.score > 0)#97
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#414 75
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#918 1299
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#407 69
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#2544 1361
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.7))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#3412
table(PromoterLike.score > 0)#2151
table(!is.na(promoter))#1089
table(!is.na(promoter) | PromoterLike.score > 0)#3180
table(!is.na(promoter) & PromoterLike.score > 0)#60
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#334 65
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#734 1018
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#312 47
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#1986 1067
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.8))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#2640
table(PromoterLike.score > 0)#1694
table(!is.na(promoter))#705
table(!is.na(promoter) | PromoterLike.score > 0)#2365
table(!is.na(promoter) & PromoterLike.score > 0)#34
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#287 50
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#574 783
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#242 35
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#1526 837
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 0.9))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#2050
table(PromoterLike.score > 0)#1352
table(!is.na(promoter))#477
table(!is.na(promoter) | PromoterLike.score > 0)#1809
table(!is.na(promoter) & PromoterLike.score > 0)#20
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#248 38
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#456 610
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#174 21
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#1215 640
detach(cREDMR.summary)

system.time(cREDMR.summary <- HighConfidentDMR(Feature.ratio = 1.0))
attach(cREDMR.summary)
table(EnhancerLike.score > 0)#1599
table(PromoterLike.score > 0)#1052
table(!is.na(promoter))#314
table(!is.na(promoter) | PromoterLike.score > 0)#1353
table(!is.na(promoter) & PromoterLike.score > 0)#13
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta > 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta > 0]))
#202 23
table(PromoterLike.T[PromoterLike.score > 0 & Median.delta < 0 ] >= 3/5*(PromoterLike.N[PromoterLike.score > 0 & Median.delta < 0]))
#356 471
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta > 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta > 0]))
#127 16
table(EnhancerLike.T[EnhancerLike.score > 0 & Median.delta < 0 ] >= 3/5*(EnhancerLike.N[EnhancerLike.score > 0 & Median.delta < 0]))
#960 490
detach(cREDMR.summary)

Cutoff = c(1:10)/10
system.time(for(i in 1:length(Cutoff))
{
  cREDMR.summary <- HighConfidentDMR(Feature.ratio = Cutoff[i])
  cREDMR.summary <- cREDMR.summary[cREDMR.summary$EnhancerLike.score > 0 | cREDMR.summary$PromoterLike.score > 0 | !is.na(cREDMR.summary$promoter),]
  write.csv(cREDMR.summary, file = paste0("cREDMR.summary_with_overlapRatio_", Cutoff[i], ".csv"), row.names=T)
})
##################################################### part 2 DMR gene regulation correlation ###############
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/HCC/Methy/Mixed/New") #51
#obtain smoothed and raw group level and perBase level methylation level for each PromoterLike pre-DMR
#smoothed methylation level
load.lbzip2("./NewSmoothed.Methylation.qc.RDataFS", n.cores=24)
load.lbzip2("PromoterLike_preDMR.summary.RDataFS", n.cores=24)
load.lbzip2("SigDMR.methy.matrix.smooth.RDataFS", n.cores=24)
preDMR.promoter <- PromoterLike.preDMR[!(rownames(PromoterLike.preDMR) %in% colnames(SigDMR.methy.matrix)), 1:3]
getDMRMethy_Smooth <- function(DMR, methy.sm=Meth.all.qc, what=c("perBase", "perRegion"))
{
  chr = DMR[1]
  start = as.integer(DMR[2])
  end = as.integer(DMR[3])
  possible_CpGs <- paste(chr, start:end, sep="_")
  avail_CpGs <- possible_CpGs[possible_CpGs %in% colnames(Meth.all.qc)]
  methy.matrix <- Meth.all.qc[,avail_CpGs]
  if(what == "perBase")
    return(methy.matrix)
  else return(rowMeans(methy.matrix))
}

#get group level MethyMatrix
library(foreach)
library(doParallel)
chunkSize = 1000
chunkNum = ceiling(nrow(preDMR.promoter)/chunkSize)-1

preDMR.promoter.methy.matrix <- foreach(j = 1:chunkNum, .combine = 'cbind')%do%
              {
                registerDoParallel(cl = 24, cores = 4)
                preDMR.promoter.methy.matrix_chunk <- foreach(i = (1+(j-1)*chunkSize):(chunkSize*j), .combine = 'cbind')%dopar%{
                  getDMRMethy_Smooth(DMR=preDMR.promoter[i,], methy.sm = Meth.all.qc, what = "perRegion")
                }
                stopImplicitCluster()
                if(ncol(preDMR.promoter.methy.matrix_chunk) == chunkNum)
                {
                  message(paste("Finished the ", j, "th chunk", sep=""))
                  preDMR.promoter.methy.matrix_chunk
                }
              }

registerDoParallel(cl = 24, cores = 4)
preDMR.promoter.methy.matrix_chunk <- foreach(i = (1+chunkNum*chunkSize):nrow(preDMR.promoter), .combine = 'cbind')%odpar%{
  getDMRMethy_Smooth(DMR=preDMR.promoter[i,], methy.sm = Meth.all.qc, what = "perRegion")
}
stopImplicitCluster()
if(ncol(preDMR.promoter.methy.matrix_chunk) == (nrow(preDMR.promoter)- chunkNum*chunkSize))
{
  preDMR.promoter.methy.matrix <- cbind(preDMR.promoter.methy.matrix, preDMR.promoter.methy.matrix_chunk)
}

colnames(preDMR.promoter.methy.matrix) <- rownames(preDMR.promoter)
rownames(preDMR.promoter.methy.matrix) <- rownames(Meth.all.qc)
save.lbzip2(preDMR.promoter.methy.matrix, file="preDMR.promoter.Smoothed_MethyMatrix.RDataFS", n.cores=24)




options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/Methy/Mixed/New") #29
#obtain smoothed and raw group level and perBase level methylation level for each PromoterLike pre-DMR
#smoothed methylation level
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores=24)
load.lbzip2("EnhancerLike_preDMR.summary.RDataFS", n.cores=24)
load.lbzip2("SigDMR.methy.matrix.smooth.RDataFS", n.cores=24)
preDMR.enhancer <- EnhancerLike.preDMR[!(rownames(EnhancerLike.preDMR) %in% colnames(SigDMR.methy.matrix)), 1:3]
getDMRMethy_Smooth <- function(DMR, methy.sm=Meth.all.qc, what=c("perBase", "perRegion"))
{
  chr = DMR[1]
  start = as.integer(DMR[2])
  end = as.integer(DMR[3])
  possible_CpGs <- paste(chr, start:end, sep="_")
  avail_CpGs <- possible_CpGs[possible_CpGs %in% colnames(Meth.all.qc)]
  methy.matrix <- Meth.all.qc[,avail_CpGs]
  if(what == "perBase")
    return(methy.matrix)
  else return(rowMeans(methy.matrix))
}

#get group level MethyMatrix
library(foreach)
library(doParallel)
chunkSize = 1000
chunkNum = ceiling(nrow(preDMR.enhancer)/chunkSize)-1

preDMR.enhancer.methy.matrix <- foreach(j = 1:chunkNum, .combine = 'cbind')%do%
  {
    registerDoParallel(cl = 24, cores = 4)
    preDMR.enhancer.methy.matrix_chunk <- foreach(i = (1+(j-1)*chunkSize):(chunkSize*j), .combine = 'cbind')%dopar%{
      getDMRMethy_Smooth(DMR=preDMR.enhancer[i,], methy.sm = Meth.all.qc, what = "perRegion")
    }
    stopImplicitCluster()
    if(ncol(preDMR.enhancer.methy.matrix_chunk) == chunkNum)
    {
      message(paste("Finished the ", j, "th chunk", sep=""))
      preDMR.enhancer.methy.matrix_chunk
    }
  }

registerDoParallel(cl = 24, cores = 4)
preDMR.enhancer.methy.matrix_chunk <- foreach(i = (1+chunkNum*chunkSize):nrow(preDMR.enhancer), .combine = 'cbind')%dopar%{
  getDMRMethy_Smooth(DMR=preDMR.enhancer[i,], methy.sm = Meth.all.qc, what = "perRegion")
}
stopImplicitCluster()
if(ncol(preDMR.enhancer.methy.matrix_chunk) == (nrow(preDMR.enhancer)- chunkNum*chunkSize))
{
  preDMR.enhancer.methy.matrix <- cbind(preDMR.enhancer.methy.matrix, preDMR.enhancer.methy.matrix_chunk)
}

colnames(preDMR.enhancer.methy.matrix) <- rownames(preDMR.enhancer)
rownames(preDMR.enhancer.methy.matrix) <- rownames(Meth.all.qc)
save.lbzip2(preDMR.enhancer.methy.matrix, file="preDMR.enhancer.Smoothed_MethyMatrix.RDataFS", n.cores=24)

############################################## Annotated Promoter/Enhancer Like DMR correlated genes #########################
#SigDMR.promoter
options(stringsAsFactors = F)
library(fastSave)
library(GenomicFeatures)
library(genomation)

setwd("/data/huangp/HCC/Methy/Mixed/New") #51
#obtain smoothed and raw group level and perBase level methylation level for each PromoterLike pre-DMR
#smoothed methylation level
load.lbzip2("PromoterLike_preDMR.summary.RDataFS", n.cores=24)
load.lbzip2("/data/huangp/HCC/Methy/Mixed/New/47.5k_DMR.promoter.methy.matrix.RDataFS", n.cores=24)
#load gene expr
load.lbzip2("/data/huangp/HCC/Methy/Mixed/New/Gene.TPM.filtered.RDataFS", n.cores=24)
Granges.gtf <- gffToGRanges("/data/huangp/HCC/Methy/Mixed/New/gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
Granges.gene <- Granges.gene[Granges.gene$gene_id %in% rownames(gene.TPM),]
range_gene_DF <- as.data.frame(Granges.gene)

library(foreach)
library(doParallel)

Distance <- function(DMR_str,DMR_end,TSS)
{
  if(TSS < DMR_str)
  {
    dis <- DMR_str - TSS
  } else if(TSS < DMR_end)
  {
    dis <- 0
  } else dis <- DMR_end-TSS
  return(dis)
}
COMBINE <- function(x, sep="|")
{
  if(length(x)!=0)
  {
    y <- x[1]
    if(length(x)>1)
    {
      for(i in 2:length(x)) y <- paste(y, x[i], sep = sep)
    }
    return(y)
  }
  else return(y=NA)
}

Promoter_GeneCorr <- function(DMR_location, distance = 2, cutoff = 0.01, AdjustMethod = "bonferroni")
{
  #DMR_location = rownames(DMR.Promoter)[1]
  DMR_chr <- strsplit(DMR_location, "_")[[1]][1]
  DMR_start <- as.integer(strsplit(DMR_location, "_")[[1]][2])
  DMR_end <- as.integer(strsplit(DMR_location, "_")[[1]][3])
  ## filter only genes in the same chromosome
  nearby_gene_DF <- range_gene_DF[range_gene_DF$seqnames==DMR_chr,]
  Distance_matrix <- foreach(j = 1:nrow(nearby_gene_DF), .combine = 'c') %do%
    {
      if(nearby_gene_DF$strand[j] == "+")
        Distance(DMR_start, DMR_end, nearby_gene_DF$start[j])
      else  Distance(DMR_start, DMR_end, nearby_gene_DF$end[j])
    }
  names(Distance_matrix) <- nearby_gene_DF$gene_id
  cisDistance <- Distance_matrix[abs(Distance_matrix) <= 1000*distance]
  if(length(cisDistance)>0)
  {
    if(length(cisDistance)==1) cisGene <- t(as.matrix(gene.TPM[names(cisDistance),])) else
      cisGene <- gene.TPM[names(cisDistance),]
    rownames(cisGene) <- names(cisDistance)
    
    expr_DMR_cor <- vector(mode="numeric", length=nrow(cisGene))
    expr_DMR_pvalue <- vector(mode="numeric", length=nrow(cisGene))
    
    for(g in 1:nrow(cisGene))
    {
      Expr = log(unlist(cisGene[g,])+0.001, 2)
      Methy = unlist(DMR.promoter.methy.matrix[,DMR_location])
      Methy = Methy[names(Expr)]
      expr_DMR_cor[g] <- round(unlist(cor.test(Expr, Methy, method = "spearman")$estimate),digits = 2)
      expr_DMR_pvalue[g] <- signif(unlist(cor.test(Expr, Methy, method = "spearman")$p.value),digits = 2)
    }
    expr_DMR_correlation <- data.frame(rho = expr_DMR_cor, pvalue = expr_DMR_pvalue, padj = signif(p.adjust(expr_DMR_pvalue, method = AdjustMethod), digits =2))
    rownames(expr_DMR_correlation) <- rownames(cisGene)
    expr_DMR_correlation_sorted <- expr_DMR_correlation[order(-abs(expr_DMR_correlation$rho)), ]
    
    expr_DMR_correlation_sum <- data.frame(DMR.ID = rep(DMR_location, nrow(expr_DMR_correlation_sorted)),
                                           gene_name = nearby_gene_DF$gene_name[match(rownames(expr_DMR_correlation_sorted), nearby_gene_DF$gene_id)], 
                                           distance = cisDistance[rownames(expr_DMR_correlation_sorted)], expr_DMR_correlation_sorted)
    expr_DMR_correlation_sum <- expr_DMR_correlation_sum[order(abs(expr_DMR_correlation_sum$distance)),]
    
    if(sum(expr_DMR_correlation_sum$padj < cutoff) >0) return(expr_DMR_correlation_sum[expr_DMR_correlation_sum$padj < cutoff,]) 
  } 
}
library(foreach)
library(doParallel)

  registerDoParallel(cl = 16, cores=2)
  system.time(DMR.promoter_Gene_pairs <- foreach(i = 1:ncol(DMR.promoter.methy.matrix), .combine = rbind) %dopar%
                {
                  Promoter_GeneCorr(DMR_location = colnames(DMR.promoter.methy.matrix)[i])
                })
  stopImplicitCluster()
  str(DMR.promoter_Gene_pairs)
  save.lbzip2(DMR.promoter_Gene_pairs, file="DMR.promoter_Gene_pairs.RDataFS", n.cores=24)


#included DMR information
load.lbzip2("PromoterLike_preDMR.summary.RDataFS", n.cores=24)
DMR.promoter_Gene_pairs <- cbind(PromoterLike.preDMR[DMR.promoter_Gene_pairs$DMR.ID, c("Pvalue", "Median.delta", "CpGIsland", "Promoter", "PromoterLike.score")], 
                             DMR.promoter_Gene_pairs)
write.csv(DMR.promoter_Gene_pairs, file="DMR.promoterLike_DMR_Gene_pairs_info.csv", row.names=F, col.names=T, quote=T)
table(Hypo = DMR.promoter_Gene_pairs$Median.delta < 0,  Classic = DMR.promoter_Gene_pairs$rho <0)
#           Classic
#Hypo    FALSE TRUE
#FALSE   190  125
#TRUE    745 3082


options(stringsAsFactors = F)
library(fastSave)
library(GenomicFeatures)
library(genomation)

setwd("/data/huangp/HCC/Methy/Mixed/New") #51
#obtain smoothed and raw group level and perBase level methylation level for each EnhancerLike pre-DMR
#smoothed methylation level
load.lbzip2("EnhancerLike_preDMR.summary.RDataFS", n.cores=24)
load.lbzip2("/data/huangp/HCC/Methy/Mixed/New23.8k_DMR.enhancer.methy.matrix.RDataFS", n.cores=24)
#load gene expr
load.lbzip2("Gene.TPM.filtered.RDataFS", n.cores=24)
Granges.gtf <- gffToGRanges("gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
Granges.gene <- Granges.gene[Granges.gene$gene_id %in% rownames(gene.TPM),]
range_gene_DF <- as.data.frame(Granges.gene)

library(foreach)
library(doParallel)

Distance <- function(DMR_str,DMR_end,TSS)
{
  if(TSS < DMR_str)
  {
    dis <- DMR_str - TSS
  } else if(TSS < DMR_end)
  {
    dis <- 0
  } else dis <- DMR_end-TSS
  return(dis)
}
COMBINE <- function(x, sep="|")
{
  if(length(x)!=0)
  {
    y <- x[1]
    if(length(x)>1)
    {
      for(i in 2:length(x)) y <- paste(y, x[i], sep = sep)
    }
    return(y)
  }
  else return(y=NA)
}

Enhancer_GeneCorr <- function(DMR_location, distance = 100, cutoff = 0.01, AdjustMethod = "bonferroni")
{
  #DMR_location = rownames(DMR.enhancer)[1]
  DMR_chr <- strsplit(DMR_location, "_")[[1]][1]
  DMR_start <- as.integer(strsplit(DMR_location, "_")[[1]][2])
  DMR_end <- as.integer(strsplit(DMR_location, "_")[[1]][3])
  ## filter only genes in the same chromosome
  nearby_gene_DF <- range_gene_DF[range_gene_DF$seqnames==DMR_chr,]
  Distance_matrix <- foreach(j = 1:nrow(nearby_gene_DF), .combine = 'c') %do%
    {
      if(nearby_gene_DF$strand[j] == "+")
        Distance(DMR_start, DMR_end, nearby_gene_DF$start[j])
      else  Distance(DMR_start, DMR_end, nearby_gene_DF$end[j])
    }
  names(Distance_matrix) <- nearby_gene_DF$gene_id
  cisDistance <- Distance_matrix[abs(Distance_matrix) <= 1000*distance]
  if(length(cisDistance)>0)
  {
    if(length(cisDistance)==1) cisGene <- t(as.matrix(gene.TPM[names(cisDistance),])) else
      cisGene <- gene.TPM[names(cisDistance),]
    rownames(cisGene) <- names(cisDistance)
    
    expr_DMR_cor <- vector(mode="numeric", length=nrow(cisGene))
    expr_DMR_pvalue <- vector(mode="numeric", length=nrow(cisGene))
    
    for(g in 1:nrow(cisGene))
    {
      Expr = log(unlist(cisGene[g,])+0.001, 2)
      Methy = unlist(DMR.enhancer.methy.matrix[,DMR_location])
      Methy = Methy[names(Expr)]
      expr_DMR_cor[g] <- round(unlist(cor.test(Expr, Methy, method = "spearman")$estimate),digits = 2)
      expr_DMR_pvalue[g] <- signif(unlist(cor.test(Expr, Methy, method = "spearman")$p.value),digits = 2)
    }
    expr_DMR_correlation <- data.frame(rho = expr_DMR_cor, pvalue = expr_DMR_pvalue, padj = signif(p.adjust(expr_DMR_pvalue, method = AdjustMethod), digits =2))
    rownames(expr_DMR_correlation) <- rownames(cisGene)
    expr_DMR_correlation_sorted <- expr_DMR_correlation[order(-abs(expr_DMR_correlation$rho)), ]
    
    expr_DMR_correlation_sum <- data.frame(DMR.ID = rep(DMR_location, nrow(expr_DMR_correlation_sorted)),
                                           gene_name = nearby_gene_DF$gene_name[match(rownames(expr_DMR_correlation_sorted), nearby_gene_DF$gene_id)], 
                                           distance = cisDistance[rownames(expr_DMR_correlation_sorted)], expr_DMR_correlation_sorted)
    expr_DMR_correlation_sum <- expr_DMR_correlation_sum[order(abs(expr_DMR_correlation_sum$distance)),]
    
    if(sum(expr_DMR_correlation_sum$padj < cutoff) >0) return(expr_DMR_correlation_sum[expr_DMR_correlation_sum$padj < cutoff,]) 
  } 
}

library(foreach)
library(doParallel)

registerDoParallel(cl = 16, cores=2)
system.time(DMR.enhancer_Gene_pairs <- foreach(i = 1:ncol(DMR.enhancer.methy.matrix), .combine = rbind) %dopar%
              {
                Enhancer_GeneCorr(DMR_location = colnames(DMR.enhancer.methy.matrix)[i])
              })
stopImplicitCluster()
str(DMR.enhancer_Gene_pairs)
save.lbzip2(DMR.enhancer_Gene_pairs, file="DMR.enhancer_Gene_pairs.RDataFS", n.cores=24)

#included DMR information
load.lbzip2("EnhancerLike_preDMR.summary.RDataFS", n.cores=24)
DMR.enhancer_Gene_pairs <- cbind(EnhancerLike.preDMR[DMR.enhancer_Gene_pairs$DMR.ID, c("Pvalue", "Median.delta", "EnhancerLike.score")], 
                                 DMR.enhancer_Gene_pairs)
write.csv(DMR.enhancer_Gene_pairs, file="DMR.enhancerLike_DMR_Gene_pairs_info.csv", row.names=F, col.names=T, quote=T)
table(Hypo = DMR.enhancer_Gene_pairs$Median.delta < 0,  Classic = DMR.enhancer_Gene_pairs$rho <0)
#         Classic
#Hypo    FALSE  TRUE
#FALSE  2537   550
#TRUE   2945 13621

############################################################### Novel active intergenic enhancer identification ######################################################
# 1. bedfiles construction
#load 
setwd("/data/huangp/HCC/Methy/Mixed/New")
library(fastSave)
load.lbzip2("preDMR.summary.RDataFS", n.cores=24)
intergenicDMR <- preDMR.summary[preDMR.summary$Intergenic,]
rm(preDMR.summary)
intergenicDMR.bed <- intergenicDMR[,c("chr", "str", "end")]
intergenicDMR.bed$DMRID <- rownames(intergenicDMR.bed)
write.table(intergenicDMR.bed, file="Intergenic.preDMR.hg38.bed", sep="\t", row.names=F, col.names=F, quote=F)
# build the shell command for eRNA profling
samplename <- substr(x = list.files(path="/data/huangp/HCC/STAR/", pattern = "RAligned.out.bam"),
                     start = 1, stop = 9)
#bedtools coverage -mean -a Intergenic.preDMR.hg38.bed -b /data/huangp/HCC/STAR/WGC066465RAligned.out.bam >meandepth.DMR.intergenic.WGC066465

for(i in 1:length(samplename))
{
  message(paste("bedtools coverage -mean -a Intergenic.preDMR.hg38.bed -b /data/huangp/HCC/STAR/", 
                samplename[i], "RAligned.out.bam >meandepth.DMR.intergenic.", samplename[i], " &",  sep=""))
}

## load the eRNA expression of each intergenic DMR
library(fastSave)
library(foreach)
fileNames <- list.files(path = "./eRNA", pattern = "meandepth.DMR.intergenic.")
sampleNames <- substr(fileNames, start = 26, stop=34)
depth.intergenicDMR.matrix <- foreach(i = 1:length(fileNames), .combine = 'rbind')%do%{
  read.table(paste("./eRNA/", fileNames[i],sep=""), header=F, stringsAsFactors = F)[[5]]
}
load("/data/huangp/HCC/Methy/Mixed/New/eRNA/Patient_Info.RData")
rownames(depth.intergenicDMR.matrix) <- Patient_ID$Original_Sample_Name[match(sampleNames, Patient_ID$WGC_ID)]
colnames(depth.intergenicDMR.matrix) <- read.table(paste("./eRNA/", fileNames[1], sep=""), header=F, stringsAsFactors = F)[[4]]
save.lbzip2(depth.intergenicDMR.matrix, file="eRNA.expression.matrix.intergenic.preDMR.RDataFS", n.cores = 32)
### filter the depth < 3 in over 50% samples in both tumor and adjacent
NA_count <- function(X, cutoff = 3)
{
  Y <- sum(X < 3)
  return(Y)
}
DMR.depth.tumor.NA_count <- apply(depth.intergenicDMR.matrix[grep("T", rownames(depth.intergenicDMR.matrix), fixed = T),], 2, NA_count)
DMR.depth.adjacent.NA_count <- apply(depth.intergenicDMR.matrix[grep("N", rownames(depth.intergenicDMR.matrix), fixed = T),], 2, NA_count)
#DMR.depth.matrix <- depth.intergenicDMR.matrix[, !(DMR.depth.tumor.NA_count > 16 & DMR.depth.adjacent.NA_count > 16)]
DMR.depth.matrix <- depth.intergenicDMR.matrix[, !(DMR.depth.tumor.NA_count > 22 & DMR.depth.adjacent.NA_count > 22)]
# depth normalization to FPM for downstream correlation test
depth <- Patient_ID$MappedDepth[match(rownames(DMR.depth.matrix), Patient_ID$Original_Sample_Name)]
library(foreach)
library(doParallel)
DMR.FPM.matrix <- foreach(i = 1:ncol(DMR.depth.matrix), .combine = 'cbind')%do%
  {
    DMR.depth.matrix[,i]/depth
  }
str(DMR.FPM.matrix)
colnames(DMR.FPM.matrix) <- colnames(DMR.depth.matrix)
rownames(DMR.FPM.matrix) <- rownames(DMR.depth.matrix)

## log2 transform of FPM
DMR.logFPM.matrix <- foreach(i = 1:nrow(DMR.FPM.matrix), .combine = 'rbind')%do%
  {
    log(DMR.FPM.matrix[i,]+0.0001, 2)
  }
rownames(DMR.logFPM.matrix) <- rownames(DMR.FPM.matrix)
colnames(DMR.logFPM.matrix) <- colnames(DMR.FPM.matrix)
save.lbzip2(DMR.logFPM.matrix, file="37k_active_Intergenic.preDMR.eRNA.logFPM+0.0001.matrix.RDataFS", n.cores=8)

load.lbzip2("EnhancerLike_preDMR.summary.RDataFS", n.cores=24)
table(colnames(DMR.logFPM.matrix) %in% rownames(EnhancerLike.preDMR))
#FALSE  TRUE
#34582  2378

################## Enhancer target gene identification of 36,960 active intergenic enhancer like DMRs #########################################
#### methylation level estimation for those 36,960 intergenic active enhancer-like DMRs
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/HCC/Methy/Mixed/New") #51
load.lbzip2("./NewSmoothed.Methylation.qc.RDataFS", n.cores=24)
load.lbzip2("37k_active_Intergenic.preDMR.eRNA.logFPM+0.0001.matrix.RDataFS", n.cores=24)
IntergenicDMR.eRNA <- data.frame(chr = sapply(strsplit(colnames(DMR.logFPM.matrix), split="_", fixed=T), "[", 1),
                                 start = as.integer(sapply(strsplit(colnames(DMR.logFPM.matrix), split="_", fixed=T), "[", 2)),
                                 end = as.integer(sapply(strsplit(colnames(DMR.logFPM.matrix), split="_", fixed=T), "[", 3)))
getDMRMethy_Smooth <- function(DMR, methy.sm=Meth.all.qc, what=c("perBase", "perRegion"))
{
  chr = DMR[1]
  start = as.integer(DMR[2])
  end = as.integer(DMR[3])
  possible_CpGs <- paste(chr, start:end, sep="_")
  avail_CpGs <- possible_CpGs[possible_CpGs %in% colnames(Meth.all.qc)]
  methy.matrix <- Meth.all.qc[,avail_CpGs]
  if(what == "perBase")
    return(methy.matrix)
  else return(rowMeans(methy.matrix))
}

#get group level MethyMatrix
library(foreach)
library(doParallel)
chunkSize = 1000
chunkNum = ceiling(nrow(IntergenicDMR.eRNA)/chunkSize)-1
for(j in 1:chunkNum)
{
  registerDoParallel(cl = 16, cores = 2)
  system.time(IntergenicDMR.eRNA.methy.matrix_chunk <- foreach(i = (1+(j-1)*chunkSize):(chunkSize*j), .combine = 'cbind')%dopar%{
    getDMRMethy_Smooth(DMR=IntergenicDMR.eRNA[i,], methy.sm = Meth.all.qc, what = "perRegion")
  })
  stopImplicitCluster()
  
  if(ncol(IntergenicDMR.eRNA.methy.matrix_chunk) == chunkSize)
  {
    message(paste("Finished the ", j, "th chunk", sep=""))
    save.lbzip2(IntergenicDMR.eRNA.methy.matrix_chunk, file=paste("IntergenicDMR.eRNA.methy.matrix_chunk_", j, ".RDataFS", sep=""), n.cores=24)
  }
}
registerDoParallel(cl = 16, cores = 2)

IntergenicDMR.eRNA.methy.matrix_chunk <- foreach(i = (1+chunkNum*chunkSize):nrow(IntergenicDMR.eRNA), .combine = 'cbind')%dopar%{
  getDMRMethy_Smooth(DMR=IntergenicDMR.eRNA[i,], methy.sm = Meth.all.qc, what = "perRegion")
}
stopImplicitCluster()
if(ncol(IntergenicDMR.eRNA.methy.matrix_chunk) == (nrow(IntergenicDMR.eRNA)- chunkNum*chunkSize))
{
  message("Finished the last chunk")
  save.lbzip2(IntergenicDMR.eRNA.methy.matrix_chunk, file="IntergenicDMR.eRNA.methy.matrix_chunk_last.RDataFS", n.cores=24)
}
rm(list=ls())
q(save="no")

### calculation of correlation between DMR methylation and eRNA expression
setwd("/data/huangp/HCC/Methy/Mixed/New")
options(stringsAsFactors = F)
library(fastSave)
#load intergenic DMR methylation and eRNA expression
load.lbzip2("/data/huangp/HCC/Methy/Mixed/New/IntergenicDMR.eRNA.methy.matrix.RDataFS", n.cores=24)
load.lbzip2("37k_active_Intergenic.preDMR.eRNA.logFPM+0.0001.matrix.RDataFS", n.cores=24)
library(foreach)
library(doParallel)
registerDoParallel(cl = 16, cores = 2)
system.time(IntergenicDMR.eRNA.methy.correlation <- foreach(i = 1:ncol(DMR.logFPM.matrix), .combine = 'rbind')%dopar%
  {
    methy <- unlist(IntergenicDMR.eRNA.methy.matrix[,i])
    eRNA <- unlist(DMR.logFPM.matrix[,i])
    eRNA <- eRNA[names(methy)]
    rho <- round(unlist(cor.test(methy, eRNA, method = "spearman")$estimate),digits = 2)
    pvalue <- signif(unlist(cor.test(methy, eRNA, method = "spearman")$p.value),digits = 2)
    c(rho, pvalue)
  })
stopImplicitCluster()
str(IntergenicDMR.eRNA.methy.correlation)
colnames(IntergenicDMR.eRNA.methy.correlation) <- c("eRNAandMethy.rho", "eRNAandMethy.pvalue")
rownames(IntergenicDMR.eRNA.methy.correlation) <- colnames(DMR.logFPM.matrix) 
IntergenicDMR.eRNA.methy.correlation <- as.data.frame(IntergenicDMR.eRNA.methy.correlation)
IntergenicDMR.eRNA.methy.correlation$eRNAandMethy.padj <- p.adjust(IntergenicDMR.eRNA.methy.correlation$eRNAandMethy.pvalue,
                                                                   method = "fdr")
#combine the DMR annotation information
load.lbzip2("preDMRs.summary.RDataFS", n.cores=24)
IntergenicDMR.eRNA.methy.correlation <- data.frame(preDMR.summary[rownames(IntergenicDMR.eRNA.methy.correlation),c("Median.delta", "Pvalue", "EnhancerLike.N", "EnhancerLike.T", "EnhancerLike.score")],
                                                   IntergenicDMR.eRNA.methy.correlation)
save.lbzip2(IntergenicDMR.eRNA.methy.correlation, file="IntergenicDMR.eRNAandMethy.correlation.RDataFS", n.cores=24)
table(Sig = IntergenicDMR.eRNA.methy.correlation$eRNAandMethy.padj < 0.05, 
      Classic = IntergenicDMR.eRNA.methy.correlation$eRNAandMethy.rho < 0)
#         Classic
#Sig     FALSE  TRUE
#FALSE  9215 20711
#TRUE   2178  4856
################## eRNA target gene identification of 37k active intergenic enhancer like DMRs #############################################
setwd("/data/huangp/Methy/Mixed/New")
options(stringsAsFactors = F)
library(fastSave)
library(foreach)
library(doParallel)
library(genomation)
library(GenomicFeatures)
options(stringsAsFactors = F)
#load gene expression
load.lbzip2("/data/huangp/Methy/Mixed/Gene.TPM.filtered.RDataFS", n.cores=32)
Granges.gtf <- gffToGRanges("../gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
Granges.gene <- Granges.gene[match(rownames(gene.TPM), Granges.gene$gene_id),]
range_gene_DF <- as.data.frame(Granges.gene)

#load eRNA normalized expression logFPM 
load.lbzip2("37k_active_Intergenic.preDMR.eRNA.logFPM+0.0001.matrix.RDataFS", n.cores=24)

Distance <- function(DMR_str,DMR_end,TSS)
{
  if(TSS < DMR_str)
  {
    dis <- DMR_str - TSS
  } else if(TSS < DMR_end)
  {
    dis <- 0
  } else dis <- DMR_end-TSS
  return(dis)
}
COMBINE <- function(x, sep="|")
{
  if(length(x)!=0)
  {
    y <- x[1]
    if(length(x)>1)
    {
      for(i in 2:length(x)) y <- paste(y, x[i], sep = sep)
    }
    return(y)
  }
  else return(y=NA)
}
eRNA_GeneCorr <- function(DMR_location, distance = 100, cutoff = 0.01, AdjustMethod = "bonferroni")
{
  #DMR_location = rownames(DMR.enhancer)[1]
  DMR_chr <- strsplit(DMR_location, "_")[[1]][1]
  DMR_start <- as.integer(strsplit(DMR_location, "_")[[1]][2])
  DMR_end <- as.integer(strsplit(DMR_location, "_")[[1]][3])
  ## filter only genes in the same chromosome
  nearby_gene_DF <- range_gene_DF[range_gene_DF$seqnames==DMR_chr,]
  Distance_matrix <- foreach(j = 1:nrow(nearby_gene_DF), .combine = 'c') %do%
    {
      if(nearby_gene_DF$strand[j] == "+")
        Distance(DMR_start, DMR_end, nearby_gene_DF$start[j])
      else  Distance(DMR_start, DMR_end, nearby_gene_DF$end[j])
    }
  names(Distance_matrix) <- nearby_gene_DF$gene_id
  cisDistance <- Distance_matrix[abs(Distance_matrix) <= 1000*distance]
  if(length(cisDistance)>0)
  {
    if(length(cisDistance)==1) cisGene <- t(as.matrix(gene.TPM[names(cisDistance),])) else
      cisGene <- gene.TPM[names(cisDistance),]
    rownames(cisGene) <- names(cisDistance)
    
    expr_eRNA_cor <- vector(mode="numeric", length=nrow(cisGene))
    expr_eRNA_pvalue <- vector(mode="numeric", length=nrow(cisGene))
    
    for(g in 1:nrow(cisGene))
    {
      Expr = log(unlist(cisGene[g,])+0.001, 2)
      eRNA = unlist(DMR.logFPM.matrix[,DMR_location])
      eRNA = eRNA[names(Expr)]
      expr_eRNA_cor[g] <- round(unlist(cor.test(Expr, eRNA, method = "spearman")$estimate),digits = 2)
      expr_eRNA_pvalue[g] <- signif(unlist(cor.test(Expr, eRNA, method = "spearman")$p.value),digits = 2)
    }
    expr_eRNA_correlation <- data.frame(rho = expr_eRNA_cor, pvalue = expr_eRNA_pvalue, padj = signif(p.adjust(expr_eRNA_pvalue, method = AdjustMethod), digits =2))
    rownames(expr_eRNA_correlation) <- rownames(cisGene)
    expr_eRNA_correlation_sorted <- expr_eRNA_correlation[order(-abs(expr_eRNA_correlation$rho)), ]
    
    expr_eRNA_correlation_sum <- data.frame(DMR.ID = rep(DMR_location, nrow(expr_eRNA_correlation_sorted)),
                                           gene_name = nearby_gene_DF$gene_name[match(rownames(expr_eRNA_correlation_sorted), nearby_gene_DF$gene_id)], 
                                           distance = cisDistance[rownames(expr_eRNA_correlation_sorted)], expr_eRNA_correlation_sorted)
    expr_eRNA_correlation_sum <- expr_eRNA_correlation_sum[order(abs(expr_eRNA_correlation_sum$distance)),]
    
    if(sum(expr_eRNA_correlation_sum$padj < cutoff) >0) return(expr_eRNA_correlation_sum[expr_eRNA_correlation_sum$padj < cutoff,]) 
  } 
}

registerDoParallel(cl = 16, cores=2)
system.time(DMReRNA_Gene_pairs <- foreach(i = 1:ncol(DMR.logFPM.matrix), .combine = rbind) %dopar%
                {
                  eRNA_GeneCorr(DMR_location = colnames(DMR.logFPM.matrix)[i])
                })
stopImplicitCluster()

str(DMReRNA_Gene_pairs)

DMReRNA_Gene_pairs <- cbind(DMR.chr = sapply(strsplit(DMReRNA_Gene_pairs$eRNAID, split = "_", fixed=T), "[", 1),
                            DMR.start = as.integer(sapply(strsplit(DMReRNA_Gene_pairs$eRNAID, split = "_", fixed=T), "[", 2)),
                            DMR.end = as.integer(sapply(strsplit(DMReRNA_Gene_pairs$eRNAID, split = "_", fixed=T), "[", 3)),
                            DMReRNA_Gene_pairs)
str(DMReRNA_Gene_pairs)
save.lbzip2(DMReRNA_Gene_pairs, file="IntergenicDMR.eRNA.Gene.pairs.RDataFS", n.cores=32)

################## overlap of eRNA regulated by DMRs and eRNA target genes ###############################
options(stringsAsFactors = F)
setwd("/data/huangp/HCC/Methy/Mixed/New")
library(fastSave)
# load methy and eRNA correlation
load.lbzip2("IntergenicDMR.eRNAandMethy.correlation.RDataFS", n.cores=24)
eRNARegulatedDMRs <- IntergenicDMR.eRNA.methy.correlation[IntergenicDMR.eRNA.methy.correlation$eRNAandMethy.padj < 0.05,]
#load eRNA correlated nearby genes
load.lbzip2("IntergenicDMR.eRNA.Gene.pairs.RDataFS", n.cores=24)
eRNARegulatedDMRs_Gene_pairs <- DMReRNA_Gene_pairs[DMReRNA_Gene_pairs$DMR.ID %in% rownames(eRNARegulatedDMRs),]
eRNARegulatedDMRs_Gene_pairs <- cbind(eRNARegulatedDMRs[eRNARegulatedDMRs_Gene_pairs$DMR.ID,],
                                      eRNARegulatedDMRs_Gene_pairs)
colnames(eRNARegulatedDMRs_Gene_pairs) <- c("DM.Median", "DM.Pvalue", "EnhancerScore.N", "EnhancerScore.T", "EnhancerScore",
                                            "eRNAandMethy.rho", "eRNAandMethy.pvalue", "eRNAandMethy.padj", "DMR.chr", "DMR.start", "DMR.end", "DMR.ID",
                                            "gene_name", "distanceToDMR", "eRNAandGene.rho", "eRNAandGene.pvalue", "eRNAandGene.padj")
eRNARegulatedDMRs_Gene_pairs <- eRNARegulatedDMRs_Gene_pairs[order(eRNARegulatedDMRs_Gene_pairs$DMR.chr, eRNARegulatedDMRs_Gene_pairs$DMR.start, eRNARegulatedDMRs_Gene_pairs$DMR.end, decreasing = F),]
save.lbzip2(eRNARegulatedDMRs_Gene_pairs, file="eRNARegulatedDMRs_Gene_pairs.RDataFS", n.cores=24)
### test of correlation between DMR and Gene
options(stringsAsFactors = F)
library(fastSave)
library(genomation)
library(GenomicFeatures)
setwd("/data/huangp/HCC/Methy/Mixed/New")
#load eRNA regulated DMRs
load.lbzip2("eRNARegulatedDMRs_Gene_pairs.RDataFS", n.cores=24)
eRNARegulatedDMRs <- unique(eRNARegulatedDMRs_Gene_pairs$DMR.ID)
#load gene expression
#load gene expression
load.lbzip2("./Gene.TPM.filtered.RDataFS", n.cores=32)
Granges.gtf <- gffToGRanges("./gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
Granges.gene <- Granges.gene[match(rownames(gene.TPM), Granges.gene$gene_id),]
range_gene_DF <- as.data.frame(Granges.gene)
#load DMR methylation level 
load.lbzip2("IntergenicDMR.eRNA.methy.matrix.RDataFS", n.cores=24)

library(foreach)
library(doParallel)

Distance <- function(DMR_str,DMR_end,TSS)
{
  if(TSS < DMR_str)
  {
    dis <- DMR_str - TSS
  } else if(TSS < DMR_end)
  {
    dis <- 0
  } else dis <- DMR_end-TSS
  return(dis)
}
COMBINE <- function(x, sep="|")
{
  if(length(x)!=0)
  {
    y <- x[1]
    if(length(x)>1)
    {
      for(i in 2:length(x)) y <- paste(y, x[i], sep = sep)
    }
    return(y)
  }
  else return(y=NA)
}

Enhancer_GeneCorr <- function(DMR_location, distance = 100, cutoff = 0.01, AdjustMethod = "bonferroni")
{
  #DMR_location = rownames(DMR.enhancer)[1]
  DMR_chr <- strsplit(DMR_location, "_")[[1]][1]
  DMR_start <- as.integer(strsplit(DMR_location, "_")[[1]][2])
  DMR_end <- as.integer(strsplit(DMR_location, "_")[[1]][3])
  ## filter only genes in the same chromosome
  nearby_gene_DF <- range_gene_DF[range_gene_DF$seqnames==DMR_chr,]
  Distance_matrix <- foreach(j = 1:nrow(nearby_gene_DF), .combine = 'c') %do%
    {
      if(nearby_gene_DF$strand[j] == "+")
        Distance(DMR_start, DMR_end, nearby_gene_DF$start[j])
      else  Distance(DMR_start, DMR_end, nearby_gene_DF$end[j])
    }
  names(Distance_matrix) <- nearby_gene_DF$gene_id
  cisDistance <- Distance_matrix[abs(Distance_matrix) <= 1000*distance]
  if(length(cisDistance)>0)
  {
    if(length(cisDistance)==1) cisGene <- t(as.matrix(gene.TPM[names(cisDistance),])) else
      cisGene <- gene.TPM[names(cisDistance),]
    rownames(cisGene) <- names(cisDistance)
    
    expr_DMR_cor <- vector(mode="numeric", length=nrow(cisGene))
    expr_DMR_pvalue <- vector(mode="numeric", length=nrow(cisGene))
    
    for(g in 1:nrow(cisGene))
    {
      Expr = log(unlist(cisGene[g,])+0.001, 2)
      Methy = unlist(IntergenicDMR.eRNA.methy.matrix[,DMR_location])
      Methy = Methy[names(Expr)]
      expr_DMR_cor[g] <- round(unlist(cor.test(Expr, Methy, method = "spearman")$estimate),digits = 2)
      expr_DMR_pvalue[g] <- signif(unlist(cor.test(Expr, Methy, method = "spearman")$p.value),digits = 2)
    }
    expr_DMR_correlation <- data.frame(rho = expr_DMR_cor, pvalue = expr_DMR_pvalue, padj = signif(p.adjust(expr_DMR_pvalue, method = AdjustMethod), digits =2))
    rownames(expr_DMR_correlation) <- rownames(cisGene)
    expr_DMR_correlation_sorted <- expr_DMR_correlation[order(-abs(expr_DMR_correlation$rho)), ]
    
    expr_DMR_correlation_sum <- data.frame(DMR.ID = rep(DMR_location, nrow(expr_DMR_correlation_sorted)),
                                           gene_name = nearby_gene_DF$gene_name[match(rownames(expr_DMR_correlation_sorted), nearby_gene_DF$gene_id)], 
                                           distance = cisDistance[rownames(expr_DMR_correlation_sorted)], expr_DMR_correlation_sorted)
    expr_DMR_correlation_sum <- expr_DMR_correlation_sum[order(abs(expr_DMR_correlation_sum$distance)),]
    
    if(sum(expr_DMR_correlation_sum$padj < cutoff) >0) return(expr_DMR_correlation_sum[expr_DMR_correlation_sum$padj < cutoff,]) 
  } 
}

library(foreach)
library(doParallel)

registerDoParallel(cl = 16, cores=2)
system.time(eRNARegulatedDMRs_correlated_Gene_pairs <- foreach(i = 1:length(eRNARegulatedDMRs), .combine = rbind) %dopar%
              {
                Enhancer_GeneCorr(DMR_location = eRNARegulatedDMRs[i])
              })
stopImplicitCluster()
str(eRNARegulatedDMRs_correlated_Gene_pairs)
eRNARegulatedDMRs_correlated_Gene_pairs$pairID <- paste(eRNARegulatedDMRs_correlated_Gene_pairs$DMR.ID, 
                                                        eRNARegulatedDMRs_correlated_Gene_pairs$gene_name, sep="|")
eRNARegulatedDMRs_Gene_pairs$pairID <- paste(eRNARegulatedDMRs_Gene_pairs$DMR.ID, eRNARegulatedDMRs_Gene_pairs$gene_name, sep="|")
table(eRNARegulatedDMRs_Gene_pairs$pairID %in% eRNARegulatedDMRs_correlated_Gene_pairs$pairID)
eRNARegulatedDMRs_Gene_pairs <- eRNARegulatedDMRs_Gene_pairs[eRNARegulatedDMRs_Gene_pairs$pairID %in% eRNARegulatedDMRs_correlated_Gene_pairs$pairID,]
eRNARegulatedDMRs_Gene_pairs <- cbind(eRNARegulatedDMRs_Gene_pairs, 
                                      eRNARegulatedDMRs_correlated_Gene_pairs[match(eRNARegulatedDMRs_Gene_pairs$pairID, eRNARegulatedDMRs_correlated_Gene_pairs$pairID),c("rho", "pvalue", "padj")])
colnames(eRNARegulatedDMRs_Gene_pairs)[19:21] <- c("MethyandGene.rho", "MethyandGene.pvalue", "MethyandGene.padj")
save.lbzip2(eRNARegulatedDMRs_Gene_pairs, file="Consistent_eRNARegulatedDMRs_Gene_pairs.RDataFS", n.cores=24)
write.csv(eRNARegulatedDMRs_Gene_pairs, file="Consistent_eRNARegulatedDMRs_Gene_pairs.csv", row.names = F, col.names = T, quote = T)





### correlation graph of gene and DMR methylation
DMRGeneScatterplot <- function(DMRID, GeneName)
{
  methy = unlist(DMR_methy[,DMRID])
  GeneID = rownames(range_gene_DF)[range_gene_DF$gene_name == GeneName]
  expr = log(unlist(expr_qc[GeneID,]), 2)
  cor(methy, expr, method = "spearman", use = "pairwise.complete.obs")
  library(ggpubr)
  library(ggplot2)
  data <- data.frame(Methylation = methy, GeneExpression = expr, Condition = c(rep(c("Tumor", "Adjacent"), each = 23), rep(c("Tumor", "Adjacent"), each=7)))
  colnames(data)[1] <- paste("Raw_Methylation_of_DMR_", DMRID, sep="")
  #colnames(data)[2] <-paste("Expression_of_Gene_", "CDKN2BAS1", sep="")
  colnames(data)[2] <-paste("Expression_of_Gene_", GeneName, sep="")
  ggscatterhist(data = data, x = colnames(data)[1], y = colnames(data)[2], fill = 'Condition', color = 'Condition',
                margin.plot = "boxplot", legend="top", ggtheme = theme_pubr())
  
  ggsave(paste("DMRGeneScatterplot of ", DMRID, " and ", GeneName, ".pdf",sep=""))
}



##### 5-aza treatment GEO dataset
setwd("/data/huangp/Methy/Mixed/5azaGEO/")
################# overlap of genic annotated enhancer-like DMRs and intergenic enhancer-like DMRs ###################
options(stringsAsFactors = F)
library(fastSave)
setwd("/data/huangp/HCC/Methy/Mixed/New")
#load annotated enhancer-like DMR_Gene_pairs
DMR.enhancerLike_Gene_pairs <- read.csv("DMR.enhancerLike_DMR_Gene_pairs_info.csv", header=T)
#load preDMRs.summary
load.lbzip2("preDMRs.summary.RDataFS", n.cores=24)
# filter out intergenic annotated enhancer-like DMR_Gene_pairs
DMR.enhancerLike_Gene_pairs <- cbind(preDMR.summary[DMR.enhancerLike_Gene_pairs$DMR.ID,c("Median.delta","Pvalue", "Intergenic","EnhancerLike.T", "EnhancerLike.N", "EnhancerLike.score")],
                                     DMR.enhancerLike_Gene_pairs[,-c(1:3)])
DMR.enhancerLike_Gene_pairs <- DMR.enhancerLike_Gene_pairs[!DMR.enhancerLike_Gene_pairs$Intergenic,]
DMR.enhancerLike_Gene_pairs$pairID <- paste(DMR.enhancerLike_Gene_pairs$DMR.ID, DMR.enhancerLike_Gene_pairs$gene_name, sep="|")
write.csv(DMR.enhancerLike_Gene_pairs[,-3], file="GenicDMR.enhancerLike_Gene_pairs.csv", col.names = T, row.names = F,quote=T)
table(DMR.enhancerLike_Gene_pairs$EnhancerLike.score)
#1    2    3    4    5    6    7    8
#6085 2931 2190 1430 1030  509  431  259
table(Hypo = DMR.enhancerLike_Gene_pairs$Median.delta < 0, Classic = DMR.enhancerLike_Gene_pairs$rho < 0)
#           Classic
#Hypo    FALSE  TRUE
#FALSE  2355   500
#TRUE   1996 10014
#load intergenic enhancer-like DMR_DEeRNA_DEG pairs
eRNARegulatedDMRs_Gene_pairs <- read.csv("Consistent_eRNARegulatedDMRs_Gene_pairs.csv", header=T)
table(Hypo = eRNARegulatedDMRs_Gene_pairs$DM.Median < 0, Classic = eRNARegulatedDMRs_Gene_pairs$MethyandGene.rho < 0)
#         Classic
#Hypo    FALSE TRUE
#FALSE    67   18
#TRUE   4642 4813
################# comparison of promoter-like DMR_Gene pairs and enhancer-like DMR_Gene pairs ###########
#rule out enhancer gene pairs regulated by promoter-like DMRs
#Promoter pairs > Enhancer pairs
setwd("/data/huangp/Methy/Mixed/New")
promoter_genes <- read.csv("DMR.promoterLike_Gene_pairs.csv", header=T)
table(Hyper = promoter_genes$Median.delta >0, Classic = promoter_genes$rho < 0)
#         Classic
#Hyper   FALSE TRUE
#FALSE   745 3082
#TRUE    190  125
hyper_down <- unique(promoter_genes$gene_name[promoter_genes$Median.delta > 0 & promoter_genes$rho < 0])#108
hypo_up <- unique(promoter_genes$gene_name[promoter_genes$Median.delta < 0 & promoter_genes$rho < 0])#2617
hyper_up <- unique(promoter_genes$gene_name[promoter_genes$Median.delta > 0 & promoter_genes$rho > 0])#182
hypo_down <- unique(promoter_genes$gene_name[promoter_genes$Median.delta < 0 & promoter_genes$rho > 0])#547
#filtered out promoter like DMR associated DEGs with delta methylation lower than 15%
promoter_genes <- promoter_genes[abs(promoter_genes$Median.delta) >=0.15,]
hyper_down <- unique(promoter_genes$gene_name[promoter_genes$Median.delta > 0 & promoter_genes$rho < 0])#91
hypo_up <- unique(promoter_genes$gene_name[promoter_genes$Median.delta < 0 & promoter_genes$rho < 0])#1860
hyper_up <- unique(promoter_genes$gene_name[promoter_genes$Median.delta > 0 & promoter_genes$rho > 0])#140
hypo_down <- unique(promoter_genes$gene_name[promoter_genes$Median.delta < 0 & promoter_genes$rho > 0])#478
write.csv(promoter_genes, file="FilteredDMR.promoterLike_Gene_pairs.csv", row.names=F)
genic_enhancer_genes <- read.csv("GenicDMR.enhancerLike_Gene_pairs.csv", header=T)
intergenic_enhancer_genes <- read.csv("Consistent_eRNARegulatedDMRs_Gene_pairs.csv", header=T)
table(genic_enhancer_genes$DMR.ID %in% promoter_genes$DMR.ID)#11898 2967
table(genic_enhancer_genes$gene_name %in% promoter_genes$gene_name)#9360  5505
genic_enhancer_genes_1 <- genic_enhancer_genes[!(genic_enhancer_genes$DMR.ID %in% promoter_genes$DMR.ID),]
table(genic_enhancer_genes_1$gene_name %in% promoter_genes$gene_name)#8102 3796
genic_enhancer_genes_2 <- genic_enhancer_genes_1[!(genic_enhancer_genes_1$gene_name %in% promoter_genes$gene_name),]
table(Hyper = genic_enhancer_genes_2$Median.delta >0, Classic = genic_enhancer_genes_2$rho < 0)
#         Classic
#Hyper   FALSE TRUE
#FALSE  1034 5082
#TRUE   1619  367
genic.hyper_down <- unique(genic_enhancer_genes_2$gene_name[genic_enhancer_genes_2$Median.delta > 0 & genic_enhancer_genes_2$rho < 0])#247
genic.hypo_up <- unique(genic_enhancer_genes_2$gene_name[genic_enhancer_genes_2$Median.delta < 0 & genic_enhancer_genes_2$rho < 0])#2824
genic.hyper_up <- unique(genic_enhancer_genes_2$gene_name[genic_enhancer_genes_2$Median.delta > 0 & genic_enhancer_genes_2$rho > 0])#1141
genic.hypo_down <- unique(genic_enhancer_genes_2$gene_name[genic_enhancer_genes_2$Median.delta < 0 & genic_enhancer_genes_2$rho > 0])#503

#genic_enhancer_genes_3 <- genic_enhancer_genes_2[genic_enhancer_genes_2$EnhancerLike.score>=3,]
#genic.hyper_down <- unique(genic_enhancer_genes_3$gene_name[genic_enhancer_genes_3$Median.delta > 0 & genic_enhancer_genes_3$rho < 0])#128
#genic.hypo_up <- unique(genic_enhancer_genes_3$gene_name[genic_enhancer_genes_3$Median.delta < 0 & genic_enhancer_genes_3$rho < 0])#1561
#genic.hyper_up <- unique(genic_enhancer_genes_3$gene_name[genic_enhancer_genes_3$Median.delta > 0 & genic_enhancer_genes_3$rho > 0])#506
#genic.hypo_down <- unique(genic_enhancer_genes_3$gene_name[genic_enhancer_genes_3$Median.delta < 0 & genic_enhancer_genes_3$rho > 0])#171
write.csv(genic_enhancer_genes_3, file="FilteredGenicDMR.enhancerLike_Gene_pairs.csv", row.names=F)

table(intergenic_enhancer_genes$DMR.ID %in% promoter_genes$DMR.ID)#9521 19
table(intergenic_enhancer_genes$gene_name %in% promoter_genes$gene_name)#4652 4888
intergenic_enhancer_genes_1 <- intergenic_enhancer_genes[!(intergenic_enhancer_genes$DMR.ID %in% promoter_genes$DMR.ID),]
table(intergenic_enhancer_genes_1$gene_name %in% promoter_genes$gene_name)#4643 4878
intergenic_enhancer_genes_2 <- intergenic_enhancer_genes_1[!(intergenic_enhancer_genes_1$gene_name %in% promoter_genes$gene_name),]
table(Hyper = intergenic_enhancer_genes_2$DM.Median >0, Classic = intergenic_enhancer_genes_2$MethyandGene.rho < 0)
write.csv(intergenic_enhancer_genes_2, file="FilteredConsistent_eRNARegulatedDMRs_Gene_pairs.csv", row.names=F)
#
############################ LIHC.TCGA Validation of Promoter-like DMR-Genes pairs ##################
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/LIHC")
load("LIHC_TCGA.phenotypes.RData")#load LIHC.data
methy.data <- LIHC.data
methy.data$sampleID <- as.character(methy.data$sampleID)
if(1 < 0){LIHC.survival <- read.table("LIHC.TCGA.survival.txt", header=T, sep="\t")
LIHC.survival$Sample <- sapply(strsplit(LIHC.survival$Sample, split="-01", fixed=T), "[", 1)
LIHC.survival <- LIHC.survival[,-3]
LIHC.survival$Description[223] <- "NA, female, white, NA, alive, 30 days"
LIHC.survival$Description[278] <- "88 years, female, NA, NA, dead, 1147 days"
library(foreach)
survival_phenotypes <- foreach(i = 1:nrow(LIHC.survival), .combine='rbind')%do%
  {
    des <- unlist(strsplit(LIHC.survival$Description[i], split=", ", fixed=T))
    if(length(des)==5)
    {
      if(length(grep("years", des, fixed=T))==0)
        des <- c("NA", des)
      else if(length(grep("male", des, fixed=T))==0)
        des <- c(des[1], "NA", des[2:5])
      else if(length(grep("stage", des, fixed=T))==0)
        des <- c(des[1:3], "NA", des[4:5])
      else if(length(grep("alive", des, fixed=T))==0 & length(grep("dead", des, fixed=T))==0)
        des <- c(des[1:4], "NA", des[5])
      else if(length(grep("days", des, fixed=T))==0)
        des <- c(des, "NA")
      else des <- c(des[1:2], "NA", des[3:5])
    }
    #des[1] <- strsplit(des[1], split=" years", fixed=T)[[1]][1]
    #des[4] <- strsplit(des[4], split=":", fixed=T)[[1]][2]
    #des[6] <- as.integer(strsplit(des[6], split=" days", fixed=T)[[1]][1])
    des
  }
colnames(survival_phenotypes) <- c("Age", "Gender", "Ethenic", "Tumor_stage", "survival_status", "survival_time")
survival_phenotypes <- as.data.frame(survival_phenotypes)
survival_phenotypes$Age[grep("NA", survival_phenotypes$Age, fixed=T)] <- NA
survival_phenotypes$Gender[grep("NA", survival_phenotypes$Gender, fixed=T)] <- NA
survival_phenotypes$Ethenic[grep("NA", survival_phenotypes$Ethenic, fixed=T)] <- NA
survival_phenotypes$Tumor_stage[grep("NA", survival_phenotypes$Tumor_stage, fixed=T)] <- NA
survival_phenotypes$survival_status[grep("NA", survival_phenotypes$survival_status, fixed=T)] <- NA
survival_phenotypes$survival_time[grep("NA", survival_phenotypes$survival_time, fixed=T)] <- NA

survival_phenotypes$Age <- as.integer(sapply(strsplit(survival_phenotypes$Age, split=" years", fixed=T), "[", 1))
survival_phenotypes$survival_time <- as.integer(sapply(strsplit(survival_phenotypes$survival_time, split=" days", fixed=T), "[", 1))
survival_phenotypes$Stages <- factor(survival_phenotypes$Tumor_stage, levels = c("stage:i", "stage:ii", "stage:iii", 
                                                                                 "stage:iiia", "stage:iiib", "stage:iiic", "stage:iv", "stage:iva", "stage:ivb"),
                                     labels = c("I", "II", "III", "III", "III", "III", "IV", "IV", "IV"))
survival_phenotypes <- cbind(SampleID = LIHC.survival$Sample, survival_phenotypes)
methy.data <- cbind(methy.data, survival_phenotypes[match(methy.data$barcode, survival_phenotypes$SampleID), c("survival_status", "survival_time")])
}
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
LIHC.expr <- z #20531Genes Gene_name+Gene_ID+423Sample
### for each promoter-like DMR-Gene pair, calculate their rho 
### load promoter-like DMR-Gene pair 
#Promoter_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/DMR.promoterLike_Gene_pairs.csv", header=T)
Promoter_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/FilteredDMR.promoterLike_Gene_pairs.csv")
#cREDMR.summary <- read.csv("/data/huangp/Methy/Mixed/New/cREDMR.summary_with_overlapRatio_1.csv", header=T)
#cREDMR.summary <- read.csv("/data/huangp/Methy/Mixed/New/cREDMR.summary_with_overlapRatio_0.5.csv", header=T)
##load DMR information for high confident cRE-like DMRs
cREDMR.summary <- read.csv("/data/huangp/Methy/Mixed/New/cREDMR.summary_with_overlapRatio_0.8.csv", header=T)
cREDMR.summary <- cREDMR.summary[!is.na(cREDMR.summary$promoter)|cREDMR.summary$PromoterLike.score>0,]
rownames(cREDMR.summary) <- paste(cREDMR.summary$chr, cREDMR.summary$str, cREDMR.summary$end, sep="_")
Promoter_Gene_pairs <- Promoter_Gene_pairs[Promoter_Gene_pairs$DMR.ID %in% rownames(cREDMR.summary),]#4142,2857,1263,577, 387
Promoter_Gene_pairs[,c("Promoter", "PromoterLike.T", "PromoterLike.N", "PromoterLike.score")] <- cREDMR.summary[match(Promoter_Gene_pairs$DMR.ID, rownames(cREDMR.summary)), 
                                                                                                                c("promoter", "PromoterLike.T", "PromoterLike.N", "PromoterLike.score")]
Promoter_Gene_pairs$pairID <- paste(Promoter_Gene_pairs$DMR.ID, Promoter_Gene_pairs$gene_name, sep="|")
length(unique(Promoter_Gene_pairs$gene_name[Promoter_Gene_pairs$Median.delta > 0 & Promoter_Gene_pairs$rho < 0]))#18#17
length(unique(Promoter_Gene_pairs$gene_name[Promoter_Gene_pairs$Median.delta < 0 & Promoter_Gene_pairs$rho < 0]))#305#247
length(unique(Promoter_Gene_pairs$gene_name[Promoter_Gene_pairs$Median.delta > 0 & Promoter_Gene_pairs$rho > 0]))#36#35
length(unique(Promoter_Gene_pairs$gene_name[Promoter_Gene_pairs$Median.delta < 0 & Promoter_Gene_pairs$rho > 0]))#38#36
Promoter_Gene_pairs$pairType <- rep("NA", nrow(Promoter_Gene_pairs))
Promoter_Gene_pairs$pairType[Promoter_Gene_pairs$Median.delta > 0 & Promoter_Gene_pairs$rho < 0] <- "HyperDown"
Promoter_Gene_pairs$pairType[Promoter_Gene_pairs$Median.delta < 0 & Promoter_Gene_pairs$rho < 0] <- "HypoUp"
Promoter_Gene_pairs$pairType[Promoter_Gene_pairs$Median.delta > 0 & Promoter_Gene_pairs$rho > 0] <- "HyperUp"
Promoter_Gene_pairs$pairType[Promoter_Gene_pairs$Median.delta < 0 & Promoter_Gene_pairs$rho > 0] <- "HypoDown"

## Alternative distribution of four type of validation (successful validation, type I, II, and III validation failure)
#Type I: No significant DEG;  
#Type II: Significant DEG but no 450k CpG available;
#Type III: Significant DEG and but non significant correlation between CpG methylation and gene expression.
library(foreach)
if(1<0){
  validation_distribution <- foreach(i = 1:nrow(Promoter_Gene_pairs), .combine = 'c') %do%
    {
      DMR = Promoter_Gene_pairs$DMR.ID[i]
      Gene = Promoter_Gene_pairs$gene_name[i]
      expr <- unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name), -c(1:2)])
      if(sum(is.na(expr)) > length(expr)/2 | sum(expr==0) > length(expr)/2) V_type <- "TypeI"#if low expression or no expression, type I validation failure
      else
      {
        expr <- log2(0.001+expr)
        phenotype.expr <- methy.data[match(names(expr), methy.data$sampleID),]
        DE.P <- signif(wilcox.test(x = expr[phenotype.expr$condition=="Tumor"], y = expr[phenotype.expr$condition=="Normal"])$p.value, digits = 2)
        if(DE.P > 0.05) V_type <- "TypeI"#if no significant DEG, type I validation failure
        else{
          DMR.chr = strsplit(DMR, split="_", fixed=T)[[1]][1]
          DMR.start = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][2])
          DMR.end = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][3])
          AllCpGs <- paste(DMR.chr, (DMR.start-0):(0+DMR.end), sep="_")
          if(sum(rownames(LIHC.methy) %in% AllCpGs) == 0) V_type <- "TypeII" #If there was no CpG, return typeII validation failure
          else{
            methy_matrix <- LIHC.methy[rownames(LIHC.methy) %in% AllCpGs,]
            # if there are only one CpG, methy_matrix need to be transformed into a matrix
            if(is.null(nrow(methy_matrix))) 
            {methy_matrix <- t(as.matrix(methy_matrix))
            rownames(methy_matrix) = intersect(AllCpGs, rownames(LIHC.methy))}
            #filtered out CpGs with too much NAs
            keep.CpGs <- rownames(methy_matrix)[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2]
            methy_matrix <- methy_matrix[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2, ]
            if(length(keep.CpGs) == 0) V_type <- "TypeII"#If there was no CpG after filtering, return typeII validation failure 
            else{
              # if there are only one CpG, methy_matrix need to be transformed into a matrix
              if(is.null(nrow(methy_matrix))) 
              {methy_matrix <- t(as.matrix(methy_matrix))
              rownames(methy_matrix) = keep.CpGs}
              #for each CpG, examine their correlation with gene expression
              matched_sample <- intersect(names(expr), colnames(methy_matrix))
              expr <- expr[matched_sample]
              methy_matrix <- methy_matrix[,matched_sample]#If nrow(methy_matrix)==1, it bacame a vector again!
              if(is.null(nrow(methy_matrix))) 
              {methy_matrix <- t(as.matrix(methy_matrix))
              rownames(methy_matrix) = keep.CpGs}
              phenotypes <- methy.data[match(matched_sample, methy.data$sampleID),]
              cor.p <- vector(mode="numeric", length=nrow(methy_matrix))
              DM.p <- vector(mode="numeric", length=nrow(methy_matrix))
              for(p in 1:nrow(methy_matrix))
              {
                DM.p[p] <- signif(wilcox.test(x = methy_matrix[p,phenotypes$condition == "Tumor"], y = methy_matrix[p,phenotypes$condition == "Normal"])$p.value, digits = 2)
                cor.p[p] = signif(cor.test(x = expr, y = methy_matrix[p,], method = "spearman")$p.value, digits = 2)
              }
              if(sum(cor.p <= 0.05 & DM.p <= 0.05) ==0) V_type <- "TypeIII"#if sig DEG, but no correlation between methy and expr, return typeIII
              else V_type <- "Validation"
            }
          }
        }
      }
      V_type
    }
  str(validation_distribution)
  table(validation_distribution)
}
#validation_distribution
#TypeI     TypeII    TypeIII Validation
#177        109         14        109

###validation classification
#Type I failure: No CpG
#Type II failure: There was at least one CpG and significant DEG, but no significant correlation between them
#Type III: Successful Validation 
library(foreach)
validation_distribution <- foreach(i = 1:nrow(Promoter_Gene_pairs), .combine = 'c') %do%
  {
    DMR = Promoter_Gene_pairs$DMR.ID[i]
    Gene = Promoter_Gene_pairs$gene_name[i]
    expr <- unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name), -c(1:2)])
    DMR.chr = strsplit(DMR, split="_", fixed=T)[[1]][1]
    DMR.start = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][2])
    DMR.end = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][3])
    AllCpGs <- paste(DMR.chr, (DMR.start-0):(0+DMR.end), sep="_")
    methy_matrix <- LIHC.methy[rownames(LIHC.methy) %in% AllCpGs,]
    #if only one CpG
    if(is.null(nrow(methy_matrix))){
      methy_matrix <- t(as.matrix(methy_matrix))
      rownames(methy_matrix) <- intersect(rownames(LIHC.methy), "chr3_158672843")
    }
    #filtered CpG with too much NAs
    keep.CpGs <- rownames(methy_matrix)[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2]
    methy_matrix <- methy_matrix[keep.CpGs,]
    if(length(keep.CpGs)==0) V_type <- "TypeI" #If there was no CpG, return typeI validation failure
    else
    {
      if(length(keep.CpGs)==1) 
        {methy_matrix <- t(as.matrix(methy_matrix))
         rownames(methy_matrix) <- keep.CpGs}
      #for each CpG, exam their correlation with expr
      matchedSample <- intersect(colnames(methy_matrix), names(expr))
      expr <- expr[matchedSample]
      #filtered out gene expresssion with too much 0 or NA
      if(sum(is.na(expr)) > length(expr)/2 | sum(expr == 0) > length(expr)/2) V_type <- "TypeII"
      else
      {
        #filtered out gene without significant DEG
        expr <- log2(0.001+expr)
        phenotype.expr <- methy.data[match(names(expr), methy.data$sampleID),]
        DE.P <- signif(wilcox.test(x = expr[phenotype.expr$condition=="Tumor"], y = expr[phenotype.expr$condition=="Normal"])$p.value, digits = 2)
        if(DE.P > 0.05) V_type <- "TypeII"
        else
        {
          #filtered methy_matrix for matched samples with expr data
          methy_matrix <- methy_matrix[,matchedSample]
          if(is.null(nrow(methy_matrix))){
            methy_matrix <- t(as.matrix(methy_matrix))
            rownames(methy_matrix) <- keep.CpGs
          }
          #filtered out gene without significant DM CpG
          #filtered out insignificant correlation between gene and expression
          cor.p <- vector(mode="numeric", length=nrow(methy_matrix))
          DM.p <- vector(mode="numeric", length=nrow(methy_matrix))
          for(p in 1:nrow(methy_matrix))
          {
            DM.p[p] <- signif(wilcox.test(x = methy_matrix[p,phenotype.expr$condition == "Tumor"], y = methy_matrix[p,phenotype.expr$condition == "Normal"])$p.value, digits = 2)
            cor.p[p] = signif(cor.test(x = expr, y = methy_matrix[p,], method = "spearman")$p.value, digits = 2)
          }
          if(sum(cor.p <= 0.05 & DM.p <= 0.05) ==0) V_type <- "TypeII"#if sig DEG, but no correlation between methy and expr, return typeIII
          else V_type <- "Validation"#still ignore the direction of correlation
        }
      }
    }
    V_type
  }
str(validation_distribution)
table(validation_distribution)
#validation_distribution
#TypeI     TypeII Validation
#169         84        92
Promoter_Gene_pairs$LIHC.validation <- validation_distribution

#Validation > TypeI >TypeII
attach(Promoter_Gene_pairs)
X <- Promoter_Gene_pairs[pairType=="HyperDown",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#4         2          11

X <- Promoter_Gene_pairs[pairType=="HypoUp",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#143        53        51

X <- Promoter_Gene_pairs[pairType=="HyperUp",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#3          19        13

X <- Promoter_Gene_pairs[pairType=="HypoDown",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#14         8         14
detach(Promoter_Gene_pairs)

library(foreach)
system.time(Validated_pairs <- foreach(i = 1:nrow(Promoter_Gene_pairs), .combine = 'rbind') %do%
              {
                DMR = Promoter_Gene_pairs$DMR.ID[i]
                Gene = Promoter_Gene_pairs$gene_name[i]
                expr <- unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name), -c(1:2)])
                if(sum(is.na(expr)) <= length(expr)/2 & sum(expr == 0) <= length(expr)/2)
                {
                  expr <- log2(0.001+expr)
                  DMR.chr = strsplit(DMR, split="_", fixed=T)[[1]][1]
                  DMR.start = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][2])
                  DMR.end = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][3])
                  AllCpGs <- paste(DMR.chr, (DMR.start-0):(0+DMR.end), sep="_")
                  methy_matrix <- LIHC.methy[rownames(LIHC.methy) %in% AllCpGs,]
                  # if there are only one CpG, methy_matrix need to be transformed into a matrix
                  if(is.null(nrow(methy_matrix))) 
                  {methy_matrix <- t(as.matrix(methy_matrix))
                  rownames(methy_matrix) = intersect(AllCpGs, rownames(LIHC.methy))}
                  #filtered CpGs with much NAs
                  keep.CpGs <- rownames(methy_matrix)[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2]
                  methy_matrix <- methy_matrix[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2, ]
                  # if there are only one CpG, methy_matrix need to be transformed into a matrix
                  if(is.null(nrow(methy_matrix))) 
                  {methy_matrix <- t(as.matrix(methy_matrix))
                  rownames(methy_matrix) = keep.CpGs}
                  if(nrow(methy_matrix) == 1)
                  {
                    matched_sample <- intersect(names(expr), colnames(methy_matrix))
                    expr <- expr[matched_sample]
                    methy_matrix <- methy_matrix[,matched_sample]#bacame a vector again!
                    phenotypes <- methy.data[match(matched_sample, methy.data$sampleID),]
                    ### wilcox test for expression
                    MeanExpr.T <- mean(expr[phenotypes$condition == "Tumor"], na.rm = T)
                    MeanExpr.N <- mean(expr[phenotypes$condition == "Normal"], na.rm=T)
                    DE.P <- signif(wilcox.test(x = expr[phenotypes$condition=="Tumor"], y = expr[phenotypes$condition=="Normal"])$p.value, digits = 2)
                    if(DE.P < 0.05)
                    {
                      ### wilcox test for methylation
                      MeanMethy.T <- mean(methy_matrix[phenotypes$condition=="Tumor"], na.rm=T)
                      MeanMethy.N <- mean(methy_matrix[phenotypes$condition=="Normal"], na.rm = T)
                      DM.P <- signif(wilcox.test(x = methy_matrix[phenotypes$condition=="Tumor"], 
                                                 y = methy_matrix[phenotypes$condition=="Normal"])$p.value, digits = 2)
                      ### spearman correlation test
                      rho = signif(cor.test(x = expr, y = methy_matrix, method = "spearman")$estimate, digits = 2)
                      pvalue = signif(cor.test(x = expr, y = methy_matrix, method = "spearman")$p.value, digits = 2)
                      if(pvalue < 0.05 & DM.P < 0.05)
                      {
                        signif_cor <- c(DMR, Gene, Promoter_Gene_pairs$Median.delta[i], 
                                        Promoter_Gene_pairs$Pvalue[i], Promoter_Gene_pairs$rho[i],
                                        keep.CpGs, MeanExpr.N, MeanExpr.T, DE.P,
                                        MeanMethy.N, MeanMethy.T, DM.P, rho, pvalue)
                        names(signif_cor) <- c("DMRID", "Gene_name","DM.median", "DM.Pvalue", "MethyandGene.rho", "CpG_pos", "MeanExpr.N", "MeanExpr.T", "DE.wilcoxP", 
                                               "MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP", "rho", "rho.pvalue")
                        signif_cor
                      }
                    }
                  } else if(nrow(methy_matrix)>1)
                  {
                    matched_sample <- intersect(names(expr), colnames(methy_matrix))
                    expr <- expr[matched_sample]
                    methy_matrix <- methy_matrix[,matched_sample]
                    phenotypes <- methy.data[match(matched_sample, methy.data$sampleID),]
                    ### wilcox test for expression
                    MeanExpr.T <- mean(expr[phenotypes$condition == "Tumor"], na.rm = T)
                    MeanExpr.N <- mean(expr[phenotypes$condition == "Normal"], na.rm=T)
                    DE.P <- signif(wilcox.test(x = expr[phenotypes$condition=="Tumor"], y = expr[phenotypes$condition=="Normal"])$p.value, digits = 2)
                    if(DE.P < 0.05)
                    {
                      ### wilcox test for methylation
                      methylation_diff <- foreach(p = 1:nrow(methy_matrix), .combine = 'rbind')%do%
                        {
                          MeanMethy.T <- mean(methy_matrix[p,phenotypes$condition=="Tumor"], na.rm=T)
                          MeanMethy.N <- mean(methy_matrix[p,phenotypes$condition=="Normal"], na.rm = T)
                          DM.P <- signif(wilcox.test(x = methy_matrix[p,phenotypes$condition=="Tumor"], 
                                                     y = methy_matrix[p,phenotypes$condition=="Normal"])$p.value, digits = 2)
                          c(MeanMethy.N, MeanMethy.T, DM.P)
                        }
                      rownames(methylation_diff) <- rownames(methy_matrix)
                      colnames(methylation_diff) <- c("MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP")
                      methylation_diff <- as.data.frame(methylation_diff)
                      ### spearman correlation test
                      spearman_cor <- foreach(q = 1:nrow(methy_matrix), .combine = 'rbind')%do%
                        {
                          rho = signif(cor.test(x = expr, y = unlist(methy_matrix[q,]), method = "spearman")$estimate, digits = 2)
                          pvalue = signif(cor.test(x = expr, y = unlist(methy_matrix[q,]), method = "spearman")$p.value, digits = 2)
                          c(rho, pvalue)
                        }
                      colnames(spearman_cor) <- c("rho", "pvalue")
                      rownames(spearman_cor) <- rownames(methy_matrix)
                      spearman_cor <- as.data.frame(spearman_cor)
                      signif_cor <- foreach(m = 1:nrow(methy_matrix), .combine='rbind')%do%
                        {
                          if(spearman_cor$pvalue[m] < 0.05 & methylation_diff$Methy.wilcoxP[m]<0.05)
                            c(DMR, Gene, Promoter_Gene_pairs$Median.delta[i], 
                              Promoter_Gene_pairs$Pvalue[i], Promoter_Gene_pairs$rho[i],
                              rownames(methy_matrix)[m], MeanExpr.N, MeanExpr.T, DE.P, unlist(methylation_diff[m,]), unlist(spearman_cor[m,]))
                        }
                      if(!is.null(signif_cor))
                      {
                        signif_cor <- as.data.frame(signif_cor)
                        colnames(signif_cor) <- c("DMRID", "Gene_name","DM.median", "DM.Pvalue", "MethyandGene.rho", "CpG_pos", "MeanExpr.N", "MeanExpr.T", "DE.wilcoxP", 
                                                  "MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP", "rho", "rho.pvalue")
                        signif_cor
                      }
                    } 
                  }
                }
              })
Validated_pairs <- as.data.frame(Validated_pairs)
Validated_pairs$pairID <- paste(Validated_pairs$DMRID, Validated_pairs$Gene_name, sep="|")
Validated_pairs$DM.median <- round(as.numeric(Validated_pairs$DM.median), digits=2)
Validated_pairs$DM.Pvalue <- signif(as.numeric(Validated_pairs$DM.Pvalue), digits=2)
Validated_pairs$MethyandGene.rho <- round(as.numeric(Validated_pairs$MethyandGene.rho), digits = 2)
Validated_pairs$MeanExpr.N <- round(as.numeric(Validated_pairs$MeanExpr.N), digits = 2)
Validated_pairs$MeanExpr.T <- round(as.numeric(Validated_pairs$MeanExpr.T), digits = 2)
Validated_pairs$DE.wilcoxP <- signif(as.numeric(Validated_pairs$DE.wilcoxP), digits = 2)
Validated_pairs$MeanMethy.N <- round(as.numeric(Validated_pairs$MeanMethy.N), digits = 2)
Validated_pairs$MeanMethy.T <- round(as.numeric(Validated_pairs$MeanMethy.T), digits = 2)
Validated_pairs$Methy.wilcoxP <- signif(as.numeric(Validated_pairs$Methy.wilcoxP), digits = 2)
Validated_pairs$rho <- round(as.numeric(Validated_pairs$rho), digits = 2)
Validated_pairs$rho.pvalue <- signif(as.numeric(Validated_pairs$rho.pvalue), digits = 2)
Validated_pairs$LFC.LIHC <- Validated_pairs$MeanExpr.T - Validated_pairs$MeanExpr.N
Validated_pairs$MethyDiff.LIHC <- Validated_pairs$MeanMethy.T - Validated_pairs$MeanMethy.N
#rule out inconsistent Methy-Expression correlation
#Validated_pairs <- Validated_pairs[Validated_pairs$MethyandGene.rho*Validated_pairs$rho>0,]
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median > 0 & Validated_pairs$MethyandGene.rho < 0]))#
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median < 0 & Validated_pairs$MethyandGene.rho < 0]))#
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median > 0 & Validated_pairs$MethyandGene.rho > 0]))#
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median < 0 & Validated_pairs$MethyandGene.rho > 0]))#

#Validated_pairs_promoter <- Validated_pairs
write.csv(Validated_pairs, file="/data/huangp/Methy/Mixed/New/LIHC_Validated_Promoter_Like_DMR_Gene_pairs_0.8.csv", row.names=F, quote=T)
#write.csv(Validated_pairs, file="/data/huangp/Methy/Mixed/New/LIHC_Validated_Promoter_Like_DMR_Gene_pairs_1.0.csv", row.names=F, quote=T)

library(fastSave)
load.lbzip2("/data/huangp/Methy/Mixed/New/Delta.gene.logTPM.RDataFS", n.cores = 32)
library(foreach)
library(genomation)
library(GenomicFeatures)
Granges.gtf <- gffToGRanges("/data/huangp/Methy/Mixed/gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
#Granges.gene <- Granges.gene[match(rownames(gene.TPM), Granges.gene$gene_id),]
range_gene_DF <- as.data.frame(Granges.gene)

LFC <- foreach(i = 1:nrow(Promoter_Gene_pairs), .combine = 'c')%do%{
  gene.id <- range_gene_DF$gene_id[range_gene_DF$gene_name==Promoter_Gene_pairs$gene_name[i]]
  LFC.median <- median(delta.expr[gene.id,])
  LFC.median
}
Promoter_Gene_pairs$LFC.median <- round(LFC,2)
write.csv(Promoter_Gene_pairs, file="/data/huangp/Methy/Mixed/New/HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", row.names=F)

library(ggplot2)
library(RColorBrewer)
data = data.frame(mRNA = expr, CpG = unlist(methy_matrix[1,]), Condition = phenotypes$condition)
ggplot(data=data, aes(x = CpG, y = mRNA, fill = Condition))+
  geom_point(shape=21,size=3, color = "white")+
  geom_smooth(method="lm", color = "black", alpha=0.8, se=F, inherit.aes = F, aes(x=CpG, y = mRNA), data=data[,-3], show.legend = F)+
  scale_fill_manual(values = brewer.pal(n = 8, name="Set1")[c(2,5)])+
  scale_color_manual(values = brewer.pal(n = 8, name="Set1")[c(2,5)])+
  labs(x="methylation level", y = "mRNA level (logTPM)")+
  guides(fill = guide_legend(title = ""))+
  theme_classic()
dev.off()

############################ LIHC.TCGA Validation of Enhancer-like DMR-Genes pairs ###########
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/LIHC")
load("LIHC_TCGA.phenotypes.RData")#load LIHC.data
methy.data <- LIHC.data
methy.data$sampleID <- as.character(methy.data$sampleID)
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
LIHC.expr <- z #20531Genes Gene_name+Gene_ID+423Sample
### for each Enhancer-like DMR-Gene pair, calculate their rho 
### load Enhancer-like DMR-Gene pair 
#Enhancer_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/GenicDMR.enhancerLike_Gene_pairs.csv", header=T)#14865
Enhancer_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/FilteredGenicDMR.enhancerLike_Gene_pairs.csv", header=T)#3414
if(1<0){
#filtered enhancer_gene pairs
table(genic_enhancer_genes$DMR.ID %in% promoter_genes$DMR.ID)#11898 2967
table(genic_enhancer_genes$gene_name %in% promoter_genes$gene_name)#9360  5505
genic_enhancer_genes_1 <- genic_enhancer_genes[!(genic_enhancer_genes$DMR.ID %in% promoter_genes$DMR.ID),]
table(genic_enhancer_genes_1$gene_name %in% promoter_genes$gene_name)#8102 3796
genic_enhancer_genes_2 <- genic_enhancer_genes_1[!(genic_enhancer_genes_1$gene_name %in% promoter_genes$gene_name),]
table(Hyper = genic_enhancer_genes_2$Median.delta >0, Classic = genic_enhancer_genes_2$rho < 0)
}
#cREDMR.summary <- read.csv("/data/huangp/Methy/Mixed/New/cREDMR.summary_with_overlapRatio_1.csv", header=T)
#cREDMR.summary <- read.csv("/data/huangp/Methy/Mixed/New/cREDMR.summary_with_overlapRatio_0.5.csv", header=T)
cREDMR.summary <- read.csv("/data/huangp/Methy/Mixed/New/cREDMR.summary_with_overlapRatio_0.8.csv", header=T)
#cREDMR.summary <- read.csv("/data/huangp/Methy/Mixed/New/cREDMR.summary_with_overlapRatio_1.csv", header=T)
cREDMR.summary <- cREDMR.summary[cREDMR.summary$EnhancerLike.score > 0,]
rownames(cREDMR.summary) <- paste(cREDMR.summary$chr, cREDMR.summary$str, cREDMR.summary$end, sep="_")
Enhancer_Gene_pairs <- Enhancer_Gene_pairs[Enhancer_Gene_pairs$DMR.ID %in% rownames(cREDMR.summary),]#
Enhancer_Gene_pairs$pairID <- paste(Enhancer_Gene_pairs$DMR.ID, Enhancer_Gene_pairs$gene_name, sep="|")
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$Median.delta > 0 & Enhancer_Gene_pairs$rho < 0]))#22
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$Median.delta < 0 & Enhancer_Gene_pairs$rho < 0]))#591
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$Median.delta > 0 & Enhancer_Gene_pairs$rho > 0]))#150
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$Median.delta < 0 & Enhancer_Gene_pairs$rho > 0]))#83
Enhancer_Gene_pairs$pairType <- rep("NA", nrow(Enhancer_Gene_pairs))
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$Median.delta > 0 & Enhancer_Gene_pairs$rho < 0] <- "HyperDown"
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$Median.delta < 0 & Enhancer_Gene_pairs$rho < 0] <- "HypoUp"
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$Median.delta > 0 & Enhancer_Gene_pairs$rho > 0] <- "HyperUp"
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$Median.delta < 0 & Enhancer_Gene_pairs$rho > 0] <- "HypoDown"

###validation classification
#Type I failure: No CpG
#Type II failure: There was at least one CpG and significant DEG, but no significant correlation between them
#Type III: Successful Validation 
library(foreach)
validation_distribution <- foreach(i = 1:nrow(Enhancer_Gene_pairs), .combine = 'c') %do%
  {
    DMR = Enhancer_Gene_pairs$DMR.ID[i]
    Gene = Enhancer_Gene_pairs$gene_name[i]
    expr <- unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name), -c(1:2)])
    DMR.chr = strsplit(DMR, split="_", fixed=T)[[1]][1]
    DMR.start = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][2])
    DMR.end = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][3])
    AllCpGs <- paste(DMR.chr, (DMR.start-0):(0+DMR.end), sep="_")
    methy_matrix <- LIHC.methy[rownames(LIHC.methy) %in% AllCpGs,]
    #if only one CpG
    if(is.null(nrow(methy_matrix))){
      methy_matrix <- t(as.matrix(methy_matrix))
      rownames(methy_matrix) <- intersect(rownames(LIHC.methy), "chr3_158672843")
    }
    #filtered CpG with too much NAs
    keep.CpGs <- rownames(methy_matrix)[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2]
    methy_matrix <- methy_matrix[keep.CpGs,]
    if(length(keep.CpGs)==0) V_type <- "TypeI" #If there was no CpG, return typeI validation failure
    else
    {
      if(length(keep.CpGs)==1) 
      {methy_matrix <- t(as.matrix(methy_matrix))
      rownames(methy_matrix) <- keep.CpGs}
      #for each CpG, exam their correlation with expr
      matchedSample <- intersect(colnames(methy_matrix), names(expr))
      expr <- expr[matchedSample]
      #filtered out gene expresssion with too much 0 or NA
      if(sum(is.na(expr)) > length(expr)/2 | sum(expr == 0) > length(expr)/2) V_type <- "TypeII"
      else
      {
        #filtered out gene without significant DEG
        expr <- log2(0.001+expr)
        phenotype.expr <- methy.data[match(names(expr), methy.data$sampleID),]
        DE.P <- signif(wilcox.test(x = expr[phenotype.expr$condition=="Tumor"], y = expr[phenotype.expr$condition=="Normal"])$p.value, digits = 2)
        if(DE.P > 0.05) V_type <- "TypeII"
        else
        {
          #filtered methy_matrix for matched samples with expr data
          methy_matrix <- methy_matrix[,matchedSample]
          if(is.null(nrow(methy_matrix))){
            methy_matrix <- t(as.matrix(methy_matrix))
            rownames(methy_matrix) <- keep.CpGs
          }
          #filtered out gene without significant DM CpG
          #filtered out insignificant correlation between gene and expression
          cor.p <- vector(mode="numeric", length=nrow(methy_matrix))
          DM.p <- vector(mode="numeric", length=nrow(methy_matrix))
          for(p in 1:nrow(methy_matrix))
          {
            DM.p[p] <- signif(wilcox.test(x = methy_matrix[p,phenotype.expr$condition == "Tumor"], y = methy_matrix[p,phenotype.expr$condition == "Normal"])$p.value, digits = 2)
            cor.p[p] = signif(cor.test(x = expr, y = methy_matrix[p,], method = "spearman")$p.value, digits = 2)
          }
          if(sum(cor.p <= 0.05 & DM.p <= 0.05) ==0) V_type <- "TypeII"#if sig DEG, but no correlation between methy and expr, return typeIII
          else V_type <- "Validation"
        }
      }
    }
    V_type
  }
str(validation_distribution)
table(validation_distribution)
#validation_distribution
#TypeI     TypeII Validation
#700        139        132
Enhancer_Gene_pairs$LIHC.validation <- validation_distribution
#Validation > TypeI >TypeII
attach(Enhancer_Gene_pairs)
X <- Enhancer_Gene_pairs[pairType=="HyperDown",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#18          2          2

X <- Enhancer_Gene_pairs[pairType=="HypoUp",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#440         73         78

X <- Enhancer_Gene_pairs[pairType=="HyperUp",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#93        35         22

X <- Enhancer_Gene_pairs[pairType=="HypoDown",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#62         6         15
detach(Enhancer_Gene_pairs)

library(fastSave)
load.lbzip2("/data/huangp/Methy/Mixed/New/Delta.gene.logTPM.RDataFS", n.cores = 32)

library(foreach)
library(genomation)
library(GenomicFeatures)
Granges.gtf <- gffToGRanges("/data/huangp/Methy/Mixed/gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
#Granges.gene <- Granges.gene[match(rownames(gene.TPM), Granges.gene$gene_id),]
range_gene_DF <- as.data.frame(Granges.gene)

LFC <- foreach(i = 1:nrow(Enhancer_Gene_pairs), .combine = 'c')%do%{
  gene.id <- range_gene_DF$gene_id[match(Enhancer_Gene_pairs$gene_name[i], range_gene_DF$gene_name)]
  LFC.median <- median(delta.expr[gene.id,])
  LFC.median
}
Enhancer_Gene_pairs$LFC.median <- round(LFC,2)
write.csv(Enhancer_Gene_pairs, file="/data/huangp/Methy/Mixed/New/HighConfidence_Genic_Enhancer_Like_DMR_Gene_pairs_0.8.csv", row.names=F)


library(foreach)
system.time(Validated_pairs <- foreach(i = 1:nrow(Enhancer_Gene_pairs), .combine = 'rbind') %do%
              {
                DMR = Enhancer_Gene_pairs$DMR.ID[i]
                Gene = Enhancer_Gene_pairs$gene_name[i]
                expr <- unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name), -c(1:2)])
                if(sum(is.na(expr)) <= length(expr)/2 & sum(expr==0) <= length(expr)/2)
                {
                  expr <- log2(0.001+expr)
                  DMR.chr = strsplit(DMR, split="_", fixed=T)[[1]][1]
                  DMR.start = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][2])
                  DMR.end = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][3])
                  AllCpGs <- paste(DMR.chr, (DMR.start-0):(0+DMR.end), sep="_")
                  methy_matrix <- LIHC.methy[rownames(LIHC.methy) %in% AllCpGs,]
                  # if there are only one CpG, methy_matrix need to be transformed into a matrix
                  if(is.null(nrow(methy_matrix))) 
                  {methy_matrix <- t(as.matrix(methy_matrix))
                  rownames(methy_matrix) = intersect(AllCpGs, rownames(LIHC.methy))}
                  #filtered CpGs with much NAs
                  keep.CpGs <- rownames(methy_matrix)[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2]
                  methy_matrix <- methy_matrix[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2, ]
                  # if there are only one CpG, methy_matrix need to be transformed into a matrix
                  if(is.null(nrow(methy_matrix))) 
                  {methy_matrix <- t(as.matrix(methy_matrix))
                  rownames(methy_matrix) = keep.CpGs}
                  if(nrow(methy_matrix) == 1)
                  {
                    matched_sample <- intersect(names(expr), colnames(methy_matrix))
                    expr <- expr[matched_sample]
                    methy_matrix <- methy_matrix[,matched_sample]#bacame a vector again!
                    phenotypes <- methy.data[match(matched_sample, methy.data$sampleID),]
                    ### wilcox test for expression
                    MeanExpr.T <- mean(expr[phenotypes$condition == "Tumor"], na.rm = T)
                    MeanExpr.N <- mean(expr[phenotypes$condition == "Normal"], na.rm=T)
                    DE.P <- signif(wilcox.test(x = expr[phenotypes$condition=="Tumor"], y = expr[phenotypes$condition=="Normal"])$p.value, digits = 2)
                    if(DE.P < 0.05)
                    {
                      ### wilcox test for methylation
                      MeanMethy.T <- mean(methy_matrix[phenotypes$condition=="Tumor"], na.rm=T)
                      MeanMethy.N <- mean(methy_matrix[phenotypes$condition=="Normal"], na.rm = T)
                      DM.P <- signif(wilcox.test(x = methy_matrix[phenotypes$condition=="Tumor"], 
                                                 y = methy_matrix[phenotypes$condition=="Normal"])$p.value, digits = 2)
                      ### spearman correlation test
                      rho = signif(cor.test(x = expr, y = methy_matrix, method = "spearman")$estimate, digits = 2)
                      pvalue = signif(cor.test(x = expr, y = methy_matrix, method = "spearman")$p.value, digits = 2)
                      if(pvalue < 0.05 & DM.P < 0.05)
                      {
                        signif_cor <- c(DMR, Gene, Enhancer_Gene_pairs$Median.delta[i], 
                                        Enhancer_Gene_pairs$Pvalue[i], Enhancer_Gene_pairs$rho[i],
                                        keep.CpGs, MeanExpr.N, MeanExpr.T, DE.P,
                                        MeanMethy.N, MeanMethy.T, DM.P, rho, pvalue)
                        names(signif_cor) <- c("DMRID", "Gene_name","DM.median", "DM.Pvalue", "MethyandGene.rho", "CpG_pos", "MeanExpr.N", "MeanExpr.T", "DE.wilcoxP", 
                                               "MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP", "rho", "rho.pvalue")
                        signif_cor
                      }
                    }
                  } else if(nrow(methy_matrix)>1)
                  {
                    matched_sample <- intersect(names(expr), colnames(methy_matrix))
                    expr <- expr[matched_sample]
                    methy_matrix <- methy_matrix[,matched_sample]
                    phenotypes <- methy.data[match(matched_sample, methy.data$sampleID),]
                    ### wilcox test for expression
                    MeanExpr.T <- mean(expr[phenotypes$condition == "Tumor"], na.rm = T)
                    MeanExpr.N <- mean(expr[phenotypes$condition == "Normal"], na.rm=T)
                    DE.P <- signif(wilcox.test(x = expr[phenotypes$condition=="Tumor"], y = expr[phenotypes$condition=="Normal"])$p.value, digits = 2)
                    if(DE.P < 0.05)
                    {
                      ### wilcox test for methylation
                      methylation_diff <- foreach(p = 1:nrow(methy_matrix), .combine = 'rbind')%do%
                        {
                          MeanMethy.T <- mean(methy_matrix[p,phenotypes$condition=="Tumor"], na.rm=T)
                          MeanMethy.N <- mean(methy_matrix[p,phenotypes$condition=="Normal"], na.rm = T)
                          DM.P <- signif(wilcox.test(x = methy_matrix[p,phenotypes$condition=="Tumor"], 
                                                     y = methy_matrix[p,phenotypes$condition=="Normal"])$p.value, digits = 2)
                          c(MeanMethy.N, MeanMethy.T, DM.P)
                        }
                      rownames(methylation_diff) <- rownames(methy_matrix)
                      colnames(methylation_diff) <- c("MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP")
                      methylation_diff <- as.data.frame(methylation_diff)
                      ### spearman correlation test
                      spearman_cor <- foreach(q = 1:nrow(methy_matrix), .combine = 'rbind')%do%
                        {
                          rho = signif(cor.test(x = expr, y = unlist(methy_matrix[q,]), method = "spearman")$estimate, digits = 2)
                          pvalue = signif(cor.test(x = expr, y = unlist(methy_matrix[q,]), method = "spearman")$p.value, digits = 2)
                          c(rho, pvalue)
                        }
                      colnames(spearman_cor) <- c("rho", "pvalue")
                      rownames(spearman_cor) <- rownames(methy_matrix)
                      spearman_cor <- as.data.frame(spearman_cor)
                      signif_cor <- foreach(m = 1:nrow(methy_matrix), .combine='rbind')%do%
                        {
                          if(spearman_cor$pvalue[m] < 0.05 & methylation_diff$Methy.wilcoxP[m]<0.05)
                            c(DMR, Gene, Enhancer_Gene_pairs$Median.delta[i], 
                              Enhancer_Gene_pairs$Pvalue[i], Enhancer_Gene_pairs$rho[i],
                              rownames(methy_matrix)[m], MeanExpr.N, MeanExpr.T, DE.P, unlist(methylation_diff[m,]), unlist(spearman_cor[m,]))
                        }
                      if(!is.null(signif_cor))
                      {
                        signif_cor <- as.data.frame(signif_cor)
                        colnames(signif_cor) <- c("DMRID", "Gene_name","DM.median", "DM.Pvalue", "MethyandGene.rho", "CpG_pos", "MeanExpr.N", "MeanExpr.T", "DE.wilcoxP", 
                                                  "MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP", "rho", "rho.pvalue")
                        signif_cor
                      }
                    } 
                  } 
                }
              })
Validated_pairs <- as.data.frame(Validated_pairs)
Validated_pairs$pairID <- paste(Validated_pairs$DMRID, Validated_pairs$Gene_name, sep="|")
Validated_pairs$DM.median <- round(as.numeric(Validated_pairs$DM.median), digits=2)
Validated_pairs$DM.Pvalue <- signif(as.numeric(Validated_pairs$DM.Pvalue), digits=2)
Validated_pairs$MethyandGene.rho <- round(as.numeric(Validated_pairs$MethyandGene.rho), digits = 2)
Validated_pairs$MeanExpr.N <- round(as.numeric(Validated_pairs$MeanExpr.N), digits = 2)
Validated_pairs$MeanExpr.T <- round(as.numeric(Validated_pairs$MeanExpr.T), digits = 2)
Validated_pairs$DE.wilcoxP <- signif(as.numeric(Validated_pairs$DE.wilcoxP), digits = 2)
Validated_pairs$MeanMethy.N <- round(as.numeric(Validated_pairs$MeanMethy.N), digits = 2)
Validated_pairs$MeanMethy.T <- round(as.numeric(Validated_pairs$MeanMethy.T), digits = 2)
Validated_pairs$Methy.wilcoxP <- signif(as.numeric(Validated_pairs$Methy.wilcoxP), digits = 2)
Validated_pairs$rho <- round(as.numeric(Validated_pairs$rho), digits = 2)
Validated_pairs$rho.pvalue <- signif(as.numeric(Validated_pairs$rho.pvalue), digits = 2)
Validated_pairs$LFC.LIHC <- Validated_pairs$MeanExpr.T - Validated_pairs$MeanExpr.N
Validated_pairs$MethyDiff.LIHC <- Validated_pairs$MeanMethy.T - Validated_pairs$MeanMethy.N
#rule out inconsistent Methy-Expression correlation
#Validated_pairs <- Validated_pairs[Validated_pairs$MethyandGene.rho*Validated_pairs$rho>0,]
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median > 0 & Validated_pairs$MethyandGene.rho < 0]))#52, 43, 25, 12, 8
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median < 0 & Validated_pairs$MethyandGene.rho < 0]))#350, 300, 162, 87, 56
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median > 0 & Validated_pairs$MethyandGene.rho > 0]))#44, 32, 18, 11, 9
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median < 0 & Validated_pairs$MethyandGene.rho > 0]))

#Validated_pairs <- cbind(Validated_pairs, Enhancer_Gene_pairs[match(Validated_pairs$pairID, paste(Enhancer_Gene_pairs$DMR.ID, Enhancer_Gene_pairs$gene_name, sep="|")), 
 #                                                             c("distance", "lfc.GSE5230", "lfc.GSE35311", "lfc.GSE112788", "lfc.GSE67318", "lfc.up", "lfc.down")],
#                         cREDMR.summary[match(Validated_pairs$DMRID, rownames(cREDMR.summary)), c("EnhancerLike.T", "EnhancerLike.N", "EnhancerLike.score")])
write.csv(Validated_pairs, file="/data/huangp/Methy/Mixed/New/LIHC_Validated_Enhancer_Like_GenicDMR_Gene_pairs_0.8.csv",  row.names=F, quote=T)

Enhancer_Gene_pairs$LIHCValidation <- Enhancer_Gene_pairs$pairID %in% Validated_pairs$pairID
Prognostic <- read.csv("/data/huangp/Methy/Mixed/New/The185PrognosticGenes.csv", header=T)[[2]]
Enhancer_Gene_pairs$Prognositc <- Enhancer_Gene_pairs$gene_name %in% Prognostic
Enhancer_Gene_pairs[,c("EnhancerLike.score", "EnhancerLike.T", "EnhancerLike.N")] <- cREDMR.summary[match(Enhancer_Gene_pairs$DMR.ID, rownames(cREDMR.summary)),c("EnhancerLike.score", "EnhancerLike.T", "EnhancerLike.N")]
write.csv(Enhancer_Gene_pairs, file="/data/huangp/Methy/Mixed/New/HighConfidence_GenicEnhancer_Like_DMR_Gene_pairs_0.8.csv", row.names = F)

write.csv(Validated_pairs, file="/data/huangp/Methy/Mixed/New/LIHC_Validated_Enhancer_Like_GenicDMR_Gene_pairs_1.0.csv",  row.names=F, quote=T)
Enhancer_Gene_pairs$LIHCValidation <- Enhancer_Gene_pairs$pairID %in% Validated_pairs$pairID
Prognostic <- read.csv("/data/huangp/Methy/Mixed/New/The185PrognosticGenes.csv", header=T)[[2]]
Enhancer_Gene_pairs$Prognositc <- Enhancer_Gene_pairs$gene_name %in% Prognostic
Enhancer_Gene_pairs[,c("EnhancerLike.score", "EnhancerLike.T", "EnhancerLike.N")] <- cREDMR.summary[match(Enhancer_Gene_pairs$DMR.ID, rownames(cREDMR.summary)),c("EnhancerLike.score", "EnhancerLike.T", "EnhancerLike.N")]
write.csv(Enhancer_Gene_pairs, file="/data/huangp/Methy/Mixed/New/HighConfidence_GenicEnhancer_Like_DMR_Gene_pairs_1.0.csv", row.names = F)
############################ LIHC Validation of Intergenic enhancer-like DMR-Gene pairs #####################
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/LIHC")
load("LIHC_TCGA.phenotypes.RData")#load LIHC.data
methy.data <- LIHC.data
methy.data$sampleID <- as.character(methy.data$sampleID)
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
LIHC.expr <- z  #20531Genes Gene_name+Gene_ID+423Sample
### for each Enhancer-like DMR-Gene pair, calculate their rho 
### load Enhancer-like DMR-Gene pair 
#Enhancer_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/GenicDMR.enhancerLike_Gene_pairs.csv", header=T)#14865
Enhancer_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/FilteredConsistent_eRNARegulatedDMRs_Gene_pairs.csv", header=T)
Enhancer_Gene_pairs <- Enhancer_Gene_pairs[Enhancer_Gene_pairs$eRNAandGene.rho >= 0.7,]
#Enhancer_Gene_pairs <- Enhancer_Gene_pairs[Enhancer_Gene_pairs$eRNAandGene.rho >= 0.8,]
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$DM.Median > 0 & Enhancer_Gene_pairs$MethyandGene.rho < 0]))#7
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$DM.Median < 0 & Enhancer_Gene_pairs$MethyandGene.rho < 0]))#80
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$DM.Median > 0 & Enhancer_Gene_pairs$MethyandGene.rho > 0]))#0
length(unique(Enhancer_Gene_pairs$gene_name[Enhancer_Gene_pairs$DM.Median < 0 & Enhancer_Gene_pairs$MethyandGene.rho > 0]))#103

Enhancer_Gene_pairs$pairType <- rep("NA", nrow(Enhancer_Gene_pairs))
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$DM.Median > 0 & Enhancer_Gene_pairs$MethyandGene.rho < 0] <- "HyperDown"
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$DM.Median < 0 & Enhancer_Gene_pairs$MethyandGene.rho < 0] <- "HypoUp"
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$DM.Median > 0 & Enhancer_Gene_pairs$MethyandGene.rho > 0] <- "HyperUp"
Enhancer_Gene_pairs$pairType[Enhancer_Gene_pairs$DM.Median < 0 & Enhancer_Gene_pairs$MethyandGene.rho > 0] <- "HypoDown"

###validation classification
#Type I failure: No CpG
#Type II failure: There was at least one CpG and significant DEG, but no significant correlation between them
#Type III: Successful Validation 
library(foreach)
validation_distribution <- foreach(i = 1:nrow(Enhancer_Gene_pairs), .combine = 'c') %do%
  {
    DMR = Enhancer_Gene_pairs$DMR.ID[i]
    Gene = Enhancer_Gene_pairs$gene_name[i]
    expr <- unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name), -c(1:2)])
    DMR.chr = strsplit(DMR, split="_", fixed=T)[[1]][1]
    DMR.start = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][2])
    DMR.end = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][3])
    AllCpGs <- paste(DMR.chr, (DMR.start-0):(0+DMR.end), sep="_")
    methy_matrix <- LIHC.methy[rownames(LIHC.methy) %in% AllCpGs,]
    #if only one CpG
    if(is.null(nrow(methy_matrix))){
      methy_matrix <- t(as.matrix(methy_matrix))
      rownames(methy_matrix) <- intersect(rownames(LIHC.methy), "chr3_158672843")
    }
    #filtered CpG with too much NAs
    keep.CpGs <- rownames(methy_matrix)[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2]
    methy_matrix <- methy_matrix[keep.CpGs,]
    if(length(keep.CpGs)==0) V_type <- "TypeI" #If there was no CpG, return typeI validation failure
    else
    {
      if(length(keep.CpGs)==1) 
      {methy_matrix <- t(as.matrix(methy_matrix))
      rownames(methy_matrix) <- keep.CpGs}
      #for each CpG, exam their correlation with expr
      matchedSample <- intersect(colnames(methy_matrix), names(expr))
      expr <- expr[matchedSample]
      #filtered out gene expresssion with too much 0 or NA
      if(sum(is.na(expr)) > length(expr)/2 | sum(expr == 0) > length(expr)/2) V_type <- "TypeII"
      else
      {
        #filtered out gene without significant DEG
        expr <- log2(0.001+expr)
        phenotype.expr <- methy.data[match(names(expr), methy.data$sampleID),]
        DE.P <- signif(wilcox.test(x = expr[phenotype.expr$condition=="Tumor"], y = expr[phenotype.expr$condition=="Normal"])$p.value, digits = 2)
        if(DE.P > 0.05) V_type <- "TypeII"
        else
        {
          #filtered methy_matrix for matched samples with expr data
          methy_matrix <- methy_matrix[,matchedSample]
          if(is.null(nrow(methy_matrix))){
            methy_matrix <- t(as.matrix(methy_matrix))
            rownames(methy_matrix) <- keep.CpGs
          }
          #filtered out gene without significant DM CpG
          #filtered out insignificant correlation between gene and expression
          cor.p <- vector(mode="numeric", length=nrow(methy_matrix))
          DM.p <- vector(mode="numeric", length=nrow(methy_matrix))
          for(p in 1:nrow(methy_matrix))
          {
            DM.p[p] <- signif(wilcox.test(x = methy_matrix[p,phenotype.expr$condition == "Tumor"], y = methy_matrix[p,phenotype.expr$condition == "Normal"])$p.value, digits = 2)
            cor.p[p] = signif(cor.test(x = expr, y = methy_matrix[p,], method = "spearman")$p.value, digits = 2)
          }
          if(sum(cor.p <= 0.05 & DM.p <= 0.05) ==0) V_type <- "TypeII"#if sig DEG, but no correlation between methy and expr, return typeIII
          else V_type <- "Validation"
        }
      }
    }
    V_type
  }
str(validation_distribution)
table(validation_distribution)
#validation_distribution
#TypeI     TypeII Validation
#700        139        132
Enhancer_Gene_pairs$LIHC.validation <- validation_distribution
#Validation > TypeI >TypeII
attach(Enhancer_Gene_pairs)
X <- Enhancer_Gene_pairs[pairType=="HyperDown",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#7         0          0

X <- Enhancer_Gene_pairs[pairType=="HypoUp",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#70         3         7

X <- Enhancer_Gene_pairs[pairType=="HyperUp",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#1        0         0

X <- Enhancer_Gene_pairs[pairType=="HypoDown",]
Gene <- unique(X$gene_name)
Gene.validationType <- vector(mode="character", length=length(Gene))
for(i in 1:length(Gene)){
  validation.type <- X$LIHC.validation[X$gene_name==Gene[i]]
  if("Validation" %in% validation.type) Gene.validationType[i] <- "Validation"
  else if("TypeI" %in% validation.type) Gene.validationType[i] <- "TypeI"
  else Gene.validationType[i] <- "TypeII"
}
table(Gene.validationType)
#Gene.validationType
#TypeI     TypeII Validation
#96        0        7
detach(Enhancer_Gene_pairs)

library(fastSave)
load.lbzip2("/data/huangp/Methy/Mixed/New/Delta.gene.logTPM.RDataFS", n.cores = 32)

library(foreach)
library(genomation)
library(GenomicFeatures)
Granges.gtf <- gffToGRanges("/data/huangp/Methy/Mixed/gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
#Granges.gene <- Granges.gene[match(rownames(gene.TPM), Granges.gene$gene_id),]
range_gene_DF <- as.data.frame(Granges.gene)

LFC <- foreach(i = 1:nrow(Enhancer_Gene_pairs), .combine = 'c')%do%{
  gene.id <- range_gene_DF$gene_id[match(Enhancer_Gene_pairs$gene_name[i], range_gene_DF$gene_name)]
  LFC.median <- median(delta.expr[gene.id,])
  LFC.median
}
Enhancer_Gene_pairs$LFC.median <- round(LFC,2)
write.csv(Enhancer_Gene_pairs, file="/data/huangp/Methy/Mixed/New/HighConfidence_Intergenic_Enhancer_Like_DMR_Gene_pairs_0.7.csv", row.names=F)


library(foreach)
system.time(Validated_pairs <- foreach(i = 1:nrow(Enhancer_Gene_pairs), .combine = 'rbind') %do%
              {
                DMR = Enhancer_Gene_pairs$DMR.ID[i]
                Gene = Enhancer_Gene_pairs$gene_name[i]
                expr <- unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name), -c(1:2)])
                if(sum(is.na(expr)) <= length(expr)/2 & sum(expr ==0) <= length(expr)/2)
                {
                  expr <- log2(0.001+expr)
                  DMR.chr = strsplit(DMR, split="_", fixed=T)[[1]][1]
                  DMR.start = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][2])
                  DMR.end = as.integer(strsplit(DMR, split="_", fixed=T)[[1]][3])
                  AllCpGs <- paste(DMR.chr, (DMR.start-0):(0+DMR.end), sep="_")
                  methy_matrix <- LIHC.methy[rownames(LIHC.methy) %in% AllCpGs,]
                  # if there are only one CpG, methy_matrix need to be transformed into a matrix
                  if(is.null(nrow(methy_matrix))) 
                  {methy_matrix <- t(as.matrix(methy_matrix))
                  rownames(methy_matrix) = intersect(AllCpGs, rownames(LIHC.methy))}
                  #filtered CpGs with much NAs
                  keep.CpGs <- rownames(methy_matrix)[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2]
                  methy_matrix <- methy_matrix[rowSums(is.na(methy_matrix)) <= ncol(methy_matrix)/2, ]
                  # if there are only one CpG, methy_matrix need to be transformed into a matrix
                  if(is.null(nrow(methy_matrix))) 
                  {methy_matrix <- t(as.matrix(methy_matrix))
                  rownames(methy_matrix) = keep.CpGs}
                  if(nrow(methy_matrix) == 1)
                  {
                    matched_sample <- intersect(names(expr), colnames(methy_matrix))
                    expr <- expr[matched_sample]
                    methy_matrix <- methy_matrix[,matched_sample]#bacame a vector again!
                    phenotypes <- methy.data[match(matched_sample, methy.data$sampleID),]
                    ### wilcox test for expression
                    MeanExpr.T <- mean(expr[phenotypes$condition == "Tumor"], na.rm = T)
                    MeanExpr.N <- mean(expr[phenotypes$condition == "Normal"], na.rm=T)
                    DE.P <- signif(wilcox.test(x = expr[phenotypes$condition=="Tumor"], y = expr[phenotypes$condition=="Normal"])$p.value, digits = 2)
                    if(DE.P < 0.05)
                    {
                      ### wilcox test for methylation
                      MeanMethy.T <- mean(methy_matrix[phenotypes$condition=="Tumor"], na.rm=T)
                      MeanMethy.N <- mean(methy_matrix[phenotypes$condition=="Normal"], na.rm = T)
                      DM.P <- signif(wilcox.test(x = methy_matrix[phenotypes$condition=="Tumor"], 
                                                 y = methy_matrix[phenotypes$condition=="Normal"])$p.value, digits = 2)
                      ### spearman correlation test
                      rho = signif(cor.test(x = expr, y = methy_matrix, method = "spearman")$estimate, digits = 2)
                      pvalue = signif(cor.test(x = expr, y = methy_matrix, method = "spearman")$p.value, digits = 2)
                      if(pvalue < 0.05 & DM.P < 0.05)
                      {
                        signif_cor <- c(DMR, Gene, Enhancer_Gene_pairs$DM.Median[i], 
                                        Enhancer_Gene_pairs$DM.Pvalue[i], Enhancer_Gene_pairs$MethyandGene.rho[i],
                                        keep.CpGs, MeanExpr.N, MeanExpr.T, DE.P,
                                        MeanMethy.N, MeanMethy.T, DM.P, rho, pvalue)
                        names(signif_cor) <- c("DMRID", "Gene_name","DM.median", "DM.Pvalue", "MethyandGene.rho", "CpG_pos", "MeanExpr.N", "MeanExpr.T", "DE.wilcoxP", 
                                               "MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP", "rho", "rho.pvalue")
                        signif_cor
                      }
                    }
                  } else if(nrow(methy_matrix)>1)
                  {
                    matched_sample <- intersect(names(expr), colnames(methy_matrix))
                    expr <- expr[matched_sample]
                    methy_matrix <- methy_matrix[,matched_sample]
                    phenotypes <- methy.data[match(matched_sample, methy.data$sampleID),]
                    ### wilcox test for expression
                    MeanExpr.T <- mean(expr[phenotypes$condition == "Tumor"], na.rm = T)
                    MeanExpr.N <- mean(expr[phenotypes$condition == "Normal"], na.rm=T)
                    DE.P <- signif(wilcox.test(x = expr[phenotypes$condition=="Tumor"], y = expr[phenotypes$condition=="Normal"])$p.value, digits = 2)
                    if(DE.P < 0.05)
                    {
                      ### wilcox test for methylation
                      methylation_diff <- foreach(p = 1:nrow(methy_matrix), .combine = 'rbind')%do%
                        {
                          MeanMethy.T <- mean(methy_matrix[p,phenotypes$condition=="Tumor"], na.rm=T)
                          MeanMethy.N <- mean(methy_matrix[p,phenotypes$condition=="Normal"], na.rm = T)
                          DM.P <- signif(wilcox.test(x = methy_matrix[p,phenotypes$condition=="Tumor"], 
                                                     y = methy_matrix[p,phenotypes$condition=="Normal"])$p.value, digits = 2)
                          c(MeanMethy.N, MeanMethy.T, DM.P)
                        }
                      rownames(methylation_diff) <- rownames(methy_matrix)
                      colnames(methylation_diff) <- c("MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP")
                      methylation_diff <- as.data.frame(methylation_diff)
                      ### spearman correlation test
                      spearman_cor <- foreach(q = 1:nrow(methy_matrix), .combine = 'rbind')%do%
                        {
                          rho = signif(cor.test(x = expr, y = unlist(methy_matrix[q,]), method = "spearman")$estimate, digits = 2)
                          pvalue = signif(cor.test(x = expr, y = unlist(methy_matrix[q,]), method = "spearman")$p.value, digits = 2)
                          c(rho, pvalue)
                        }
                      colnames(spearman_cor) <- c("rho", "pvalue")
                      rownames(spearman_cor) <- rownames(methy_matrix)
                      spearman_cor <- as.data.frame(spearman_cor)
                      signif_cor <- foreach(m = 1:nrow(methy_matrix), .combine='rbind')%do%
                        {
                          if(spearman_cor$pvalue[m] < 0.05 & methylation_diff$Methy.wilcoxP[m]<0.05)
                            c(DMR, Gene, Enhancer_Gene_pairs$DM.Median[i], 
                              Enhancer_Gene_pairs$DM.Pvalue[i], Enhancer_Gene_pairs$MethyandGene.rho[i],
                              rownames(methy_matrix)[m], MeanExpr.N, MeanExpr.T, DE.P, unlist(methylation_diff[m,]), unlist(spearman_cor[m,]))
                        }
                      if(!is.null(signif_cor))
                      {
                        signif_cor <- as.data.frame(signif_cor)
                        colnames(signif_cor) <- c("DMRID", "Gene_name","DM.median", "DM.Pvalue", "MethyandGene.rho", "CpG_pos", "MeanExpr.N", "MeanExpr.T", "DE.wilcoxP", 
                                                  "MeanMethy.N", "MeanMethy.T", "Methy.wilcoxP", "rho", "rho.pvalue")
                        signif_cor
                      }
                    } 
                  } 
                }
              })
Validated_pairs <- as.data.frame(Validated_pairs)
Validated_pairs$pairID <- paste(Validated_pairs$DMRID, Validated_pairs$Gene_name, sep="|")
Validated_pairs$DM.median <- round(as.numeric(Validated_pairs$DM.median), digits=2)
Validated_pairs$DM.Pvalue <- signif(as.numeric(Validated_pairs$DM.Pvalue), digits=2)
Validated_pairs$MethyandGene.rho <- round(as.numeric(Validated_pairs$MethyandGene.rho), digits = 2)
Validated_pairs$MeanExpr.N <- round(as.numeric(Validated_pairs$MeanExpr.N), digits = 2)
Validated_pairs$MeanExpr.T <- round(as.numeric(Validated_pairs$MeanExpr.T), digits = 2)
Validated_pairs$DE.wilcoxP <- signif(as.numeric(Validated_pairs$DE.wilcoxP), digits = 2)
Validated_pairs$MeanMethy.N <- round(as.numeric(Validated_pairs$MeanMethy.N), digits = 2)
Validated_pairs$MeanMethy.T <- round(as.numeric(Validated_pairs$MeanMethy.T), digits = 2)
Validated_pairs$Methy.wilcoxP <- signif(as.numeric(Validated_pairs$Methy.wilcoxP), digits = 2)
Validated_pairs$rho <- round(as.numeric(Validated_pairs$rho), digits = 2)
Validated_pairs$rho.pvalue <- signif(as.numeric(Validated_pairs$rho.pvalue), digits = 2)
Validated_pairs$LFC.LIHC <- Validated_pairs$MeanExpr.T - Validated_pairs$MeanExpr.N
Validated_pairs$MethyDiff.LIHC <- Validated_pairs$MeanMethy.T - Validated_pairs$MeanMethy.N
#rule out inconsistent Methy-Expression correlation
#Validated_pairs <- Validated_pairs[Validated_pairs$MethyandGene.rho*Validated_pairs$rho>0,]

length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median > 0 & Validated_pairs$MethyandGene.rho < 0]))#52, 43, 25, 12, 8
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median < 0 & Validated_pairs$MethyandGene.rho < 0]))#350, 300, 162, 87, 56
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median > 0 & Validated_pairs$MethyandGene.rho > 0]))#44, 32, 18, 11, 9
length(unique(Validated_pairs$Gene_name[Validated_pairs$DM.median < 0 & Validated_pairs$MethyandGene.rho > 0]))

#Validated_pairs <- cbind(Validated_pairs, Enhancer_Gene_pairs[match(Validated_pairs$DMRID, Enhancer_Gene_pairs$DMR.ID), 
#                                                              c("distanceToDMR", "lfc.GSE5230", "lfc.GSE35311", "lfc.GSE112788", "lfc.GSE67318", "lfc.up", "lfc.down", "EnhancerScore.T", "EnhancerScore.N", "EnhancerScore")])
#
write.csv(Validated_pairs, file="/data/huangp/Methy/Mixed/New/LIHC_Validated_IntergenicEnhancer_Like_DMR_Gene_pairs_0.7.csv",  row.names=F, quote=T)

Enhancer_Gene_pairs$LIHCValidation <- Enhancer_Gene_pairs$gene_name %in% Validated_pairs$Gene_name
Prognostic <- read.csv("/data/huangp/Methy/Mixed/New/The185PrognosticGenes.csv", header=T)[[2]]
Enhancer_Gene_pairs$Prognositc <- Enhancer_Gene_pairs$gene_name %in% Prognostic
write.csv(Enhancer_Gene_pairs, file="/data/huangp/Methy/Mixed/New/HighConfidence_InterGenicEnhancer_Like_DMR_Gene_pairs_0.7.csv", row.names = F)
########################## Clinical relevance of high confident cRE-like DMRs in discovery cohort ######################
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/New")#29
library(fastSave)
load("../All66HCC.phenotypes.RData")#load clinical phenotypes
#load group level methylation
load.lbzip2("23.8k_DMR.enhancer.methy.matrix.RDataFS", n.cores=24)
load.lbzip2("47.5k_DMR.promoter.methy.matrix.RDataFS", n.cores=24)
load.lbzip2("IntergenicDMR.eRNA.methy.matrix.RDataFS", n.cores = 24)
load.lbzip2("../NewSmoothed.Methylation.qc.RDataFS", n.cores =64)#load smoothed methylation
Meth.sm <- t(Meth.all.qc)
rm(Meth.all.qc)
load.lbzip2("NewRaw.Methylation.masked.RDataFS", n.cores=32)#load raw methylation
Meth.raw <- Meth.raw[rownames(Meth.sm), colnames(Meth.sm)]
#load gene expression
load.lbzip2("../Gene.TPM.filtered.RDataFS", n.cores=32)
library(genomation)
library(GenomicFeatures)
Granges.gtf <- gffToGRanges("../gencode.v29.annotation.gtf")
Granges.gene <- Granges.gtf[Granges.gtf$type=="gene",]
Granges.gene <- Granges.gene[match(rownames(gene.TPM), Granges.gene$gene_id),]
range_gene_DF <- as.data.frame(Granges.gene)
#load high confident promoter-like DMRs and enhancer-like DMRs that show the potential to regulate gene transcription
# promoter-like DMR
promoterDMR <- read.csv("HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", header=T)
library(foreach)
ClinicalRevelance.promoterDMR <- foreach(i = 1:nrow(promoterDMR), .combine = 'rbind')%do%
  {
    DMR.ID <- promoterDMR$DMR.ID[i]
    Gene <- promoterDMR$gene_name[i]
    Gene.ID <- range_gene_DF$gene_id[match(Gene, range_gene_DF$gene_name)]
    methy <- DMR.promoter.methy.matrix[,DMR.ID]
    expr <- gene.TPM[Gene.ID,] 
    expr <- expr[names(methy)]
    df <- data.frame(methy = methy, expr = expr, differentiation = HCC.phenotypes$differentiation[match(names(methy), HCC.phenotypes$SampleID)], 
                     tumorSize = HCC.phenotypes$tumor_size[match(names(methy), HCC.phenotypes$SampleID)], 
                     Condition = HCC.phenotypes$Condition[match(names(methy), HCC.phenotypes$SampleID)])
    df <- df[df$Condition == "T",]
    df$expr <- log(df$expr+0.0001, 2)
    # ANOVA for tumor differentialtion
    differentiation.methy.p <- signif(summary(aov(methy~differentiation, data = df))[[1]][[5]][1], digits = 2)
    differentiation.expr.p <- signif(summary(aov(expr~differentiation, data = df))[[1]][[5]][1], digits = 2)
    Mean.methy.low <- round(mean(df$methy[df$differentiation=="low"]), digits = 2)
    Mean.methy.lowtomoderate <- round(mean(df$methy[df$differentiation == "low-to-moderate"]),2)
    Mean.methy.moderate <- round(mean(df$methy[df$differentiation == "moderate"]),2)
    Mean.expr.low <- round(mean(df$expr[df$differentiation=="low"]), digits = 2)
    Mean.expr.lowtomoderate <- round(mean(df$expr[df$differentiation == "low-to-moderate"]),2)
    Mean.expr.moderate <- round(mean(df$expr[df$differentiation == "moderate"]),2)
    # Spearman correlation for tumor size
    tumorSize.methy.cor <- round(cor.test(df$methy, df$tumorSize, method = "spearman")$estimate, 2)
    tumorSize.methy.p <- signif(cor.test(df$methy, df$tumorSize, method = "spearman")$p.value, 2)
    tumorSize.expr.cor <- round(cor.test(df$expr, df$tumorSize, method = "spearman")$estimate, 2)
    tumorSize.expr.p <- signif(cor.test(df$expr, df$tumorSize, method = "spearman")$p.value, 2)
    c(Mean.methy.low, Mean.methy.lowtomoderate, Mean.methy.moderate, differentiation.methy.p, Mean.expr.low, Mean.expr.lowtomoderate, Mean.expr.moderate, differentiation.expr.p,
      tumorSize.methy.cor, tumorSize.methy.p, tumorSize.expr.cor, tumorSize.expr.p)
  }
colnames(ClinicalRevelance.promoterDMR) <- c("Mean.methy.low", "Mean.methy.lowtomoderate", "Mean.methy.moderate", "differentiation.methy.pvalue", 
                                             "Mean.expr.low", "Mean.expr.lowtomoderate", "Mean.expr.moderate", "differentiation.expr.pvalue",
                                             "tumorSize.methy.cor", "tumorSize.methy.pvalue", "tumorSize.expr.cor", "tumorSize.expr.pvalue")
ClinicalRevelance.promoterDMR <- cbind(promoterDMR[, c("DMR.ID", "Median.delta", "gene_name", "rho")], ClinicalRevelance.promoterDMR)
ClinicalRevelance.promoterDMR$Median.delta <- round(ClinicalRevelance.promoterDMR$Median.delta, 2)
# enhancer-like DMR
GenicEnhancerDMR <- read.csv("HighConfidence_GenicEnhancer_Like_DMR_Gene_pairs_0.8.csv", header=T)
library(foreach)
ClinicalRevelance.GenicEnhancerDMR <- foreach(i = 1:nrow(GenicEnhancerDMR), .combine = 'rbind')%do%
  {
    DMR.ID <- GenicEnhancerDMR$DMR.ID[i]
    Gene <- GenicEnhancerDMR$gene_name[i]
    Gene.ID <- range_gene_DF$gene_id[match(Gene, range_gene_DF$gene_name)]
    methy <- DMR.enhancer.methy.matrix[,DMR.ID]
    expr <- gene.TPM[Gene.ID,] 
    expr <- expr[names(methy)]
    df <- data.frame(methy = methy, expr = expr, differentiation = HCC.phenotypes$differentiation[match(names(methy), HCC.phenotypes$SampleID)], 
                     tumorSize = HCC.phenotypes$tumor_size[match(names(methy), HCC.phenotypes$SampleID)], 
                     Condition = HCC.phenotypes$Condition[match(names(methy), HCC.phenotypes$SampleID)])
    df <- df[df$Condition == "T",]
    df$expr <- log(df$expr+0.0001, 2)
    # ANOVA for tumor differentialtion
    differentiation.methy.p <- signif(summary(aov(methy~differentiation, data = df))[[1]][[5]][1], digits = 2)
    differentiation.expr.p <- signif(summary(aov(expr~differentiation, data = df))[[1]][[5]][1], digits = 2)
    Mean.methy.low <- round(mean(df$methy[df$differentiation=="low"]), digits = 2)
    Mean.methy.lowtomoderate <- round(mean(df$methy[df$differentiation == "low-to-moderate"]),2)
    Mean.methy.moderate <- round(mean(df$methy[df$differentiation == "moderate"]),2)
    Mean.expr.low <- round(mean(df$expr[df$differentiation=="low"]), digits = 2)
    Mean.expr.lowtomoderate <- round(mean(df$expr[df$differentiation == "low-to-moderate"]),2)
    Mean.expr.moderate <- round(mean(df$expr[df$differentiation == "moderate"]),2)
    # Spearman correlation for tumor size
    tumorSize.methy.cor <- round(cor.test(df$methy, df$tumorSize, method = "spearman")$estimate, 2)
    tumorSize.methy.p <- signif(cor.test(df$methy, df$tumorSize, method = "spearman")$p.value, 2)
    tumorSize.expr.cor <- round(cor.test(df$expr, df$tumorSize, method = "spearman")$estimate, 2)
    tumorSize.expr.p <- signif(cor.test(df$expr, df$tumorSize, method = "spearman")$p.value, 2)
    c(Mean.methy.low, Mean.methy.lowtomoderate, Mean.methy.moderate, differentiation.methy.p, Mean.expr.low, Mean.expr.lowtomoderate, Mean.expr.moderate, differentiation.expr.p,
      tumorSize.methy.cor, tumorSize.methy.p, tumorSize.expr.cor, tumorSize.expr.p)
  }
colnames(ClinicalRevelance.GenicEnhancerDMR) <- c("Mean.methy.low", "Mean.methy.lowtomoderate", "Mean.methy.moderate", "differentiation.methy.pvalue", 
                                                  "Mean.expr.low", "Mean.expr.lowtomoderate", "Mean.expr.moderate", "differentiation.expr.pvalue",
                                                  "tumorSize.methy.cor", "tumorSize.methy.pvalue", "tumorSize.expr.cor", "tumorSize.expr.pvalue")
ClinicalRevelance.GenicEnhancerDMR <- cbind(GenicEnhancerDMR[, c("DMR.ID", "Median.delta", "gene_name", "rho")], ClinicalRevelance.GenicEnhancerDMR)
ClinicalRevelance.GenicEnhancerDMR$Median.delta <- round(ClinicalRevelance.GenicEnhancerDMR$Median.delta, 2)
#intergenic-like DMR
InterGenicEnhancerDMR <- read.csv("HighConfidence_InterGenicEnhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)
library(foreach)
ClinicalRevelance.InterGenicEnhancerDMR <- foreach(i = 1:nrow(InterGenicEnhancerDMR), .combine = 'rbind')%do%
  {
    DMR.ID <- InterGenicEnhancerDMR$DMR.ID[i]
    Gene <- InterGenicEnhancerDMR$gene_name[i]
    Gene.ID <- range_gene_DF$gene_id[match(Gene, range_gene_DF$gene_name)]
    methy <- IntergenicDMR.eRNA.methy.matrix[,DMR.ID]
    expr <- gene.TPM[Gene.ID,] 
    expr <- expr[names(methy)]
    df <- data.frame(methy = methy, expr = expr, differentiation = HCC.phenotypes$differentiation[match(names(methy), HCC.phenotypes$SampleID)], 
                     tumorSize = HCC.phenotypes$tumor_size[match(names(methy), HCC.phenotypes$SampleID)], 
                     Condition = HCC.phenotypes$Condition[match(names(methy), HCC.phenotypes$SampleID)])
    df <- df[df$Condition == "T",]
    df$expr <- log(df$expr+0.0001, 2)
    # ANOVA for tumor differentialtion
    differentiation.methy.p <- signif(summary(aov(methy~differentiation, data = df))[[1]][[5]][1], digits = 2)
    differentiation.expr.p <- signif(summary(aov(expr~differentiation, data = df))[[1]][[5]][1], digits = 2)
    Mean.methy.low <- round(mean(df$methy[df$differentiation=="low"]), digits = 2)
    Mean.methy.lowtomoderate <- round(mean(df$methy[df$differentiation == "low-to-moderate"]),2)
    Mean.methy.moderate <- round(mean(df$methy[df$differentiation == "moderate"]),2)
    Mean.expr.low <- round(mean(df$expr[df$differentiation=="low"]), digits = 2)
    Mean.expr.lowtomoderate <- round(mean(df$expr[df$differentiation == "low-to-moderate"]),2)
    Mean.expr.moderate <- round(mean(df$expr[df$differentiation == "moderate"]),2)
    # Spearman correlation for tumor size
    tumorSize.methy.cor <- round(cor.test(df$methy, df$tumorSize, method = "spearman")$estimate, 2)
    tumorSize.methy.p <- signif(cor.test(df$methy, df$tumorSize, method = "spearman")$p.value, 2)
    tumorSize.expr.cor <- round(cor.test(df$expr, df$tumorSize, method = "spearman")$estimate, 2)
    tumorSize.expr.p <- signif(cor.test(df$expr, df$tumorSize, method = "spearman")$p.value, 2)
    c(Mean.methy.low, Mean.methy.lowtomoderate, Mean.methy.moderate, differentiation.methy.p, Mean.expr.low, Mean.expr.lowtomoderate, Mean.expr.moderate, differentiation.expr.p,
      tumorSize.methy.cor, tumorSize.methy.p, tumorSize.expr.cor, tumorSize.expr.p)
  }
colnames(ClinicalRevelance.InterGenicEnhancerDMR) <- c("Mean.methy.low", "Mean.methy.lowtomoderate", "Mean.methy.moderate", "differentiation.methy.pvalue", 
                                                       "Mean.expr.low", "Mean.expr.lowtomoderate", "Mean.expr.moderate", "differentiation.expr.pvalue",
                                                       "tumorSize.methy.cor", "tumorSize.methy.pvalue", "tumorSize.expr.cor", "tumorSize.expr.pvalue")
ClinicalRevelance.InterGenicEnhancerDMR <- cbind(InterGenicEnhancerDMR[, c("DMR.ID", "DM.Median", "gene_name", "MethyandGene.rho")], ClinicalRevelance.InterGenicEnhancerDMR)
colnames(ClinicalRevelance.InterGenicEnhancerDMR)[c(2,4)] <- c("Median.delta", "rho") 
ClinicalRevelance.InterGenicEnhancerDMR$Median.delta <- round(ClinicalRevelance.InterGenicEnhancerDMR$Median.delta, 2)

# combined 
ClinicalRevelance.cREDMR <- rbind(ClinicalRevelance.promoterDMR, ClinicalRevelance.GenicEnhancerDMR, ClinicalRevelance.InterGenicEnhancerDMR)
ClinicalRevelance.cREDMR$type = c(rep("promoter", nrow(promoterDMR)), rep("GenicEnhancer", nrow(GenicEnhancerDMR)), rep("InterGenicEnhancer", nrow(InterGenicEnhancerDMR)))
attach(ClinicalRevelance.cREDMR)
table(differentiation.methy.pvalue < 0.05 & differentiation.expr.pvalue < 0.05)
table(tumorSize.methy.pvalue < 0.05 & tumorSize.expr.pvalue < 0.05)
detach(ClinicalRevelance.cREDMR)
write.csv(ClinicalRevelance.cREDMR, file = "Clinical_Relevance_of_cRE_DMR_in_discovery_cohort.csv", row.names = F)
############################# Clinical relevance of LIHC validated high-confident cRE-like DMRs in TCGA ###############
options(stringsAsFactors = F)
setwd("/data/huangp/Methy/Mixed/LIHC")
#load Overall Survival, Progression-free Survival phenotypes of TCGA
if(1<0){survival.phenotypes <- read.csv("SurvivalPhenotypes.csv", header=T)
survival.phenotypes$OS.time[survival.phenotypes$OS.time == "#N/A"] <- NA
survival.phenotypes$PFI.time[survival.phenotypes$PFI.time == "#N/A"] <- NA
survival.phenotypes$OS.time <- as.integer(survival.phenotypes$OS.time)
survival.phenotypes$PFI.time <- as.integer(survival.phenotypes$PFI.time)
load("LIHC_TCGA.phenotypes.RData")
LIHC.data <- methy.data
rm(methy.data)
LIHC.data <- cbind(LIHC.data, survival.phenotypes[match(LIHC.data$barcode, survival.phenotypes$bcr_patient_barcode),c("OS", "OS.time", "PFI", "PFI.time")])
LIHC.data$stage <- survival.phenotypes$ajcc_pathologic_tumor_stage[match(LIHC.data$barcode, survival.phenotypes$bcr_patient_barcode)]
LIHC.data$stage[LIHC.data$stage == "[Discrepancy]"] <-  NA
LIHC.data$stage[LIHC.data$stage == "[Not Available]"] <- NA
LIHC.data$stage2 <- LIHC.data$stage
LIHC.data$stage2[LIHC.data$stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
LIHC.data$stage2[LIHC.data$stage %in% c("Stage IVA", "Stage IVB")] <- "Stage IV"
LIHC.data$stage2 <- factor(LIHC.data$stage2, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
save(LIHC.data, file="LIHC_TCGA.phenotypes.RData")
}
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
LIHC.expr <- z #20531Genes Gene_name+Gene_ID+423Sample

#for 661 high confident methylation driven genes, examine their prognostic significance and progression association at expression level
#load all 661 high confident methylation driven genes
Promoter_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/HighConfidence_Promoter_Like_DMR_Gene_pairs_0.8.csv", header=T)
InterGenicEnhancer_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/HighConfidence_Intergenic_Enhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)
GenicEnhancer_Gene_pairs <- read.csv("/data/huangp/Methy/Mixed/New/XXXHighConfidence_genic_Enhancer_Like_DMR_Gene_pairs_0.8.csv", header=T)
x <- unique(Promoter_Gene_pairs$gene_name)#335
y <- unique(GenicEnhancer_Gene_pairs$gene_name)#834
z <- unique(InterGenicEnhancer_Gene_pairs$gene_name)#186
MethyGenes <- union(x, union(y,z))#591

library(foreach)
ClinicalRelevance.expr <- foreach(i = 1:length(MethyGenes), .combine = 'rbind')%do%{
  Gene = MethyGenes[i]
  if(Gene %in% LIHC.expr$Gene_name)
  {
    expr = unlist(LIHC.expr[match(Gene, LIHC.expr$Gene_name),-(1:2)])
    #filtered genes with too much NAs or 0s
    if(sum(is.na(expr)) <= length(expr)/2 | sum(expr == 0) <= length(expr)/2)
    {
      expr <- log2(expr+0.0001)
      library(survival)
      library(survminer)
      # cox for expr
      df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
      df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
      expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
      df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
      df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
      df.expr <- df.expr[df.expr$condition == "Tumor" & df.expr$stage2 != "Stage IV",]#keep only tumor samples and stage 1-3 (only 5 stage 4 samples)
      df.expr <- df.expr[!is.na(df.expr$age) & !is.na(df.expr$race) & !is.na(df.expr$gender) & !is.na(df.expr$stage2) & !is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                           !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr) & !is.na(df.expr$expr.scaled),]
      OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
      OS.pvalue.expr <- signif(summary(OS.cox)[[7]][8,5],2)
      OS.HR.expr <- round(summary(OS.cox)[[8]][8,1],2)
      PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
      PFS.pvalue.expr <- signif(summary(PFS.cox)[[7]][8,5], 2)
      PFS.HR.expr <- round(summary(PFS.cox)[[8]][8,1], 2)
      expr.survival <- c(OS.HR.expr, OS.pvalue.expr, PFS.HR.expr, PFS.pvalue.expr)
      names(expr.survival) <- c("OS.HR.expr", "OS.pvalue.expr", "PFS.HR.expr", "PFS.pvalue.expr")
      #  ANOVA of expression for tumor stage
      df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
      df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
      expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
      df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
      df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
      df.expr <- df.expr[df.expr$condition == "Tumor" & df.expr$stage2 != "Stage IV",]#keep only tumor samples
      stage.expr.p <- signif(summary(aov(expr ~ stage2, data = df.expr))[[1]][[5]][1], digits = 2)
      stageI.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage I"], na.rm=T),2)
      stageII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage II"], na.rm=T),2)
      stageIII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage III"], na.rm=T),2)
      #stageIV.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage IV"], na.rm=T),2)
      stage.expr <- c(stageI.expr, stageII.expr, stageIII.expr, stage.expr.p)
      names(stage.expr) <- c("mean.expr.stageI", "mean.expr.stageII", "mean.expr.stageIII", "stage.pvalue.expr")
      CR <- c(Gene, expr.survival,  stage.expr)
      CR
    }
  }
}
colnames(ClinicalRelevance.expr)[1] <- "Gene"
ClinicalRelevance.expr <- as.data.frame(ClinicalRelevance.expr)
for(j in 2:ncol(ClinicalRelevance.expr))
{
  ClinicalRelevance.expr[,j] <- as.numeric(ClinicalRelevance.expr[,j])
}
str(ClinicalRelevance.expr)
rm(list=c("Gene", "OS.HR.expr", "OS.pvalue.expr", "PFS.HR.expr", "PFS.pvalue.expr"))

attach(ClinicalRelevance.expr)
table(OS.pvalue.expr <= 0.05)#172
table(PFS.pvalue.expr <= 0.05)#140
table(stage.pvalue.expr <= 0.05)#276
table(OS.pvalue.expr <= 0.05 | PFS.pvalue.expr <= 0.05)#194
table((OS.pvalue.expr <= 0.05 | PFS.pvalue.expr <= 0.05) & stage.pvalue.expr <= 0.05)#104
table(OS.pvalue.expr <= 0.05 & PFS.pvalue.expr <= 0.05)#118
table((OS.pvalue.expr <= 0.05 & PFS.pvalue.expr <= 0.05) & stage.pvalue.expr <= 0.05)#68
detach(ClinicalRelevance.expr)
#describe the direction of each prognostic biomarker
#obtain the LFC and pairType for each gene
DMRType <- vector(mode="character", length=nrow(ClinicalRelevance.expr))
pairType <- vector(mode="character", length=nrow(ClinicalRelevance.expr))
LFC.median <- vector(mode = "numeric", length=nrow(ClinicalRelevance.expr))
for(i in 1:nrow(ClinicalRelevance.expr))
{
  Gene = ClinicalRelevance.expr$Gene[i]
  if(Gene %in% Promoter_Gene_pairs$gene_name)
  {
    LFC.median[i] = Promoter_Gene_pairs$LFC.median[match(Gene, Promoter_Gene_pairs$gene_name)] 
    pairType[i] = Promoter_Gene_pairs$pairType[match(Gene, Promoter_Gene_pairs$gene_name)]
    DMRType[i] = "Promoter"
  }
  else if(Gene %in% InterGenicEnhancer_Gene_pairs$gene_name)
  {
    LFC.median[i] = InterGenicEnhancer_Gene_pairs$LFC.median[match(Gene, InterGenicEnhancer_Gene_pairs$gene_name)]
    pairType[i] = InterGenicEnhancer_Gene_pairs$pairType[match(Gene, InterGenicEnhancer_Gene_pairs$gene_name)]
    DMRType[i] = "IntergenicEnhancer"
  }
  else {
    LFC.median[i] = GenicEnhancer_Gene_pairs$LFC.median[match(Gene, GenicEnhancer_Gene_pairs$gene_name)]
    pairType[i] = GenicEnhancer_Gene_pairs$pairType[match(Gene, GenicEnhancer_Gene_pairs$gene_name)]
    DMRType[i] = "GenicEnhancer"
  }
}
ClinicalRelevance.expr <- cbind(Gene = ClinicalRelevance.expr[[1]], DMRType = DMRType,
                                pairType = pairType, LFC.median = LFC.median, ClinicalRelevance.expr[,-1])
rm(list=c("DMRType", "Gene", "LFC.median", "pairType"))
OS.effect <- vector(mode="character", length = nrow(ClinicalRelevance.expr))
PFS.effect <- vector(mode = "character", length = nrow(ClinicalRelevance.expr))
Stage.association <- vector(mode="character", length = nrow(ClinicalRelevance.expr))
attach(ClinicalRelevance.expr)
OS.effect[OS.pvalue.expr > 0.05] <- "N.S."
OS.effect[OS.pvalue.expr <= 0.05 & OS.HR.expr > 1] <- "Unfavorable"
OS.effect[OS.pvalue.expr <= 0.05 & OS.HR.expr < 1] <- "Favorable"
PFS.effect[PFS.pvalue.expr > 0.05] <- "N.S."
PFS.effect[PFS.pvalue.expr <= 0.05 & PFS.HR.expr > 1] <- "Unfavorable"
PFS.effect[PFS.pvalue.expr <= 0.05 & PFS.HR.expr < 1] <- "Favorable"
Stage.association[stage.pvalue.expr > 0.05] <- "N.S."
Stage.association[stage.pvalue.expr <= 0.05 & (mean.expr.stageI - mean.expr.stageIII) > 0] <- "Inverse"
Stage.association[stage.pvalue.expr <= 0.05 & (mean.expr.stageI - mean.expr.stageIII) < 0] <- "Positive"
#Stage.association[stage.pvalue.expr <= 0.05 & (mean.expr.stageI - mean.expr.stageII) >= 0 & (mean.expr.stageII - mean.expr.stageIII) >= 0] <- "Suppressor"
#Stage.association[stage.pvalue.expr <= 0.05 & (mean.expr.stageI - mean.expr.stageII) <= 0 & (mean.expr.stageII - mean.expr.stageIII) <= 0] <- "Oncogene"
detach(ClinicalRelevance.expr)
ClinicalRelevance.expr <- cbind(ClinicalRelevance.expr[,1:6], OS.effect = OS.effect, 
                                ClinicalRelevance.expr[,7:8], PFS.effect = PFS.effect,
                                ClinicalRelevance.expr[,9:12], Stage.effect = Stage.association)
rm(list=c("OS.effect", "PFS.effect", "Stage.association"))
attach(ClinicalRelevance.expr)
table(pairType[OS.effect == "Favorable"])
#HyperDown  HypoDown    HypoUp
#14        12         4
table(pairType[OS.effect == "Unfavorable"])
#HyperUp HypoDown   HypoUp
#20        1      121
table(pairType[PFS.effect == "Favorable"])
#HyperDown   HyperUp  HypoDown    HypoUp
#8         3        16         1
table(pairType[PFS.effect == "Unfavorable"])
#HyperUp  HypoUp
#16      96
table(pairType[Stage.effect == "Inverse"])
#HyperDown   HyperUp  HypoDown    HypoUp
#21         6        41        50
table(pairType[Stage.effect == "Positive"])
#HyperDown   HyperUp  HypoDown    HypoUp
#1        17         4       136
### adjust the inconsistent OS effect 
ClinicalRelevance.expr$OS.effect[OS.effect == "Favorable" & LFC.median > 0] <- "Inconsistent"
ClinicalRelevance.expr$OS.effect[OS.effect == "Unfavorable" & LFC.median < 0] <- "Inconsistent"
### adjust the inconsistent PFS effect 
ClinicalRelevance.expr$PFS.effect[PFS.effect == "Favorable" & LFC.median > 0] <- "Inconsistent"
ClinicalRelevance.expr$PFS.effect[PFS.effect == "Unfavorable" & LFC.median < 0] <- "Inconsistent"
### adjust the inconsistent tumor progression effect 
ClinicalRelevance.expr$Stage.effect[Stage.effect == "Inverse" & LFC.median > 0] <- "Inconsistent"
ClinicalRelevance.expr$Stage.effect[Stage.effect == "Positive" & LFC.median < 0] <- "Inconsistent"
detach(ClinicalRelevance.expr)
OS.biomarker <- ClinicalRelevance.expr$Gene[ClinicalRelevance.expr$OS.effect %in% c("Favorable", "Unfavorable")]
#168
PFS.biomarker <-  ClinicalRelevance.expr$Gene[ClinicalRelevance.expr$PFS.effect %in% c("Favorable", "Unfavorable")]
#135
Stage.biomarker <- ClinicalRelevance.expr$Gene[ClinicalRelevance.expr$Stage.effect %in% c("Inverse", "Positive")]
#217
length(union(OS.biomarker, PFS.biomarker))# 186
length(union(union(OS.biomarker, PFS.biomarker),Stage.biomarker))#303
length(intersect(union(OS.biomarker, PFS.biomarker),Stage.biomarker))#100
length(intersect(intersect(OS.biomarker, PFS.biomarker),Stage.biomarker))#66
#save ClinicalRelevance.expr
write.csv(ClinicalRelevance.expr, file="/data/huangp/Methy/Mixed/New/ClinicalRelevance_highconfident_methylation_driven_genes_at_expression_level.csv", row.names=F, quote=F)

#for 1335 high confident methylation driven genes, examine their prognostic significance and progression association at methylation level
Validated_promoter_pairs <- read.csv("/data/huangp/Methy/Mixed/New/LIHC_Validated_Promoter_Like_DMR_Gene_pairs_0.8.csv", header=T)
Validated_enhancer_pairs <- read.csv("/data/huangp/Methy/Mixed/New/LIHC_Validated_Enhancer_Like_GenicDMR_Gene_pairs_0.8.csv", header=T)
Validated_intergenic_enhancer_pairs <- read.csv("/data/huangp/Methy/Mixed/New/LIHC_Validated_IntergenicEnhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)
X <- Validated_promoter_pairs[,c("DMRID", "Gene_name", "pairID", "DM.median", "MethyandGene.rho", "CpG_pos")]
Y <- Validated_enhancer_pairs[,c("DMRID", "Gene_name", "pairID", "DM.median", "MethyandGene.rho", "CpG_pos")]
Z <- Validated_intergenic_enhancer_pairs[,c("DMRID", "Gene_name", "pairID", "DM.median", "MethyandGene.rho", "CpG_pos")]
Validated_pairs <- rbind(X, Y, Z)
Validated_pairs$DMRType <- c(rep("Promoter", nrow(Validated_promoter_pairs)), rep("GenicEnhancer", nrow(Validated_enhancer_pairs)),
                             rep("IntergenicEnhancer", nrow(Validated_intergenic_enhancer_pairs)))
MethyCpGs <- unique(Validated_pairs$CpG_pos)

library(foreach)
ClinicalRelevance.methy <- foreach(i = 1:length(MethyCpGs), .combine = 'rbind')%do%{
  CpG = MethyCpGs[i]
  if(CpG %in% rownames(LIHC.methy))
  {
    methy <- LIHC.methy[CpG,]
    #filtered CpG with too much NAs 
    if(sum(is.na(methy)) <= length(methy)/2)
    {
      
      library(survival)
      library(survminer)
      # cox for methy
      df.methy <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
      df.methy <- df.methy[df.methy$sampleID %in% names(methy),]#keep the sample with both expression and phenotypes
      methy.scaled <- scale(methy)[,1]# transform the expression data to Z score
      df.methy$methy <- round(methy[match(df.methy$sampleID, names(methy))],2)
      df.methy$methy.scaled <- round(methy.scaled[match(df.methy$sampleID, names(methy.scaled))],2)
      df.methy <- df.methy[df.methy$condition == "Tumor" & df.methy$stage2 != "Stage IV",]#keep only tumor samples and stage 1-3 (only 5 stage 4 samples)
      df.methy <- df.methy[!is.na(df.methy$age) & !is.na(df.methy$race) & !is.na(df.methy$gender) & !is.na(df.methy$stage2) & !is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                           !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy) & !is.na(df.methy$methy.scaled),]
      OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
      OS.pvalue.methy <- signif(summary(OS.cox)[[7]][8,5],2)
      OS.HR.methy <- round(summary(OS.cox)[[8]][8,1],2)
      PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
      PFS.pvalue.methy <- signif(summary(PFS.cox)[[7]][8,5], 2)
      PFS.HR.methy <- round(summary(PFS.cox)[[8]][8,1], 2)
      methy.survival <- c(OS.HR.methy, OS.pvalue.methy, PFS.HR.methy, PFS.pvalue.methy)
      names(methy.survival) <- c("OS.HR.methy", "OS.pvalue.methy", "PFS.HR.methy", "PFS.pvalue.methy")
      #  ANOVA of methy for tumor stage
      df.methy <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
      df.methy <- df.methy[df.methy$sampleID %in% names(methy),]#keep the sample with both methyession and phenotypes
      methy.scaled <- scale(methy)[,1]# transform the methyession data to Z score
      df.methy$methy <- round(methy[match(df.methy$sampleID, names(methy))],2)
      df.methy$methy.scaled <- round(methy.scaled[match(df.methy$sampleID, names(methy.scaled))],2)
      df.methy <- df.methy[df.methy$condition == "Tumor" & df.methy$stage2 != "Stage IV",]#keep only tumor samples
      stage.methy.p <- signif(summary(aov(methy ~ stage2, data = df.methy))[[1]][[5]][1], digits = 2)
      stageI.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage I"], na.rm=T),2)
      stageII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage II"], na.rm=T),2)
      stageIII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage III"], na.rm=T),2)
      #stageIV.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage IV"], na.rm=T),2)
      stage.methy <- c(stageI.methy, stageII.methy, stageIII.methy, stage.methy.p)
      names(stage.methy) <- c("mean.methy.stageI", "mean.methy.stageII", "mean.methy.stageIII", "stage.pvalue.methy")
      CR <- c(CpG, methy.survival,  stage.methy)
      CR
    }
  }
}
colnames(ClinicalRelevance.methy)[1] <- "CpG"
ClinicalRelevance.methy <- as.data.frame(ClinicalRelevance.methy)
for(j in 2:ncol(ClinicalRelevance.methy))
{
  ClinicalRelevance.methy[,j] <- as.numeric(ClinicalRelevance.methy[,j])
}
str(ClinicalRelevance.methy)
rm(list=c("CpG", "OS.HR.methy", "OS.pvalue.methy", "PFS.HR.methy", "PFS.pvalue.methy"))

attach(ClinicalRelevance.methy)
table(OS.pvalue.methy <= 0.05)#48
table(PFS.pvalue.methy <= 0.05)#61
table(stage.pvalue.methy <= 0.05)#41
table(OS.pvalue.methy <= 0.05 | PFS.pvalue.methy <= 0.05)#72
table((OS.pvalue.methy <= 0.05 | PFS.pvalue.methy <= 0.05) & stage.pvalue.methy <= 0.05)#10
table(OS.pvalue.methy <= 0.05 & PFS.pvalue.methy <= 0.05)#37
table((OS.pvalue.methy <= 0.05 & PFS.pvalue.methy <= 0.05) & stage.pvalue.methy <= 0.05)#8
detach(ClinicalRelevance.methy)
#describe the direction of each prognostic biomarker
#obtain the LFC and pairType for each gene
DMRType <- vector(mode="character", length=nrow(ClinicalRelevance.methy))
DM.median <- vector(mode = "numeric", length=nrow(ClinicalRelevance.methy))
for(i in 1:nrow(ClinicalRelevance.methy))
{
  CpG = ClinicalRelevance.methy$CpG[i]
  if(CpG %in% Validated_promoter_pairs$CpG_pos)
  {
    DM.median[i] = Validated_promoter_pairs$DM.median[match(CpG, Validated_promoter_pairs$CpG_pos)] 
    DMRType[i] = "Promoter"
  }
  else if(CpG %in% Validated_intergenic_enhancer_pairs$CpG_pos)
  {
    DM.median[i] = Validated_intergenic_enhancer_pairs$DM.median[match(CpG, Validated_intergenic_enhancer_pairs$CpG_pos)]
    DMRType[i] = "IntergenicEnhancer"
  }
  else {
    DM.median[i] = Validated_enhancer_pairs$DM.median[match(CpG, Validated_enhancer_pairs$CpG_pos)]
    DMRType[i] = "GenicEnhancer"
  }
}
ClinicalRelevance.methy <- cbind(CpG = ClinicalRelevance.methy[[1]], DMRType = DMRType,
                                 DM.median = DM.median, ClinicalRelevance.methy[,-1])
rm(list=c("DMRType", "CpG", "DM.median"))
OS.effect <- vector(mode="character", length = nrow(ClinicalRelevance.methy))
PFS.effect <- vector(mode = "character", length = nrow(ClinicalRelevance.methy))
Stage.association <- vector(mode="character", length = nrow(ClinicalRelevance.methy))
attach(ClinicalRelevance.methy)
OS.effect[OS.pvalue.methy > 0.05] <- "N.S."
OS.effect[OS.pvalue.methy <= 0.05 & OS.HR.methy > 1] <- "Unfavorable"
OS.effect[OS.pvalue.methy <= 0.05 & OS.HR.methy < 1] <- "Favorable"
PFS.effect[PFS.pvalue.methy > 0.05] <- "N.S."
PFS.effect[PFS.pvalue.methy <= 0.05 & PFS.HR.methy > 1] <- "Unfavorable"
PFS.effect[PFS.pvalue.methy <= 0.05 & PFS.HR.methy < 1] <- "Favorable"
Stage.association[stage.pvalue.methy > 0.05] <- "N.S."
Stage.association[stage.pvalue.methy <= 0.05 & (mean.methy.stageI - mean.methy.stageIII) > 0] <- "Inverse"
Stage.association[stage.pvalue.methy <= 0.05 & (mean.methy.stageI - mean.methy.stageIII) < 0] <- "Positive"
#Stage.association[stage.pvalue.methy <= 0.05 & (mean.methy.stageI - mean.methy.stageII) >= 0 & (mean.methy.stageII - mean.methy.stageIII) >= 0] <- "Suppressor"
#Stage.association[stage.pvalue.methy <= 0.05 & (mean.methy.stageI - mean.methy.stageII) <= 0 & (mean.methy.stageII - mean.methy.stageIII) <= 0] <- "Oncogene"
detach(ClinicalRelevance.methy)
ClinicalRelevance.methy <- cbind(ClinicalRelevance.methy[,1:5], OS.effect = OS.effect, 
                                ClinicalRelevance.methy[,6:7], PFS.effect = PFS.effect,
                                ClinicalRelevance.methy[,8:11], Stage.effect = Stage.association)
rm(list=c("OS.effect", "PFS.effect", "Stage.association"))
attach(ClinicalRelevance.methy)
table(DM.median[OS.effect == "Favorable"] < 0)
#TRUE
#38
table(DM.median[OS.effect == "Unfavorable"] < 0)
#FALSE
#10
table(DM.median[PFS.effect == "Favorable"] < 0)
#TRUE
#41
table(DM.median[PFS.effect == "Unfavorable"] < 0)
#FALSE
#20
table(DM.median[Stage.effect == "Inverse"] < 0)
#FALSE TRUE
#4   20
table(DM.median[Stage.effect == "Positive"] < 0)
#FALSE  TRUE
#10     6
### adjust the inconsistent OS effect 
ClinicalRelevance.methy$OS.effect[OS.effect == "Favorable" & DM.median > 0] <- "Inconsistent"
ClinicalRelevance.methy$OS.effect[OS.effect == "Unfavorable" & DM.median < 0] <- "Inconsistent"
### adjust the inconsistent PFS effect 
ClinicalRelevance.methy$PFS.effect[PFS.effect == "Favorable" & DM.median > 0] <- "Inconsistent"
ClinicalRelevance.methy$PFS.effect[PFS.effect == "Unfavorable" & DM.median < 0] <- "Inconsistent"
### adjust the inconsistent tumor progression effect 
ClinicalRelevance.methy$Stage.effect[Stage.effect == "Inverse" & DM.median > 0] <- "Inconsistent"
ClinicalRelevance.methy$Stage.effect[Stage.effect == "Positive" & DM.median < 0] <- "Inconsistent"
detach(ClinicalRelevance.methy)
OS.biomarker <- ClinicalRelevance.methy$CpG[ClinicalRelevance.methy$OS.effect %in% c("Favorable", "Unfavorable")]
#48
PFS.biomarker <-  ClinicalRelevance.methy$CpG[ClinicalRelevance.methy$PFS.effect %in% c("Favorable", "Unfavorable")]
#61
Stage.biomarker <- ClinicalRelevance.methy$CpG[ClinicalRelevance.methy$Stage.effect %in% c("Inverse", "Positive")]
#30
OS.gene <- unique(Validated_pairs$Gene_name[Validated_pairs$CpG_pos %in% OS.biomarker])
#42
PFS.gene <- unique(Validated_pairs$Gene_name[Validated_pairs$CpG_pos %in% PFS.biomarker])
#52
Stage.gene <- unique(Validated_pairs$Gene_name[Validated_pairs$CpG_pos %in% Stage.biomarker])
#32

length(union(OS.gene, PFS.gene))#59
length(union(union(OS.gene, PFS.gene),Stage.gene))#75
length(intersect(union(OS.gene, PFS.gene),Stage.gene))#16
length(intersect(intersect(OS.gene, PFS.gene),Stage.gene))#13

Validated_pairs <- cbind(CpG = Validated_pairs$CpG_pos, Gene = Validated_pairs$Gene_name, DMRType = Validated_pairs$DMRType, DM.median = Validated_pairs$DM.median,
                         ClinicalRelevance.methy[match(Validated_pairs$CpG_pos, ClinicalRelevance.methy$CpG),c(4:14)])
#save ClinicalRelevance.methy
ClinicalRelevance.methy <- Validated_pairs
write.csv(ClinicalRelevance.methy, file="/data/huangp/Methy/Mixed/New/ClinicalRelevance_highconfident_methylation_driven_genes_at_methylation_level.csv", row.names=F, quote=F)


#test each CpG/gene pair independently via Cox regression
#LIHC validated promoter CpG-gene pairs
promoterCpGs_Genes <- read.csv("../New/LIHC_Validated_Promoter_Like_DMR_Gene_pairs_0.8.csv", header=T)
library(foreach)
ClinicalRelevance.promoter <- foreach(i = 1:nrow(promoterCpGs_Genes), .combine = 'rbind')%do%{
  CpGID = promoterCpGs_Genes$CpG_pos[i]
  Gene = promoterCpGs_Genes$Gene_name[i]
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$age) & !is.na(df.methy$race) & !is.na(df.methy$gender) & !is.na(df.methy$stage2) & !is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy) & !is.na(df.methy$methy.scaled),]
  OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
  OS.pvalue.methy <- signif(summary(OS.cox)[[7]][8,5],2)
  OS.HR.methy <- round(summary(OS.cox)[[8]][8,1],2)
  PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
  PFS.pvalue.methy <- signif(summary(PFS.cox)[[7]][8,5], 2)
  PFS.HR.methy <- round(summary(PFS.cox)[[8]][8,1], 2)
  methy.survival <- c(OS.HR.methy, OS.pvalue.methy, PFS.HR.methy, PFS.pvalue.methy)
  names(methy.survival) <- c("OS.HR.methy", "OS.pvalue.methy", "PFS.HR.methy", "PFS.pvalue.methy")
  # cox for expr
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$age) & !is.na(df.expr$race) & !is.na(df.expr$gender) & !is.na(df.expr$stage2) & !is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr) & !is.na(df.expr$expr.scaled),]
  OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
  OS.pvalue.expr <- signif(summary(OS.cox)[[7]][8,5],2)
  OS.HR.expr <- round(summary(OS.cox)[[8]][8,1],2)
  PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
  PFS.pvalue.expr <- signif(summary(PFS.cox)[[7]][8,5], 2)
  PFS.HR.expr <- round(summary(PFS.cox)[[8]][8,1], 2)
  expr.survival <- c(OS.HR.expr, OS.pvalue.expr, PFS.HR.expr, PFS.pvalue.expr)
  names(expr.survival) <- c("OS.HR.expr", "OS.pvalue.expr", "PFS.HR.expr", "PFS.pvalue.expr")
  #  ANOVA of methylation for tumor stage
  df.methy <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.methy <- df.methy[df.methy$sampleID %in% names(methy),]#keep the sample with both methylation and phenotypes
  methy.scaled <- scale(methy)[,1]# transform the methylation data to Z score
  df.methy$methy <- round(methy[match(df.methy$sampleID, names(methy))],2)
  df.methy$methy.scaled <- round(methy.scaled[match(df.methy$sampleID, names(methy.scaled))],2)
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  stage.methy.p <- signif(summary(aov(methy ~ stage2, data = df.methy))[[1]][[5]][1], digits = 2)
  stageI.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage I"], na.rm=T),2)
  stageII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage II"], na.rm=T),2)
  stageIII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage III"], na.rm=T),2)
  stageIV.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage IV"], na.rm=T),2)
  stage.methy <- c(stageI.methy, stageII.methy, stageIII.methy, stageIV.methy, stage.methy.p)
  names(stage.methy) <- c("mean.methy.stageI", "mean.methy.stageII", "mean.methy.stageIII", "mean.methy.stageIV", "stage.pvalue.methy")
  #  ANOVA of expression for tumor stage
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  stage.expr.p <- signif(summary(aov(expr ~ stage2, data = df.expr))[[1]][[5]][1], digits = 2)
  stageI.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage I"], na.rm=T),2)
  stageII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage II"], na.rm=T),2)
  stageIII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage III"], na.rm=T),2)
  stageIV.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage IV"], na.rm=T),2)
  stage.expr <- c(stageI.expr, stageII.expr, stageIII.expr, stageIV.expr, stage.expr.p)
  names(stage.expr) <- c("mean.expr.stageI", "mean.expr.stageII", "mean.expr.stageIII", "mean.expr.stageIV", "stage.pvalue.expr")
  CR <- c(methy.survival, expr.survival, stage.methy, stage.expr)
  CR
}
ClinicalRelevance.promoter <- cbind(promoterCpGs_Genes[,c("Gene_name", "CpG_pos", "DMRID","DM.median", "MethyandGene.rho", "LFC.LIHC", "MethyDiff.LIHC")], ClinicalRelevance.promoter)
rm(OS.HR.expr)
rm(OS.HR.methy)
rm(OS.pvalue.expr)
rm(OS.pvalue.methy)
rm(PFS.cox)
rm(PFS.HR.expr)
rm(PFS.HR.methy)
rm(PFS.pvalue.expr)
rm(PFS.pvalue.methy)
attach(ClinicalRelevance.promoter)
table(OS.pvalue.methy < 0.05 & OS.pvalue.expr < 0.05)#6
table(PFS.pvalue.methy < 0.05 & PFS.pvalue.expr < 0.05)#7
table(stage.pvalue.methy < 0.05 & stage.pvalue.expr < 0.05)#6
table((OS.pvalue.expr < 0.05 & PFS.pvalue.expr < 0.05) & (OS.pvalue.methy < 0.05 & PFS.pvalue.methy < 0.05))#3
table((OS.pvalue.expr < 0.05 & PFS.pvalue.expr < 0.05) & (OS.pvalue.methy < 0.05 & PFS.pvalue.methy < 0.05) & (stage.pvalue.methy < 0.05 & stage.pvalue.expr < 0.05))#1
detach(ClinicalRelevance.promoter)
write.csv(ClinicalRelevance.promoter, file="../New/ClinicalRelevance.LIHC_validated_promoterCpG_Genes.csv", row.names = F)
#LIHC validated genic enhancer CpG-gene pairs
GenicEnhancerCpGs_Genes <- read.csv("../New/LIHC_Validated_Enhancer_Like_GenicDMR_Gene_pairs_0.8.csv", header=T)
library(foreach)
ClinicalRelevance.GenicEnhancer <- foreach(i = 1:nrow(GenicEnhancerCpGs_Genes), .combine = 'rbind')%do%{
  CpGID = GenicEnhancerCpGs_Genes$CpG_pos[i]
  Gene = GenicEnhancerCpGs_Genes$Gene_name[i]
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$age) & !is.na(df.methy$race) & !is.na(df.methy$gender) & !is.na(df.methy$stage2) & !is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy) & !is.na(df.methy$methy.scaled),]
  OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
  OS.pvalue.methy <- signif(summary(OS.cox)[[7]][8,5],2)
  OS.HR.methy <- round(summary(OS.cox)[[8]][8,1],2)
  PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
  PFS.pvalue.methy <- signif(summary(PFS.cox)[[7]][8,5], 2)
  PFS.HR.methy <- round(summary(PFS.cox)[[8]][8,1], 2)
  methy.survival <- c(OS.HR.methy, OS.pvalue.methy, PFS.HR.methy, PFS.pvalue.methy)
  names(methy.survival) <- c("OS.HR.methy", "OS.pvalue.methy", "PFS.HR.methy", "PFS.pvalue.methy")
  # cox for expr
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$age) & !is.na(df.expr$race) & !is.na(df.expr$gender) & !is.na(df.expr$stage2) & !is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr) & !is.na(df.expr$expr.scaled),]
  OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
  OS.pvalue.expr <- signif(summary(OS.cox)[[7]][8,5],2)
  OS.HR.expr <- round(summary(OS.cox)[[8]][8,1],2)
  PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
  PFS.pvalue.expr <- signif(summary(PFS.cox)[[7]][8,5], 2)
  PFS.HR.expr <- round(summary(PFS.cox)[[8]][8,1], 2)
  expr.survival <- c(OS.HR.expr, OS.pvalue.expr, PFS.HR.expr, PFS.pvalue.expr)
  names(expr.survival) <- c("OS.HR.expr", "OS.pvalue.expr", "PFS.HR.expr", "PFS.pvalue.expr")
  #  ANOVA of methylation for tumor stage
  df.methy <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.methy <- df.methy[df.methy$sampleID %in% names(methy),]#keep the sample with both methylation and phenotypes
  methy.scaled <- scale(methy)[,1]# transform the methylation data to Z score
  df.methy$methy <- round(methy[match(df.methy$sampleID, names(methy))],2)
  df.methy$methy.scaled <- round(methy.scaled[match(df.methy$sampleID, names(methy.scaled))],2)
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  stage.methy.p <- signif(summary(aov(methy ~ stage2, data = df.methy))[[1]][[5]][1], digits = 2)
  stageI.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage I"], na.rm=T),2)
  stageII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage II"], na.rm=T),2)
  stageIII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage III"], na.rm=T),2)
  stageIV.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage IV"], na.rm=T),2)
  stage.methy <- c(stageI.methy, stageII.methy, stageIII.methy, stageIV.methy, stage.methy.p)
  names(stage.methy) <- c("mean.methy.stageI", "mean.methy.stageII", "mean.methy.stageIII", "mean.methy.stageIV", "stage.pvalue.methy")
  #  ANOVA of expression for tumor stage
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  stage.expr.p <- signif(summary(aov(expr ~ stage2, data = df.expr))[[1]][[5]][1], digits = 2)
  stageI.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage I"], na.rm=T),2)
  stageII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage II"], na.rm=T),2)
  stageIII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage III"], na.rm=T),2)
  stageIV.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage IV"], na.rm=T),2)
  stage.expr <- c(stageI.expr, stageII.expr, stageIII.expr, stageIV.expr, stage.expr.p)
  names(stage.expr) <- c("mean.expr.stageI", "mean.expr.stageII", "mean.expr.stageIII", "mean.expr.stageIV", "stage.pvalue.expr")
  CR <- c(methy.survival, expr.survival, stage.methy, stage.expr)
  CR
}
ClinicalRelevance.GenicEnhancer <- cbind(GenicEnhancerCpGs_Genes[,c("Gene_name","CpG_pos","DMRID", "DM.median", "MethyandGene.rho", "LFC.LIHC", "MethyDiff.LIHC")], ClinicalRelevance.GenicEnhancer)
rm(OS.HR.expr)
rm(OS.HR.methy)
rm(OS.pvalue.expr)
rm(OS.pvalue.methy)
rm(PFS.cox)
rm(PFS.HR.expr)
rm(PFS.HR.methy)
rm(PFS.pvalue.expr)
rm(PFS.pvalue.methy)
attach(ClinicalRelevance.GenicEnhancer)
table(OS.pvalue.methy < 0.05 & OS.pvalue.expr < 0.05)#7
table(PFS.pvalue.methy < 0.05 & PFS.pvalue.expr < 0.05)#7
table(stage.pvalue.methy < 0.05 & stage.pvalue.expr < 0.05)#17
table((OS.pvalue.expr < 0.05 & PFS.pvalue.expr < 0.05) & (OS.pvalue.methy < 0.05 & PFS.pvalue.methy < 0.05))#6
table((OS.pvalue.expr < 0.05 & PFS.pvalue.expr < 0.05) & (OS.pvalue.methy < 0.05 & PFS.pvalue.methy < 0.05) & (stage.pvalue.methy < 0.05 & stage.pvalue.expr < 0.05))#1
detach(ClinicalRelevance.GenicEnhancer)
write.csv(ClinicalRelevance.GenicEnhancer, file="../New/ClinicalRelevance.LIHC_validated_GenicEnhancerCpG_Genes.csv", row.names = F)

#LIHC validated intergenic enhancer CpG-gene pairs
InterGenicEnhancerCpGs_Genes <- read.csv("../New/LIHC_Validated_IntergenicEnhancer_Like_DMR_Gene_pairs_0.7.csv", header=T)
library(foreach)
ClinicalRelevance.InterGenicEnhancer <- foreach(i = 1:nrow(InterGenicEnhancerCpGs_Genes), .combine = 'rbind')%do%{
  CpGID = InterGenicEnhancerCpGs_Genes$CpG_pos[i]
  Gene = InterGenicEnhancerCpGs_Genes$Gene_name[i]
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$age) & !is.na(df.methy$race) & !is.na(df.methy$gender) & !is.na(df.methy$stage2) & !is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy) & !is.na(df.methy$methy.scaled),]
  OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
  OS.pvalue.methy <- signif(summary(OS.cox)[[7]][8,5],2)
  OS.HR.methy <- round(summary(OS.cox)[[8]][8,1],2)
  PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + methy.scaled, data = df.methy)
  PFS.pvalue.methy <- signif(summary(PFS.cox)[[7]][8,5], 2)
  PFS.HR.methy <- round(summary(PFS.cox)[[8]][8,1], 2)
  methy.survival <- c(OS.HR.methy, OS.pvalue.methy, PFS.HR.methy, PFS.pvalue.methy)
  names(methy.survival) <- c("OS.HR.methy", "OS.pvalue.methy", "PFS.HR.methy", "PFS.pvalue.methy")
  # cox for expr
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$age) & !is.na(df.expr$race) & !is.na(df.expr$gender) & !is.na(df.expr$stage2) & !is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr) & !is.na(df.expr$expr.scaled),]
  OS.cox <- coxph(Surv(OS.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
  OS.pvalue.expr <- signif(summary(OS.cox)[[7]][8,5],2)
  OS.HR.expr <- round(summary(OS.cox)[[8]][8,1],2)
  PFS.cox <- coxph(Surv(PFI.time, OS) ~ age + race + gender + stage2 + expr.scaled, data = df.expr)
  PFS.pvalue.expr <- signif(summary(PFS.cox)[[7]][8,5], 2)
  PFS.HR.expr <- round(summary(PFS.cox)[[8]][8,1], 2)
  expr.survival <- c(OS.HR.expr, OS.pvalue.expr, PFS.HR.expr, PFS.pvalue.expr)
  names(expr.survival) <- c("OS.HR.expr", "OS.pvalue.expr", "PFS.HR.expr", "PFS.pvalue.expr")
  #  ANOVA of methylation for tumor stage
  df.methy <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.methy <- df.methy[df.methy$sampleID %in% names(methy),]#keep the sample with both methylation and phenotypes
  methy.scaled <- scale(methy)[,1]# transform the methylation data to Z score
  df.methy$methy <- round(methy[match(df.methy$sampleID, names(methy))],2)
  df.methy$methy.scaled <- round(methy.scaled[match(df.methy$sampleID, names(methy.scaled))],2)
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  stage.methy.p <- signif(summary(aov(methy ~ stage2, data = df.methy))[[1]][[5]][1], digits = 2)
  stageI.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage I"], na.rm=T),2)
  stageII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage II"], na.rm=T),2)
  stageIII.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage III"], na.rm=T),2)
  stageIV.methy <- round(mean(df.methy$methy[df.methy$stage2 == "Stage IV"], na.rm=T),2)
  stage.methy <- c(stageI.methy, stageII.methy, stageIII.methy, stageIV.methy, stage.methy.p)
  names(stage.methy) <- c("mean.methy.stageI", "mean.methy.stageII", "mean.methy.stageIII", "mean.methy.stageIV", "stage.pvalue.methy")
  #  ANOVA of expression for tumor stage
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both expression and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the expression data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  stage.expr.p <- signif(summary(aov(expr ~ stage2, data = df.expr))[[1]][[5]][1], digits = 2)
  stageI.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage I"], na.rm=T),2)
  stageII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage II"], na.rm=T),2)
  stageIII.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage III"], na.rm=T),2)
  stageIV.expr <- round(mean(df.expr$expr[df.expr$stage2 == "Stage IV"], na.rm=T),2)
  stage.expr <- c(stageI.expr, stageII.expr, stageIII.expr, stageIV.expr, stage.expr.p)
  names(stage.expr) <- c("mean.expr.stageI", "mean.expr.stageII", "mean.expr.stageIII", "mean.expr.stageIV", "stage.pvalue.expr")
  CR <- c(methy.survival, expr.survival, stage.methy, stage.expr)
  CR
}
ClinicalRelevance.InterGenicEnhancer <- cbind(InterGenicEnhancerCpGs_Genes[,c("Gene_name","CpG_pos", "DMRID", "DM.median", "MethyandGene.rho", "LFC.LIHC", "MethyDiff.LIHC")], ClinicalRelevance.InterGenicEnhancer)
rm(OS.HR.expr)
rm(OS.HR.methy)
rm(OS.pvalue.expr)
rm(OS.pvalue.methy)
rm(PFS.cox)
rm(PFS.HR.expr)
rm(PFS.HR.methy)
rm(PFS.pvalue.expr)
rm(PFS.pvalue.methy)
attach(ClinicalRelevance.InterGenicEnhancer)
table(OS.pvalue.methy < 0.05 & OS.pvalue.expr < 0.05)#0
table(PFS.pvalue.methy < 0.05 & PFS.pvalue.expr < 0.05)#1
table(stage.pvalue.methy < 0.05 & stage.pvalue.expr < 0.05)#1
table((OS.pvalue.expr < 0.05 & PFS.pvalue.expr < 0.05) & (OS.pvalue.methy < 0.05 & PFS.pvalue.methy < 0.05))#0
table((OS.pvalue.expr < 0.05 & PFS.pvalue.expr < 0.05) & (OS.pvalue.methy < 0.05 & PFS.pvalue.methy < 0.05) & (stage.pvalue.methy < 0.05 & stage.pvalue.expr < 0.05))#0
detach(ClinicalRelevance.InterGenicEnhancer)
write.csv(ClinicalRelevance.InterGenicEnhancer, file="../New/ClinicalRelevance.LIHC_validated_InterGenicEnhancerCpG_Genes.csv", row.names = F)

############################# Survival plot ###################################################
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
df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                       !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                     !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])

###determined the optimal cutoff 
cutoffPick.Expr <- function(DF.both=df.both)
{
  cutoffList <- unique(DF.both$expr)
  cutoffList <- cutoffList[cutoffList >= quantile(DF.both$expr, 0.1) & cutoffList <= quantile(DF.both$expr, 0.9)]
  cutoffList <- cutoffList[order(cutoffList, decreasing=F)]
  library(foreach)
  ChisqList <- foreach(i = 1:length(cutoffList),.combine = 'rbind')%do%
    {
      DF.both$Group <- DF.both$expr >= cutoffList[i]
      DF.both$Group <- factor(DF.both$Group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.both)[[5]]
      Chisq.PFI <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.both)[[5]]
      c(cutoffList[i], Chisq.OS, Chisq.PFI)
    }
  ChisqList <- as.data.frame(ChisqList)
  colnames(ChisqList) <- c("Cutoff", "Chisq.OS", "Chisq.PFS")
  Cutoff.OS <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Cutoff"][1]
  MaxChisq.OS <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Chisq.OS"][1]
  Cutoff.PFS <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T), "Cutoff"][1]
  MaxChisq.PFS <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T),"Chisq.PFS"][1]
  return(c(Cutoff.OS, MaxChisq.OS, Cutoff.PFS, MaxChisq.PFS))
}
cutoffPick.Methy <- function(DF.both = df.both)
{
  cutoffList <- unique(DF.both$methy)
  cutoffList <- cutoffList[cutoffList >= quantile(DF.both$methy, 0.1) & cutoffList <= quantile(DF.both$methy, 0.9)]
  cutoffList <- cutoffList[order(cutoffList, decreasing=F)]
  library(foreach)
  ChisqList <- foreach(i = 1:length(cutoffList),.combine = 'rbind')%do%
    {
      DF.both$Group <- DF.both$methy >= cutoffList[i]
      DF.both$Group <- factor(DF.both$Group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.both)[[5]]
      Chisq.PFI <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.both)[[5]]
      c(cutoffList[i], Chisq.OS, Chisq.PFI)
    }
  ChisqList <- as.data.frame(ChisqList)
  colnames(ChisqList) <- c("Cutoff", "Chisq.OS", "Chisq.PFS")
  Cutoff.OS <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Cutoff"][1]
  MaxChisq.OS <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Chisq.OS"][1]
  Cutoff.PFS <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T), "Cutoff"][1]
  MaxChisq.PFS <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T),"Chisq.PFS"][1]
  return(c(Cutoff.OS, MaxChisq.OS, Cutoff.PFS, MaxChisq.PFS))
}

cutoffPick.Combined <- function(DF.both = df.both)
{
  #methy top 30 cutoff
  cutoffList <- unique(DF.both$methy)
  cutoffList <- cutoffList[cutoffList >= quantile(DF.both$methy, 0.1) & cutoffList <= quantile(DF.both$methy, 0.9)]
  cutoffList <- cutoffList[order(cutoffList, decreasing=F)]
  library(foreach)
  ChisqList <- foreach(i = 1:length(cutoffList),.combine = 'rbind')%do%
    {
      DF.both$Group <- DF.both$methy >= cutoffList[i]
      DF.both$Group <- factor(DF.both$Group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.both)[[5]]
      Chisq.PFI <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.both)[[5]]
      c(cutoffList[i], Chisq.OS, Chisq.PFI)
    }
  ChisqList <- as.data.frame(ChisqList)
  colnames(ChisqList) <- c("Cutoff", "Chisq.OS", "Chisq.PFS")
  Cutoff.OS.methy30 <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Cutoff"][1:30]
  Cutoff.PFS.methy30 <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T), "Cutoff"][1:30]
  #expr top 30 cutoff
  cutoffList <- unique(DF.both$expr)
  cutoffList <- cutoffList[cutoffList >= quantile(DF.both$expr, 0.1) & cutoffList <= quantile(DF.both$expr, 0.9)]
  cutoffList <- cutoffList[order(cutoffList, decreasing=F)]
  library(foreach)
  ChisqList <- foreach(i = 1:length(cutoffList),.combine = 'rbind')%do%
    {
      DF.both$Group <- DF.both$expr >= cutoffList[i]
      DF.both$Group <- factor(DF.both$Group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.both)[[5]]
      Chisq.PFI <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.both)[[5]]
      c(cutoffList[i], Chisq.OS, Chisq.PFI)
    }
  ChisqList <- as.data.frame(ChisqList)
  colnames(ChisqList) <- c("Cutoff", "Chisq.OS", "Chisq.PFS")
  Cutoff.OS.expr30 <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Cutoff"][1:30]
  Cutoff.PFS.expr30 <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T), "Cutoff"][1:30]
  # combined optimal cutoff
  #OS
  cutoffList.OS <- data.frame(ExprCutoff = rep(Cutoff.OS.expr30, each =30), MethyCutoff = rep(Cutoff.OS.methy30, 30))
  library(foreach)
  ChisqList.OS <- foreach(i = 1:nrow(cutoffList.OS),.combine = 'rbind')%do%
    {
      ExprGroup <- DF.both$expr >= cutoffList.OS[i,1] 
      ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
      MethyGroup <- DF.both$methy >= cutoffList.OS[i,2]
      MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
      DF.both$Group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
      DF.combined <- DF.both[DF.both$Group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
      DF.combined$Group <- factor(DF.combined$Group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.combined)[[5]]
      c(unlist(cutoffList.OS[i,]), Chisq.OS)
    }
  ChisqList.OS <- as.data.frame(ChisqList.OS)
  colnames(ChisqList.OS) <- c("ExprCutoff", "MethyCutoff", "Chisq")
  ChisqList.OS <- ChisqList.OS[order(ChisqList.OS$Chisq, decreasing = T),]
  Cutoff.OS.combined <- unlist(ChisqList.OS[1,1:2])
  MaxChisq.OS.combined <- ChisqList.OS[1,3]
  #PFS
  cutoffList.PFS <- data.frame(ExprCutoff = rep(Cutoff.PFS.expr30, each =30), MethyCutoff = rep(Cutoff.PFS.methy30, 30))
  library(foreach)
  ChisqList.PFS <- foreach(i = 1:nrow(cutoffList.PFS),.combine = 'rbind')%do%
    {
      ExprGroup <- DF.both$expr >= cutoffList.PFS[i,1] 
      ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
      MethyGroup <- DF.both$methy >= cutoffList.PFS[i,2]
      MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
      DF.both$Group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
      DF.combined <- DF.both[DF.both$Group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
      DF.combined$Group <- factor(DF.combined$Group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
      Chisq.PFS <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.combined)[[5]]
      c(unlist(cutoffList.PFS[i,]), Chisq.PFS)
    }
  ChisqList.PFS <- as.data.frame(ChisqList.PFS)
  colnames(ChisqList.PFS) <- c("ExprCutoff", "MethyCutoff", "Chisq")
  ChisqList.PFS <- ChisqList.PFS[order(ChisqList.PFS$Chisq, decreasing = T),]
  Cutoff.PFS.combined <- unlist(ChisqList.PFS[1,1:2])
  MaxChisq.PFS.combined <- ChisqList.PFS[1,3]
  
  return(c(Cutoff.OS.combined, MaxChisq.OS.combined, Cutoff.PFS.combined, MaxChisq.PFS.combined))
}

cutoffPick.Combined.Pos <- function(DF.both = df.both)
{
  #methy top 30 cutoff
  cutoffList <- unique(DF.both$methy)
  cutoffList <- cutoffList[cutoffList >= quantile(DF.both$methy, 0.1) & cutoffList <= quantile(DF.both$methy, 0.9)]
  cutoffList <- cutoffList[order(cutoffList, decreasing=F)]
  library(foreach)
  ChisqList <- foreach(i = 1:length(cutoffList),.combine = 'rbind')%do%
    {
      DF.both$Group <- DF.both$methy >= cutoffList[i]
      DF.both$Group <- factor(DF.both$Group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.both)[[5]]
      Chisq.PFI <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.both)[[5]]
      c(cutoffList[i], Chisq.OS, Chisq.PFI)
    }
  ChisqList <- as.data.frame(ChisqList)
  colnames(ChisqList) <- c("Cutoff", "Chisq.OS", "Chisq.PFS")
  Cutoff.OS.methy30 <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Cutoff"][1:30]
  Cutoff.PFS.methy30 <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T), "Cutoff"][1:30]
  #expr top 30 cutoff
  cutoffList <- unique(DF.both$expr)
  cutoffList <- cutoffList[cutoffList >= quantile(DF.both$expr, 0.1) & cutoffList <= quantile(DF.both$expr, 0.9)]
  cutoffList <- cutoffList[order(cutoffList, decreasing=F)]
  library(foreach)
  ChisqList <- foreach(i = 1:length(cutoffList),.combine = 'rbind')%do%
    {
      DF.both$Group <- DF.both$expr >= cutoffList[i]
      DF.both$Group <- factor(DF.both$Group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.both)[[5]]
      Chisq.PFI <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.both)[[5]]
      c(cutoffList[i], Chisq.OS, Chisq.PFI)
    }
  ChisqList <- as.data.frame(ChisqList)
  colnames(ChisqList) <- c("Cutoff", "Chisq.OS", "Chisq.PFS")
  Cutoff.OS.expr30 <- ChisqList[order(ChisqList$Chisq.OS, decreasing = T),"Cutoff"][1:30]
  Cutoff.PFS.expr30 <- ChisqList[order(ChisqList$Chisq.PFS, decreasing = T), "Cutoff"][1:30]
  # combined optimal cutoff
  #OS
  cutoffList.OS <- data.frame(ExprCutoff = rep(Cutoff.OS.expr30, each =30), MethyCutoff = rep(Cutoff.OS.methy30, 30))
  library(foreach)
  ChisqList.OS <- foreach(i = 1:nrow(cutoffList.OS),.combine = 'rbind')%do%
    {
      ExprGroup <- DF.both$expr >= cutoffList.OS[i,1] 
      ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
      MethyGroup <- DF.both$methy >= cutoffList.OS[i,2]
      MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
      DF.both$Group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
      DF.combined <- DF.both[DF.both$Group %in% c("HyperMethy & HighExpr", "HypoMethy & LowExpr"),]
      DF.combined$Group <- factor(DF.combined$Group, levels = c("HyperMethy & HighExpr", "HypoMethy & LowExpr"))
      Chisq.OS <- survdiff(Surv(OS.time, OS)~Group, data=DF.combined)[[5]]
      c(unlist(cutoffList.OS[i,]), Chisq.OS)
    }
  ChisqList.OS <- as.data.frame(ChisqList.OS)
  colnames(ChisqList.OS) <- c("ExprCutoff", "MethyCutoff", "Chisq")
  ChisqList.OS <- ChisqList.OS[order(ChisqList.OS$Chisq, decreasing = T),]
  Cutoff.OS.combined <- unlist(ChisqList.OS[1,1:2])
  MaxChisq.OS.combined <- ChisqList.OS[1,3]
  #PFS
  cutoffList.PFS <- data.frame(ExprCutoff = rep(Cutoff.PFS.expr30, each =30), MethyCutoff = rep(Cutoff.PFS.methy30, 30))
  library(foreach)
  ChisqList.PFS <- foreach(i = 1:nrow(cutoffList.PFS),.combine = 'rbind')%do%
    {
      ExprGroup <- DF.both$expr >= cutoffList.PFS[i,1] 
      ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
      MethyGroup <- DF.both$methy >= cutoffList.PFS[i,2]
      MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
      DF.both$Group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
      DF.combined <- DF.both[DF.both$Group %in% c("HyperMethy & HighExpr", "HypoMethy & LowExpr"),]
      DF.combined$Group <- factor(DF.combined$Group, levels = c("HyperMethy & HighExpr", "HypoMethy & LowExpr"))
      Chisq.PFS <- survdiff(Surv(PFI.time, PFI)~Group, data=DF.combined)[[5]]
      c(unlist(cutoffList.PFS[i,]), Chisq.PFS)
    }
  ChisqList.PFS <- as.data.frame(ChisqList.PFS)
  colnames(ChisqList.PFS) <- c("ExprCutoff", "MethyCutoff", "Chisq")
  ChisqList.PFS <- ChisqList.PFS[order(ChisqList.PFS$Chisq, decreasing = T),]
  Cutoff.PFS.combined <- unlist(ChisqList.PFS[1,1:2])
  MaxChisq.PFS.combined <- ChisqList.PFS[1,3]
  
  return(c(Cutoff.OS.combined, MaxChisq.OS.combined, Cutoff.PFS.combined, MaxChisq.PFS.combined))
}
cutoffPick.Methy(df.both)
# 0.440000 15.506634  0.440000  9.959812
cutoffPick.Expr(df.both)
#8.24000 28.49248  8.21000 16.21972
cutoffPick.Combined(df.both)
#8.24000     0.40000    40.71206     7.98000     0.44000    19.13941
####OS
df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[1]
df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[1]
df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))

ExprGroup <- DF.both$expr >= cutoffPick.Combined(df.both)[1] 
ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
MethyGroup <- DF.both$methy >= cutoffPick.Combined(df.both)[2]
MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))

survdiff(Surv(OS.time, OS)~ expr.group, data=df.both)
survdiff(Surv(OS.time, OS)~ methy.group, data=df.both)
survdiff(Surv(OS.time, OS)~group, data=df.combined)

library(survminer)
library(ggpubr)

OS.combined <- ggsurvplot(survfit(Surv(OS.time, OS) ~ group,
                                  data = df.combined),
                          data = df.combined,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("CDC20 OS.combined.pdf", width=11.69/2, height = 8.268/2)
OS.combined
dev.off()

OS.expr <- ggsurvplot(survfit(Surv(OS.time, OS) ~ expr.group,
                                  data = df.both),
                          data = df.both,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("CDC20 OS.expr.pdf", width=11.69/2, height = 8.268/2)
OS.expr
dev.off()

OS.methy <- ggsurvplot(survfit(Surv(OS.time, OS) ~ methy.group,
                              data = df.both),
                      data = df.both,
                      risk.table = F,
                      pval = TRUE,
                      break.time.by = 500,
                      ggtheme = theme_pubr(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("CDC20 OS.methy.pdf", width=11.69/2, height = 8.268/2)
OS.methy
dev.off()

####PFS
df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))

ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[4] 
ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[5]
MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))

survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
survdiff(Surv(PFI.time, PFI)~group, data=df.combined)

library(survminer)
library(ggpubr)

PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                  data = df.combined),
                          data = df.combined,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
  ylab("Progression-free survival probability")
pdf("CDC20 PFS.combined.pdf", width=11.69/2, height = 8.268/2)
PFI.combined
dev.off()

PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                              data = df.both),
                      data = df.both,
                      risk.table = F,
                      pval = TRUE,
                      break.time.by = 500,
                      ggtheme = theme_pubr(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)+
  ylab("Progression-free survival probability")
pdf("CDC20 PFS.expr.pdf", width=11.69/2, height = 8.268/2)
PFI.expr
dev.off()

PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                               data = df.both),
                       data = df.both,
                       risk.table = F,
                       pval = TRUE,
                       break.time.by = 500,
                       ggtheme = theme_pubr(),
                       risk.table.y.text.col = TRUE,
                       risk.table.y.text = FALSE)+
  ylab("Progression-free survival probability")
pdf("CDC20 PFS.methy.pdf", width=11.69/2, height = 8.268/2)
PFI.methy
dev.off()

CpGID = "chr1_43360233"
Gene = "CDC20"
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
df.methy <- df.methy[!is.na(df.methy$stage2) & !is.na(df.methy$methy) & df.methy$stage2 != "Stage IV",]
df.methy$group <- rep("Adjacent", nrow(df.methy))
df.methy$group[df.methy$condition == "Tumor" & df.methy$stage2 == "Stage I"] <- "Stage I"
df.methy$group[df.methy$condition == "Tumor" & df.methy$stage2 == "Stage II"] <- "Stage II"
df.methy$group[df.methy$condition == "Tumor" & df.methy$stage2 == "Stage III"] <- "Stage III"


df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
df.expr <- df.expr[!is.na(df.expr$stage2) & !is.na(df.expr$expr) & df.expr$stage2 != "Stage IV",]
df.expr$group <- rep("Adjacent", nrow(df.expr))
df.expr$group[df.expr$condition == "Tumor" & df.expr$stage2 == "Stage I"] <- "Stage I"
df.expr$group[df.expr$condition == "Tumor" & df.expr$stage2 == "Stage II"] <- "Stage II"
df.expr$group[df.expr$condition == "Tumor" & df.expr$stage2 == "Stage III"] <- "Stage III"

# boxplot of methy in different tumor stage
library(ggplot2)
library(ggsci)
library(ggpubr)
stage.methy <- ggplot(data=df.methy, aes(x = group, y = methy))+
  geom_jitter(aes(fill = group), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_lancet(guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("")+
  ylab("Methylation of CDC20")+
  theme_pubr()+ggtitle("p = 0.0096")
postscript("CDC20 methy tumor stage.eps", width=11.69/3, height = 8.268/3)
stage.methy
dev.off()

# boxplot of expr in different tumor stage
stage.expr <- ggplot(data=df.expr, aes(x = group, y = expr))+
  geom_jitter(aes(fill = group), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_lancet(guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("")+
  ylab("mRNA expression of CDC20")+
  theme_pubr()+ggtitle("p = 0.00013")
postscript("CDC20 expr tumor stage.eps", width=11.69/3, height = 8.268/3)
stage.expr  
dev.off()

# UCK2 chr1_165857033
# OS
CpGID = "chr1_165857033"
Gene = "UCK2"
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
df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                       !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                     !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])

cutoffPick.Methy(df.both)
# 0.40000 26.10802  0.34000 18.28246
cutoffPick.Expr(df.both)
# 9.52000 54.75569  9.43000 14.61707
cutoffPick.Combined(df.both)
# 9.52000     0.70000    54.13664     8.46000     0.34000    23.99151
####OS
df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[1]
df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[1]
df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))

ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[1] 
ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[2]
MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))

survdiff(Surv(OS.time, OS)~ expr.group, data=df.both)
survdiff(Surv(OS.time, OS)~ methy.group, data=df.both)
survdiff(Surv(OS.time, OS)~group, data=df.combined)

library(survminer)
library(ggpubr)

OS.combined <- ggsurvplot(survfit(Surv(OS.time, OS) ~ group,
                                  data = df.combined),
                          data = df.combined,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("UCK2 OS.combined.pdf", width=11.69/2, height = 8.268/2)
OS.combined
dev.off()

OS.expr <- ggsurvplot(survfit(Surv(OS.time, OS) ~ expr.group,
                              data = df.both),
                      data = df.both,
                      risk.table = F,
                      pval = TRUE,
                      break.time.by = 500,
                      ggtheme = theme_pubr(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("UCK2 OS.expr.pdf", width=11.69/2, height = 8.268/2)
OS.expr
dev.off()

OS.methy <- ggsurvplot(survfit(Surv(OS.time, OS) ~ methy.group,
                               data = df.both),
                       data = df.both,
                       risk.table = F,
                       pval = TRUE,
                       break.time.by = 500,
                       ggtheme = theme_pubr(),
                       risk.table.y.text.col = TRUE,
                       risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("UCK2 OS.methy.pdf", width=11.69/2, height = 8.268/2)
OS.methy
dev.off()

####PFS
####PFS
df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))

ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[4] 
ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[5]
MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))

survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
survdiff(Surv(PFI.time, PFI)~group, data=df.combined)

library(survminer)
library(ggpubr)

PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                   data = df.combined),
                           data = df.combined,
                           risk.table = F,
                           pval = TRUE,
                           break.time.by = 500,
                           ggtheme = theme_pubr(),
                           risk.table.y.text.col = TRUE,
                           risk.table.y.text = FALSE)+
  ylab("Progression-free survival probability")
pdf("UCK2 PFS.combined.pdf", width=11.69/2, height = 8.268/2)
PFI.combined
dev.off()

PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                               data = df.both),
                       data = df.both,
                       risk.table = F,
                       pval = TRUE,
                       break.time.by = 500,
                       ggtheme = theme_pubr(),
                       risk.table.y.text.col = TRUE,
                       risk.table.y.text = FALSE)+
  ylab("Progression-free survival probability")
pdf("UCK2 PFS.expr.pdf", width=11.69/2, height = 8.268/2)
PFI.expr
dev.off()

PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                                data = df.both),
                        data = df.both,
                        risk.table = F,
                        pval = TRUE,
                        break.time.by = 500,
                        ggtheme = theme_pubr(),
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE)+
  ylab("Progression-free survival probability")
pdf("UCK2 PFS.methy.pdf", width=11.69/2, height = 8.268/2)
PFI.methy
dev.off()

# boxplot of methy in different tumor stage
stage.methy <- ggplot(data=df.methy, aes(x = stage2, y = methy))+
  geom_jitter(aes(fill = stage2), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_brewer(palette = "Reds", guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("Tumor stage")+
  ylab("Methylation of UCK2")+
  theme_pubr()+ggtitle("p = 0.033")
stage.methy
# boxplot of expr in different tumor stage
stage.expr <- ggplot(data=df.expr, aes(x = stage2, y = expr))+
  geom_jitter(aes(fill = stage2), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_brewer(palette = "Reds", guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("Tumor stage")+
  ylab("mRNA expression of UCK2")+
  theme_pubr()+ggtitle("p = 0.013")
stage.expr 


# OS SLC22A15 chr1:115978918
if(1<0){
CpGID = "chr1_115978918"
Gene = "SLC22A15"
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
df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                       !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                     !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])

cutoffPick.Methy(df.both)
#  0.47000 12.16883  0.43000  3.74578
cutoffPick.Expr(df.both)
# 6.930000 14.315366  6.890000  8.258751
cutoffPick.Combined(df.both)
# 6.93000     0.66000    18.86492     6.35000     0.36000    14.28142
####OS
df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[1]
df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[1]
df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))

ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[1] 
ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[2]
MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))

survdiff(Surv(OS.time, OS)~ expr.group, data=df.both)
survdiff(Surv(OS.time, OS)~ methy.group, data=df.both)
survdiff(Surv(OS.time, OS)~group, data=df.combined)

library(survminer)
library(ggpubr)

OS.combined <- ggsurvplot(survfit(Surv(OS.time, OS) ~ group,
                                  data = df.combined),
                          data = df.combined,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("SLC22A15 OS.combined.pdf", width=11.69/2, height = 8.268/2)
OS.combined
dev.off()

OS.expr <- ggsurvplot(survfit(Surv(OS.time, OS) ~ expr.group,
                              data = df.both),
                      data = df.both,
                      risk.table = F,
                      pval = TRUE,
                      break.time.by = 500,
                      ggtheme = theme_pubr(),
                      risk.table.y.text.col = TRUE,
                      risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("SLC22A15 OS.expr.pdf", width=11.69/2, height = 8.268/2)
OS.expr
dev.off()

OS.methy <- ggsurvplot(survfit(Surv(OS.time, OS) ~ methy.group,
                               data = df.both),
                       data = df.both,
                       risk.table = F,
                       pval = TRUE,
                       break.time.by = 500,
                       ggtheme = theme_pubr(),
                       risk.table.y.text.col = TRUE,
                       risk.table.y.text = FALSE)+
  ylab("Overall survival probability")
pdf("SLC22A15 OS.methy.pdf", width=11.69/2, height = 8.268/2)
OS.methy
dev.off()
}
# OS CACYBP chr1:175073161
if(1<0){
  CpGID = "chr1_175073161"
  Gene = "CACYBP"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.320000  5.922966  0.170000 14.821506
  cutoffPick.Expr(df.both)
  # 11.140000 18.968206 11.140000  5.597727
  cutoffPick.Combined(df.both)
  # 11.39000     0.18000    30.65307    11.37000     0.17000    30.01720
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[1]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[1]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[1] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[2]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
  
  survdiff(Surv(OS.time, OS)~ expr.group, data=df.both)
  survdiff(Surv(OS.time, OS)~ methy.group, data=df.both)
  survdiff(Surv(OS.time, OS)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  OS.combined <- ggsurvplot(survfit(Surv(OS.time, OS) ~ group,
                                    data = df.combined),
                            data = df.combined,
                            risk.table = F,
                            pval = TRUE,
                            break.time.by = 500,
                            ggtheme = theme_pubr(),
                            risk.table.y.text.col = TRUE,
                            risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("CACYBP OS.combined.pdf", width=11.69/2, height = 8.268/2)
  OS.combined
  dev.off()
  
  OS.expr <- ggsurvplot(survfit(Surv(OS.time, OS) ~ expr.group,
                                data = df.both),
                        data = df.both,
                        risk.table = F,
                        pval = TRUE,
                        break.time.by = 500,
                        ggtheme = theme_pubr(),
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("CACYBP OS.expr.pdf", width=11.69/2, height = 8.268/2)
  OS.expr
  dev.off()
  
  OS.methy <- ggsurvplot(survfit(Surv(OS.time, OS) ~ methy.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("CACYBP OS.methy.pdf", width=11.69/2, height = 8.268/2)
  OS.methy
  dev.off()
}
# OS MTA3 chr2:42492645
if(1<0){
  CpGID = "chr2_42492645"
  Gene = "MTA3"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.680000  5.271556  0.920000 13.320991
  cutoffPick.Expr(df.both)
  # 9.310000 7.147405 9.310000 2.887169
  cutoffPick.Combined.Pos(df.both)
  # 9.25000     0.92000    11.15357     9.10000     0.92000    19.77474
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[1]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[1]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined.Pos(df.both)[1] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined.Pos(df.both)[2]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & HighExpr", "HypoMethy & LowExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & HighExpr", "HypoMethy & LowExpr"))
  
  survdiff(Surv(OS.time, OS)~ expr.group, data=df.both)
  survdiff(Surv(OS.time, OS)~ methy.group, data=df.both)
  survdiff(Surv(OS.time, OS)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  OS.combined <- ggsurvplot(survfit(Surv(OS.time, OS) ~ group,
                                    data = df.combined),
                            data = df.combined,
                            risk.table = F,
                            pval = TRUE,
                            break.time.by = 500,
                            ggtheme = theme_pubr(),
                            risk.table.y.text.col = TRUE,
                            risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("MTA3 OS.combined.pdf", width=11.69/2, height = 8.268/2)
  OS.combined
  dev.off()
  
  OS.expr <- ggsurvplot(survfit(Surv(OS.time, OS) ~ expr.group,
                                data = df.both),
                        data = df.both,
                        risk.table = F,
                        pval = TRUE,
                        break.time.by = 500,
                        ggtheme = theme_pubr(),
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("MTA3 OS.expr.pdf", width=11.69/2, height = 8.268/2)
  OS.expr
  dev.off()
  
  OS.methy <- ggsurvplot(survfit(Surv(OS.time, OS) ~ methy.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("MTA3 OS.methy.pdf", width=11.69/2, height = 8.268/2)
  OS.methy
  dev.off()
}
# OS SLC9A3R2 chr16:2090869
if(1<0){
  CpGID = "chr16_2090869"
  Gene = "SLC9A3R2"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.48000 10.54369  0.86000  1.63693
  cutoffPick.Expr(df.both)
  # 10.71000 14.52505 11.86000 15.55355
  cutoffPick.Combined(df.both)
  # 10.82000     0.66000    30.37569    11.86000     0.82000    12.10629
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[1]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[1]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[1] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[2]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
  
  survdiff(Surv(OS.time, OS)~ expr.group, data=df.both)
  survdiff(Surv(OS.time, OS)~ methy.group, data=df.both)
  survdiff(Surv(OS.time, OS)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  OS.combined <- ggsurvplot(survfit(Surv(OS.time, OS) ~ group,
                                    data = df.combined),
                            data = df.combined,
                            risk.table = F,
                            pval = TRUE,
                            break.time.by = 500,
                            ggtheme = theme_pubr(),
                            risk.table.y.text.col = TRUE,
                            risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("SLC9A3R2 OS.combined.pdf", width=11.69/2, height = 8.268/2)
  OS.combined
  dev.off()
  
  OS.expr <- ggsurvplot(survfit(Surv(OS.time, OS) ~ expr.group,
                                data = df.both),
                        data = df.both,
                        risk.table = F,
                        pval = TRUE,
                        break.time.by = 500,
                        ggtheme = theme_pubr(),
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("SLC9A3R2 OS.expr.pdf", width=11.69/2, height = 8.268/2)
  OS.expr
  dev.off()
  
  OS.methy <- ggsurvplot(survfit(Surv(OS.time, OS) ~ methy.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Overall survival probability")
  pdf("SLC9A3R2 OS.methy.pdf", width=11.69/2, height = 8.268/2)
  OS.methy
  dev.off()
}
#PFS LRRC4 chr7:128032099
if(1<0){
  CpGID = "chr7_128032099"
  Gene = "LRRC4"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.460000  6.433256  0.480000 20.303422
  cutoffPick.Expr(df.both)
  # 2.370000  8.168262  2.940000 10.541918
  cutoffPick.Combined(df.both)
  # 2.37000     0.40000    11.80686     3.38000     0.48000    21.98796
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[4] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[5]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
  
  survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                    data = df.combined),
                            data = df.combined,
                            risk.table = F,
                            pval = TRUE,
                            break.time.by = 500,
                            ggtheme = theme_pubr(),
                            risk.table.y.text.col = TRUE,
                            risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("LRRC4 PFS.combined.pdf", width=11.69/2, height = 8.268/2)
  PFI.combined
  dev.off()
  
  PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                                data = df.both),
                        data = df.both,
                        risk.table = F,
                        pval = TRUE,
                        break.time.by = 500,
                        ggtheme = theme_pubr(),
                        risk.table.y.text.col = TRUE,
                        risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("LRRC4 PFS.expr.pdf", width=11.69/2, height = 8.268/2)
  PFI.expr
  dev.off()
  
  PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("LRRC4 PFS.methy.pdf", width=11.69/2, height = 8.268/2)
  PFI.methy
  dev.off()
}
#PFS MTA3 chr2:42492645
if(1<0){
  CpGID = "chr2_42492645"
  Gene = "MTA3"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.680000  5.271556  0.920000 13.320991
  cutoffPick.Expr(df.both)
  #  9.310000 7.147405 9.310000 2.887169
  cutoffPick.Combined.Pos(df.both)
  # 9.25000     0.92000    11.15357     9.10000     0.92000    19.77474
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined.Pos(df.both)[4] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined.Pos(df.both)[5]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & HighExpr", "HypoMethy & LowExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & HighExpr", "HypoMethy & LowExpr"))
  
  survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                     data = df.combined),
                             data = df.combined,
                             risk.table = F,
                             pval = TRUE,
                             break.time.by = 500,
                             ggtheme = theme_pubr(),
                             risk.table.y.text.col = TRUE,
                             risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("MTA3 PFS.combined.pdf", width=11.69/2, height = 8.268/2)
  PFI.combined
  dev.off()
  
  PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("MTA3 PFS.expr.pdf", width=11.69/2, height = 8.268/2)
  PFI.expr
  dev.off()
  
  PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                                  data = df.both),
                          data = df.both,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("MTA3 PFS.methy.pdf", width=11.69/2, height = 8.268/2)
  PFI.methy
  dev.off()
}
#PFS HEATR6 chr17:60079852
if(1<0){
  CpGID = "chr17_60079852"
  Gene = "HEATR6"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.400000  8.594469  0.400000 14.549966
  cutoffPick.Expr(df.both)
  # 8.60000 19.22969  8.43000 11.12039
  cutoffPick.Combined(df.both)
  # 8.30000     0.40000    17.78959     8.19000     0.40000    26.59083
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[4] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[5]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
  
  survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                     data = df.combined),
                             data = df.combined,
                             risk.table = F,
                             pval = TRUE,
                             break.time.by = 500,
                             ggtheme = theme_pubr(),
                             risk.table.y.text.col = TRUE,
                             risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("HEATR6 PFS.combined.pdf", width=11.69/2, height = 8.268/2)
  PFI.combined
  dev.off()
  
  PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("HEATR6 PFS.expr.pdf", width=11.69/2, height = 8.268/2)
  PFI.expr
  dev.off()
  
  PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                                  data = df.both),
                          data = df.both,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("HEATR6 PFS.methy.pdf", width=11.69/2, height = 8.268/2)
  PFI.methy
  dev.off()
}
#PFS CACYBP chr1:175073161
if(1<0){
  CpGID = "chr1_175073161"
  Gene = "CACYBP"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #   0.320000  5.922966  0.170000 14.821506
  cutoffPick.Expr(df.both)
  # 11.140000 18.968206 11.140000  5.597727
  cutoffPick.Combined(df.both)
  # 11.39000     0.18000    30.65307    11.37000     0.17000    30.01720
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[4] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[5]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
  
  survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                     data = df.combined),
                             data = df.combined,
                             risk.table = F,
                             pval = TRUE,
                             break.time.by = 500,
                             ggtheme = theme_pubr(),
                             risk.table.y.text.col = TRUE,
                             risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("CACYBP PFS.combined.pdf", width=11.69/2, height = 8.268/2)
  PFI.combined
  dev.off()
  
  PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("CACYBP PFS.expr.pdf", width=11.69/2, height = 8.268/2)
  PFI.expr
  dev.off()
  
  PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                                  data = df.both),
                          data = df.both,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("CACYBP PFS.methy.pdf", width=11.69/2, height = 8.268/2)
  PFI.methy
  dev.off()
}
#PFS SLC9A3R2 chr16:2090869
if(1<0){
  CpGID = "chr16_2090869"
  Gene = "SLC9A3R2"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.48000 10.54369  0.86000  1.63693
  cutoffPick.Expr(df.both)
  # 10.71000 14.52505 11.86000 15.55355
  cutoffPick.Combined(df.both)
  # 10.82000     0.66000    30.37569    11.86000     0.82000    12.10629
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined(df.both)[4] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined(df.both)[5]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & LowExpr", "HypoMethy & HighExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & LowExpr", "HypoMethy & HighExpr"))
  
  survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                     data = df.combined),
                             data = df.combined,
                             risk.table = F,
                             pval = TRUE,
                             break.time.by = 500,
                             ggtheme = theme_pubr(),
                             risk.table.y.text.col = TRUE,
                             risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("SLC9A3R2 PFS.combined.pdf", width=11.69/2, height = 8.268/2)
  PFI.combined
  dev.off()
  
  PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("SLC9A3R2 PFS.expr.pdf", width=11.69/2, height = 8.268/2)
  PFI.expr
  dev.off()
  
  PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                                  data = df.both),
                          data = df.both,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("SLC9A3R2 PFS.methy.pdf", width=11.69/2, height = 8.268/2)
  PFI.methy
  dev.off()
}
#PFS ITK chr5:157058637
if(1<0){
  CpGID = "chr5_157058637"
  Gene = "ITK"
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
  df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
  df.methy <- df.methy[!is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                         !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy),]
  df.expr <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
  df.expr <- df.expr[df.expr$sampleID %in% names(expr),]#keep the sample with both exprlation and phenotypes
  expr.scaled <- scale(expr)[,1]# transform the exprlation data to Z score
  df.expr$expr <- round(expr[match(df.expr$sampleID, names(expr))],2)
  df.expr$expr.scaled <- round(expr.scaled[match(df.expr$sampleID, names(expr.scaled))],2)
  df.expr <- df.expr[df.expr$condition == "Tumor",]#keep only tumor samples
  df.expr <- df.expr[!is.na(df.expr$OS) & !is.na(df.expr$OS.time) &
                       !is.na(df.expr$PFI) & !is.na(df.expr$PFI.time) & !is.na(df.expr$expr),]
  samples.both <- intersect(df.methy$sampleID, df.expr$sampleID)
  df.both <- cbind(df.methy[match(samples.both, df.methy$sampleID),], df.expr[match(samples.both, df.expr$sampleID),c("expr", "expr.scaled")])
  
  cutoffPick.Methy(df.both)
  #  0.360000 16.945211  0.680000  2.027449
  cutoffPick.Expr(df.both)
  #  5.35000 11.12727  4.09000 13.88016
  cutoffPick.Combined.Pos(df.both)
  # 4.41000     0.37000    36.60506     4.09000     0.75000    14.80800
  ####OS
  df.both$methy.group <- df.both$methy >= cutoffPick.Methy(df.both)[3]
  df.both$methy.group <- factor(df.both$methy.group, levels = c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$expr.group <- df.both$expr >= cutoffPick.Expr(df.both)[3]
  df.both$expr.group <- factor(df.both$expr.group, levels = c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  
  ExprGroup <- df.both$expr >= cutoffPick.Combined.Pos(df.both)[4] 
  ExprGroup <- factor(ExprGroup, levels=c(TRUE, FALSE), labels = c("HighExpr", "LowExpr"))
  MethyGroup <- df.both$methy >= cutoffPick.Combined.Pos(df.both)[5]
  MethyGroup <- factor(MethyGroup, levels=c(TRUE, FALSE), labels = c("HyperMethy", "HypoMethy"))
  df.both$group  <- paste(MethyGroup, " & ", ExprGroup, sep="")
  df.combined <- df.both[df.both$group %in% c("HyperMethy & HighExpr", "HypoMethy & LowExpr"),]
  df.combined$group <- factor(df.combined$group, levels = c("HyperMethy & HighExpr", "HypoMethy & LowExpr"))
  
  survdiff(Surv(PFI.time, PFI)~ expr.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~ methy.group, data=df.both)
  survdiff(Surv(PFI.time, PFI)~group, data=df.combined)
  
  library(survminer)
  library(ggpubr)
  
  PFI.combined <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group,
                                     data = df.combined),
                             data = df.combined,
                             risk.table = F,
                             pval = TRUE,
                             break.time.by = 500,
                             ggtheme = theme_pubr(),
                             risk.table.y.text.col = TRUE,
                             risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("ITK PFS.combined.pdf", width=11.69/2, height = 8.268/2)
  PFI.combined
  dev.off()
  
  PFI.expr <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ expr.group,
                                 data = df.both),
                         data = df.both,
                         risk.table = F,
                         pval = TRUE,
                         break.time.by = 500,
                         ggtheme = theme_pubr(),
                         risk.table.y.text.col = TRUE,
                         risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("ITK PFS.expr.pdf", width=11.69/2, height = 8.268/2)
  PFI.expr
  dev.off()
  
  PFI.methy <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ methy.group,
                                  data = df.both),
                          data = df.both,
                          risk.table = F,
                          pval = TRUE,
                          break.time.by = 500,
                          ggtheme = theme_pubr(),
                          risk.table.y.text.col = TRUE,
                          risk.table.y.text = FALSE)+
    ylab("Progression-free survival probability")
  pdf("ITK PFS.methy.pdf", width=11.69/2, height = 8.268/2)
  PFI.methy
  dev.off()
}
# boxplot of methy in different tumor stage
CpGID = "chr16_2037931"
methy = unlist(LIHC.methy[CpGID,])
library(survival)
library(survminer)
#cox for methy
df.methy <- LIHC.data[,c("sampleID","condition", "age", "race", "gender", "stage2", "OS", "OS.time", "PFI", "PFI.time")]
df.methy <- df.methy[df.methy$sampleID %in% names(methy),]#keep the sample with both methylation and phenotypes
methy.scaled <- scale(methy)[,1]# transform the methylation data to Z score
df.methy$methy <- round(methy[match(df.methy$sampleID, names(methy))],2)
df.methy$methy.scaled <- round(methy.scaled[match(df.methy$sampleID, names(methy.scaled))],2)
df.methy <- df.methy[df.methy$condition == "Tumor",]#keep only tumor samples
df.methy <- df.methy[!is.na(df.methy$age) & !is.na(df.methy$race) & !is.na(df.methy$gender) & !is.na(df.methy$stage2) & !is.na(df.methy$OS) & !is.na(df.methy$OS.time) &
                       !is.na(df.methy$PFI) & !is.na(df.methy$PFI.time) & !is.na(df.methy$methy) & !is.na(df.methy$methy.scaled),]

stage.methy <- ggplot(data=df.methy, aes(x = stage2, y = methy))+
  geom_jitter(aes(fill = stage2), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_brewer(palette = "Reds", guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("Tumor stage")+
  ylab("Methylation of UCK2")+
  theme_pubr()+ggtitle("p = 0.016")
stage.methy
# boxplot of expr in different tumor stage
stage.expr <- ggplot(data=df.expr, aes(x = stage2, y = expr))+
  geom_jitter(aes(fill = stage2), position = position_jitter(0.3), shape=21, size=1.5, color = "black")+
  scale_fill_brewer(palette = "Reds", guide=F)+
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom = "pointrange", color = "black", size= 1.2)+
  stat_summary(fun.y = "mean", fun.args = list(mult=1), geom="point", color = "white", size= 4)+
  xlab("Tumor stage")+
  ylab("mRNA expression of UCK2")+
  theme_pubr()+ggtitle("p = 0.016")
stage.expr


