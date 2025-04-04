library(lme4)
setwd("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID/heatmap_trend_Barplot_Boxplot")

#### 
dt.orig <- read.csv("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_subset_data_merged_w_soc_expansion1_fix_sampleID.csv",stringsAsFactors = FALSE)#56284
table(dt.orig[,"cell_type"])

table(dt.orig[,"treatment"])
dt.orig[dt.orig[,"patient"]=="s2"|dt.orig[,"patient"]=="s3","treatment"] <- "pre"
dt.orig[dt.orig[,"patient"]=="s1","treatment"] <- NA

dt.orig[,"Response"] <- NA
dt.orig[dt.orig[,"patient"]=="a2"|dt.orig[,"patient"]=="a4","Response"] <- "Y"
dt.orig[dt.orig[,"patient"]=="a3"|dt.orig[,"patient"]=="a5","Response"] <- "N"
dt.orig[grep("s",dt.orig[,"patient"]),"patient_type"] <- "soc"
dt.orig[grep("a",dt.orig[,"patient"]),"patient_type"] <- "ateo"

#dt.patient <- aggregate(dt.orig[,"cell_type"], by=list(patient=dt.orig$patient), FUN=function(x) sum(!is.na(x)))
#dt.cell_type <- aggregate(dt.orig[,"cell_type"], by=list(cell_type=dt.orig$cell_type), FUN=function(x) sum(!is.na(x)))

idx <- dt.orig$TCR_size_gt1=="TCR_cluster_gt1" # we only looked the clusters with cluster size larger than 1
dt.orig <- dt.orig[idx,]#16597

#######
dt.TCR_cluster.bycelltype <- aggregate(dt.orig[,"cell_type"], 
                                      by=list(cell_type=dt.orig$cell_type,
                                      cluster=dt.orig$TCR_cluster_id,
                                      comp=dt.orig$comp,
                                      patient=dt.orig$patient,
                                      tissue=dt.orig$tissue,
                                      treatment=dt.orig$treatment), 
                                      FUN=function(x) sum(!is.na(x))) # 4865

dt.TCR_cluster.ncluster.bycelltype <- aggregate(dt.TCR_cluster.bycelltype[,"comp"], 
                                      by=list(cell_type=dt.TCR_cluster.bycelltype$cell_type,
                                             comp=dt.TCR_cluster.bycelltype$comp,
                                             patient=dt.TCR_cluster.bycelltype$patient,
                                             tissue=dt.TCR_cluster.bycelltype$tissue,
                                             treatment=dt.TCR_cluster.bycelltype$treatment), 
                                     FUN=function(x) sum(!is.na(x))) # x is the # of clusters in each combination, # 251

dt.TCR_cluster.ncluster.sum <- aggregate(dt.TCR_cluster.bycelltype[,"comp"], 
                                                by=list(comp=dt.TCR_cluster.bycelltype$comp,
                                                        patient=dt.TCR_cluster.bycelltype$patient,
                                                        tissue=dt.TCR_cluster.bycelltype$tissue,
                                                        treatment=dt.TCR_cluster.bycelltype$treatment), 
                                                FUN=function(x) sum(!is.na(x))) # x is the # of clusters in each combination, # 251

match.ls <- match(paste(dt.TCR_cluster.ncluster.bycelltype[,"comp"],dt.TCR_cluster.ncluster.bycelltype[,"patient"],dt.TCR_cluster.ncluster.bycelltype[,"tissue"],dt.TCR_cluster.ncluster.bycelltype[,"treatment"]),
                  paste(dt.TCR_cluster.ncluster.sum[,"comp"],dt.TCR_cluster.ncluster.sum[,"patient"],dt.TCR_cluster.ncluster.sum[,"tissue"],dt.TCR_cluster.ncluster.sum[,"treatment"]))
dt.TCR_cluster.ncluster.bycelltype[,"total"] <- dt.TCR_cluster.ncluster.sum[match.ls,"x"]
dt.TCR_cluster.ncluster.bycelltype[,"log10x"] <- log10(dt.TCR_cluster.ncluster.bycelltype[,"x"]/dt.TCR_cluster.ncluster.bycelltype[,"total"])
dt.TCR_cluster.ncluster.bycelltype[,"Response"] <- NA
dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"patient"]=="a3"|dt.TCR_cluster.ncluster.bycelltype[,"patient"]=="a5","Response"] <- "N"
dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"patient"]=="a2"|dt.TCR_cluster.ncluster.bycelltype[,"patient"]=="a4","Response"] <- "Y"
dt.TCR_cluster.ncluster.bycelltype[grep("a",dt.TCR_cluster.ncluster.bycelltype[,"patient"]),"patient_type"] <- "ateo"
dt.TCR_cluster.ncluster.bycelltype[grep("s",dt.TCR_cluster.ncluster.bycelltype[,"patient"]),"patient_type"] <- "soc"


dt.TCR_cluster.clustersize.sum <- aggregate(dt.orig[,"TCR_cluster_id"], 
                                           by=list(cluster=dt.orig$TCR_cluster_id,
                                                   comp=dt.orig$comp,
                                                   patient=dt.orig$patient,
                                                   tissue=dt.orig$tissue,
                                                   treatment=dt.orig$treatment), 
                                                   FUN=function(x) length(x)) # cluster size for each cluster. #2683

match.ls <- match(paste(dt.TCR_cluster.bycelltype[,"cluster"],dt.TCR_cluster.bycelltype[,"treatment"]),
                  paste(dt.TCR_cluster.clustersize.sum[,"cluster"],dt.TCR_cluster.clustersize.sum[,"treatment"]))
dt.TCR_cluster.bycelltype[,"total"] <- dt.TCR_cluster.clustersize.sum[match.ls,"x"]
dt.TCR_cluster.bycelltype[,"log10x"] <- log10(dt.TCR_cluster.bycelltype[,"x"]/dt.TCR_cluster.bycelltype[,"total"])

dt.TCR_cluster.clustersize.bycelltype <- aggregate(dt.TCR_cluster.bycelltype[,"log10x"], 
                                                by=list(cell_type=dt.TCR_cluster.bycelltype$cell_type,
                                                        comp=dt.TCR_cluster.bycelltype$comp,
                                                        patient=dt.TCR_cluster.bycelltype$patient,
                                                        tissue=dt.TCR_cluster.bycelltype$tissue,
                                                        treatment=dt.TCR_cluster.bycelltype$treatment), 
                                                FUN=function(x) max(x)) # x is the # of clusters in each combination, # 251

dt.TCR_cluster.clustersize.bycelltype[,"Response"] <- NA
dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"patient"]=="a3"|dt.TCR_cluster.clustersize.bycelltype[,"patient"]=="a5","Response"] <- "N"
dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"patient"]=="a2"|dt.TCR_cluster.clustersize.bycelltype[,"patient"]=="a4","Response"] <- "Y"
dt.TCR_cluster.clustersize.bycelltype[grep("a",dt.TCR_cluster.clustersize.bycelltype[,"patient"]),"patient_type"] <- "ateo"
dt.TCR_cluster.clustersize.bycelltype[grep("s",dt.TCR_cluster.clustersize.bycelltype[,"patient"]),"patient_type"] <- "soc"

### include shared information
dt.TCR_cluster.shared.bycelltype <- aggregate(dt.orig[,"cell_type"], 
                                              by=list(cell_type=dt.orig$cell_type,
                                                      cluster=dt.orig$TCR_cluster_id,
                                                      comp=dt.orig$comp,
                                                      patient=dt.orig$patient,
                                                      tissue=dt.orig$tissue,
                                                      treatment=dt.orig$treatment,
                                                      blood_tumor=dt.orig$blood_tumor), 
                                              FUN=function(x) sum(!is.na(x))) # 4865

dt.TCR_cluster.shared.ncluster.bycelltype <- aggregate(dt.TCR_cluster.shared.bycelltype[,"comp"], 
                                                       by=list(cell_type=dt.TCR_cluster.shared.bycelltype$cell_type,
                                                               comp=dt.TCR_cluster.shared.bycelltype$comp,
                                                               patient=dt.TCR_cluster.shared.bycelltype$patient,
                                                               tissue=dt.TCR_cluster.shared.bycelltype$tissue,
                                                               treatment=dt.TCR_cluster.shared.bycelltype$treatment,
                                                               blood_tumor=dt.TCR_cluster.shared.bycelltype$blood_tumor), 
                                                       FUN=function(x) sum(!is.na(x))) # x is the # of clusters in each combination, # 251

dt.TCR_cluster.shared.ncluster.sum <- aggregate(dt.TCR_cluster.shared.bycelltype[,"comp"], 
                                                by=list(comp=dt.TCR_cluster.shared.bycelltype$comp,
                                                        patient=dt.TCR_cluster.shared.bycelltype$patient,
                                                        tissue=dt.TCR_cluster.shared.bycelltype$tissue,
                                                        treatment=dt.TCR_cluster.shared.bycelltype$treatment,
                                                        blood_tumor=dt.TCR_cluster.shared.bycelltype$blood_tumor), 
                                                FUN=function(x) sum(!is.na(x))) # x is the # of clusters in each combination, # 251

match.ls <- match(paste(dt.TCR_cluster.shared.ncluster.bycelltype[,"comp"],dt.TCR_cluster.shared.ncluster.bycelltype[,"patient"],dt.TCR_cluster.shared.ncluster.bycelltype[,"tissue"],dt.TCR_cluster.shared.ncluster.bycelltype[,"treatment"],dt.TCR_cluster.shared.ncluster.bycelltype[,"blood_tumor"]),
                  paste(dt.TCR_cluster.shared.ncluster.sum[,"comp"],dt.TCR_cluster.shared.ncluster.sum[,"patient"],dt.TCR_cluster.shared.ncluster.sum[,"tissue"],dt.TCR_cluster.shared.ncluster.sum[,"treatment"],dt.TCR_cluster.shared.ncluster.sum[,"blood_tumor"]))
dt.TCR_cluster.shared.ncluster.bycelltype[,"total"] <- dt.TCR_cluster.shared.ncluster.sum[match.ls,"x"]
dt.TCR_cluster.shared.ncluster.bycelltype[,"log10x"] <- log10(dt.TCR_cluster.shared.ncluster.bycelltype[,"x"]/dt.TCR_cluster.shared.ncluster.bycelltype[,"total"])
dt.TCR_cluster.shared.ncluster.bycelltype[,"Response"] <- NA
dt.TCR_cluster.shared.ncluster.bycelltype[dt.TCR_cluster.shared.ncluster.bycelltype[,"patient"]=="a3"|dt.TCR_cluster.shared.ncluster.bycelltype[,"patient"]=="a5","Response"] <- "N"
dt.TCR_cluster.shared.ncluster.bycelltype[dt.TCR_cluster.shared.ncluster.bycelltype[,"patient"]=="a2"|dt.TCR_cluster.shared.ncluster.bycelltype[,"patient"]=="a4","Response"] <- "Y"
dt.TCR_cluster.shared.ncluster.bycelltype[grep("a",dt.TCR_cluster.shared.ncluster.bycelltype[,"patient"]),"patient_type"] <- "ateo"
dt.TCR_cluster.shared.ncluster.bycelltype[grep("s",dt.TCR_cluster.shared.ncluster.bycelltype[,"patient"]),"patient_type"] <- "soc"


dt.TCR_cluster.shared.clustersize.sum <- aggregate(dt.orig[,"TCR_cluster_id"], 
                                                   by=list(cluster=dt.orig$TCR_cluster_id,
                                                           comp=dt.orig$comp,
                                                           patient=dt.orig$patient,
                                                           tissue=dt.orig$tissue,
                                                           treatment=dt.orig$treatment,
                                                           blood_tumor=dt.orig$blood_tumor), 
                                                   FUN=function(x) length(x)) # cluster size for each cluster. #2683

match.ls <- match(paste(dt.TCR_cluster.shared.bycelltype[,"cluster"],dt.TCR_cluster.shared.bycelltype[,"treatment"],dt.TCR_cluster.shared.bycelltype[,"blood_tumor"]),
                  paste(dt.TCR_cluster.shared.clustersize.sum[,"cluster"],dt.TCR_cluster.shared.clustersize.sum[,"treatment"],dt.TCR_cluster.shared.clustersize.sum[,"blood_tumor"]))
dt.TCR_cluster.shared.bycelltype[,"total"] <- dt.TCR_cluster.shared.clustersize.sum[match.ls,"x"]
dt.TCR_cluster.shared.bycelltype[,"log10x"] <- log10(dt.TCR_cluster.shared.bycelltype[,"x"]/dt.TCR_cluster.shared.bycelltype[,"total"])

dt.TCR_cluster.shared.clustersize.bycelltype <- aggregate(dt.TCR_cluster.shared.bycelltype[,"log10x"], 
                                                          by=list(cell_type=dt.TCR_cluster.shared.bycelltype$cell_type,
                                                                  comp=dt.TCR_cluster.shared.bycelltype$comp,
                                                                  patient=dt.TCR_cluster.shared.bycelltype$patient,
                                                                  tissue=dt.TCR_cluster.shared.bycelltype$tissue,
                                                                  treatment=dt.TCR_cluster.shared.bycelltype$treatment,
                                                                  blood_tumor=dt.TCR_cluster.shared.bycelltype$blood_tumor), 
                                                          FUN=function(x) max(x)) # x is the # of clusters in each combination, # 251

dt.TCR_cluster.shared.clustersize.bycelltype[,"Response"] <- NA
dt.TCR_cluster.shared.clustersize.bycelltype[dt.TCR_cluster.shared.clustersize.bycelltype[,"patient"]=="a3"|dt.TCR_cluster.shared.clustersize.bycelltype[,"patient"]=="a5","Response"] <- "N"
dt.TCR_cluster.shared.clustersize.bycelltype[dt.TCR_cluster.shared.clustersize.bycelltype[,"patient"]=="a2"|dt.TCR_cluster.shared.clustersize.bycelltype[,"patient"]=="a4","Response"] <- "Y"
dt.TCR_cluster.shared.clustersize.bycelltype[grep("a",dt.TCR_cluster.shared.clustersize.bycelltype[,"patient"]),"patient_type"] <- "ateo"
dt.TCR_cluster.shared.clustersize.bycelltype[grep("s",dt.TCR_cluster.shared.clustersize.bycelltype[,"patient"]),"patient_type"] <- "soc"


## calucated pvalue based on anova test for each cell
cell_type <- names(table(dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]))
## CD4 only, ncluster
rt.pvalue.ncluster.bycell.CD4 <- matrix(NA,nrow=length(cell_type)-3,ncol=7)
k <- 1
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD4",]
  fm00 <- lmer(log10x ~ (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm0 <- lmer(log10x ~ patient_type+(1 | patient), dt)
  fm.txt <- lmer(log10x ~ treatment + patient_type+(1 | patient), dt)
  fm.response <- lmer(log10x ~ Response + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  
  fm.tissue <- lmer(log10x ~ tissue + patient_type+(1 | patient), dt)
  fm.txt.tissue <- lmer(log10x ~ treatment +tissue+ patient_type+(1 | patient), dt)
  
  fm.tissue0 <- lmer(log10x ~ tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm.response.tissue <- lmer(log10x ~ Response+tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  ##
  dt.shared <- dt.TCR_cluster.shared.ncluster.bycelltype[dt.TCR_cluster.shared.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.ncluster.bycelltype[,"comp"]=="CD4",]
  fm.blood_tumor0 <- lmer(log10x ~ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor <- lmer(log10x ~ blood_tumor + patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue0 <- lmer(log10x ~ tissue+ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue <- lmer(log10x ~ blood_tumor +tissue+ patient_type+(1 | patient), dt.shared)
  
  rt.pvalue.ncluster.bycell.CD4[k,1:7] <- c(anova(fm.txt,fm0)[2,"Pr(>Chisq)"],# treatment effect
                                   anova(fm.response,fm00)[2,"Pr(>Chisq)"],# response
                                   anova(fm.blood_tumor,fm.blood_tumor0)[2,"Pr(>Chisq)"],# overlap vs. nonoverlap
                                   anova(fm.tissue,fm0)[2,"Pr(>Chisq)"], # tissue or blood effect
                                   anova(fm.txt.tissue,fm.tissue)[2,"Pr(>Chisq)"],# treatment effect after accounting speciment type
                                   anova(fm.response.tissue,fm.tissue0)[2,"Pr(>Chisq)"],# response after accounting speciment type
                                   anova(fm.blood_tumor.tissue,fm.blood_tumor.tissue0)[2,"Pr(>Chisq)"])# overlap vs. nonoverlap after accounting speciment type
  
  k <- k+1
}
rt.pvalue.ncluster.bycell.CD4 <- round(rt.pvalue.ncluster.bycell.CD4,digits=3)
colnames(rt.pvalue.ncluster.bycell.CD4) <- c("treatment","response","blood_tumor","tissue",
                                    "treatment.tissue","response.tissue","blood_tumor.tissue")
rownames(rt.pvalue.ncluster.bycell.CD4) <- cell_type[-c(1,6,7)]
write.csv(rt.pvalue.ncluster.bycell.CD4, file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/rt.pvalue.ncluster.bycell.CD4.csv")

## CD4 only, clustersize
rt.pvalue.clustersize.bycell.CD4 <- matrix(NA,nrow=length(cell_type)-3,ncol=7)
k <- 1
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD4",]
  fm00 <- lmer(x ~ (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm0 <- lmer(x ~ patient_type+(1 | patient), dt)
  fm.txt <- lmer(x ~ treatment + patient_type+(1 | patient), dt)
  fm.response <- lmer(x ~ Response + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  
  fm.tissue <- lmer(x ~ tissue + patient_type+(1 | patient), dt)
  fm.txt.tissue <- lmer(x ~ treatment +tissue+ patient_type+(1 | patient), dt)
  
  fm.tissue0 <- lmer(x ~ tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm.response.tissue <- lmer(x ~ Response+tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  ##
  dt.shared <- dt.TCR_cluster.shared.clustersize.bycelltype[dt.TCR_cluster.shared.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.clustersize.bycelltype[,"comp"]=="CD4",]
  fm.blood_tumor0 <- lmer(x ~ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor <- lmer(x ~ blood_tumor + patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue0 <- lmer(x ~ tissue+ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue <- lmer(x ~ blood_tumor +tissue+ patient_type+(1 | patient), dt.shared)
  
  rt.pvalue.clustersize.bycell.CD4[k,1:7] <- c(anova(fm.txt,fm0)[2,"Pr(>Chisq)"],# treatment effect
                                               anova(fm.response,fm00)[2,"Pr(>Chisq)"],# response
                                               anova(fm.blood_tumor,fm.blood_tumor0)[2,"Pr(>Chisq)"],# overlap vs. nonoverlap
                                               anova(fm.tissue,fm0)[2,"Pr(>Chisq)"], # tissue or blood effect
                                               anova(fm.txt.tissue,fm.tissue)[2,"Pr(>Chisq)"],# treatment effect after accounting speciment type
                                               anova(fm.response.tissue,fm.tissue0)[2,"Pr(>Chisq)"],# response after accounting speciment type
                                               anova(fm.blood_tumor.tissue,fm.blood_tumor.tissue0)[2,"Pr(>Chisq)"])# overlap vs. nonoverlap after accounting speciment type
  
  k <- k+1
}
rt.pvalue.clustersize.bycell.CD4 <- round(rt.pvalue.clustersize.bycell.CD4,digits=3)
colnames(rt.pvalue.clustersize.bycell.CD4) <- c("treatment","response","blood_tumor","tissue",
                                                "treatment.tissue","response.tissue","blood_tumor.tissue")
rownames(rt.pvalue.clustersize.bycell.CD4) <- cell_type[-c(1,6,7)]
write.csv(rt.pvalue.clustersize.bycell.CD4, file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/rt.pvalue.clustersize.bycell.CD4.csv")


## CD8 only, ncluster
rt.pvalue.ncluster.bycell.CD8 <- matrix(NA,nrow=length(cell_type)-2,ncol=7)
k <- 1
for (i in cell_type[-c(1,6)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD8",]
  fm00 <- lmer(log10x ~ (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm0 <- lmer(log10x ~ patient_type+(1 | patient), dt)
  fm.txt <- lmer(log10x ~ treatment + patient_type+(1 | patient), dt)
  fm.response <- lmer(log10x ~ Response + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  
  fm.tissue <- lmer(log10x ~ tissue + patient_type+(1 | patient), dt)
  fm.txt.tissue <- lmer(log10x ~ treatment +tissue+ patient_type+(1 | patient), dt)
  
  fm.tissue0 <- lmer(log10x ~ tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm.response.tissue <- lmer(log10x ~ Response+tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  ##
  dt.shared <- dt.TCR_cluster.shared.ncluster.bycelltype[dt.TCR_cluster.shared.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.ncluster.bycelltype[,"comp"]=="CD8",]
  fm.blood_tumor0 <- lmer(log10x ~ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor <- lmer(log10x ~ blood_tumor + patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue0 <- lmer(log10x ~ tissue+ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue <- lmer(log10x ~ blood_tumor +tissue+ patient_type+(1 | patient), dt.shared)
  
  rt.pvalue.ncluster.bycell.CD8[k,1:7] <- c(anova(fm.txt,fm0)[2,"Pr(>Chisq)"],# treatment effect
                                            anova(fm.response,fm00)[2,"Pr(>Chisq)"],# response
                                            anova(fm.blood_tumor,fm.blood_tumor0)[2,"Pr(>Chisq)"],# overlap vs. nonoverlap
                                            anova(fm.tissue,fm0)[2,"Pr(>Chisq)"], # tissue or blood effect
                                            anova(fm.txt.tissue,fm.tissue)[2,"Pr(>Chisq)"],# treatment effect after accounting speciment type
                                            anova(fm.response.tissue,fm.tissue0)[2,"Pr(>Chisq)"],# response after accounting speciment type
                                            anova(fm.blood_tumor.tissue,fm.blood_tumor.tissue0)[2,"Pr(>Chisq)"])# overlap vs. nonoverlap after accounting speciment type
  
  k <- k+1
}
rt.pvalue.ncluster.bycell.CD8 <- round(rt.pvalue.ncluster.bycell.CD8,digits=3)
colnames(rt.pvalue.ncluster.bycell.CD8) <- c("treatment","response","blood_tumor","tissue",
                                             "treatment.tissue","response.tissue","blood_tumor.tissue")
rownames(rt.pvalue.ncluster.bycell.CD8) <- cell_type[-c(1,6)]
write.csv(rt.pvalue.ncluster.bycell.CD8, file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/rt.pvalue.ncluster.bycell.CD8.csv")

## CD8 only, clustersize
rt.pvalue.clustersize.bycell.CD8 <- matrix(NA,nrow=length(cell_type)-2,ncol=7)
k <- 1
for (i in cell_type[-c(1,6)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD8",]
  fm00 <- lmer(x ~ (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm0 <- lmer(x ~ patient_type+(1 | patient), dt)
  fm.txt <- lmer(x ~ treatment + patient_type+(1 | patient), dt)
  fm.response <- lmer(x ~ Response + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  
  fm.tissue <- lmer(x ~ tissue + patient_type+(1 | patient), dt)
  fm.txt.tissue <- lmer(x ~ treatment +tissue+ patient_type+(1 | patient), dt)
  
  fm.tissue0 <- lmer(x ~ tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  fm.response.tissue <- lmer(x ~ Response+tissue + (1 | patient), dt[dt[,"patient_type"]=="ateo",])
  ##
  dt.shared <- dt.TCR_cluster.shared.clustersize.bycelltype[dt.TCR_cluster.shared.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.clustersize.bycelltype[,"comp"]=="CD8",]
  fm.blood_tumor0 <- lmer(x ~ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor <- lmer(x ~ blood_tumor + patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue0 <- lmer(x ~ tissue+ patient_type+(1 | patient), dt.shared)
  fm.blood_tumor.tissue <- lmer(x ~ blood_tumor +tissue+ patient_type+(1 | patient), dt.shared)
  
  rt.pvalue.clustersize.bycell.CD8[k,1:7] <- c(anova(fm.txt,fm0)[2,"Pr(>Chisq)"],# treatment effect
                                               anova(fm.response,fm00)[2,"Pr(>Chisq)"],# response
                                               anova(fm.blood_tumor,fm.blood_tumor0)[2,"Pr(>Chisq)"],# overlap vs. nonoverlap
                                               anova(fm.tissue,fm0)[2,"Pr(>Chisq)"], # tissue or blood effect
                                               anova(fm.txt.tissue,fm.tissue)[2,"Pr(>Chisq)"],# treatment effect after accounting speciment type
                                               anova(fm.response.tissue,fm.tissue0)[2,"Pr(>Chisq)"],# response after accounting speciment type
                                               anova(fm.blood_tumor.tissue,fm.blood_tumor.tissue0)[2,"Pr(>Chisq)"])# overlap vs. nonoverlap after accounting speciment type
  
  k <- k+1
}
rt.pvalue.clustersize.bycell.CD8 <- round(rt.pvalue.clustersize.bycell.CD8,digits=3)
colnames(rt.pvalue.clustersize.bycell.CD8) <- c("treatment","response","blood_tumor","tissue",
                                                "treatment.tissue","response.tissue","blood_tumor.tissue")
rownames(rt.pvalue.clustersize.bycell.CD8) <- cell_type[-c(1,6)]
write.csv(rt.pvalue.clustersize.bycell.CD8, file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/rt.pvalue.clustersize.bycell.CD8.csv")



## get details for each significant findings
### function
coef.fn <- function(fm)
{
  coefs <- data.frame(coef(summary(fm)))
  coefs[,"pvalue"] <- round(2*(1 - pnorm(abs(coefs$t.value))),digits=3)
  coefs[,"LCI"] <- round(coefs[,"Estimate"]-1.96*coefs[,"Std..Error"],digits=2)
  coefs[,"UCI"] <- round(coefs[,"Estimate"]+1.96*coefs[,"Std..Error"],digits=2)
  coefs[,"Estimate"] <- round(coefs[,"Estimate"],digits=2)
  coefs[,"rt"] <- paste(coefs[,"Estimate"],"(",coefs[,"LCI"],",",coefs[,"UCI"],"p=",coefs[,"pvalue"],")",sep=" ")
  return(coefs)
}

### CD4 and ncluster
coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(log10x ~ tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.tissue.CD4.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(log10x ~ treatment + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.treatment.CD4.csv")

# no signficant
coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(log10x ~ Response + tissue+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.response.CD4.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.shared.ncluster.bycelltype[dt.TCR_cluster.shared.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.ncluster.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(log10x ~ blood_tumor + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.shared.CD4.csv")
### CD4 and clustersize
coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(x ~ tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.tissue.CD4.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(x ~ treatment + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.treatment.CD4.csv")

# no signficant
coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(x ~ Response + tissue+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.response.CD4.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,6,7)])
{
  dt <- dt.TCR_cluster.shared.clustersize.bycelltype[dt.TCR_cluster.shared.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.clustersize.bycelltype[,"comp"]=="CD4",]
  fm <- lmer(x ~ blood_tumor + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.shared.CD4.csv")

### CD8 and ncluster
coef.rt <- NULL
for (i in cell_type[-c(1,6)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(log10x ~ tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.tissue.CD8.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,6)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(log10x ~ treatment + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.treatment.CD8.csv")

# no signficant
coef.rt <- NULL
for (i in cell_type[-c(1,6)])
{
  dt <- dt.TCR_cluster.ncluster.bycelltype[dt.TCR_cluster.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.ncluster.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(log10x ~ Response + tissue+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.response.CD8.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,6)])
{
  dt <- dt.TCR_cluster.shared.ncluster.bycelltype[dt.TCR_cluster.shared.ncluster.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.ncluster.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(log10x ~ blood_tumor + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.ncluster.shared.CD8.csv")
### CD8 and clustersize
coef.rt <- NULL
for (i in cell_type[-c(1,5,6)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(x ~ tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.tissue.CD8.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,5,6)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(x ~ treatment + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.treatment.CD8.csv")

# no signficant
coef.rt <- NULL
for (i in cell_type[-c(1,5,6)])
{
  dt <- dt.TCR_cluster.clustersize.bycelltype[dt.TCR_cluster.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.clustersize.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(x ~ Response + tissue+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.response.CD8.csv")

coef.rt <- NULL
for (i in cell_type[-c(1,5,6)])
{
  dt <- dt.TCR_cluster.shared.clustersize.bycelltype[dt.TCR_cluster.shared.clustersize.bycelltype[,"cell_type"]==i & dt.TCR_cluster.shared.clustersize.bycelltype[,"comp"]=="CD8",]
  fm <- lmer(x ~ blood_tumor + tissue+patient_type+(1 | patient), dt)
  coef.rt <- rbind(coef.rt,cbind(i,coef.fn(fm)))
}
coef.rt <- coef.rt[-grep("Intercept",rownames(coef.rt)),]
coef.rt <- coef.rt[coef.rt[,"pvalue"]< 0.05,c("i","rt")]
write.csv(coef.rt,file="./Results/rt_cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype/coef.clustersize.shared.CD8.csv")

save.image("./Results/cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype.Rd")

