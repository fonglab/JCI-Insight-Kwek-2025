# convert raw output of TCR to match format with Tony's paired TCR format
library(plyr)

setwd("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/Data")
Atezo_Bladder_Total_paired_tcr <- read.csv("Atezo_Bladder_Total_paired_tcr_fix_sampleID.csv")
#### fix the issue of a2 tumor and blood from two pations
# Get levels and add "a6"
df <- Atezo_Bladder_Total_paired_tcr
levels <- levels(df$patient)
levels[length(levels) + 1] <- "a6"
df$patient <- factor(df$patient, levels = levels)
df$patient[df$patient == "a2" & df$tissue == "blood"] <- 'a6'
Atezo_Bladder_Total_paired_tcr <- df
write.csv(Atezo_Bladder_Total_paired_tcr, file = "Atezo_Bladder_Total_paired_tcr_fix_sampleID.csv", row.names = FALSE)


cell_type_info <- read.csv("bladder_total_adata_multiplex_Tcellonly_obs_data_w_TCR_info_fix_sampleID.csv")
dim(cell_type_info)
#### fix the issue of a2 tumor and blood from two pations
df <- cell_type_info
levels <- levels(df$patient)
levels[length(levels) + 1] <- "a6"
df$patient <- factor(df$patient, levels = levels)
df$patient[df$patient == "a2" & df$tissue == "blood"] <- 'a6'
cell_type_info <- df
write.csv(cell_type_info, file = "bladder_total_adata_multiplex_Tcellonly_obs_data_w_TCR_info_fix_sampleID.csv", row.names = FALSE)


input6 <- read.csv("./TCR_raw_files/1203CD4-731_TCR-clones_CBsummary.csv")
dim(input6)
# 4876   13
head(input6$cell.barcode)

length(unique(input6$cell.barcode))
# 3154

barcode_table <- as.data.frame(table(input6$cell.barcode))

dim(barcode_table[barcode_table$Freq ==2,])
# 1722    2
barcode_paired <- barcode_table[barcode_table$Freq ==2,]$Var1
length(barcode_paired) # 1722

TCRa <- input6[input6$Chains == "TRA" & input6$cell.barcode %in% barcode_paired, ]
drops <- c("Chains", "cloneFraction", "clonalSequenceQuality", "bestVHitScore", "bestJHitScore", "bestDHitScore" )
TCRa <-  TCRa[ , !(names(TCRa) %in% drops)]
colnames(TCRa) <- c("alphacloneCount","alphaClonalSeq", "alphaCDR3", "alphaVGene", "alphaJGene", "alphaDGene", "cell.barcode")

TCRb <- input6[input6$Chains == "TRB" & input6$cell.barcode %in% barcode_paired, ]
TCRb <-  TCRb[ , !(names(TCRb) %in% drops)]
colnames(TCRb) <- c("betacloneCount","betaClonalSeq", "betaCDR3", "betaVGene", "betaJGene", "betaDGene", "cell.barcode")

paired_TCRab <- merge(TCRa, TCRb, by = "cell.barcode")
dim(paired_TCRab)
# [1] 1722   13

paired_TCRab <- merge(paired_TCRab, cell_type_info[,c("cell.barcode", "cell_calls")], by = "cell.barcode",  all.x = TRUE)
colnames(paired_TCRab)[length(paired_TCRab)] <- "cell_type" 
table(paired_TCRab$cell_type)
# Activated           CM      CXCL13+        GZMB+        GZMK+         IFN+         MAIT         Mito        Naive Prolif GZMK+        Tregs 
#         0          197          213          213           21            0            2          375          587           38           68 
# compare with all a5 CD4 blood atezo on cell info
sub_all_cell_info <- cell_type_info[cell_type_info$patient == "a5" & cell_type_info$tissue == "blood" & cell_type_info$comp == "CD4" & cell_type_info$treatment == "atezo",]
table(sub_all_cell_info$cell_calls)
# Activated           CM      CXCL13+        GZMB+        GZMK+         IFN+         MAIT         Mito        Naive Prolif GZMK+        Tregs 
#         0          398          600          614           48            2            3         1142         1582           52          143
nrow(sub_all_cell_info) # 4584
nrow(paired_TCRab) # 1722
# 1722/4584 = 37.6%

paired_TCRab$Sample_ID <- "Anti-PD-L1 D blood CD4 post"
paired_TCRab$cloneCount <- paired_TCRab$alphacloneCount + paired_TCRab$betacloneCount
paired_TCRab$patient <- "a5"
paired_TCRab$comp <- "CD4"
paired_TCRab$tissue <- "blood"
paired_TCRab$treatment <- "atezo"
paired_TCRab$ab <- paste0(paired_TCRab$alphaClonalSeq, ";", paired_TCRab$betaClonalSeq)
dim(paired_TCRab) # 1722   21

# merge the sample into pervious paired TCR file
update_Total_paired_tcr <- rbind.fill(Atezo_Bladder_Total_paired_tcr, paired_TCRab)
dim(Atezo_Bladder_Total_paired_tcr) # 76084    26
dim(update_Total_paired_tcr) # 77806    28
# 1722 more paired cells

write.csv(update_Total_paired_tcr, "Atezo_Bladder_Total_paired_tcr_updated_fix_sampleID.csv", row.names=FALSE)
# no row order change in pervious cells.
