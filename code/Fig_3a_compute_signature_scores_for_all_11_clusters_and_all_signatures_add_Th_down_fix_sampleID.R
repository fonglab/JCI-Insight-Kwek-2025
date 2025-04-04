library(ggplot2)
library(reshape2)
library(readxl)

setwd("/Users/hai/Documents/workstation/bladder_blood_analysis/scRNA/output_for_paper/A_gene_signatures/")
excel_gene_list_dir <- "/Users/hai/Documents/workstation/bladder_blood_analysis/metadata/Hollbacher\ 2020_Supplemental_Table_2.xlsx"
excel_sheets(path = excel_gene_list_dir)


#### gene signatures ####
gene_list_name <- "Kumar_2017_Trm"
Kumar_2017_Trm_up <- c('CA10', 'ITGA1', 'ITGAE', 'IL2', 'IL10', 'CXCR6', 'CXCL13', 'KCNK5', 'RGS1', 'CRTAM', 'DUSP6', 'PDCD1', 'IL23R' )
Kumar_2017_Trm_down <- c('STK38', 'TTC16', 'SELL', 'KLF3', 'KLF2', 'D4S234E', 'SBK1', 'FAM65B', 'TTYH2', 'NPDC1', 'KRT72', 'S1PR1', 'SOX13', 'KRT73', 'TSPAN18', 'PTGDS', 'PRAP1GAP2', 'CX3CR1')

gene_list_name <- "Crawford_Ex"
Crawford_Ex_up <- c("Nap1l3", "Ifi44" , "Lag3", "Eomes", "Mki67", "Rsad2", "Arap2", "Ifih1" , "Naip5", "Kif15", "Mx1", "Serpina3f", "Usp18", "Irf4", "Nab1" , "Tnfsf4", 'Cd86', "Plod2", "Raet1b", 'Diap3' , "Ccnb2" , "Tpx2" , "Sh2d1a", "Ccna2", "fi205", "Depdc1a", "Cenpe", "Gpr141", "Dtl"  )
Crawford_Ex_down <- c("Taf1a" , "Entpd5", "Idh2", "Dph5", "Nudt15", "Nlk" , "Snord35b", "Txnip" , "Klhl6" , "Grap2" , "Satb1")

gene_list_name <- "Miller_Ex_progenitors"
Miller_Ex_progenitors_genes <- c('XCL1', 'MS4A4C', 'ID3', 'CXCL10', 'SLAMF6', 'TCF7', 'IL7R', 'TNFSF8', 'CTLA2A', 'RPS20', 'LTB', 'CD9', 'RPL36A', 'SOCS3', 'RPL35', 'RPS2', 'RPS19', 'JUN', 'RPS17', 'GPR183', 'RPLP0', 'RPL13', 'TRAF1', 'RPL31', 'ITGB1', 'MT-CYTB', 'RPS18', 'RPL10-PS3', 'RPL12', 'GM8730', 'RPS28', 'EMB', 'RPL21', 'RPS29', 'RPS26', 'RPSA', 'RPS7', 'GM10073', 'MT-CO3'  )

gene_list_name <- "Miller_Ex_terminal"
Miller_Ex_terminal_genes <- c('CD7', 'CD160', 'RGS1', 'CXCR6', 'GZMA', 'CCL3', 'CCL4', 'NR4A2', 'LAG3', 'KIAA1671', 'PLAC8', 'ABI3', 'SH2D2A', 'PTGER4', 'AW112010', 'ISG15', 'CD3G', 'SERPINA3G', 'LAX1', 'PDCD1', 'GZMK', 'GZMB', 'MBNL1', 'GLRX', 'ID2', 'RGS3', 'ARL6IP1', 'CD8A', 'EFHD2', 'FASL', 'PTPN22', 'VMP1', 'DUSP2', 'IFI47', 'GIMAP7', 'ITPKB', 'TAPBPL', 'STAT1', 'SHISA5'  )

gene_list_name <- "Cano_effector"
Cano_effector_genes <- c('TNFRSF18', 'LMNA', 'GNLY', 'IQCG', 'HOPX', 'GZMA', 'TMEM173', 'HLA-DRA', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DPB1', 'ACTB', 'CTSW', 'C10orf128', 'IFNG', 'GZMB', 'NFKBIA', 'LGMN', 'CCL5', 'CCL4', 'NKG7',  'EPAS1',  'CSF2', 'CCL3')
  
gene_list_name <- "Hollbacher_Th1"
Hollbacher_Th1up_genes <- c('ABHD17A', 'ABI3', 'ACRC', 'ACSL3', 'ADAMTS10', 'AGPAT5', 'AKNA', 'APMAP', 'APPBP2', 'ARAP2', 'ARHGAP26', 'ARPC5L', 'ASAH2B', 'ATP10A', 'BACH2', 'BBC3', 'BCR', 'BTBD11', 'BTN3A1', 'C9orf91', 'CARHSP1', 'CBLB', 'CCDC64', 'CCDC78', 'CCNDBP1', 'CDCA7L', 'CDK5R1', 'CLHC1', 'CLIC6', 'CLN8', 'CMC1', 'CMPK2', 'CNN3', 'COL6A2', 'CRTAM', 'CTSW', 'CXCR3', 'DFNB31', 'DIP2A', 'DKK3', 'DOCK3', 'DRAM1', 'DTHD1', 'E2F3', 'ENC1', 'ENG', 'ENTPD6', 'EOMES', 'EPB41L5', 'EXPH5', 'F2R', 'FCRL6', 'FGFBP2', 'FLT4', 'FOXJ1', 'FYN', 'GBP4', 'GBP5', 'GFOD1', 'GGT7', 'GNB3', 'GSE1', 'GZMK', 'GZMM', 'HAVCR1', 'HBEGF', 'HIC1', 'IFI27', 'IFNG', 'IGF1R', 'IGFBP3', 'IGFBP4', 'IL12RB2', 'INPP5A', 'IRF8', 'IVD', 'JAK2', 'KIAA1522', 'KIAA1671', 'L3MBTL4', 'LDLRAP1', 'LGR6', 'LITAF', 'LMO7', 'LPCAT1', 'LUZP1', 'LYAR', 'MAN2A1', 'MANBA', 'MAPRE3', 'MATK', 'MCTP2', 'ME3', 'MPP1', 'MZB1', 'NCALD', 'P3H3', 'PACSIN1', 'PDE9A', 'PDSS1', 'PDZD8', 'PHLPP1', 'PMAIP1', 'PMEPA1', 'PRF1', 'PRMT2', 'PSME4', 'PTK2B', 'PTPRE', 'PTPRJ', 'PWP2', 'QSOX2', 'RAB24', 'RASGEF1B', 'RASGRF2', 'RHBDF2', 'ROBO3', 'ROR2', 'RRAGD', 'RSAD2', 'RTP5', 'RUFY4', 'RUSC2', 'S100PBP', 'S1PR5', 'SAMD3', 'SAV1', 'SCD', 'SH2B3', 'SH3RF3', 'SIPA1', 'SIRT3', 'SLA2', 'SLAMF7', 'SLC4A4', 'SMAD3', 'SMAD7', 'SNTB1', 'SOS2', 'SPATA20', 'SRPK3', 'SYNM', 'TBKBP1', 'TBX21', 'TMEM109', 'TMEM62', 'TNS4', 'TRAPPC10', 'TRPS1', 'TSPAN33', 'TTC23', 'TULP4', 'UPP1', 'XCL1', 'YARS', 'YES1', 'ZBTB49', 'ZC4H2', 'ZCCHC18', 'ZIK1', 'ZNF135', 'ZNF541', 'ZNF618', 'ZNF629', 'ZNF711'  )
Th1_down <- read_excel(path = excel_gene_list_dir, sheet = "Th1 down",col_names =F)
Hollbacher_Th1down_genes <- Th1_down$...1


gene_list_name <- "Hollbacher_Th2"
Hollbacher_Th2up_genes <- c('ADGRE1', 'AEBP1', 'AKAP12', 'ALS2CL', 'ANKRD26', 'ARMCX4', 'ASAP1', 'AXIN2', 'BCAS4', 'BEND5', 'C1orf228', 'CACNA1I', 'CADM1', 'CCDC141', 'CCND1', 'CD55', 'CD8B', 'CDHR3', 'CDO1', 'CLECL1', 'CLIP3', 'CNKSR2', 'CRLF2', 'CYSLTR1', 'DBH', 'DCP2', 'DENND5A', 'DLX2', 'EDEM3', 'ENDOD1', 'EVI2A', 'FHIT', 'FNIP2', 'FXYD5', 'GALT', 'GATA3', 'GCNT4', 'GCOM1', 'GNA11', 'GNG8', 'GOLGA7', 'HEATR5B', 'HS3ST3B1', 'ICOS', 'IKZF2', 'IL11RA', 'IL16', 'IL4R', 'IL6ST', 'JAKMIP3', 'KCNQ1', 'KIF5B', 'KLF5', 'KLF7', 'KRT2', 'LAPTM4A', 'LEF1', 'LETM2', 'LPAR2', 'MAN1C1', 'MAP4K4', 'MASTL', 'MICU3', 'MRC2', 'NDFIP1', 'NDFIP2', 'NDUFAF7', 'NET1', 'NOL11', 'PASK', 'PKP2', 'PLXDC1', 'PPP1R9B', 'PRKCA', 'PTGDR2', 'PTK2', 'PUM1', 'R3HDM4', 'RBL1', 'RHBDD2', 'RILPL2', 'RIN3', 'S1PR1', 'SARAF', 'SBNO2', 'SELL', 'SESN3', 'SGSM3', 'SLC22A23', 'SLC25A33', 'SLC40A1', 'SLFN5', 'SNED1', 'SPINT2', 'STK4', 'TAB3', 'TESPA1', 'TMEM201', 'TNFSF11', 'TNNC1', 'TRIM16', 'UBXN2A', 'VIPR1', 'ZBED3'  )
Th2_down <- read_excel(path = excel_gene_list_dir, sheet = "Th2 down",col_names =F)
Hollbacher_Th2down_genes <- Th2_down$...1

gene_list_name <- "Hollbacher_Th17"
Hollbacher_Th17up_genes <- c('AATK', 'ABCA1', 'ABCB1', 'ADAM10', 'AHCY', 'ALG1', 'AMMECR1', 'ANK1', 'ANO6', 'AP3S2', 'ATF5', 'BCCIP', 'C16orf74', 'C2CD4B', 'C4orf46', 'CASR', 'CCDC6', 'CEBPG', 'CEP97', 'CERK', 'CFL2', 'CMTM6', 'COMMD3', 'CR1', 'CTLA4', 'DAB1', 'DNAJC10', 'DST', 'DUSP16', 'EXT1', 'FAM126A', 'FAM63B', 'FBL', 'FRMD4A', 'FURIN', 'FUT10', 'GDF7', 'GMDS', 'GRAMD3', 'GTF3C4', 'HADHA', 'HERC6', 'HRH4', 'IKZF3', 'IL12RB1', 'IL1R2', 'INPP5F', 'JMY', 'KLHL42', 'KLHL5', 'LDLRAD4', 'LSM10', 'LTA4H', 'LTB', 'MAN1A1', 'MAP3K4', 'MGAT4A', 'MIEN1', 'MRPL34', 'NPAS2', 'PALB2', 'PANK4', 'PDDC1', 'PDE4D', 'PHF21A', 'PLXNC1', 'PPA1', 'PPIF', 'PRKRIR', 'PSD3', 'PSMG3', 'PTPN4', 'RAD54B', 'RANBP9', 'RPLP0', 'SBK1', 'SEMA6C', 'SESN1', 'SKA3', 'SLAMF1', 'SLC35F2', 'SLC35G1', 'SOCS2', 'SRD5A3', 'SYNGAP1', 'TBC1D4', 'TLR2', 'TMED8', 'TMEM156', 'TMEM60', 'TP53INP1', 'TRAF3IP2', 'UBAP1', 'UFC1', 'UROS', 'USP10', 'YEATS4', 'YWHAH', 'ZC3H12D'  )
Th17_down <- read_excel(path = excel_gene_list_dir, sheet = "Th17 down",col_names =F)
Hollbacher_Th17down_genes <- Th17_down$...1

gene_list_name <- "Hollbacher_Th22"
Hollbacher_Th22up_genes <- c('ACAP3', 'ACER3', 'ACTN4', 'ADAM8', 'ADAMTSL4', 'ADCY4', 'ANKLE1', 'ANKRD28', 'ANKRD32', 'ARHGEF4', 'ATG4A', 'ATP10D', 'AZIN2', 'BCL7A', 'BMPR2', 'BRSK1', 'BTBD7', 'C10orf67', 'C11orf80', 'C12orf60', 'C15orf53', 'C20orf194', 'C9orf89', 'CACNB1', 'CACNB3', 'CAMK2G', 'CAMSAP2', 'CAMTA1', 'CAPN2', 'CCDC88A', 'CCR10', 'CD101', 'CD109', 'CD63', 'CD9', 'CDC42BPB', 'CDS1', 'CEL', 'CFAP45', 'CIART', 'CLU', 'CMIP', 'CNTNAP1', 'CORO2A', 'CTBP2', 'DAGLA', 'DAPK2', 'DBP', 'DDIT4', 'DGCR8', 'DIAPH1', 'DIP2C', 'DLGAP4', 'DUSP7', 'ELL', 'EML6', 'ETS2', 'FAM102B', 'FAM65C', 'FGF18', 'FHL3', 'FHOD1', 'FLNA', 'FMO5', 'FSCN1', 'GAK', 'GALNT6', 'GARNL3', 'GLB1L2', 'GLG1', 'GSG2', 'GTF2IRD1', 'GYLTL1B', 'HEATR3', 'HIPK2', 'HIST1H3D', 'HOMER2', 'HPGD', 'ITGAE', 'ITPR3', 'KIAA0232', 'KIAA0355', 'KIAA0825', 'KIFC2', 'L1CAM', 'LCP1', 'LGALS1', 'LMNA', 'LMNTD2', 'LMO4', 'LONRF2', 'LOXL1', 'LPAR6', 'LSS', 'MAP4', 'MAP6', 'MARCH9', 'MB21D1', 'MC1R', 'MCTP1', 'MECP2', 'MFHAS1', 'MGLL', 'MICAL1', 'MMRN2', 'MPP6', 'MTM1', 'MTSS1', 'MVB12B', 'MXRA7', 'MYH9', 'MYO18A', 'MYO1C', 'MYO1G', 'NBEAL1', 'NBEAL2', 'NEDD4L', 'OTUD1', 'PAFAH1B1', 'PAG1', 'PAK1', 'PAQR6', 'PDLIM2', 'PDZD2', 'PER1', 'PER3', 'PFKL', 'PHLDA1', 'PI16', 'PICALM', 'PIK3R3', 'PLCB3', 'PLEC', 'PLIN2', 'PLOD3', 'PLP2', 'PLXNA3', 'PLXNB1', 'POGZ', 'PPFIBP1', 'PPP4R1', 'PROSER2', 'PRRC2B', 'PRRT2', 'PTTG1IP', 'RAB11FIP5', 'RABAC1', 'RAD9A', 'RALGPS1', 'RBMS1', 'RCBTB2', 'RECQL', 'REEP4', 'RFX7', 'ROGDI', 'RPGR', 'S100A10', 'S100A4', 'S100A6', 'SCO2', 'SCPEP1', 'SEPT11', 'SGOL1', 'SLC16A2', 'SLC16A6', 'SLC20A1', 'SLC22A3', 'SLC22A5', 'SLC25A30', 'SLC4A8', 'SLC9A6', 'SMARCD3', 'SOS1', 'SOWAHC', 'SPATS2L', 'SPTBN5', 'STARD9', 'STX1A', 'STX4', 'STXBP1', 'SUN1', 'SUSD6', 'SVIL', 'SYCP2', 'TAOK2', 'TCF25', 'TEF', 'TENM4', 'TESK1', 'TMEM173', 'TNFSF10', 'TOB1', 'TPM4', 'TRAM2', 'TRERF1', 'TRIM2', 'TUBG2', 'UBAP1L', 'UBL3', 'UCKL1', 'USP5', 'VAV2', 'VIM', 'VWA7', 'WHSC1', 'ZDHHC11', 'ZFYVE28', 'ZNF438', 'ZYG11A'  )
Th22_down <- read_excel(path = excel_gene_list_dir, sheet = "Th22 down",col_names =F)
Hollbacher_Th22down_genes <- Th22_down$...1

gene_list_name_list <- c("Kumar_2017_Trm", "Crawford_Ex", "Miller_Ex_progenitors", "Miller_Ex_terminal", "Cano_effector", "Hollbacher_Th1", "Hollbacher_Th2", "Hollbacher_Th17", "Hollbacher_Th22")
# gene_list_name_list <- c( "Hollbacher_Th1", "Hollbacher_Th2", "Hollbacher_Th17", "Hollbacher_Th22")



#### blood & tumor, CD4  all cell types ####

input_data <- "/Users/hai/Documents/workstation/bladder_blood_analysis/metadata/CD4_matrix_and_metadata/adata_Blood_and_tumor_CD4_all_cell_types_raw_matrix.rds"
p2 <- read.csv("/Users/hai/Documents/workstation/bladder_blood_analysis/metadata/CD4_matrix_and_metadata/adata_Blood_and_tumor_CD4_all_cell_types_metadata.csv",header = T, row.names = 1)
# fix sample ID mis label issue
# a2 blood --> a6 blood
# Get levels and add "a6"
df <- p2
levels <- levels(df$patient)
levels[length(levels) + 1] <- "a6"
df$patient <- factor(df$patient, levels = levels)
df$patient[df$patient == "a2" & df$tissue == "blood"] <- 'a6'

levels <- levels(df$sample)
levels[length(levels) + 1] <- "atezo6_pre"
levels[length(levels) + 1] <- "atezo6_post"
df$sample <- factor(df$sample, levels = levels)
df$sample[df$sample == "atezo2_pre" ] <- 'atezo6_pre'
df$sample[df$sample == "atezo2_post" ] <- 'atezo6_post'

# if(pat_id %in% c("a2","a6")) {input_raw_data_select[input_raw_data_select$cluster_id==i,]$blood_tumor <- "NA"}
levels <- levels(df$blood_tumor)
levels[length(levels) + 1] <- "NA"
df$blood_tumor <- factor(df$blood_tumor, levels = levels)
df$blood_tumor[df$patient == "a2" ] <- 'NA'
df$blood_tumor[df$patient == "a6" ] <- 'NA'

# refine cell type order of plots
levels(df$cell_calls) <- c(levels(df$cell_calls), "EA") 
df$cell_calls[df$cell_calls == "Activated"]  <- "EA" 

levels(df$cell_calls) <- c(levels(df$cell_calls), "Prolif") 
df$cell_calls[df$cell_calls == "Prolif GZMK+"]  <- "Prolif" 

levels(df$cell_calls) <- c(levels(df$cell_calls), "IFN") 
df$cell_calls[df$cell_calls == "IFN+"]  <- "IFN" 

df$cell_calls <- factor(df$cell_calls, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))


p2 <- df
write.csv(p2, file = "/Users/hai/Documents/workstation/bladder_blood_analysis/metadata/CD4_matrix_and_metadata/adata_Blood_and_tumor_CD4_all_cell_types_metadata_fix_sampleID.csv", row.names = FALSE)
# p2 <- read.csv("/Users/hai/Documents/workstation/bladder_blood_analysis/metadata/CD4_matrix_and_metadata/adata_Blood_and_tumor_CD4_all_cell_types_metadata_fix_sampleID.csv",header = T)



####
MyData <- readRDS(input_data)
MyData[is.na(MyData)] <- 0
dim(MyData)

expr <- data.matrix(MyData)
expr.medians <- apply(expr, 1, median, na.rm=TRUE)
expr.sd <- apply(expr, 1, sd, na.rm=TRUE)
expr.zscores <- (expr - expr.medians) / expr.sd


for (gene_list_name in gene_list_name_list) {
  print(gene_list_name)
  
  #### match gene list name with gene list ####
  if(gene_list_name ==  "Kumar_2017_Trm") {
    pos.genes <- Kumar_2017_Trm_up
    neg.genes <- Kumar_2017_Trm_down
  }
  
  if(gene_list_name ==  "Crawford_Ex") {
    pos.genes <- Crawford_Ex_up
    neg.genes <- Crawford_Ex_down
  }

  if(gene_list_name ==  "Miller_Ex_progenitors") {
    pos.genes <- Miller_Ex_progenitors_genes
  }
  
  if(gene_list_name ==  "Miller_Ex_terminal") {
    pos.genes <- Miller_Ex_terminal_genes
  }
  
  if(gene_list_name ==  "Cano_effector") {
    pos.genes <- Cano_effector_genes
  }
  
  if(gene_list_name ==  "Hollbacher_Th1") {
    pos.genes <- Hollbacher_Th1up_genes
    neg.genes <- Hollbacher_Th1down_genes
  }
  
  if(gene_list_name ==  "Hollbacher_Th2") {
    pos.genes <- Hollbacher_Th2up_genes
    neg.genes <- Hollbacher_Th2down_genes
    
  }
  
  if(gene_list_name ==  "Hollbacher_Th17") {
    pos.genes <- Hollbacher_Th17up_genes
    neg.genes <- Hollbacher_Th17down_genes
    
  }
  
  if(gene_list_name ==  "Hollbacher_Th22") {
    pos.genes <- Hollbacher_Th22up_genes
    neg.genes <- Hollbacher_Th22down_genes
    
  }
  
  #### ####
  
  pos.genes <- lapply(pos.genes, toupper)
  pos.genes.measured <- intersect(rownames(expr), pos.genes)  # 4
  pos.mat <- expr.zscores[unlist(pos.genes.measured),]
  pos.mat <- pos.mat[expr.sd[names(expr.sd) %in% pos.genes.measured] > 0,]
  pos.values <- colMeans(pos.mat)
  p2$sig_score <- pos.values
  print(paste0("In ", gene_list_name, " pos.genes list, ", length(unlist(pos.genes.measured)), " out of ", length(pos.genes), " are in our data."))
  print(unlist(pos.genes.measured))
  
  if (gene_list_name %in% c("Kumar_2017_Trm", "Crawford_Ex", "Hollbacher_Th1", "Hollbacher_Th2", "Hollbacher_Th17", "Hollbacher_Th22")) {
  neg.genes <- lapply(neg.genes, toupper)
  neg.genes.measured <- intersect(rownames(expr), neg.genes)  # 3
  neg.mat <- expr.zscores[unlist(neg.genes.measured),]
  neg.mat <- neg.mat[expr.sd[names(expr.sd) %in% neg.genes.measured] > 0,]
  neg.values <- colMeans(neg.mat) 
  print(paste0("In ", gene_list_name, " neg.genes list, ", length(unlist(neg.genes.measured)), " out of ", length(neg.genes), " are in our data."))
  print(unlist(neg.genes.measured))
  
  p2$sig_score <- pos.values - neg.values
  }
  
  
  # mac.pos.info.ordered
  # mac.neg.info.ordered
  
  plot_input <- p2[,c("tissue", "blood_tumor","cell_calls", "sig_score")]
  # head(plot_input)
  
  table(plot_input[,c("blood_tumor",  "cell_calls","tissue")])

  
  options(repr.plot.width = 12, repr.plot.height = 12, repr.plot.res = 100) 
  
  barplot0 <- ggplot(data = plot_input, aes(x = cell_calls, y = sig_score, fill = cell_calls)) +
    #         geom_bar(stat = "identity", position=position_dodge()) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean") +
    # facet_wrap(~tissue,  ncol=2) +
    labs(x = "Phenotype", y = gene_list_name) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("Data: Blood and Tumor CD4")
  barplot0
  
  
  barplot1 <- ggplot(data = plot_input, aes(x = cell_calls, y = sig_score, fill = tissue)) +
    #         geom_bar(stat = "identity", position=position_dodge()) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean") +
    # facet_wrap(~tissue,  ncol=2) +
    labs(x = "Phenotype", y = gene_list_name) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("Data: Blood and Tumor CD4")
  barplot1
  
  # mean(plot_input[plot_input$tissue == "tumor" & plot_input$blood_tumor == "shared" & plot_input$cell_call == "GZMK+",]$sig_score)
  
  # related GZMB+ --> GZMK+?
  barplot2 <- ggplot(data = plot_input, aes(x =  cell_calls, y = sig_score, fill = blood_tumor)) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
    # facet_wrap(~ blood_tumor,  ncol=2) +
    labs(x = "Phenotype", y = gene_list_name) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("Data: Blood and Tumor CD4")
  barplot2
  # from this plot, we can not get this conclusion
  # for GZMK+, cells from blood direct to tumor, the ext level from low to high;
  
  
  barplot3 <- ggplot(data = plot_input, aes(x = cell_calls, y = sig_score, fill =  blood_tumor)) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
    facet_wrap(~tissue,  ncol=2) +
    labs(x = "Phenotype", y = gene_list_name) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("Data: Blood and Tumor CD4")
  barplot3
  
  barplot3_1 <- ggplot(data = plot_input[plot_input$blood_tumor != "NA",], aes(x = cell_calls, y = sig_score, fill =  blood_tumor)) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
    facet_wrap(~tissue,  ncol=2) +
    labs(x = "Phenotype", y = gene_list_name) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("Data: Blood and Tumor CD4")
  barplot3_1
  
  barplot4 <- ggplot(data = plot_input, aes(x = cell_calls, y = sig_score, fill =  tissue)) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean") + 
    facet_wrap(~  blood_tumor ,  ncol=2) +
    labs(x = "Phenotype", y = gene_list_name) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("Data: Blood and Tumor CD4")
  barplot4
  
  pdf(paste0(gene_list_name, " gene signature barplot panels on data CD4 all cell types comparing blood vs tumor, share vs non-shared.pdf"))
  # print(plot_celltypes)
  print(barplot0)
        print(barplot1)
              print(barplot2)
                    print(barplot3)
                    print(barplot3_1)
                          print(barplot4)
  dev.off()
  write.csv(p2, paste0(gene_list_name, " gene signature barplot panels on data CD4 all cell types comparing blood vs tumor, share vs non-shared.csv"))
  
}
