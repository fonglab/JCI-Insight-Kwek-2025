# bar plot for paper
library(ggplot2)
library(stringr)
library(dplyr)

setwd("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID")

# input_data <- read.csv("./Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_subset_data_merged_w_soc_fix_sampleID.csv")
input_data <- read.csv("./Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_subset_data_merged_w_soc_expansion1_fix_sampleID.csv")

# reorder TCR_size_gt1 factor
# input_data$cell_type<-factor(input_data$cell_type, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN+'))
input_data$cell_type<-factor(input_data$cell_type, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))

# add shared into the TCR_size_gt1
levels(input_data$TCR_size_gt1) <- c(levels(input_data$TCR_size_gt1), "shared_cluster") 
input_data$TCR_size_gt1[input_data$blood_tumor == "shared" & input_data$TCR_size_gt1 == "TCR_cluster_gt1"] <- "shared_cluster"
table(input_data$TCR_size_gt1)

input_data$TCR_size_gt1 <- factor(input_data$TCR_size_gt1, levels = c("shared_cluster","TCR_cluster_gt1", "one_cell_cluster"))


CD_list <- c("CD4", "CD8")
tissue_list <- c("blood", "tumor")

pdf("Barplot on percentage_expansion1 add shared cluster.pdf")
for(CD_index in CD_list) {
  gg_data <- input_data[input_data$comp == CD_index,]
  melt_dat <- gg_data %>% 
    group_by(cell_type, TCR_size_gt1) %>% 
    summarise(count = n()) %>% 
    mutate(percent = count/sum(count)*100)
  
  p <-   ggplot(melt_dat, aes(x = cell_type, y = percent, fill = TCR_size_gt1)) +
    geom_bar(stat="identity", width = 0.7) +
    # labs(x = "Groupchange", y = "percent", fill = "Symscore") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = CD_index)
  print(p)
}
  

  
  
  for(tissue_index in tissue_list) {
    gg_data <- input_data[ input_data$tissue == tissue_index, ]
    melt_dat <- gg_data %>% 
      group_by(cell_type, TCR_size_gt1) %>% 
      summarise(count = n()) %>% 
      mutate(percent = count/sum(count)*100)
    
    p <- ggplot(melt_dat, aes(x = cell_type, y = percent, fill = TCR_size_gt1)) +
      geom_bar(stat="identity", width = 0.7) +
      # labs(x = "Groupchange", y = "percent", fill = "Symscore") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title = paste(tissue_index))
    print(p)
  }
  
  
  
for(CD_index in CD_list) {
  for(tissue_index in tissue_list) {
    gg_data <- input_data[input_data$comp == CD_index & input_data$tissue == tissue_index, ]
    melt_dat <- gg_data %>% 
      group_by(cell_type, TCR_size_gt1) %>% 
      summarise(count = n()) %>% 
      mutate(percent = count/sum(count)*100)
    
    p <- ggplot(melt_dat, aes(x = cell_type, y = percent, fill = TCR_size_gt1)) +
      geom_bar(stat="identity", width = 0.7) +
      # labs(x = "Groupchange", y = "percent", fill = "Symscore") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title = paste(CD_index, tissue_index))
    print(p)
    
  }
}

dev.off()






# the number of expanded and non-expanded, shared and non shared for patient a4's CD4 blood, CD8 blood, CD4 tumor, CD8 tumor only

input_data_sel <- input_data[input_data$patient == "Anti-PD-L1 C" & input_data$comp == "CD4" & input_data$tissue == "blood",]
table(input_data_sel$blood_tumor)
# non-shared     shared 
#       5994        203
table(input_data_sel$TCR_size_gt1)
# TCR_cluster_gt1 one_cell_cluster 
#             402             5795 

input_data_sel <- input_data[input_data$patient == "Anti-PD-L1 C" & input_data$comp == "CD8" & input_data$tissue == "blood",]
table(input_data_sel$blood_tumor)
# non-shared     shared 
#       3891        502 
table(input_data_sel$TCR_size_gt1)
# TCR_cluster_gt1 one_cell_cluster 
#            2111             2282 

input_data_sel <- input_data[input_data$patient == "Anti-PD-L1 C" & input_data$comp == "CD4" & input_data$tissue == "tumor",]
table(input_data_sel$blood_tumor)
# non-shared     shared 
#       1696         49
table(input_data_sel$TCR_size_gt1)
# TCR_cluster_gt1 one_cell_cluster 
#             223             1522  

input_data_sel <- input_data[input_data$patient == "Anti-PD-L1 C" & input_data$comp == "CD8" & input_data$tissue == "tumor",]
table(input_data_sel$blood_tumor)
# non-shared     shared 
#       474        100 
table(input_data_sel$TCR_size_gt1)
# TCR_cluster_gt1 one_cell_cluster 
#             210              364 
