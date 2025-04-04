# create a heatmap on percentage of cell types

# data from cell level data 
# get percentage from each subset data
# create a heatmap 

# four plots: CD4 vs CD8; Tumor vs Blood
#  Overall vs expanded clones vs shared clones


setwd("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/")
#### merege data ####
file_list <- list.files(pattern = "Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_subset_data_on_")

merge_file <- NULL
for( file_index in file_list) {
  read_file <- read.csv(file_index)
  print(dim(read_file))
  merge_file <- rbind(merge_file, read_file)
}
dim(merge_file)

# read merge file
merge_file <- read.csv("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/backup/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_subset_data_merged_w_soc_fix_sampleID.csv")
table(merge_file$cell_type)

# create TCR_cluster by (cluster size >3)
merge_file$TCR_cluster_id <- paste0(merge_file$patient, merge_file$comp, merge_file$tissue, merge_file$cluster_id)
TCR_cluster_table <- as.data.frame(table(merge_file$TCR_cluster_id))
TCR_cluster_gt1_list <- TCR_cluster_table[TCR_cluster_table$Freq > 1,]$Var1
merge_file$TCR_size_gt1 <- ifelse(merge_file$TCR_cluster_id %in% TCR_cluster_gt1_list, "TCR_cluster_gt1", "one_cell_cluster")

# refine cell type order of plots
levels(merge_file$cell_type) <- c(levels(merge_file$cell_type), "EA") 
merge_file$cell_type[merge_file$cell_type == "Activated"]  <- "EA" 

levels(merge_file$cell_type) <- c(levels(merge_file$cell_type), "Prolif") 
merge_file$cell_type[merge_file$cell_type == "Prolif GZMK+"]  <- "Prolif" 

levels(merge_file$cell_type) <- c(levels(merge_file$cell_type), "IFN") 
merge_file$cell_type[merge_file$cell_type == "IFN+"]  <- "IFN" 

merge_file$cell_type <- factor(merge_file$cell_type)
# levels(merge_file$cell_type) <- c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN+') # wrong!!!
merge_file$cell_type <- factor(merge_file$cell_type, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))


# refine patient id for plots
levels(merge_file$patient) <- c(levels(merge_file$patient), "Anti-PD-L1 A") 
merge_file$patient[merge_file$patient == "a2"]  <- "Anti-PD-L1 A" 

levels(merge_file$patient) <- c(levels(merge_file$patient), "Anti-PD-L1 B") 
merge_file$patient[merge_file$patient == "a3"]  <- "Anti-PD-L1 B" 

levels(merge_file$patient) <- c(levels(merge_file$patient), "Anti-PD-L1 C") 
merge_file$patient[merge_file$patient == "a4"]  <- "Anti-PD-L1 C" 

levels(merge_file$patient) <- c(levels(merge_file$patient), "Anti-PD-L1 D") 
merge_file$patient[merge_file$patient == "a5"]  <- "Anti-PD-L1 D" 

levels(merge_file$patient) <- c(levels(merge_file$patient), "Anti-PD-L1 E") 
merge_file$patient[merge_file$patient == "a6"]  <- "Anti-PD-L1 E" 

levels(merge_file$patient) <- c(levels(merge_file$patient), "Chemo") 
merge_file$patient[merge_file$patient == "s1"]  <- "Chemo" 

levels(merge_file$patient) <- c(levels(merge_file$patient), "Untreated A") 
merge_file$patient[merge_file$patient == "s2"]  <- "Untreated A" 

levels(merge_file$patient) <- c(levels(merge_file$patient), "Untreated B") 
merge_file$patient[merge_file$patient == "s3"]  <- "Untreated B" 

merge_file$patient <- factor(merge_file$patient)
table(merge_file$patient)
levels(merge_file$patient)


write.csv(merge_file, "../Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_subset_data_merged_w_soc_expansion1_fix_sampleID.csv")



#### heatmap plot ####
# pat_list <- c("a2", "a3", "a4", "a5", "a6", "s1", "s2", "s3")
pat_list <- levels(merge_file$patient)
# subset_group <- c("All", "expanded","shared", "TCR_cluster_gt1", "pre_blood", "atezo_blood")
subset_group <- c("All", "expanded","shared", "TCR_cluster_gt1", "pre_blood_all","pre_blood_expanded","pre_blood_shared", "atezo_blood_all","atezo_blood_expanded","atezo_blood_shared")
CD_list <- c("CD4", "CD8")
tissue_list <- c("blood", "tumor")
# cell_type_table_all <- NULL


# for(subset_index in subset_group) {
#   print(subset_index)
#   if(subset_index == "All") { input_data <- cell_level_data }
#   if(subset_index == "expanded") { input_data <- cell_level_data[cell_level_data$expansion == 'expanded',] }
#   if(subset_index == "shared") { input_data <- cell_level_data[cell_level_data$blood_tumor == 'shared',] }
# }

# CD_index <- "CD8"
# tissue_index <- "blood"
# subset_index <- "shared"
# pat_id <- "a4"
pdf("../heatmap_on_percentage_of_cell_type_w_soc_expansion1.pdf")

for(CD_index in CD_list) {
  for(tissue_index in tissue_list) {
    for(subset_index in subset_group) {
      # for(pat_id in pat_list) {
        # print(pat_id)
        # cell_level_data <- read.csv(paste0("Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_subset_data_on_", CD_index, "_patient_ID_", pat_id,"_tissue_", tissue_index,".csv"))
        cell_level_data <- merge_file[merge_file$comp == CD_index & merge_file$tissue == tissue_index,]
        dim(cell_level_data)
        # input_data <- cell_level_data
        print(subset_index)
        if(subset_index == "All") { input_data <- cell_level_data }
        if(subset_index == "expanded") { input_data <- cell_level_data[cell_level_data$expansion == 'expanded',] }
        if(subset_index == "shared") { input_data <- cell_level_data[cell_level_data$blood_tumor == 'shared',] }
        if(subset_index == "TCR_cluster_gt1") { input_data <- cell_level_data[cell_level_data$TCR_size_gt1 == 'TCR_cluster_gt1',] }
        if(subset_index == "pre_blood_all") { input_data <- cell_level_data[cell_level_data$tissue == 'blood' & cell_level_data$treatment == "pre",] }
        if(subset_index == "pre_blood_expanded") { input_data <- cell_level_data[cell_level_data$tissue == 'blood' & cell_level_data$treatment == "pre" & cell_level_data$expansion == 'expanded',] }
        if(subset_index == "pre_blood_shared") { input_data <- cell_level_data[cell_level_data$tissue == 'blood' & cell_level_data$treatment == "pre" & cell_level_data$blood_tumor == 'shared',] }
        
        if(subset_index == "atezo_blood_all") { input_data <- cell_level_data[cell_level_data$tissue == 'blood' & cell_level_data$treatment == "atezo",] }
        if(subset_index == "atezo_blood_expanded") { input_data <- cell_level_data[cell_level_data$tissue == 'blood' & cell_level_data$treatment == "atezo" & cell_level_data$expansion == 'expanded',] }
        if(subset_index == "atezo_blood_shared") { input_data <- cell_level_data[cell_level_data$tissue == 'blood' & cell_level_data$treatment == "atezo" & cell_level_data$blood_tumor == 'shared',] }
        
        
        # cell_type_table <- as.data.frame(table(input_data$cell_type))
        cell_type_table <- as.data.frame(table(input_data$patient,input_data$cell_type))
        # cell_type_table$percent <- round(cell_type_table$Freq/nrow(input_data)*100,1)
        cell_type_table$percent <- 0
        for(pat_id in pat_list) {
          cell_type_table[cell_type_table$Var1 == pat_id,]$percent <- round(cell_type_table[cell_type_table$Var1 == pat_id,]$Freq/sum(cell_type_table[cell_type_table$Var1 == pat_id,]$Freq)*100,1)
        }
        
        colnames(cell_type_table) <- c("patient", "cell_types",  "Freq", "percent")
        # cell_type_table
        cell_type_table$subset_data <- paste0(CD_index, "_", tissue_index, "_", subset_index)
        write.csv(cell_type_table, paste0("../cell_type_table_pencentage_and_cell_counts/cell_type_table_pencentage_and_cell_counts_",  CD_index, "_", tissue_index, "_", subset_index, "_expansion1.csv"))
        
      

      cell_type_table_pct <- cell_type_table[,c("cell_types", "patient", "percent")]

      # switch row and column
      # library(reshape2)
      # cell_type_table_all_t <- as.data.frame(t(cell_type_table_pct))
      # colnames(cell_type_table_all_t) <- as.character(unlist(cell_type_table_all_t[1,]))
      # cell_type_table_all_t <- cell_type_table_all_t[-1,]
      # colnames(cell_type_table_all_t)[1] <- NA
      # melted_cell_type_table_all <- melt(as.matrix(cell_type_table_all_t))
      # melted_cell_type_table_all$value <- as.numeric(as.character(melted_cell_type_table_all$value))
      
      
      library(ggplot2)
      
      p <- ggplot(cell_type_table_pct, aes( patient,cell_types)) +
        geom_tile(aes(fill = percent)) + 
        geom_text(aes(label = round(percent, 1))) +
        scale_fill_gradient(low = "white", high = "red", na.value = 'white') +
        ggtitle(paste0(subset_index, " ", CD_index, " ", tissue_index)) + 
        labs(x = "Patient", y = "Cell type") +
        theme(text = element_text(size=6), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      print(p)
    }
  }
}
dev.off()



#### create trend plots ####

pdf("../trend_plots_of_pencentage_from_heatmap_by_cell_types_expansion1.pdf")

# select_subset_list <- c("All", "expanded", "TCR_cluster_gt1", "shared")
select_subset_list <- c("All", "TCR_cluster_gt1", "shared")

cell_type_list <- c("GZMB+", "GZMK+", "Prolif GZMK+", "CXCL13+")
CD_list <- c("CD4", "CD8")
tissue_list <- c("blood", "tumor")

# CD_index <- "CD8"
# tissue_index <- "blood"
# cell_type_index <- "GZMB+"

for(CD_index in CD_list) {
  print(CD_index)
  for(tissue_index in tissue_list) {
    print(tissue_index)
    merge_select_data <- NULL
    for(select_subset_index in select_subset_list) {
      print(select_subset_index)
      read_data <- read.csv(paste0("../cell_type_table_pencentage_and_cell_counts/cell_type_table_pencentage_and_cell_counts_",  CD_index, "_", tissue_index, "_", select_subset_index, "_expansion1.csv"))
      print(dim(read_data))
      merge_select_data <- rbind(merge_select_data, read_data)
      print(dim(merge_select_data))
    }
    write.csv(merge_select_data, paste0("../cell_type_table_pencentage_and_cell_counts/merge_cell_type_table_pencentage_and_cell_counts_",  CD_index, "_", tissue_index, "_expansion1.csv"))
    
    
    # x_axis_label <- c("All", "Expanded", "TCR cluster > 3", "Shared")
    x_axis_label <- c("All", "TCR cluster > 1", "Shared")
    
    for(cell_type_index in cell_type_list ) {
      subtable <- merge_select_data[merge_select_data$cell_types == cell_type_index,]
      trend_plot <- ggplot(subtable, aes(x = subset_data, y = percent, group = patient)) + 
        geom_point(aes(color = patient)) + 
        geom_line(aes(color = patient)) +
        ggtitle(paste0(CD_index, "_", tissue_index, "_", cell_type_index)) +
        scale_x_discrete(labels= x_axis_label) +
        theme(panel.background = element_blank())
      print(trend_plot)
    }
  }
}
dev.off()
