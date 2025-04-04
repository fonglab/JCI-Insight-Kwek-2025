# box plot for paper
library(ggplot2)
library(stringr)

setwd("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID")

input_data1 <- read.csv("./cell_type_table_pencentage_and_cell_counts/merge_cell_type_table_pencentage_and_cell_counts_CD4_blood_expansion1.csv")
input_data2 <- read.csv("./cell_type_table_pencentage_and_cell_counts/merge_cell_type_table_pencentage_and_cell_counts_CD8_tumor_expansion1.csv")
input_data3 <- read.csv("./cell_type_table_pencentage_and_cell_counts/merge_cell_type_table_pencentage_and_cell_counts_CD4_tumor_expansion1.csv")
input_data4 <- read.csv("./cell_type_table_pencentage_and_cell_counts/merge_cell_type_table_pencentage_and_cell_counts_CD8_blood_expansion1.csv")

input_data <- rbind(input_data1,input_data2, input_data3, input_data4)
input_data$comp <- str_split_fixed(input_data$subset_data, "_", 3)[,1]
input_data$tissue <- str_split_fixed(input_data$subset_data, "_", 3)[,2]
input_data$sub_info <- str_split_fixed(input_data$subset_data, "_", 3)[,3]

levels(input_data$cell_types) <- c(levels(input_data$cell_types), "IFN") 
input_data$cell_types[input_data$cell_types == "IFN+"]  <- "IFN" 

head(input_data)
summary(input_data$subset_data)
input_data$sub_info <- factor(input_data$sub_info, levels = c("All", "expanded", "TCR_cluster_gt1", "shared"))
input_data$cell_types <-factor(input_data$cell_types, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))

CD_list <- c("CD4", "CD8")
tissue_list <- c("blood", "tumor")
# subset_group <- c("All", "expanded","shared", "TCR_cluster_gt1")
subset_group <- c("All","shared", "TCR_cluster_gt1")


# pdf("Boxplot on percentage on comparsion_expansion1 reorder cell types.pdf")
# for(CD_index in CD_list) {
#   gg_data <- input_data[input_data$comp == CD_index,]
#   p <- ggplot(gg_data, aes(x=cell_types, y=percent, fill=tissue)) + 
#     geom_boxplot() +
#     facet_wrap(vars(sub_info)) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     labs(title = CD_index) +
#     theme(panel.background = element_blank())
#   
#   print(p)
#   
#   p <- ggplot(gg_data, aes(x=cell_types, y=percent, fill=sub_info )) + 
#     geom_boxplot() +
#     facet_wrap(vars(tissue)) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     labs(title = CD_index) +
#     theme(panel.background = element_blank())
#   print(p)
#   
#   
#   for(tissue_index in tissue_list) {
#     gg_data <- input_data[input_data$comp == CD_index & input_data$tissue == tissue_index, ]
#     p <- ggplot(gg_data, aes(x=cell_types, y=percent, fill=sub_info)) + 
#       geom_boxplot() +
#       # facet_wrap(vars(sub_info)) +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       labs(title = paste(CD_index, tissue_index)) +
#       theme(panel.background = element_blank())
#     print(p)
#   }
#   
#   
#   for(subset_index in subset_group) {
#     gg_data <- input_data[input_data$comp == CD_index  & input_data$sub_info == subset_index, ]
#     p <- ggplot(gg_data, aes(x=cell_types, y=percent, fill=tissue)) + 
#       geom_boxplot() +
#       # facet_wrap(vars(sub_info)) +
#       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       labs(title = paste(CD_index, subset_index)) +
#       theme(panel.background = element_blank())
#     print(p)
#   }
#   
# }
# 
# dev.off()

pdf(paste0("./heatmap_trend_Barplot_Boxplot/Boxplot_on_percentage_on_subset_data_expansion1_reorder_cell_types_output_csv/","Boxplot on percentage on subset data_expansion1 reorder cell types_11.7.2023.pdf"))
for(CD_index in CD_list) {
  for(tissue_index in tissue_list) {
    for(subset_index in subset_group) {
      gg_data <- input_data[input_data$comp == CD_index & input_data$tissue == tissue_index & input_data$sub_info == subset_index, ]
      p <- ggplot(gg_data, aes(x=cell_types, y=percent)) + 
        geom_boxplot() +
        # facet_wrap(vars(sub_info)) +
        ylim(0,100) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(title = paste(CD_index, tissue_index, subset_index)) +
        theme(panel.background = element_blank())
      print(p)
      # save the corresponding data to csv file
      write.csv(gg_data, paste0("./heatmap_trend_Barplot_Boxplot/Boxplot_on_percentage_on_subset_data_expansion1_reorder_cell_types_output_csv/",CD_index, "_", tissue_index, "_", subset_index, "_", "boxplot.csv"))
    }
  }
}
dev.off()




