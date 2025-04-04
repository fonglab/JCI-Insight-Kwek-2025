

# CD4
library(ggplot2)
library(reshape2)
library(readxl)

# load data
input_dir <- "/Users/hai/Documents/workstation/bladder_blood_analysis/scRNA/output_for_paper/A_gene_signatures/"

comp_list <- c("CD4", "CD8")
gene_list_name_list <- c( "Hollbacher_Th1", "Hollbacher_Th2", "Kumar_2017_Trm")
for (comp_index in comp_list) {
    for(gene_list_name in gene_list_name_list) {

# comp_index <- "CD4"
# gene_list_name <- "Hollbacher_Th1"

p2 <- read.csv(paste0(input_dir, gene_list_name, " gene signature barplot panels on data ", comp_index, " all cell types comparing blood vs tumor, share vs non-shared.csv"  ))
p2$cell_calls <- factor(p2$cell_calls, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))

plot_input <- p2[,c("tissue", "blood_tumor","cell_calls", "sig_score")]

  barplot0 <- ggplot(data = plot_input, aes(x = cell_calls, y = sig_score, fill = cell_calls)) +
    #         geom_bar(stat = "identity", position=position_dodge()) +
    geom_bar(position = "dodge", stat = "summary", fun = "mean") +
    # facet_wrap(~tissue,  ncol=2) +
    labs(x = "Phenotype", y = gene_list_name) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle("Data: Blood and Tumor CD8")
  barplot0

# get mean of the sig_score for each cell_calls
plot_input_mean <- aggregate(sig_score ~ cell_calls, data = plot_input, FUN = mean)
print(paste0(comp_index, " ",gene_list_name))
print(plot_input_mean)

    }
}
