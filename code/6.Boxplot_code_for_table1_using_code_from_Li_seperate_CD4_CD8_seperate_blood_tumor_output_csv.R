library(ggplot2)

# data create from code: cluster_glmm11302020_normalizedwsoc_greater1_comparisonbycelltype.R (run line 1-155)

#### ncluster ####
input_data_all <- dt.TCR_cluster.shared.ncluster.bycelltype 
# input_data_all$cell_type <-factor(input_data_all$cell_type, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN+'))
table(input_data_all$cell_type)
# rename the cell types
levels(input_data_all$cell_type) <- c(levels(input_data_all$cell_type), "IFN") 
input_data_all$cell_type[input_data_all$cell_type == "IFN+"]  <- "IFN" 
levels(input_data_all$cell_type) <- c(levels(input_data_all$cell_type), "Prolif") 
input_data_all$cell_type[input_data_all$cell_type == "Prolif GZMK+"]  <- "Prolif" 
levels(input_data_all$cell_type) <- c(levels(input_data_all$cell_type), "EA") 
input_data_all$cell_type[input_data_all$cell_type == "Activated"]  <- "EA" 


input_data_all$cell_type <-factor(input_data_all$cell_type, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))
table(input_data_all$cell_type)

# for( select_tissue in c("blood","tumor")) {
#   input_data <- input_data_all[input_data_all$tissue== select_tissue,]
# 
#   head(input_data)
#   plot_CD4CD8 <- ggplot(input_data,aes(x=cell_type,y=log10x,fill=blood_tumor))+
#     geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
#     facet_wrap(~comp) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     labs(x = "Cell types") +
#     labs(y = "log10(number of clusters)")
# 
#   plot_CD4 <- ggplot(input_data[input_data$comp=="CD4",],aes(x=cell_type,y=log10x,fill=blood_tumor))+
#     geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
#     # facet_wrap(~comp) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     labs(x = "Cell types") +
#     labs(y = "log10(number of clusters)")
# 
#   plot_CD8 <- ggplot(input_data[input_data$comp=="CD8",],aes(x=cell_type,y=log10x,fill=blood_tumor))+
#     geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
#     # facet_wrap(~comp) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     labs(x = "Cell types") +
#     labs(y = "log10(number of clusters)")
# 
#   pdf(paste0("boxplot for Table 1 subset on CD4 CD8 comparing on ncluster shared vs non shared for ", select_tissue," data only.pdf"))
#   print(plot_CD4CD8)
#   print(plot_CD4)
#   print(plot_CD8)
#   dev.off()
# }



# head(input_data)
for( select_comp in c("CD4","CD8")) {
  input_data <- input_data_all[input_data_all$comp== select_comp,]
  
  plot_tissue <- ggplot(input_data,aes(x=cell_type,y=log10x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_comp)) +
    facet_wrap(~tissue) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(number of clusters)")
  
  plot_blood <- ggplot(input_data[input_data$tissue=="blood",],aes(x=cell_type,y=log10x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_comp)) +
    # facet_wrap(~comp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(number of clusters)") +
    scale_x_discrete(drop=FALSE)
  # save the corresponding data to csv file
  write.csv(input_data[input_data$tissue=="blood",], paste0("boxplot for Table 1 subset on blood tumor comparing on ncluster shared vs non shared", "_", select_comp, "_", "blood boxplot.csv"))
  
  plot_tumor <- ggplot(input_data[input_data$tissue=="tumor",],aes(x=cell_type,y=log10x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_comp)) +
    # facet_wrap(~comp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(number of clusters)") +
    scale_x_discrete(drop=FALSE)
  # save the corresponding data to csv file
  write.csv(input_data[input_data$tissue=="blood",], paste0("boxplot for Table 1 subset on blood tumor comparing on ncluster shared vs non shared", "_", select_comp, "_", "tumor boxplot.csv"))
  
  pdf(paste0("boxplot for Table 1 subset on blood tumor comparing on ncluster shared vs non shared for ",select_comp," data only_11.7.2023.pdf"))
  print(plot_tissue)
  print(plot_blood)
  print(plot_tumor)
  dev.off()
}


#### clustersize ####
input_data_all <- dt.TCR_cluster.shared.clustersize.bycelltype
input_data_all$cell_type <-factor(input_data_all$cell_type, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))

for( select_tissue in c("blood","tumor")) {
  input_data <- input_data_all[input_data_all$tissue== select_tissue,]
  head(input_data)
  
  
  plot_CD4CD8 <- ggplot(input_data,aes(x=cell_type,y=x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
    facet_wrap(~comp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(cluster size)")
  
  plot_CD4 <- ggplot(input_data[input_data$comp=="CD4",],aes(x=cell_type,y=x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
    # facet_wrap(~comp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(cluster size)")
  
  plot_CD8 <- ggplot(input_data[input_data$comp=="CD8",],aes(x=cell_type,y=x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
    # facet_wrap(~comp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(cluster size)")
  
  pdf(paste0("boxplot for Table 1 subset on CD4 CD8 comparing on clustersize shared vs non shared for ", select_tissue," data only.pdf"))
  print(plot_CD4CD8)
  print(plot_CD4)
  print(plot_CD8)
  dev.off()
}


for( select_comp in c("CD4","CD8")) {
  input_data <- input_data_all[input_data_all$comp== select_comp,]
  
  head(input_data)
  plot_tissue <- ggplot(input_data,aes(x=cell_type,y=x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
    facet_wrap(~tissue) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(cluster size)")
  
  plot_blood <- ggplot(input_data[input_data$tissue=="blood",],aes(x=cell_type,y=x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
    # facet_wrap(~comp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(cluster size)")
  
  plot_tumor <- ggplot(input_data[input_data$tissue=="tumor",],aes(x=cell_type,y=x,fill=blood_tumor))+
    geom_boxplot()+ ggtitle(paste0("input data: ", select_tissue)) +
    # facet_wrap(~comp) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Cell types") +
    labs(y = "log10(cluster size)")
  
  pdf(paste0("boxplot for Table 1 subset on blood tumor comparing on clustersize shared vs non shared for ",select_comp," data only.pdf"))
  print(plot_tissue)
  print(plot_blood)
  print(plot_tumor)
  dev.off()
}
