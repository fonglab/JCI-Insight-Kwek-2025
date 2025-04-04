# code for paper 10.13
# update on a5 cd4 blood atezo sample
# update on a2 blood sample --> a6
# subset on CD4 or CD8
# plot on blood and tumor
# plot on (pre vs post) or (atezo vs soc) node shape
# plot on cell type
# no node size
# plot on patient level


library(ggplot2)
library(ggseqlogo)
library(RecordLinkage)
library(igraph)
library(reshape)
library(RColorBrewer)
library(Matrix)
library(sparklyr)
library(dplyr)
library(ggraph)
library(graphlayouts)
library(snahelper)
library(data.table)
library(gtools)
library(repr)

# Change plot size to 12 x 12
options(repr.plot.width=12, repr.plot.height=12)

sc <- spark_connect(master = "local", version = "2.1.0")
levenshteinDistMatrixSpark_2 <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=length(x),ncol=length(x))
  cos <- as.data.frame(m)
  df_tbl <- copy_to(sc, data.frame(a=x))
  for(i in 1:length(x)) {
    print(i)
    b <- x[i]
    output <- df_tbl %>% mutate(dist = levenshtein(a, b))
    cos[,i] <- collect(select(output,dist))
  }
  return(cos)
}

setwd("/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID")
# output file directory
output_dir <- "/Users/hai/Documents/workstation/bladder_blood_analysis/TCR_analysis/output/cell_level_cluster_level_meta_data_plots_fix_sampleID/"
source("/Users/hai/Documents/workstation/Dendreon_BCR/distance_matrix/code_latest_version/SelectSubsetNodes_function.R")


#### load data #####
raw_data_with_cluster_id <- read.csv("../../Data/Atezo_Bladder_Total_paired_tcr_updated_fix_sampleID.csv")
dim(raw_data_with_cluster_id)


# refine cell type order of plots
levels(raw_data_with_cluster_id$cell_type) <- c(levels(raw_data_with_cluster_id$cell_type), "EA") 
raw_data_with_cluster_id$cell_type[raw_data_with_cluster_id$cell_type == "Activated"]  <- "EA" 

levels(raw_data_with_cluster_id$cell_type) <- c(levels(raw_data_with_cluster_id$cell_type), "Prolif") 
raw_data_with_cluster_id$cell_type[raw_data_with_cluster_id$cell_type == "Prolif GZMK+"]  <- "Prolif" 

levels(raw_data_with_cluster_id$cell_type) <- c(levels(raw_data_with_cluster_id$cell_type), "IFN") 
raw_data_with_cluster_id$cell_type[raw_data_with_cluster_id$cell_type == "IFN+"]  <- "IFN" 

raw_data_with_cluster_id$cell_type <- factor(raw_data_with_cluster_id$cell_type)
# levels(raw_data_with_cluster_id$cell_type) <- c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN+','')
raw_data_with_cluster_id$cell_type <- factor(raw_data_with_cluster_id$cell_type, levels = c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN'))



# cell type in a5, CD4, blood, atezo are missing.
raw_data_a5_CD4_blood_atezo <- raw_data_with_cluster_id[raw_data_with_cluster_id$patient == "a5" & raw_data_with_cluster_id$comp == "CD4" & raw_data_with_cluster_id$tissue == "blood" & raw_data_with_cluster_id$treatment == "atezo",]
dim(raw_data_a5_CD4_blood_atezo)


pat_id_list <- unique(raw_data_with_cluster_id$patient)
for (pat_id in pat_id_list) {
# for (pat_id in c("a5")) {
  print(pat_id)
  # Calculate distance matrix
  input_raw_data <- raw_data_with_cluster_id[raw_data_with_cluster_id$patient == pat_id,]
  dim(input_raw_data)
  input_raw_data <- input_raw_data[,!(names(input_raw_data) %in% "Pop..ID")]
  dim(input_raw_data)
  input_raw_data_unique <- input_raw_data[!duplicated(input_raw_data$cell.barcode), ]
  dim(input_raw_data_unique)
  input_raw_data <- input_raw_data_unique
  dim(input_raw_data)
  write.csv(input_raw_data, file = paste0("./distance_matrix_by_patient/raw_data_", pat_id, ".csv"))
  
  # calculate distance matrix
  distance_matrix_final_input_table_alphaCDR3_AA_unique <- levenshteinDistMatrixSpark_2(input_raw_data$alphaCDR3)
  saveRDS(distance_matrix_final_input_table_alphaCDR3_AA_unique, file =  paste0("./distance_matrix_by_patient/distance_matrix_alphaCDR3_", pat_id, ".rds"))

  distance_matrix_final_input_table_betaCDR3_AA_unique <- levenshteinDistMatrixSpark_2(input_raw_data$betaCDR3)
  saveRDS(distance_matrix_final_input_table_betaCDR3_AA_unique, file =  paste0("./distance_matrix_by_patient/distance_matrix_betaCDR3_", pat_id, ".rds"))
}



# pat_id_list <- sort(unique(raw_data_with_cluster_id$patient))
sink("./Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/soc atezo patient level summary table info_fix_sampleID.txt")
# pat_id_list <- c( 'a5')
pat_id_list <- c('a2', 'a3', 'a4', 'a5', 'a6', 's1', 's2', 's3') # no health control here.
# pat_id_list <- c('a2', 'a6') # no health control here.
comp_list <- c("CD4", "CD8")
tissue_list <- c("blood","tumor") # no normal 
# tissue_list <- c("blood","tumor","normal")
# tumor is in another data set as the plot include soc patients.
for (pat_id in pat_id_list) {
  for (comp_id in comp_list) {
    for (tissue_id in tissue_list) {
      # tissue_id
      # select patient
      # pat_id <- "a2"
      # comp_id <- "CD4"
      # tissue_id <- "tumor"
      # pat_id <- "a5"
      # comp_id <- "CD4"
      # tissue_id <- "blood"
      print(pat_id)
      print(comp_id)
      print(tissue_id)
      
      input_raw_data <- raw_data_with_cluster_id[raw_data_with_cluster_id$patient == pat_id,]
      dim(input_raw_data)
      distance_matrix_final_input_table_alphaCDR3_AA_unique <- readRDS( paste0("./distance_matrix_by_patient/distance_matrix_alphaCDR3_", pat_id, ".rds"))
      distance_matrix_final_input_table_betaCDR3_AA_unique  <- readRDS(paste0("./distance_matrix_by_patient/distance_matrix_betaCDR3_", pat_id, ".rds"))
      
      # calculate TCR cluster and expansion and share by blood and tumor info, and add to cell level
      # input_raw_data_select <- input_raw_data[input_raw_data$tissue %in% c("blood","tumor") & input_raw_data$comp %in% c(comp_id),]
      # exclued cells without cell type
      input_raw_data_select <- input_raw_data[input_raw_data$tissue %in% c("blood","tumor") & input_raw_data$comp %in% c(comp_id) & !(input_raw_data$cell_type %in% c("")) & !is.na(input_raw_data$cell_type),]
      # input_raw_data_select <- input_raw_data[input_raw_data$tissue %in% c("blood","tumor") & input_raw_data$comp %in% c(comp_id),]
      
      
      dim(input_raw_data_select)
      select_index <- rownames(input_raw_data) %in% rownames(input_raw_data_select)
      
      distance_matrix_final_input_table_alphaCDR3_AA_select <- distance_matrix_final_input_table_alphaCDR3_AA_unique[select_index,select_index]
      distance_matrix_final_input_table_betaCDR3_AA_select <- distance_matrix_final_input_table_betaCDR3_AA_unique[select_index,select_index]
      dim(distance_matrix_final_input_table_alphaCDR3_AA_select)
      dim(distance_matrix_final_input_table_betaCDR3_AA_select)
      # input_raw_data <- input_raw_data_select
      
      # for alphaCDR3 and betaCDR3 chain, make distance = 0 as 1
      input_distance_matrix_alphaCDR3 <- distance_matrix_final_input_table_alphaCDR3_AA_select
      input_distance_matrix_alphaCDR3[input_distance_matrix_alphaCDR3 == 0] <- 9999
      input_distance_matrix_alphaCDR3[input_distance_matrix_alphaCDR3 != 9999] <- 0
      input_distance_matrix_alphaCDR3[input_distance_matrix_alphaCDR3 == 9999] <- 1
      input_distance_matrix_betaCDR3 <- distance_matrix_final_input_table_betaCDR3_AA_select
      input_distance_matrix_betaCDR3[input_distance_matrix_betaCDR3 == 0] <- 9999
      input_distance_matrix_betaCDR3[input_distance_matrix_betaCDR3 != 9999] <- 0
      input_distance_matrix_betaCDR3[input_distance_matrix_betaCDR3 == 9999] <- 1
      
      
      # for combined distance matrix, exclude the case only one of alphaCDR3 / betaCDR3 chain has one distance, the other one have big distance.
      input_distance_matrix <- input_distance_matrix_alphaCDR3 + input_distance_matrix_betaCDR3
      input_distance_matrix[input_distance_matrix == 1] <- 0
      input_distance_matrix[input_distance_matrix == 2] <- 1
      
      
      input_distance_matrix <- as.data.frame(as.matrix(input_distance_matrix))       
      dim(input_distance_matrix)
      matrix_2_net <- as.matrix(input_distance_matrix)
      net <- graph_from_adjacency_matrix(matrix_2_net)
      net <- as.undirected(simplify(net, remove.multiple = T, remove.loops = T))
      
      # add cluster_id in raw data  
      cfg <- cluster_fast_greedy(as.undirected(net))
      input_raw_data_select$cluster_id <- cfg$membership
      
      # add shared by blood and tumor info 
      input_raw_data_select$blood_tumor <- ""
      
      for (i in unique(input_raw_data_select$cluster_id)) {
        index <- ("blood" %in%  unique(input_raw_data_select[input_raw_data_select$cluster_id==i,]$tissue)) &
          ("tumor" %in%  unique(input_raw_data_select[input_raw_data_select$cluster_id==i,]$tissue))
        
        input_raw_data_select[input_raw_data_select$cluster_id==i,]$blood_tumor <- ifelse(index, "shared", "non-shared")
        if(pat_id %in% c("a2","a6")) {input_raw_data_select[input_raw_data_select$cluster_id==i,]$blood_tumor <- "NA"}
      }
      
      print(table(input_raw_data_select$blood_tumor))
      
      # # add expansion info
      input_raw_data_select$expansion <- ""
      deg <- degree(net)
      index_deg <- (deg > 0)
      input_raw_data_select$expansion <- ifelse(index_deg, "expanded", "non-expanded")
      table(input_raw_data_select$expansion)
      
      
      # index_deg3 <- (deg > 3)
      # input_raw_data_select$TCR_size_lt3 <- ifelse(index_deg3, "TCR_clustr_lt3", "TCR_clustr_sm")
      # table(input_raw_data_select$TCR_size_lt3)
      # table(input_raw_data_select[input_raw_data_select$TCR_size_lt3 == "TCR_clustr_lt3" & input_raw_data_select$tissue == "blood",]$cluster_id)
      # 14 15 17 19 
      # 8  8  5  5
      
      
      
      # subset data for plot
      # select tumor and post time point sample only
      dim(input_raw_data_select)
      tmp <- input_raw_data_select
      input_raw_data_select <- tmp[!duplicated(tmp$cell.barcode), ]
      dim(input_raw_data_select)
      # input_raw_data_select_post_tumor <- input_raw_data_select[input_raw_data_select$tissue %in% c( tissue_id ) & input_raw_data_select$comp %in% c(comp_id) & input_raw_data_select$treatment %in% c("atezo","pre"),]
      input_raw_data_select_post_tumor <- input_raw_data_select[input_raw_data_select$tissue %in% c( tissue_id ) & input_raw_data_select$comp %in% c(comp_id) & input_raw_data_select$treatment %in% c("atezo","soc","pre"),]
      
      dim(input_raw_data_select_post_tumor)
      post_tumor_index <- rownames(input_raw_data) %in% rownames(input_raw_data_select_post_tumor)
      
      # distance_matrix_final_input_table_alphaCDR3_AA_unique <- distance_matrix_final_input_table_alphaCDR3_AA_unique[select_index,select_index]
      # distance_matrix_final_input_table_betaCDR3_AA_unique <- distance_matrix_final_input_table_betaCDR3_AA_unique[select_index,select_index]
      
      distance_matrix_final_input_table_alphaCDR3_AA_post_tumor <- distance_matrix_final_input_table_alphaCDR3_AA_unique[post_tumor_index,post_tumor_index]
      distance_matrix_final_input_table_betaCDR3_AA_post_tumor <- distance_matrix_final_input_table_betaCDR3_AA_unique[post_tumor_index,post_tumor_index]
      dim(distance_matrix_final_input_table_alphaCDR3_AA_post_tumor)
      dim(distance_matrix_final_input_table_betaCDR3_AA_post_tumor)
      input_raw_data_sel <- input_raw_data_select_post_tumor
      
      # for alphaCDR3 and betaCDR3 chain, make distance = 0 as 1
      input_distance_matrix_alphaCDR3 <- distance_matrix_final_input_table_alphaCDR3_AA_post_tumor
      input_distance_matrix_alphaCDR3[input_distance_matrix_alphaCDR3 == 0] <- 9999
      input_distance_matrix_alphaCDR3[input_distance_matrix_alphaCDR3 != 9999] <- 0
      input_distance_matrix_alphaCDR3[input_distance_matrix_alphaCDR3 == 9999] <- 1
      input_distance_matrix_betaCDR3 <- distance_matrix_final_input_table_betaCDR3_AA_post_tumor
      input_distance_matrix_betaCDR3[input_distance_matrix_betaCDR3 == 0] <- 9999
      input_distance_matrix_betaCDR3[input_distance_matrix_betaCDR3 != 9999] <- 0
      input_distance_matrix_betaCDR3[input_distance_matrix_betaCDR3 == 9999] <- 1
      
      # for combined distance matrix, exclude the case only one of alphaCDR3 / betaCDR3 chain has one distance, the other one have big distance.
      input_distance_matrix <- input_distance_matrix_alphaCDR3 + input_distance_matrix_betaCDR3
      input_distance_matrix[input_distance_matrix == 1] <- 0
      input_distance_matrix[input_distance_matrix == 2] <- 1
      
      
      input_distance_matrix <- as.data.frame(as.matrix(input_distance_matrix))       
      dim(input_distance_matrix)
      matrix_2_net <- as.matrix(input_distance_matrix)
      net <- graph_from_adjacency_matrix(matrix_2_net)
      net <- as.undirected(simplify(net, remove.multiple = T, remove.loops = T))
      
      #### subset network by degree > 0 ####
      deg <- degree(net)
      print(table(deg))
      # print(pat_id)
      # deg > 0
      dim(input_distance_matrix)
      #   sub_select_matrix <- input_distance_matrix[deg > 0,deg > 0]
      sub_select_matrix <- input_distance_matrix
      dim(sub_select_matrix)
      
      #   sub_input_raw_data <- input_raw_data_sel[deg > 0,]
      sub_input_raw_data <- input_raw_data_sel
      dim(sub_input_raw_data)
      
      # convert the matrix into igraph format
      net <- graph_from_adjacency_matrix(as.matrix(sub_select_matrix), weighted=TRUE)
      net <- as.undirected(simplify(net, remove.multiple = T, remove.loops = T)) 
      
      # add cluster_id in raw data  
      if(nrow(sub_input_raw_data)==0) next
      cfg <- cluster_fast_greedy(as.undirected(net))
      sub_input_raw_data$cluster_id <- cfg$membership
      
      # # refine cell type order of plots
      # levels(sub_input_raw_data$cell_type) <- c(levels(sub_input_raw_data$cell_type), "EA") 
      # sub_input_raw_data$cell_type[sub_input_raw_data$cell_type == "Activated"]  <- "EA" 
      # 
      # levels(sub_input_raw_data$cell_type) <- c(levels(sub_input_raw_data$cell_type), "Prolif") 
      # sub_input_raw_data$cell_type[sub_input_raw_data$cell_type == "Prolif GZMK+"]  <- "Prolif" 
      # 
      # sub_input_raw_data$cell_type <- factor(sub_input_raw_data$cell_type)
      # levels(sub_input_raw_data$cell_type) <- c('Naive', 'GZMB+', 'Mito', 'CM', 'GZMK+', 'CXCL13+', 'Tregs', 'Prolif', 'MAIT', 'EA', 'IFN+')
      # 
      # add expansion info
      # sub_input_raw_data$expansion <- ""
      # deg <- degree(net)
      # index_deg <- (deg > 0)
      # sub_input_raw_data$expansion <- ifelse(index_deg, "expanded", "non-expanded")
      # print(table(sub_input_raw_data$expansion))
      
      #### plot with color on different feature ####
      pdf(paste0("./social_network_plot_on_each_patient_cell_level/social network plot on each patient cell level, subset data on ", comp_id, " patient ID ",pat_id, " tissue ", tissue_id, "_fix_sampleID.pdf"))
      ### best pick layout
      set.seed(9999)
      l <- layout_components(net)
      
      # match figure a on expand vs non-expand
      expansion <- as.factor(sub_input_raw_data$expansion)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = expansion),size=0.5) + 
              # geom_node_point(aes(color = cell_type,size = cloneCount)) + scale_size(range = c(0.1,log(max(cloneCount))/2.5)) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      # color on treatment
      treatment <- as.factor(sub_input_raw_data$treatment)
      print(table(sub_input_raw_data$treatment))
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = treatment),size=0.5) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      # note: the big clusters are all pre, Not any atezo treatment effect.
      
      # match figure a on shared by blood and tumor
      blood_tumor <- as.factor(sub_input_raw_data$blood_tumor)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = blood_tumor),size=0.5) + 
              # geom_node_point(aes(color = cell_type,size = cloneCount)) + scale_size(range = c(0.1,log(max(cloneCount))/2.5)) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      # color on cell_type
      # if (pat_id %in% c("h1","h2","h3")) {next}
      # cloneCount <- sub_input_raw_data$cloneCount
      cell_type <- as.factor(sub_input_raw_data$cell_type)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = cell_type),size=0.5) + 
              # geom_node_point(aes(color = cell_type,size = cloneCount)) + scale_size(range = c(0.1,log(max(cloneCount))/2.5)) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      cluster_id_table <- as.data.frame(table(sub_input_raw_data$cluster_id))
      colnames(cluster_id_table) <- c("cluster_id","Freq")
      head(cluster_id_table[order(-cluster_id_table$Freq),])
      
      # cell count of all cell types among select data
      print(table(sub_input_raw_data$cell_type))
      
      # cell count of the largest cluster
      # get the membership id from the output of "# find the membership of big cluster"  (upper two code sections)
      print(table(sub_input_raw_data[sub_input_raw_data$cluster_id == 1,]$cell_type))
      # Note: this one cluster has 1/4 of the GZMB+ cells.
      # Note: there is more than one cell type in a TCR cluster
      
      #### color on GZMB+ ####
      sub_input_raw_data$GZMB[sub_input_raw_data$cell_type == "GZMB+"] <- 1
      sub_input_raw_data$GZMB[sub_input_raw_data$cell_type != "GZMB+"] <- 0
      GZMB <- as.factor(sub_input_raw_data$GZMB)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = GZMB),size=0.5) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      ### color on GZMK+ ####
      sub_input_raw_data$GZMK[sub_input_raw_data$cell_type == "GZMK+"] <- 1
      sub_input_raw_data$GZMK[sub_input_raw_data$cell_type != "GZMK+"] <- 0
      GZMK <- as.factor(sub_input_raw_data$GZMK)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = GZMK),size=0.5) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      ### color on Prolif GZMK+ ####
      #     sub_input_raw_data$Prolif_GZMK[sub_input_raw_data$cell_type == "Prolif GZMK+"] <- 1
      #     sub_input_raw_data$Prolif_GZMK[sub_input_raw_data$cell_type != "Prolif GZMK+"] <- 0
      #     Prolif_GZMK <- as.factor(sub_input_raw_data$Prolif_GZMK)
      #     print(ggraph(net,layout = l)+
      #             geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
      #             geom_node_point(aes(color = Prolif_GZMK)) +
      #             theme_graph(base_family="sans") +
      #             ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      
      ### color on Prolif GZMK+, GZMK+, GZMB+ ####
      sub_input_raw_data$GZM[sub_input_raw_data$cell_type %in% c("Prolif GZMK+",  "GZMK+", "GZMB+")] <- 1
      sub_input_raw_data$GZM[!(sub_input_raw_data$cell_type %in% c("Prolif GZMK+",  "GZMK+", "GZMB+"))] <- 0
      GZM <- as.factor(sub_input_raw_data$GZM)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = GZM),size=0.5) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      # Note: the three largest cluster are mostly from GZM cell types.
      
      # make sure the data is from blood
      tissue <- as.factor(sub_input_raw_data$tissue)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              # geom_node_point(aes(color = tissue,size = gene_count)) +
              geom_node_point(aes(color = tissue),size=0.5) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      ### color on Tregs ####
      sub_input_raw_data$Tregs_index[sub_input_raw_data$cell_type == "Tregs"] <- 1
      sub_input_raw_data$Tregs_index[sub_input_raw_data$cell_type != "Tregs"] <- 0
      Tregs <- as.factor(sub_input_raw_data$Tregs_index)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = Tregs),size=0.5) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      
      
      # make sure data is from CD4
      # color on comp
      comp <- as.factor(sub_input_raw_data$comp)
      print(ggraph(net,layout = l)+
              geom_edge_link0(width=0.1,colour="grey",alpha =0.8 )+
              geom_node_point(aes(color = comp),size=0.5) +
              theme_graph(base_family="sans") +
              ggtitle(paste0('Alpha and Beta chain, Patient: ', pat_id, " ", comp_id," ", tissue_id)) )
      dev.off()
      
      
      
      #### add statistics info ####
      #### add network statistics ####
      # input data: sub_input_raw_data
      print(dim(sub_input_raw_data))

      sub_input_raw_data$betaCDR3_length <- apply(sub_input_raw_data[,2,drop=FALSE],2,nchar)[,1]
      sub_input_raw_data$alphaCDR3_length <- apply(sub_input_raw_data[,3,drop=FALSE],2,nchar)[,1]

      sub_input_raw_data$deg <- degree(net)
      sub_input_raw_data$transitivity <- transitivity(net, type="local")
      sub_input_raw_data$closeness <- closeness(net, mode="all", weights=NA)
      sub_input_raw_data$centr_clo_res <- centr_clo(net, mode="all", normalized=T)$res
      sub_input_raw_data$eigen_centrality <- eigen_centrality(net, directed=T, weights=NA)$vector
      sub_input_raw_data$centr_eigen <- centr_eigen(net, directed=T, normalized=T)$vector
      sub_input_raw_data$betweenness <- betweenness(net, directed=T, weights=NA)
      sub_input_raw_data$centr_betw <- centr_betw(net, directed=T, normalized=T)$res
      sub_input_raw_data$authority_score <- authority_score(net, weights=NA)$vector
      sub_input_raw_data$coreness <- coreness(net, mode="all")
      sub_input_raw_data$page_rank <- page_rank(net)$vector

      cfg <- cluster_fast_greedy(as.undirected(net))
      sub_input_raw_data$membership_stat <- cfg$membership
      dir.create("Bladder_TCR_match_paper_cell_level_w_all_cluster_stats", showWarnings = FALSE)
      write.csv(sub_input_raw_data, file = paste0( "./Bladder_TCR_match_paper_cell_level_w_all_cluster_stats/Bladder_TCR_match_paper_cell_level_w_all_cluster_stats_","subset_data_on_", comp_id, "_patient_ID_",pat_id, "_tissue_", tissue_id, "_fix_sampleID.csv"), row.names=FALSE)





      #### create membership table for each patient ####
      dim(sub_input_raw_data)
      table(sub_input_raw_data$membership_stat)
      membership_table <- as.data.frame(table(sub_input_raw_data$membership_stat))
      colnames(membership_table) <- c("membership","node_count")



      total_membership <- length(membership_table$membership)
      membership_table$motif_top_deg_alpha <- ""
      membership_table$motif_top_deg_beta <- ""
      membership_table$motif_freq_50_alpha <- ""
      membership_table$motif_freq_50_beta <- ""
      membership_table$deg <- 0
      membership_table$betaCDR3_length <- 0
      membership_table$alphaCDR3_length <- 0

      membership_table$count <- 0
      # membership_table$Count_PRE_INFUSION <- 0
      # membership_table$Count_DOSE_2 <- 0
      # social network properites
      membership_table$deg_avg <- 0
      membership_table$diam_length <- 0
      membership_table$assortativity <- 0
      membership_table$transitivity <- 0
      membership_table$edge_density <- 0
      membership_table$centr_degree <- 0
      membership_table$centr_clo <- 0
      membership_table$eigen_centrality <- 0
      membership_table$centr_eigen <- 0


      for(membership_id in 1:total_membership) {
        # print(paste0(membership_id, " out of ", total_membership))
        # membership_id <- 1
        membership_table[membership_table$membership == membership_id,]$deg  <- round(mean(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$deg),2)
        membership_table[membership_table$membership == membership_id,]$betaCDR3_length  <- round(mean(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$betaCDR3_length),2)
        membership_table[membership_table$membership == membership_id,]$alphaCDR3_length  <- round(mean(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$alphaCDR3_length),2)

        membership_table[membership_table$membership == membership_id,]$count  <- sum(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$count)
        # membership_table[membership_table$membership == membership_id,]$Count_PRE_INFUSION  <- sum(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$Count_PRE_INFUSION)
        # membership_table[membership_table$membership == membership_id,]$Count_DOSE_2  <- sum(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$Count_DOSE_2)
        max_deg <- max(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]$deg)
        membership_table[membership_table$membership == membership_id,]$motif_top_deg_alpha  <- as.character(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$deg == max_deg,]$alphaCDR3[1])
        membership_table[membership_table$membership == membership_id,]$motif_top_deg_beta  <- as.character(sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$deg == max_deg,]$betaCDR3[1])

        # get the representative motif
        betaCDR3_length_round <- round(membership_table[membership_table$membership == membership_id,]$betaCDR3_length)
        for(i in 1:betaCDR3_length_round) {
          # print(i)
          string_list <- sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$betaCDR3_length == betaCDR3_length_round,]$betaCDR3
          freq_table <- as.data.frame(table(substring(string_list, i,i)))
          select_letter <- ifelse(dim(freq_table[freq_table$Freq > length(string_list)/2,])[1] == 1, as.character(freq_table[freq_table$Freq > length(string_list)/2,]$Var1), "*")
          # print(select_letter)
          select_letter_index <- paste0("char_", i)
          assign(select_letter_index, select_letter)
        }
        membership_table[membership_table$membership == membership_id,]$motif_freq_50_beta  <- paste(noquote(mget(mixedsort(ls(pattern= "char_")))), collapse = "")
        rm(list = ls(pattern= "char_"))

        alphaCDR3_length_round <- round(membership_table[membership_table$membership == membership_id,]$alphaCDR3_length)
        for(i in 1:alphaCDR3_length_round) {
          # print(i)
          string_list <- sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id & sub_input_raw_data$alphaCDR3_length == alphaCDR3_length_round,]$alphaCDR3
          freq_table <- as.data.frame(table(substring(string_list, i,i)))
          select_letter <- ifelse(dim(freq_table[freq_table$Freq > length(string_list)/2,])[1] == 1, as.character(freq_table[freq_table$Freq > length(string_list)/2,]$Var1), "*")
          # print(select_letter)
          select_letter_index <- paste0("char_", i)
          assign(select_letter_index, select_letter)
        }
        membership_table[membership_table$membership == membership_id,]$motif_freq_50_alpha  <- paste(noquote(mget(mixedsort(ls(pattern= "char_")))), collapse = "")
        rm(list = ls(pattern= "char_"))



        ### network properties ###
        input_matrix <- sub_select_matrix
        input_data_membership <- sub_input_raw_data[sub_input_raw_data$membership_stat == membership_id,]
        input_data_membership_index <- rownames(sub_input_raw_data) %in% rownames(input_data_membership)

        input_matrix_membership <- input_matrix[input_data_membership_index,input_data_membership_index]
        input_matrix_membership[input_matrix_membership > 1] <- 0

        matrix_2_net <- as.matrix(input_matrix_membership)
        net <- graph_from_adjacency_matrix(matrix_2_net)
        net <- as.undirected(simplify(net, remove.multiple = T, remove.loops = T))

        deg <- degree(net, mode="all")
        # table(deg) # there should be no deg== 0.
        membership_table[membership_table$membership == membership_id,]$deg_avg <- round(mean(deg),2)

        # Diameter (longest geodesic distance)
        diam <- get_diameter(net, directed=T)
        # diam
        membership_table[membership_table$membership == membership_id,]$diam_length <- length(diam)

        # Assortativity
        membership_table[membership_table$membership == membership_id,]$assortativity <- assortativity_degree(net, directed=F)

        # Transitivity
        membership_table[membership_table$membership == membership_id,]$transitivity <- transitivity(net, type="global")  # net is treated as an undirected network

        # Density
        # The proportion of present edges from all possible ties.
        membership_table[membership_table$membership == membership_id,]$edge_density <- edge_density(net, loops=F)

        # centralization on degree
        membership_table[membership_table$membership == membership_id,]$centr_degree <- centr_degree(net, mode="in", normalized=T)$centralization

        # centralization on Closeness (centrality based on distance to others in the graph)
        membership_table[membership_table$membership == membership_id,]$centr_clo <- centr_clo(net, mode="all", normalized=T)$centralization

        # centralization on Eigenvector (centrality proportional to the sum of connection centralities)
        # Values of the first eigenvector of the graph adjacency matrix
        membership_table[membership_table$membership == membership_id,]$eigen_centrality <- eigen_centrality(net, directed=T, weights=NA)$value
        membership_table[membership_table$membership == membership_id,]$centr_eigen<- centr_eigen(net, directed=T, normalized=T)$centralization
      }
        head(membership_table)  
        dir.create("cluster_level_info", showWarnings = FALSE)
        write.csv(membership_table, file = paste0("./cluster_level_info/cluster_level_info_match_paper_","subset_data_on_", comp_id, "_patient_ID_",pat_id, "_tissue_", tissue_id, ".csv"), row.names=FALSE)

      
    }   
  }
}
sink()




