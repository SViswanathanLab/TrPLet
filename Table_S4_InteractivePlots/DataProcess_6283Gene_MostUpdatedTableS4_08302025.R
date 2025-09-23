######## Read Table S4, prepare data for boxplots
library(readxl)
library(dplyr)
library(tidyverse)
library(stringr)
library(maftools)

# load Table S4 excel sheets
Table_S4 <- excel_sheets("/Users/yantong/TrPLet_RShiny/Table_S4.xlsx")
Table_S4_sheets <- lapply(Table_S4, function(x) read_excel("/Users/yantong/TrPLet_RShiny/Table_S4.xlsx", sheet=x))
Table_S4[28] <- "TCGA_sRCC"
names(Table_S4_sheets) <- Table_S4
# Remove columns "mean", "CCLE", and "diff" from each sheet in Table_S4_dfs
Table_S4_dfs <- lapply(Table_S4_sheets, as.data.frame) # 23sheets --> 38 sheets
Table_S4_dfs <- lapply(Table_S4_dfs, function(df) df[ , !(names(df) %in% c("mean", "CCLE_mean", "diff"))]) # make each sheet contain only chronos scores and gene names

# TCGA_srRCC manually modification: 1. remove duplicated columns (keep the first 50 + Gene) 2. change name to TCGA_sRCC
Table_S4_dfs[["TCGA_sRCC"]] <- Table_S4_dfs[["TCGA_sRCC"]] %>%
  select(c(1:50, 54))
colnames(Table_S4_dfs[["TCGA_sRCC"]]) <- gsub("\\.\\.\\..*$", "", colnames(Table_S4_dfs[["TCGA_sRCC"]]))

# TCGA_pRCCT2 contains 36 empty columns, remove them
Table_S4_dfs[["TCGA_pRCCT2"]] <- Table_S4_dfs[["TCGA_pRCCT2"]] %>%
  select(c(37:99))

genes <- Table_S4_dfs[[1]]$Gene
saveRDS(genes, "/Users/yantong/TrPLet_RShiny/6283gene_names.rds")

TSS <- read.delim("/Users/yantong/TrPLet_RShiny/tissueSourceSite.tsv", header=TRUE, sep="\t")
TSS$TSS.Code <- paste0("TCGA-", TSS$TSS.Code) # modify the TSS column to match with the column names in TCGA_all
TCGA_abb <- read.delim("/Users/yantong/TrPLet_RShiny/diseaseStudy.tsv", header=TRUE, sep="\t")
map_df <- merge(TSS, TCGA_abb, by="Study.Name")

# Annotate Unscreened_CCLE cellline name
cellline_map <- read.delim("/Users/yantong/TrPLet_RShiny/sample_info.csv", header=T, sep=",")
cellline_map <- setNames(cellline_map$stripped_cell_line_name, cellline_map$DepMap_ID)
cellline_map <- append(cellline_map, c("Gene"="Gene"))
colnames(Table_S4_dfs[["Unscreened_CCLE"]]) <- ifelse(colnames(Table_S4_dfs[["Unscreened_CCLE"]]) %in% names(cellline_map), 
                                                      cellline_map[colnames(Table_S4_dfs[["Unscreened_CCLE"]])],
                                                      colnames(Table_S4_dfs[["Unscreened_CCLE"]]))

# Create one dataframe for each gene across all different tumor grouping (38 sheets) 6283 genes --> 6283 dfs
plot_data <- list()
for (i in 1:length(genes)){
  df <- data.frame(matrix(ncol=3, nrow=0))
  
  gene_name <- genes[i]
  for (tumor_type in Table_S4){
    data <- Table_S4_dfs[[tumor_type]]
    
    if (gene_name %in% data$Gene){
      sub_data <- data[data$Gene == gene_name, ]
      sub_data_long <- sub_data %>%
        pivot_longer(cols = -c(Gene), 
                     names_to = "Tumor_ID",
                     values_to = "Chronos_score")
      sub_data_long$Tumor_type <- tumor_type
      sub_data_long$Tumor_Sample_ID <- ifelse(sub_data_long$Tumor_type == "TCGA_all",
                                              sub("^([A-Za-z0-9-]+-[0-9]+[A-Z]-[0-9]+)(?:R-[A-Za-z0-9]+-[0-9]+.*)?$", "\\1", sub_data_long$Tumor_ID),
                                              NA)
      df <- rbind(df, sub_data_long)
    }
    
  }
  
  df$Grouping_L1 <- ifelse(df$Tumor_type %in% c("tRCC_CCLs", "ASPS_CCLs", "Unscreened_CCLE"), 
                           "Celllines",
                           ifelse(df$Tumor_type %in% c("TCGA_tRCC", "TCGA_CIMPpRCC", "TCGA_FHdeficientRCC", "TCGA_MDchRCC",
                                                       "TCGA_Eosinophilic_chRCC", "TCGA_oncocytoma", "TCGA_ccRCC", "TCGA_pRCCT1",
                                                       "TCGA_pRCCT2", "TCGA_chRCC", "TCGA_sRCC", "Sun_tRCC", "Qu_tRCC", "Brugarolas_tRCC", 
                                                       "Wang2016_CDC", "Msaouel_RMC", "Coutinho_Wilms", "Drost_Wilms", "Drost_MRTK",
                                                       "Drost_NephrogenicRest","Drost_pediatricRCC", "Brugarolas_ccRCC", "Braun_ccRCC",
                                                       "Wang2020_PEComa_RenalAML"),
                                  "Kidney Tumor Subtypes",
                                  ifelse(df$Tumor_type == "TCGA_all",
                                         "TCGA Lineages",
                                         "Sarcoma Subtypes")))
  df$TSS <- sub("^([A-Z0-9]+-[A-Z0-9]+)-.*", "\\1", df$Tumor_ID)
  df_mapped <- merge(df, map_df[, c("TSS.Code", "Study.Abbreviation")], by.x="TSS", by.y="TSS.Code", all.x=T)
  
  df_mapped$Grouping_L2 <- ifelse(df_mapped$Tumor_type == "TCGA_all",
                                  df_mapped$Study.Abbreviation,
                                  df_mapped$Tumor_type)
  
  plot_data[[i]] <- df_mapped
}
names(plot_data) <- genes

plot_data <- lapply(plot_data, as.data.frame)
plot_data <- plot_data[order(names(plot_data))] # sort the order in alphabetical order
saveRDS(plot_data, "/Users/yantong/TrPLet_RShiny/6283dfs_merged_across_Grouping.rds")


########### mutation
library(TCGAmutations)
library(maftools)

TCGA_lineages <- tcga_available()$Study_Abbreviation[1:33]
mutation_df <- data.frame(matrix(ncol=3, nrow=0))
for (i in TCGA_lineages){
  subdf <- tcga_load(study=i)
  subdf <- subdf@data[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  subdf$Tumor_Sample_ID <- sub("^([A-Za-z0-9-]+-[0-9]+[A-Z]-[0-9]+).*", "\\1", subdf$Tumor_Sample_Barcode)
  subdf$Lineage <- i
  mutation_df <- rbind(mutation_df, subdf)
  
}
saveRDS(mutation_df, "/Users/yantong/TrPLet_RShiny/mutation.rds")



######## Prepare reference matrix (650 genes * 1095 celllines, from DepMap_23Q2)
library(tidyverse)

# load data into r (download from https://figshare.com/articles/dataset/DepMap_23Q2_Public/22765112?file=40448555)
DepMap <- read.csv("/Users/yantong/TrPLet_RShiny/CRISPR_DepMap_Internal_23Q2_Score_Chronos.csv", check.names = FALSE, header = TRUE)

Model_meta <- read.csv("/Users/yantong/TrPLet_RShiny/Model.csv", , check.names = FALSE, header = TRUE)

rownames(DepMap) <- DepMap[,1]
DepMap <- DepMap[, 2:ncol(DepMap)] 
DepMap$Model_ID <- rownames(DepMap)
DepMap <- merge(DepMap, Model_meta[,c("ModelID", "StrippedCellLineName")], by.x = "Model_ID", by.y = "ModelID", all.x=T)
DepMap <- DepMap[!is.na(DepMap$StrippedCellLineName),]
rownames(DepMap) <- DepMap$StrippedCellLineName
DepMap <- DepMap[, !colnames(DepMap) %in% c("StrippedCellLineName", "Model_ID")]

DepMap <- as.data.frame(t(DepMap)) # transpose the dataframe to have genes as rows, celllines as columns

DepMap <- DepMap %>%
  rownames_to_column(var="Gene")
DepMap$Gene <- gsub("\\s*\\(\\d+\\)", "", DepMap$Gene) # remove (number) in the gene names
DepMap_GeneList <- DepMap$Gene

DepMap_6283genes <- DepMap[DepMap$Gene %in% genes, ] # 6208 genes 


# create df for each gene from the reference matrix to create boxplot for each gene later
ref_data <- list()
for (i in 1:length(genes)){
  gene_name <- genes[i]
  data <- DepMap_6283genes[DepMap_6283genes$Gene == gene_name,]
  data_long <- data %>%
    pivot_longer(cols = -c(Gene),
                 names_to = "CellLineName",
                 values_to = "Chronos_score")
  
  ref_data[[i]] <- data_long
}
names(ref_data) <- genes
ref_data <- lapply(ref_data, as.data.frame)
saveRDS(ref_data, "/Users/yantong/TrPLet_RShiny/6283dfs_DepMap.rds")


# Save r information for a subset of 6283 genes
models_df <- read.csv("/Users/yantong/TrPLet_RShiny/summary_sheet_bestmodel_bestMpermodel_noKNN.csv", header=T)
saveRDS(models_df, "/Users/yantong/TrPLet_RShiny/AllgenesModels.rds")
Models_6283genes <- models_df[models_df$Correl_Coeffs_5CV >= 0.2, ]
saveRDS(Models_6283genes, "/Users/yantong/TrPLet_RShiny/6283Models.rds")



################## module & feature selection
feature_select_path <- "/Users/yantong/TrPLet_RShiny/model_and_feature_number_vs_performance_no_knn"
features_files <- list.files(path = feature_select_path, pattern = "*.csv", full.names = TRUE)
features_files_list <- setNames(lapply(features_files, read.csv), 
                                tools::file_path_sans_ext(basename(features_files)) %>%
                                  str_extract("(?<=_5foldCV_).*?(?=_testdata)"))
module_genes <- features_files_list[[1]]$Gene_Name

features_files_list_updated <- list()
for (i in 1:length(module_genes)){
  df <- data.frame(matrix(ncol=4, nrow=0))
  gene_name <- module_genes[i]
  
  for (module_name in names(features_files_list)){
    data <- features_files_list[[module_name]]
    if (gene_name %in% data$Gene_Name){
      sub_data <- data[data$Gene_Name == gene_name, ]
      
      sub_data$Number_of_features <- sub("_.*", "", module_name)
      sub_data$Module <- sub("^[^_]+_", "", module_name)
      df <- rbind(df, sub_data)
      
    }
    
  }
  features_files_list_updated[[i]] <- df
  
}
names(features_files_list_updated) <- module_genes
saveRDS(features_files_list_updated, "/Users/yantong/TrPLet_RShiny/16845dfs_mergedModuleSelections.rds")





############################ Process predictor weights files
# Define the folder path
folder_path <- "/Users/yantong/TrPLet_RShiny/coefficients"

# Get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# Read all CSV files into a list of data frames
df_list <- setNames(lapply(csv_files, read.csv), 
                    tools::file_path_sans_ext(basename(csv_files)) %>% 
                      gsub("_coefficients_.*", "", .))

df_list_updated <- c()
for (g in names(df_list)){
  df <- df_list[[g]]
  df$Features <- gsub("_x$", "", df$Features)
  df$abs_coef <- abs(df$Coefficients)
  df_sorted <- df[order(-df$abs_coef), ]
  df_sorted$Rank <- 1:nrow(df_sorted)
  df_list_updated[[g]] <- df_sorted
}

# Save the list as an RDS file
saveRDS(df_list_updated, file = "/Users/yantong/TrPLet_RShiny/6283featureWeights.rds")



############################# Process data source table
data_sources <- read_excel("/Users/yantong/TrPLet_RShiny/Data_sources.xlsx")
subset_data_sources <- data_sources[, c(1,2,3)]
# Add hyperlinks for the third column
subset_data_sources$`Study details` <- ifelse(
  startsWith(subset_data_sources$`Study details`, "https"),
  paste0('<a href="', subset_data_sources$`Study details`, '" target="_blank">', subset_data_sources$`Study details`, '</a>'),
  subset_data_sources$`Study details`
)


View(subset_data_sources)
saveRDS(subset_data_sources, file = "/Users/yantong/TrPLet_RShiny/dataSources.rds")
