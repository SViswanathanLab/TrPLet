######## Read Table S4, prepare data for boxplots
library(readxl)
library(dplyr)
library(tidyverse)

# load data
Table_S4 <- excel_sheets("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/Table_S4.xlsx")
Table_S4_sheets <- lapply(Table_S4, function(x) read_excel("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/Table_S4.xlsx", sheet=x))
names(Table_S4_sheets) <- Table_S4
Table_S4_dfs <- lapply(Table_S4_sheets, as.data.frame) #23 sheets in total for now

genes <- Table_S4_dfs[[1]]$Gene
saveRDS(genes, "/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657gene_names.rds")

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
                     values_to = "Chromos_score")
      sub_data_long$Tumor_type <- tumor_type
      df <- rbind(df, sub_data_long)
    }
    
  }
  
  plot_data[[i]] <- df
}

names(plot_data) <- genes
plot_data <- lapply(plot_data, as.data.frame)
saveRDS(plot_data, "/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657dfs_merged_across_23TumorTypes.rds")




######## Prepare reference matrix (650 genes * 1095 celllines, from DepMap_23Q2)
library(tidyverse)

# load data into r (download from https://figshare.com/articles/dataset/DepMap_23Q2_Public/22765112?file=40448555)
DepMap <- read.csv("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/CRISPRGeneEffect.csv", check.names = FALSE, header = TRUE)

rownames(DepMap) <- DepMap[,1]
DepMap <- DepMap[, 2:ncol(DepMap)]
DepMap <- as.data.frame(t(DepMap)) # transpose the dataframe to have genes as rows, celllines as columns

DepMap <- DepMap %>%
  rownames_to_column(var="Gene")
DepMap$Gene <- gsub("\\s*\\(\\d+\\)", "", DepMap$Gene) # remove (number) in the gene names

write.csv(DepMap, "/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/CRISPRGeneEffect_DepMap_Public_23Q2.csv")
saveRDS(DepMap, "/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/CRISPRGeneEffect_DepMap_Public_23Q2.rds")



DepMap_657genes <- DepMap[DepMap$Gene %in% genes, ] # 650 genes 
print(genes[!genes %in% intersect(DepMap$Gene, genes)]) # 7 genes of the 657 are missing: "C15orf41" "FGFR1OP"  "CASC1"    "TMEM189"  "C7orf61"  "CXorf56"  "TAZ"
# Check synonyms for these 7 genes
aliases <- c("CDIN1", "CEP43", "DNAI7", "PEDS1", "SPACDR", "STEEP1", "TAFAZZIN")
DepMap_657genes <- DepMap[DepMap$Gene %in% genes | DepMap$Gene %in% aliases, ]

DepMap_657genes$Gene[DepMap_657genes$Gene == "CDIN1"] <- "C15orf41"
DepMap_657genes$Gene[DepMap_657genes$Gene == "CEP43"] <- "FGFR1OP"
DepMap_657genes$Gene[DepMap_657genes$Gene == "DNAI7"] <- "CASC1"
DepMap_657genes$Gene[DepMap_657genes$Gene == "PEDS1"] <- "TMEM189"
DepMap_657genes$Gene[DepMap_657genes$Gene == "SPACDR"] <- "C7orf61"
DepMap_657genes$Gene[DepMap_657genes$Gene == "STEEP1"] <- "CXorf56"
DepMap_657genes$Gene[DepMap_657genes$Gene == "TAFAZZIN"] <- "TAZ"

write.csv(DepMap_657genes, "/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657genes_CRISPRGeneEffect_DepMap_Public_23Q2.csv")
saveRDS(DepMap_657genes, "/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657genes_CRISPRGeneEffect_DepMap_Public_23Q2.rds")



# create df for each gene from the reference matrix to create boxplot for each gene later
ref_data <- list()
for (i in 1:length(genes)){
  gene_name <- genes[i]
  data <- DepMap_657genes[DepMap_657genes$Gene == gene_name,]
  data_long <- data %>%
    pivot_longer(cols = -c(Gene),
                 names_to = "Model_ID",
                 values_to = "Chromos_score")
  
  ref_data[[i]] <- data_long
}
names(ref_data) <- genes
ref_data <- lapply(ref_data, as.data.frame)
saveRDS(ref_data, "/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657dfs_DepMap.rds")

