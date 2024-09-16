
Table_S4.xlsx (23 sheets)

CRISPRGeneEffect.csv (DepMap Public 23Q2, downloaded from https://figshare.com/articles/dataset/DepMap_23Q2_Public/22765112?file=40448555)


R SCRIPTS:
1. DataProcess.R - Prepare data for plotting
	Generate:
		657gene_names.rds
		657dfs_merged_across_23TumorTypes.rds
		CRISPRGeneEffect_DepMap_Public_23Q2.csv & CRISPRGeneEffect_DepMap_Public_23Q2.rds
		650genes_CRISPRGeneEffect_DepMap_Public_23Q2.csv & 650genes_CRISPRGeneEffect_DepMap_Public_23Q2.rds
		657dfs_DepMap.rds (7 empty dfs, 7 genes of the 657 are missed in CRISPRGeneEffect.csv)

		
2. Boxplot_1filter.R - Version1 Boxplot
	Apply one filter for gene names, boxplot of chromos scores of all tumors across 23 tumor types for the selected gene (from Table S4) + boxplot of DepMap data for the selected gene
	
	Run:
	> library(shiny)	> runApp('/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/Boxplot_1filter.R') 

3. Boxplot_2filters.R - Version2 Boxplot
	Apply 2 filters for gene names and tumor types, boxplot of chromos scores of all tumors for the selected gene and selected tumor type (from Table S4) + boxplot of DepMap data for the selected gene 
	
	Run:
	> library(shiny)	> runApp('/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/Boxplot_2filters.R') 

4. Boxplot_1filter_colors.R - Version3 Boxplot
	Apply one filter for gene names, boxplots of chromos scores of all tumors for the selected gene, colored by different tumor types (from Table S4) 
	
	Run:
	> library(shiny)	> runApp('/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/Boxplot_1filter_colors.R') 




