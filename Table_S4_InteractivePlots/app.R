set.seed(1234)
library(shiny)
library(plotly)
library(dplyr)
library(DT)


# load data
plot_data <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/6283dfs_merged_across_Grouping.rds")
genes <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/6283gene_names.rds")

ref_data <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/6283dfs_DepMap.rds")

models_df <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/AllgenesModels.rds")
features_df <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/6283featureWeights.rds")
module_selection <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/16845dfs_mergedModuleSelections.rds")
mutation_df <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/mutation.rds")

dataSources_df <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet/TrPLet_Table_S4/trplet_rshiny_final/data/dataSources.rds")

genes_6283 <- sort(genes) 
genes_all <- sort(models_df$Gene_Name)

# -------------------- global colors reused across renders ---------------------
TCGA_COLORS <- c(
  "#1f77b4", "#FFA07A", "#2ca02c", "#AC6A9F", "#9467bd", "#008080", "#e377c2",
  "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
  "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
  "#637939", "#8c6d31", "#5f9ea0", "#FFE119", "#5254a3", "#6b6ecf", "#9c9ede",
  "#cedb9c", "#e7cb94", "#e7ba52", "#bd9e39", "#8c6d78"
)

BASE_BAR_COLORS <- c(
  SVR = "#aec7e8",
  ridge_regression = "#c5b6d2",
  linear_SVR = "#cedb9c",
  lasso_regression = "#fffacd",
  elasticnet_regression = "#f7b6d2"
)

# ------------------------------ UI -------------------------------
ui <- fluidPage(
  
  # Centered Title
  tags$div(
    style = "text-align: center; margin-top: 20px;",
    tags$h1("Interactive Dashboard for TrPLet - Transcriptional Prediction of Lethality")
  ),
  
  tabsetPanel(
    tabPanel("TrPLet Overview", 
             # Centered Title
             tags$h3(
               style = "text-align: center; margin-top: 10px;",
               "Welcome to TrPLet - Cancer Dependency Prediction from RNA-seq Data"
             ),
             # Project Introduction Section
             fluidRow(
               column(12,
                      h3("Background"),
                      tags$div(style = "font-size: 16px;", 
                               tags$p(
                                 "Large-scale functional genetic screens have enabled the discovery of many selective cancer dependencies. 
                                 However, rare cancers have often been underrepresented in prior efforts, 
                                 and some cancer models may not be amenable to high-throughput screening."
                               ),
                               tags$p(
                                 "This project aims to predict gene dependency scores based on tumor or cell line transcriptional profile (RNA-Seq)."
                               ),
                               tags$p(
                                 "This interactive dashboard is designed to enable users to explore predicted dependencies across 11K+ TCGA tumor samples, 
                                 850+ kidney tumors across 13 different molecular subtypes, and 600+ previously unscreened cancer cell lines."
                               ),
                               tags$p(
                                 "We also provide visualizations to report performance for each gene dependency prediction model as well as to explore top features driving each model."
                               )
                      )
               )
             ),
             
             # Tab Explanation Section
             fluidRow(
               column(12,
                      h3("User Guide:"),
                      tags$div(style = "font-size: 16px;", 
                               tags$p("The application consists of multiple interactive tabs:"),
                               tags$ul(
                                 tags$li(strong("Dependency Visualizations:"), HTML(" Interactive visualization of dependency scores predicted based on RNA-seq data from the TCGA, other tumor cohorts, 
                                         or cell line models. Experimentally-derived dependency scores from the <a href='https://depmap.org/portal/' target='_blank'>Cancer DepMap</a> are also shown for comparison by lineage (CRISPR KO, 23Q2).")),
                                 tags$li(strong("Model & Feature Selection:"), " Compare model performance based on different number of features and gene models for a given gene of interest."),
                                 tags$li(strong("Feature Importance:"), " Explore the contribution of each feature to the dependency prediction for 6283 highly predictable genes."),
                                 tags$li(strong("Prediction Consistency:"), " Assess the performance (r) for 6283 highly predictable genes (all with r ≥ 0.2 between predicted and experimentally observed).")
                               )
                      )
               )
             ),
             
             # Citation Section
             fluidRow(
               column(12,
                      h4("Useful Links:"),
                      tags$p("Snakemake Pipeline: ",
                             tags$a(
                               href = "https://github.com/SViswanathanLab/TrPLet/tree/main/snakemake_TrPLet",
                               "https://github.com/SViswanathanLab/TrPLet/tree/main/snakemake_TrPLet",
                               target = "_blank"
                             )
                      ),
                      tags$p("Please cite: ",
                             tags$a(
                               href = "https://doi.org/10.1101/2024.10.24.620074",
                               "https://doi.org/10.1101/2024.10.24.620074",
                               target = "_blank"
                             )
                      )
               )
             )
    ),
    tabPanel("Dependency Visualizations", 
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 # Dropdown for gene selection
                 selectInput("selected_gene", "Select a gene:", choices = genes_6283, selected = genes_6283[1], selectize = TRUE),
                 
                 # Dropdown for tumor type selection (initially populated based on the first gene)
                 selectInput("selected_grouping_l1", "Select RNA-seq type:", 
                             choices = unique(plot_data[[genes_6283[1]]]$Grouping_L1), 
                             selected = unique(plot_data[[genes_6283[1]]]$Grouping_L1)[1],
                             selectize = TRUE),
                 # Checkbox group for Grouping_L2 options (populated dynamically)
                 uiOutput("grouping_l2_ui")
               ),
               
               mainPanel(
                 width = 9, 
                 plotlyOutput("boxplot_plot", height = "650px", width = "1000px")
               )
             ),
    ),
    tabPanel("Model & Feature Selection", 
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 # Dropdown for gene selection
                 selectInput("selected_gene_module2", "Select a gene:", choices = names(module_selection), selected = names(module_selection)[1], selectize = TRUE),
                 checkboxInput("highest_bar", "Highlight Best Performance", value = FALSE),
               ),
               
               mainPanel(
                 width = 9, 
                 plotlyOutput("bar_plot2", height = "600px", width = "1100px")
               )
             ),
    ),    
    tabPanel("Feature Importance", 
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 # Dropdown for gene selection
                 selectInput("selected_gene_feature", "Select a gene:", choices = genes_6283, selected = genes_6283[1], selectize = TRUE)
               ),
               
               mainPanel(
                 width = 10, 
                 plotlyOutput("feature_plot", height = "600px", width = "1200px")
               )
             ),
    ),
    tabPanel("Prediction Consistency", 
             sidebarLayout(
               sidebarPanel(
                 width = 2,
                 # Dropdown for gene selection
                 selectInput("selected_gene_scatter", "Select a gene:", choices = genes_all, selected = genes_all[1], selectize = TRUE)
               ),
               
               mainPanel(
                 width = 10, 
                 plotlyOutput("scatter_plot", height = "600px", width = "1200px")
               )
             ),
    ),
    tabPanel("Acknowledgement & License",
             fluidRow(
               column(12,
                      h3("Acknowledgement:"),
                      p("The interactive data visualization dashboard is powered by the following R libraries: "),
                      tags$ul(
                        tags$li(strong("shiny:"), " Used for building the web application interface"),
                        tags$li(strong("plotly:"), " Used for interactive plotting and visualizations"),
                        tags$li(strong("TCGAmutations:"), " Source for somatic mutations from TCGA cohorts"),
                        tags$li(strong("DepMap:"), "Used for comparison with predicted cancer dependency scores")
                      ),
                      p("This web application was developed by Yantong Cui at Viswanathan Lab.")
               )
             ),
             #### add another section for data sources (2 columns: 1column for subtype name and 2nd column for link of the study)
             fluidRow(
               column(9,
                      h3("Data Sources:"),
                      DTOutput("dataSources_table")
               ),
             ),
             fluidRow(
               column(12,
                      h3("Copyright and License Information:"),
                      p("©2025 Viswanathan Lab at Dana-Farber Cancer Institute")
               )
             )
    )
  )
)

# --------------------------------- server -------------------------------------
server <- function(input, output, session) {
  
  # Cache per-gene frames so we don't reconstruct on every downstream reactive
  gene_plot_data <- reactive({
    req(input$selected_gene)
    plot_data[[input$selected_gene]]
  }) %>% bindCache(input$selected_gene)
  
  ref_plot_data <- reactive({
    req(input$selected_gene)
    ref_data[[input$selected_gene]]
  }) %>% bindCache(input$selected_gene)
  
  sub_mutation <- reactive({
    if (is.null(input$gene_mutated)) return(mutation_df[0, c("Tumor_Sample_ID","Hugo_Symbol")])
    mutation_df[mutation_df$Hugo_Symbol == input$gene_mutated, c("Tumor_Sample_ID","Hugo_Symbol"), drop = FALSE]
  }) %>% bindCache(input$gene_mutated)
  
  grouping_l2_choices <- reactive({
    req(input$selected_gene, input$selected_grouping_l1)
    pdf <- gene_plot_data()
    unique(pdf$Grouping_L2[pdf$Grouping_L1 == input$selected_grouping_l1])
  }) %>% bindCache(input$selected_gene, input$selected_grouping_l1)
  
  output$grouping_l2_ui <- renderUI({
    req(input$selected_grouping_l1)
    choices <- grouping_l2_choices()
    req(length(choices) > 0)
    
    prev_selection <- isolate(input$selected_grouping_l2)
    selected <- if (!is.null(prev_selection) && prev_selection %in% choices) prev_selection else choices[1]
    
    dropdown_ui <- selectInput(
      inputId = "selected_grouping_l2",
      label = "Select subtype:",
      choices = choices,
      selected = selected,
      selectize = TRUE
    )
    
    additional_ui <- if (input$selected_grouping_l1 == "TCGA Lineages") {
      tagList(
        checkboxInput("tcga_all_lineages", "Compare All Lineages", value = FALSE),
        conditionalPanel(
          condition = "input.tcga_all_lineages",
          selectInput(
            inputId = "gene_mutated",
            label = "Show mutations - Select gene:",
            choices = genes_6283,
            selected = genes_6283[1],
            selectize = TRUE
          )
        )
      )
    } else {
      checkboxInput("compare_all_checkbox", "Compare All", value = FALSE)
    }
    
    tagList(dropdown_ui, additional_ui)
  })
  
  # Keep selected gene synced across tabs 
  observeEvent(input$selected_gene, {
    if (!identical(input$selected_gene_module2, input$selected_gene))
      updateSelectInput(session, "selected_gene_module2", selected = input$selected_gene)
    if (!identical(input$selected_gene_feature, input$selected_gene))
      updateSelectInput(session, "selected_gene_feature", selected = input$selected_gene)
    if (!identical(input$selected_gene_scatter, input$selected_gene))
      updateSelectInput(session, "selected_gene_scatter", selected = input$selected_gene)
  }, ignoreInit = TRUE)
  
  # -------------------- Dependency Visualizations -----------------------
  output$boxplot_plot <- renderPlotly({
    # Use cached frames and reuse local variables to avoid incidental copies
    gene_plot_data_local <- gene_plot_data()
    ref_plot_data_local  <- ref_plot_data()
    
    # factor only if needed 
    if (!is.factor(gene_plot_data_local$Grouping_L2))
      gene_plot_data_local$Grouping_L2 <- factor(gene_plot_data_local$Grouping_L2)
    
    if (input$selected_grouping_l1 == "TCGA Lineages") {
      if (!is.null(input$tcga_all_lineages) && input$tcga_all_lineages) {
        subset_data <- gene_plot_data_local[gene_plot_data_local$Tumor_type == "TCGA_all", , drop = FALSE]
        sm <- sub_mutation()
        if (nrow(sm)) {
          subset_data$mutation <- ifelse(subset_data$Tumor_Sample_ID %in% sm$Tumor_Sample_ID, "mutated", "not-mutated")
        } else {
          subset_data$mutation <- "not-mutated"
        }
        
        colors <- TCGA_COLORS[seq_len(min(length(TCGA_COLORS), length(unique(subset_data$Study.Abbreviation))))]
        
        # main box
        p1 <- plot_ly(
          data = subset_data,
          x = ~as.factor(Study.Abbreviation),
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 3, opacity = 0.5),
          color = ~as.factor(Study.Abbreviation),
          colors = colors,
          text = ~paste('Tumor Lineage:', Study.Abbreviation, '<br>Tumor ID:', Tumor_Sample_ID, '<br>Chronos score:', Chronos_score),
          hoverinfo = 'text',
          name = 'TCGA Tumor All Lineages'
        )
        
        if (!is.null(input$gene_mutated)) {
          p1 <- p1 %>% add_trace(
            data = subset_data[subset_data$mutation == "mutated", , drop = FALSE],
            x = ~as.factor(Study.Abbreviation),
            y = ~Chronos_score,
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 5, color = "darkred", opacity = 1),
            text = ~paste('Tumor Lineage:', Study.Abbreviation, '<br>Tumor ID:', Tumor_Sample_ID, '<br>Chronos score:', Chronos_score),
            hoverinfo = 'text',
            showlegend = TRUE,
            name = "Mutated Samples"
          )
        }
        
        p2 <- plot_ly(
          data = ref_plot_data_local,
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 3, opacity = 0.9, color = '#EF8636'),
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),
          hoverinfo = 'text',
          color = ~as.factor(Gene),
          colors = "#EF8636",
          name = 'DepMap Chronos Scores'
        )
        
        return(
          subplot(p1, p2, nrows = 1, shareY = TRUE) %>%
            layout(
              title = list(text = paste0("Boxplots of Chronos Scores for ", input$selected_gene), x = 0.5),
              yaxis = list(title = "Chronos Scores"),
              showlegend = FALSE,
              margin = list(t = 100, b = 100),
              xaxis = list(domain = c(0, 0.97), tickangle = 90, title = "TCGA Tumor Lineages"),
              xaxis2 = list(domain = c(0.98, 1), tickangle = 90)
            )
        )
      } else {
        subset_data <- gene_plot_data_local[gene_plot_data_local$Grouping_L2 == input$selected_grouping_l2, , drop = FALSE]
        p_val <- suppressWarnings(wilcox.test(subset_data$Chronos_score, ref_plot_data_local$Chronos_score, alternative = "two.sided"))$p.value
        
        gene_boxplot <- plot_ly(
          data = subset_data,
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 7, opacity = 0.5),
          text = ~paste('Type:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),
          hoverinfo = 'text',
          name = paste0(input$selected_grouping_l2, "<br>Predicted Chronos Scores")
        )
        
        ref_boxplot <- plot_ly(
          data = ref_plot_data_local,
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 7, opacity = 0.5),
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),
          hoverinfo = 'text',
          name = 'DepMap Chronos Scores'
        )
        
        return(
          subplot(gene_boxplot, ref_boxplot, nrows = 1, shareY = TRUE) %>%
            layout(
              title = list(text = paste0("Boxplots of Chronos Scores for ", input$selected_gene), x = 0.5),
              annotations = list(
                list(
                  text = paste0("p-value: ", format(p_val, digits = 3, scientific = TRUE)),
                  x = 1.15, y = 0.85, showarrow = FALSE, xref = "paper", yref = "paper",
                  xanchor = "center", yanchor = "top", font = list(size = 16)
                )
              ),
              yaxis = list(title = "Chronos Scores"),
              showlegend = TRUE,
              margin = list(t = 100, b = 100),
              xaxis = list(domain = c(0, 0.45)),
              xaxis2 = list(domain = c(0.55, 1))
            )
        )
      }
    } else {
      if (!is.null(input$compare_all_checkbox) && input$compare_all_checkbox) {
        subset_data <- gene_plot_data_local[gene_plot_data_local$Grouping_L1 == input$selected_grouping_l1, , drop = FALSE]
        subset_data <- droplevels(subset_data)
        subset_data$Grouping_L2 <- factor(subset_data$Grouping_L2, levels = unique(subset_data$Grouping_L2))
        
        if (input$selected_grouping_l1 == "Kidney Tumor Subtypes") {
          subtype_legend_map <- data.frame(
            Grouping_L2 = c("Brugarolas_ccRCC","Braun_ccRCC","Wang2016_CDC","Wang2020_PEComa_RenalAML","Drost_MRTK",
                            "Drost_NephrogenicRest","Drost_pediatricRCC","Msaouel_RMC","Sun_tRCC","TCGA_pRCCT1",
                            "TCGA_pRCCT2","TCGA_sRCC","TCGA_tRCC","TCGA_ccRCC","TCGA_chRCC","TCGA_CIMPpRCC",
                            "TCGA_MDchRCC","TCGA_FHdeficientRCC","TCGA_Eosinophilic_chRCC","TCGA_oncocytoma",
                            "Brugarolas_tRCC","Drost_Wilms","Qu_tRCC","Coutinho_Wilms"),
            color_group = c("ccRCC","ccRCC","CDC","RenalAML","MRTK","NephrogenicRest","pediatricRCC","RMC",
                            "tRCC","pRCC","pRCC","sRCC","tRCC","ccRCC","chRCC","pRCC","chRCC","FHdeficientRCC",
                            "chRCC","oncocytoma","tRCC","Wilms","tRCC","Wilms"),
            stringsAsFactors = FALSE
          )
          legend_colors <- c(
            ccRCC="#1f77b4", CDC="#2ca02c", RenalAML="#f7b6d2", MRTK="#9467bd",
            NephrogenicRest="#008080", pediatricRCC="#e377c2", RMC="#7f7f7f", tRCC="#AC6A9F",
            pRCC="#aec7e8", sRCC="#ffbb78", chRCC="#c49c94", FHdeficientRCC="#dbdb8d",
            oncocytoma="#5f9ea0", Wilms="#9edae5"
          )
          subset_data <- merge(subset_data, subtype_legend_map, by = "Grouping_L2", all.x = TRUE)
          
          p1 <- plot_ly(
            data = subset_data,
            x = ~as.factor(Grouping_L2),
            y = ~Chronos_score,
            type = 'box',
            boxpoints = 'all',
            jitter = 0.3,
            pointpos = 0,
            marker = list(size = 5, opacity = 0.5),
            color = ~color_group,
            split = ~color_group,
            colors = legend_colors,
            text = ~paste('Tumor Subtype:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),
            hoverinfo = 'text',
            name = ~color_group,
            showlegend = TRUE
          )
        } else if (input$selected_grouping_l1 == "Sarcoma Subtypes") {
          colors <- c("#1f77b4", "#ffbb78", "#2ca02c", "#AC6A9F", "#9467bd", "#008080", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
          p1 <- plot_ly(
            data = subset_data,
            x = ~as.factor(Grouping_L2),
            y = ~Chronos_score,
            type = 'box',
            boxpoints = 'all',
            jitter = 0.3,
            pointpos = 0,
            marker = list(size = 5, opacity = 0.5),
            color = ~as.factor(Grouping_L2),
            colors = colors,
            text = ~paste('Tumor Subtype:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),
            hoverinfo = 'text',
            name = paste('All', input$selected_grouping_l1)
          )
        } else {
          colors <- c("#1f77b4", "#c5b0d5", "#5f9ea0")
          p1 <- plot_ly(
            data = subset_data,
            x = ~as.factor(Grouping_L2),
            y = ~Chronos_score,
            type = 'box',
            boxpoints = 'all',
            jitter = 0.3,
            pointpos = 0,
            marker = list(size = 5, opacity = 0.5),
            color = ~as.factor(Grouping_L2),
            colors = colors,
            text = ~paste('Tumor Subtype:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),
            hoverinfo = 'text',
            name = paste('All', input$selected_grouping_l1)
          )
        }
        
        p2 <- plot_ly(
          data = ref_plot_data_local,
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 5, opacity = 0.9, color = '#EF8636'),
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),
          hoverinfo = 'text',
          color = ~as.factor(Gene),
          colors = '#EF8636',
          name = 'DepMap Chronos Scores'
        )
        
        return(
          subplot(p1, p2, nrows = 1, shareY = TRUE) %>%
            layout(
              title = list(text = paste0("Boxplots of Chronos Scores for ", input$selected_gene), x = 0.5),
              yaxis = list(title = "Chronos Scores"),
              showlegend = TRUE,
              margin = list(t = 100, b = 100),
              xaxis = list(domain = c(0, 0.94), tickangle = 90, title = paste("All", input$selected_grouping_l1)),
              xaxis2 = list(domain = c(0.95, 1), tickangle = 90)
            )
        )
      } else {
        subset_data <- gene_plot_data_local[gene_plot_data_local$Grouping_L2 == input$selected_grouping_l2, , drop = FALSE]
        p_val <- suppressWarnings(wilcox.test(subset_data$Chronos_score, ref_plot_data_local$Chronos_score, alternative = "two.sided"))$p.value
        
        gene_boxplot <- plot_ly(
          data = subset_data,
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 7, opacity = 0.5),
          text = ~paste('Type:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),
          hoverinfo = 'text',
          name = paste0(input$selected_grouping_l2, "<br>Predicted Chronos Scores")
        )
        
        ref_boxplot <- plot_ly(
          data = ref_plot_data_local,
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(size = 7, opacity = 0.5),
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),
          hoverinfo = 'text',
          name = 'DepMap Chronos Scores'
        )
        
        return(
          subplot(gene_boxplot, ref_boxplot, nrows = 1, shareY = TRUE) %>%
            layout(
              title = list(text = paste0("Boxplots of Chronos Scores for ", input$selected_gene), x = 0.5),
              annotations = list(
                list(
                  text = paste0("p-value: ", format(p_val, digits = 3, scientific = TRUE)),
                  x = 1.15, y = 0.85, showarrow = FALSE, xref = "paper", yref = "paper",
                  xanchor = "center", yanchor = "top", font = list(size = 16)
                )
              ),
              yaxis = list(title = "Chronos Scores"),
              showlegend = TRUE,
              margin = list(t = 100, b = 100),
              xaxis = list(domain = c(0, 0.45)),
              xaxis2 = list(domain = c(0.55, 1))
            )
        )
      }
    }
  }) %>% bindCache(
    input$selected_gene, input$selected_grouping_l1, input$selected_grouping_l2,
    input$tcga_all_lineages, input$gene_mutated, input$compare_all_checkbox
  )
  
  # -------------------- Module & Feature Selection ------------------------
  output$bar_plot2 <- renderPlotly({
    gene_df <- module_selection[[input$selected_gene_module2]]
    
    sub_df_sorted <- gene_df %>%
      mutate(Number_of_features = factor(
        Number_of_features,
        levels = unique(gene_df$Number_of_features[order(as.numeric(gsub("N", "", gene_df$Number_of_features)))])
      ))
    
    max_corr <- max(sub_df_sorted$Correl_Coeffs_5CV, na.rm = TRUE)
    highlight_best <- isTRUE(input$highest_bar)
    
    p <- plot_ly()
    for (mod in unique(sub_df_sorted$Module)) {
      df_mod  <- sub_df_sorted[sub_df_sorted$Module == mod, , drop = FALSE]
      is_best <- df_mod$Correl_Coeffs_5CV == max_corr
      trace_color <- if (highlight_best && any(is_best)) ifelse(is_best, "darkred", BASE_BAR_COLORS[mod]) else BASE_BAR_COLORS[mod]
      
      p <- p %>% add_trace(
        data = df_mod, x = ~Number_of_features, y = ~Correl_Coeffs_5CV, type = "bar",
        name = mod, marker = list(color = trace_color, opacity = 0.7),
        text = ~paste("Model:", Module,
                      "<br>Number of Features:", gsub("N", "", Number_of_features),
                      "<br>Correlation:", round(Correl_Coeffs_5CV, 3)),
        hoverinfo = "text", textposition = "none"
      )
    }
    
    p %>% layout(
      title = "Grouped Bar Plot: Model Performance",
      xaxis = list(title = "Number of Features",
                   categoryorder = "array",
                   categoryarray = levels(sub_df_sorted$Number_of_features)),
      yaxis = list(title = "5-fold CV Correlation of Predicted vs. Observed Chronos Score in Test Data (R)"),
      barmode = "group"
    )
  }) %>% bindCache(input$selected_gene_module2, input$highest_bar)
  
  # -------------------- Feature Importance -------------------------------
  output$feature_plot <- renderPlotly({
    feature_df <- features_df[[input$selected_gene_feature]]
    
    gray_points <- feature_df[feature_df$Rank > 10, , drop = FALSE]
    red_point   <- feature_df[feature_df$Rank <= 10, , drop = FALSE]
    
    plot_ly() %>%
      add_trace(
        data = gray_points, type = "scatter", mode = "markers",
        x = ~Rank, y = ~abs_coef,
        marker = list(color = "lightgray", size = 5),
        textposition = "outside",
        text = ~paste('Feature:', Features, '<br>Coefficient:', Coefficients, '<br>Rank:', Rank),
        hoverinfo = 'text',
        name = "Lowly correlated Features"
      ) %>%
      add_trace(
        data = red_point, type = "scatter", mode = "markers",
        x = ~Rank, y = ~abs_coef,
        marker = list(color = "darkred", size = 10),
        textposition = "outside",
        text = ~paste('Feature:', Features, '<br>Coefficient:', Coefficients, '<br>Rank:', Rank),
        hoverinfo = 'text',
        name = "Top correlated Features"
      ) %>%
      layout(
        title = paste0(input$selected_gene_feature, " Chronos Score Predictors"),
        xaxis = list(title = "Feature Rank",
                     zeroline = TRUE,
                     range = c(-(max(feature_df$Rank))*0.05, max(feature_df$Rank+(max(feature_df$Rank))*0.05))),
        yaxis = list(title = "Relative Feature Importance",
                     zeroline = TRUE,
                     range = c(-(max(feature_df$abs_coef))*0.05, max(feature_df$abs_coef+(max(feature_df$abs_coef))*0.05))),
        showlegend = TRUE
      )
  }) %>% bindCache(input$selected_gene_feature)
  
  # -------------------- Prediction Consistency --------------------------
  output$scatter_plot <- renderPlotly({
    models_df_sorted <- models_df[order(-models_df$Correl_Coeffs_5CV), , drop = FALSE]
    models_df_sorted$Rank <- seq_len(nrow(models_df_sorted))
    
    gray_points <- models_df_sorted[models_df_sorted$Gene_Name != input$selected_gene_scatter, , drop = FALSE]
    red_point   <- models_df_sorted[models_df_sorted$Gene_Name == input$selected_gene_scatter, , drop = FALSE]
    
    plot_ly() %>%
      add_trace(
        data = gray_points, type = "scatter", mode = "markers",
        x = ~Rank, y = ~Correl_Coeffs_5CV,
        marker = list(color = "lightgray", size = 5),
        textposition = "outside",
        text = ~paste('Gene:', Gene_Name, '<br>r:', Correl_Coeffs_5CV, '<br>Rank:', Rank),
        hoverinfo = 'text',
        name = "Unselected Genes"
      ) %>%
      add_trace(
        data = red_point, type = "scatter", mode = "markers",
        x = ~Rank, y = ~Correl_Coeffs_5CV,
        marker = list(color = "darkred", size = 15),
        textposition = "outside",
        text = ~paste('Gene:', Gene_Name, '<br>r:', Correl_Coeffs_5CV, '<br>Rank:', Rank),
        hoverinfo = 'text',
        name = input$selected_gene_scatter
      ) %>%
      layout(
        title = "Scatterplot of Averaged 5-Fold Cross-Validation Correlation Coefficients for the Best Model per Gene",
        xaxis = list(title = "Rank"),
        yaxis = list(title = "Correlation of Predicted vs. Observed Chronos Score in Test Data (R)"),
        showlegend = TRUE
      )
  }) %>% bindCache(input$selected_gene_scatter)
  
  # -------------------- Data sources table --------------------------
  output$dataSources_table <- renderDT({
    datatable(dataSources_df, escape = FALSE, options = list(pageLength = 5))
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

