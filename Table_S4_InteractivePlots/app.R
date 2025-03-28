set.seed(1234)
library(shiny)
library(plotly)
library(dplyr)
library(DT)


# load data
plot_data <- readRDS("6283dfs_merged_across_Grouping.rds")
genes <- readRDS("6283gene_names.rds")

ref_data <- readRDS("6283dfs_DepMap.rds")

models_df <- readRDS("6283Models.rds")
features_df <- readRDS("6283featureWeights.rds")
module_selection <- readRDS("16845dfs_mergedModuleSelections.rds")
mutation_df <- readRDS("mutation.rds")

dataSources_df <- readRDS("dataSources.rds")


# Define UI for the dashboard
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
                 selectInput("selected_gene", "Select a gene:", choices = genes, selected = genes[1], selectize = TRUE),
                 
                 # Dropdown for tumor type selection (initially populated based on the first gene)
                 selectInput("selected_grouping_l1", "Select RNA-seq type:", 
                             choices = unique(plot_data[[genes[1]]]$Grouping_L1), 
                             selected = unique(plot_data[[genes[1]]]$Grouping_L1)[1],
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
                 selectInput("selected_gene_feature", "Select a gene:", choices = genes, selected = genes[1], selectize = TRUE)
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
                 selectInput("selected_gene_scatter", "Select a gene:", choices = genes, selected = genes[1], selectize = TRUE)
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
                      p("This web application was developed by Yantong Cui at Viswanathan Lab."),
                      p("This project is funded by ...")
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


# Define server logic for generating the boxplot
server <- function(input, output, session) {
  
  grouping_l2_choices <- reactive({
    req(input$selected_gene, input$selected_grouping_l1)  # Ensure inputs are not NULL
    plot_data_subset <- plot_data[[input$selected_gene]]
    unique(plot_data_subset$Grouping_L2[plot_data_subset$Grouping_L1 == input$selected_grouping_l1])
  })
  
  
  # After users select the RNA-seq type, the level-3 selection pops up
  output$grouping_l2_ui <- renderUI({
    req(input$selected_grouping_l1, grouping_l2_choices())
    
    # Grab the previous selection directly here
    prev_selection <- isolate(input$selected_grouping_l2)
    choices <- grouping_l2_choices()
    
    selected <- if (!is.null(prev_selection) && prev_selection %in% grouping_l2_choices()) {
      prev_selection
    } else {
      grouping_l2_choices()[1]
    }
    
    # level 3 Subtype dropdown list UI
    dropdown_ui <- selectInput(
      inputId = "selected_grouping_l2",
      label = "Select subtype:",
      choices = grouping_l2_choices(),
      selected = selected,
      selectize = TRUE
    )
    
    # Checkbox UI
    additional_ui <- if (input$selected_grouping_l1 == "TCGA Lineages") {
      tagList(
        checkboxInput("tcga_all_lineages", "Compare All Lineages", value = FALSE),
        conditionalPanel(
          condition = "input.tcga_all_lineages",
          checkboxInput("mutation_button", "Show Mutations", value = FALSE)
        )
      )
    } else {
      checkboxInput("compare_all_checkbox", "Compare All", value = FALSE)
    }
    
    # Combine dropdown + conditionally displayed checkboxes
    tagList(
      dropdown_ui,
      additional_ui
    )
  })
  
  
  output$boxplot_plot <- renderPlotly({
    
    # Get the data for the selected gene
    gene_plot_data <- plot_data[[input$selected_gene]]
    sub_mutation <- mutation_df[mutation_df$Hugo_Symbol == input$selected_gene, ]
    
    gene_plot_data$Grouping_L2 <- as.factor(gene_plot_data$Grouping_L2)
    
    ref_plot_data <- ref_data[[input$selected_gene]]
    
    
    if (input$selected_grouping_l1 == "TCGA Lineages"){
      if (!is.null(input$tcga_all_lineages) && input$tcga_all_lineages){
        subset_data <- gene_plot_data[gene_plot_data$Tumor_type == "TCGA_all",]
        subset_data$mutation <- ifelse(subset_data$Tumor_Sample_ID %in% sub_mutation$Tumor_Sample_ID,
                                       "mutated",
                                       "not-mutated")
        p_val <- NA
        
        
        colors <- c(
          "#1f77b4", "#ff7f0e", "#2ca02c", "#AC6A9F", "#9467bd", "#008080", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
          "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
          "#393b79", "#637939", "#8c6d31", "#5f9ea0", "#FFE119", "#5254a3", "#6b6ecf", "#9c9ede", "#cedb9c", "#e7cb94",
          "#e7ba52", "#bd9e39", "#8c6d78"
        )
        
        study_colors <- setNames(colors[1:length(unique(subset_data$Study.Abbreviation))], unique(subset_data$Study.Abbreviation))
        
        if (!is.null(input$mutation_button) && input$mutation_button){
          p1 <- plot_ly(
            data = subset_data,
            x = ~as.factor(Study.Abbreviation),
            y = ~Chronos_score,
            type = 'box',
            boxpoints = 'all',
            jitter = 0.3,
            pointpos = 0,
            marker = list(
              size = 3, 
              opacity = 0.5),
            color = ~as.factor(Study.Abbreviation),
            colors = colors,
            text = ~paste('Tumor Lineage:', Study.Abbreviation, '<br>Tumor ID:', Tumor_Sample_ID, '<br>Chronos score:', Chronos_score),  # Tooltip info
            hoverinfo = 'text',               # Show only the text info on hover
            name = 'TCGA Tumor All Lineages'
          ) %>% add_trace(
            data = subset_data[subset_data$mutation == "mutated", ],  # Select only mutated samples
            x = ~as.factor(Study.Abbreviation),
            y = ~Chronos_score,
            type = 'scatter',
            mode = 'markers',
            marker = list(size = 5, color = "darkred", opacity = 1),  # Highlighted points in dark red
            text = ~paste('Tumor Lineage:', Study.Abbreviation, '<br>Tumor ID:', Tumor_Sample_ID, '<br>Chronos score:', Chronos_score),  
            hoverinfo = 'text',
            showlegend = TRUE,
            name = "Mutated Samples"
          )
        } else{
          p1 <- plot_ly(
            data = subset_data,
            x = ~as.factor(Study.Abbreviation),
            y = ~Chronos_score,
            type = 'box',
            boxpoints = 'all',
            jitter = 0.3,
            pointpos = 0,
            marker = list(
              size = 3, 
              opacity = 0.5),
            color = ~as.factor(Study.Abbreviation),
            colors = colors,
            text = ~paste('Tumor Lineage:', Study.Abbreviation, '<br>Tumor ID:', Tumor_Sample_ID, '<br>Chronos score:', Chronos_score),  # Tooltip info
            hoverinfo = 'text',               # Show only the text info on hover
            name = 'TCGA Tumor All Lineages'
          ) 
        }
        
        
        p2 <- plot_ly(
          data = ref_plot_data,
          y = ~Chronos_score,               # Chronos scores on the y-axis for the reference data
          type = 'box',                     # Boxplot type
          boxpoints = 'all',                # Show all data points on the same box
          jitter = 0.3,                     # Add some jitter to prevent points from overlapping
          pointpos = 0,                  # Position the points inside the box
          marker = list(size = 3, opacity = 0.9, color = '#EF8636'),  # Uniform color for the reference plot
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),  # Tooltip info
          hoverinfo = 'text',               # Show only the text info on hover
          color = ~as.factor(Gene),
          colors = "#EF8636",
          name = 'DepMap Mean Chronos Scores'
        )
        
        # Combine with subplot
        subplot(p1, p2, nrows = 1, shareY = TRUE) %>%
          layout(
            title = list(
              text = paste0("Boxplots of Chronos_scores for ", input$selected_gene),
              x = 0.5
            ),
            yaxis = list(title = "Chronos Scores"),
            showlegend = FALSE,
            margin = list(t = 100, b = 100),  # Increased bottom margin for better label spacing
            xaxis = list(
              domain = c(0, 0.97),
              tickangle = 90,  # Rotate x-axis labels
              title = "TCGA Tumor Lineages"
            ),
            xaxis2 = list(
              domain = c(0.98, 1),
              tickangle = 90  # Rotate the label for the last box
            )
          )

      } else {
        subset_data <- gene_plot_data[gene_plot_data$Grouping_L2 == input$selected_grouping_l2, ]
        p_val <- suppressWarnings(wilcox.test(subset_data$Chronos_score, ref_plot_data$Chronos_score, alternative="two.sided"))$p.value
        
        
        # Create the first boxplot for gene_data
        gene_boxplot <- plot_ly(
          data = subset_data,
          y = ~Chronos_score,               # Chronos scores on the y-axis
          type = 'box',                     # Boxplot type
          boxpoints = 'all',                # Show all data points on the same box
          jitter = 0.3,                     # Add some jitter to prevent points from overlapping
          pointpos = 0,                  # Position the points inside the box
          marker = list(size = 7, opacity = 0.5),  # Color points by Tumor_type
          text = ~paste('Type:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),  # Tooltip info
          hoverinfo = 'text',               # Show only the text info on hover
          name = paste0(input$selected_grouping_l2, "<br>Predicted Chronos Scores")
        )
        
        # Create the second boxplot for ref_plot_data
        ref_boxplot <- plot_ly(
          data = ref_plot_data,
          y = ~Chronos_score,               # Chronos scores on the y-axis for the reference data
          type = 'box',                     # Boxplot type
          boxpoints = 'all',                # Show all data points on the same box
          jitter = 0.3,                     # Add some jitter to prevent points from overlapping
          pointpos = 0,                  # Position the points inside the box
          marker = list(size = 7, opacity = 0.5),  # Uniform color for the reference plot
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),  # Tooltip info
          hoverinfo = 'text',               # Show only the text info on hover
          name = 'DepMap Mean Chronos Scores'
        )
        
        # Combine the two boxplots side by side using subplot
        subplot(gene_boxplot, ref_boxplot, nrows = 1, shareY = TRUE) %>%
          layout(
            title = list(
              text = paste0("Boxplots of Chronos Scores for ", input$selected_gene),
              x=0.5
            ),
            annotations = list(
              list(
                text = paste0("p-value: ", format(p_val, digits = 3, scientific=TRUE)),
                x=1.15,
                y=0.85,
                showarrow = FALSE,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="top",
                font=list(size=16)
              )
            ),
            yaxis = list(title = "Chronos Scores"),
            showlegend = TRUE,              # Show legend for Tumor_type colors
            margin = list(t = 100, b = 100),  # Adjust margins for better spacing
            xaxis = list(domain = c(0, 0.45)), # First boxplot occupies 45% of width
            xaxis2 = list(domain = c(0.55, 1)) # Second boxplot starts after 55%
          )
      }
    } else{
      if (!is.null(input$compare_all_checkbox) && input$compare_all_checkbox) {
        
        subset_data <- gene_plot_data[gene_plot_data$Grouping_L1 == input$selected_grouping_l1, ]
        subset_data <- droplevels(subset_data) # drop unused levels
        subset_data$Grouping_L2 <- factor(subset_data$Grouping_L2, levels = unique(subset_data$Grouping_L2))

        p_val <- NA
        
        if(input$selected_grouping_l1 == "Kidney Tumor Subtypes"){
          colors <- c(
            "#1f77b4", "#e7ba52", "#2ca02c", "#AC6A9F", "#9467bd", "#008080", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
            "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
            "#393b79", "#637939", "#8c6d31", "#5f9ea0"
          )
          
          study_colors <- setNames(colors[1:length(unique(subset_data$Study.Abbreviation))], unique(subset_data$Study.Abbreviation))
        } else if (input$selected_grouping_l1 == "Sarcoma Subtypes"){
          colors <- c(
            "#1f77b4", "#ffbb78", "#2ca02c", "#AC6A9F", "#9467bd", "#008080", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
          )
          
          study_colors <- setNames(colors[1:length(unique(subset_data$Study.Abbreviation))], unique(subset_data$Study.Abbreviation))
        } else {
          colors <- c(
            "#1f77b4", "#c5b0d5", "#5f9ea0"
          )
          study_colors <- setNames(colors[1:length(unique(subset_data$Study.Abbreviation))], unique(subset_data$Study.Abbreviation))
        }
        
        
        
        p1 <- plot_ly(
          data = subset_data,
          x = ~as.factor(Grouping_L2),
          y = ~Chronos_score,
          type = 'box',
          boxpoints = 'all',
          jitter = 0.3,
          pointpos = 0,
          marker = list(
            size = 5, 
            opacity = 0.5),
          color = ~as.factor(Grouping_L2),
          colors = colors,
          text = ~paste('Tumor Subtype:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),  # Tooltip info
          hoverinfo = 'text',               # Show only the text info on hover
          name = paste('All', input$selected_grouping_l1)
        ) 
        
        
        p2 <- plot_ly(
          data = ref_plot_data,
          y = ~Chronos_score,               # Chronos scores on the y-axis for the reference data
          type = 'box',                     # Boxplot type
          boxpoints = 'all',                # Show all data points on the same box
          jitter = 0.3,                     # Add some jitter to prevent points from overlapping
          pointpos = 0,                  # Position the points inside the box
          marker = list(size = 5, opacity = 0.9, color = '#EF8636'),  # Uniform color for the reference plot
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),  # Tooltip info
          hoverinfo = 'text',               # Show only the text info on hover
          color = ~as.factor(Gene),
          colors = "#EF8636",
          name = 'DepMap Mean Chronos Scores'
        )
        
        # Combine with subplot
        subplot(p1, p2, nrows = 1, shareY = TRUE) %>%
          layout(
            title = list(
              text = paste0("Boxplots of Chronos_scores for ", input$selected_gene),
              x = 0.5
            ),
            yaxis = list(title = "Chronos Scores"),
            showlegend = FALSE,
            margin = list(t = 100, b = 100),  # Increased bottom margin for better label spacing
            xaxis = list(
              domain = c(0, 0.94),
              tickangle = 90,  # Rotate x-axis labels
              title = paste("All", input$selected_grouping_l1)
            ),
            xaxis2 = list(
              domain = c(0.95, 1),
              tickangle = 90  # Rotate the label for the last box
            )
          )
        
      } else {
        subset_data <- gene_plot_data[gene_plot_data$Grouping_L2 == input$selected_grouping_l2, ]
        p_val <- suppressWarnings(wilcox.test(subset_data$Chronos_score, ref_plot_data$Chronos_score, alternative="two.sided"))$p.value
        
        # Create the first boxplot for gene_data
        gene_boxplot <- plot_ly(
          data = subset_data,
          y = ~Chronos_score,               # Chronos scores on the y-axis
          type = 'box',                     # Boxplot type
          boxpoints = 'all',                # Show all data points on the same box
          jitter = 0.3,                     # Add some jitter to prevent points from overlapping
          pointpos = 0,                  # Position the points inside the box
          marker = list(size = 7, opacity = 0.5),  # Color points by Tumor_type
          text = ~paste('Type:', Grouping_L2, '<br>ID:', Tumor_ID, '<br>Chronos score:', Chronos_score),  # Tooltip info
          hoverinfo = 'text',               # Show only the text info on hover
          name = paste0(input$selected_grouping_l2, "<br>Predicted Chronos Scores")
        )
        
        # Create the second boxplot for ref_plot_data
        ref_boxplot <- plot_ly(
          data = ref_plot_data,
          y = ~Chronos_score,               # Chronos scores on the y-axis for the reference data
          type = 'box',                     # Boxplot type
          boxpoints = 'all',                # Show all data points on the same box
          jitter = 0.3,                     # Add some jitter to prevent points from overlapping
          pointpos = 0,                  # Position the points inside the box
          marker = list(size = 7, opacity = 0.5),  # Uniform color for the reference plot
          text = ~paste('Cellline Name:', CellLineName, '<br>Chronos score:', Chronos_score),  # Tooltip info
          hoverinfo = 'text',               # Show only the text info on hover
          name = 'DepMap Mean Chronos Scores'
        )
        
        # Combine the two boxplots side by side using subplot
        subplot(gene_boxplot, ref_boxplot, nrows = 1, shareY = TRUE) %>%
          layout(
            title = list(
              text = paste0("Boxplots of Chronos Scores for ", input$selected_gene),
              x=0.5
            ),
            annotations = list(
              list(
                text = paste0("p-value: ", format(p_val, digits = 3, scientific=TRUE)),
                x=1.15,
                y=0.85,
                showarrow = FALSE,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="top",
                font=list(size=16)
              )
            ),
            yaxis = list(title = "Chronos Scores"),
            showlegend = TRUE,              # Show legend for Tumor_type colors
            margin = list(t = 100, b = 100),  # Adjust margins for better spacing
            xaxis = list(domain = c(0, 0.45)), # First boxplot occupies 45% of width
            xaxis2 = list(domain = c(0.55, 1)) # Second boxplot starts after 55%
          )
      }
    }
  })
  
  ################## Tab for Module & Feature Selection 
  output$bar_plot2 <- renderPlotly({
    gene_df <- module_selection[[input$selected_gene_module2]]
    
    # Convert Number_of_features to ordered factor based on numeric values
    sub_df_sorted <- gene_df %>%
      mutate(Number_of_features = factor(Number_of_features, 
                                         levels = unique(gene_df$Number_of_features[order(as.numeric(gsub("N", "", gene_df$Number_of_features)))])))
    
    sub_df_sorted$bar_colors <- ifelse(sub_df_sorted$Correl_Coeffs_5CV == max(sub_df_sorted$Correl_Coeffs_5CV),
                                       "Best Performance",
                                       sub_df_sorted$Module)
    
    
    # Define colors
    module_colors <- setNames(c("#aec7e8", "#c5b6d2", "#cedb9c", "#fffacd", "#f7b6d2", "darkred"),
                              c("SVR", "ridge_regression", "linear_SVR", "lasso_regression", "elasticnet_regression", "Best Performance"))
    
    if (!is.null(input$highest_bar) && input$highest_bar){
      plot_ly(
        data = sub_df_sorted,
        x = ~Number_of_features,
        y = ~Correl_Coeffs_5CV,
        color = ~Module,
        colors = module_colors,
        type = "bar",
        marker = list(opacity = 0.3),
        textposition = "none",
        text = ~paste("Model:", Module, "<br>Number of Features:", gsub("N", "", Number_of_features), "<br>Correlation:", round(Correl_Coeffs_5CV, 3)),
        hoverinfo = "text"
      ) %>%
        add_trace(
          data = sub_df_sorted[sub_df_sorted$Correl_Coeffs_5CV == max(sub_df_sorted$Correl_Coeffs_5CV), ],
          x = ~Number_of_features,
          y = ~Correl_Coeffs_5CV,
          type = "bar",
          mode = "markers",
          marker = list(color = "darkred", opacity = 1),
          textposition = "none",
          text = ~paste("Model:", Module, "<br>Number of Features:", gsub("N", "", Number_of_features), "<br>Correlation:", round(Correl_Coeffs_5CV, 3)),
          hoverinfo = "text",
          showlegend = TRUE,
          name = "Best Performance"
        ) %>%
        layout(
          title = "Grouped Bar Plot: Model Performance",
          xaxis = list(title = "Number of Features", categoryorder = "array", categoryarray = levels(sub_df_sorted$Number_of_features)),
          yaxis = list(title = "5-fold CV Correlation of Predicted vs. Observed Chronos Score in Test Data (R)"),
          barmode = "group"
        )
      
    } else {
      # Create grouped bar plot
      plot_ly(sub_df_sorted,
              x = ~Number_of_features,
              y = ~Correl_Coeffs_5CV,
              color = ~Module,
              colors = module_colors,
              type = "bar",
              marker = list(opacity = 0.7),
              textposition = "none",
              text = ~paste("Model:", Module, "<br>Number of Features:", gsub("N", "", Number_of_features), "<br>Correlation:", round(Correl_Coeffs_5CV, 3)),
              hoverinfo = "text") %>%
        layout(
          title = "Grouped Bar Plot: Model Performance",
          xaxis = list(title = "Number of Features", categoryorder = "array", categoryarray = levels(sub_df_sorted$Number_of_features)),
          yaxis = list(title = "5-fold CV Correlation of Predicted vs. Observed Chronos Score in Test Data (R)"),
          barmode = "group"
        )
    }
    
  })
  
  
  
  ################## Tab for feature importance
  output$feature_plot <- renderPlotly({
    
    feature_df <- features_df[[input$selected_gene_feature]]
    
    # Separate data for gray and red points
    gray_points <- feature_df[feature_df$Rank > 10, ]
    red_point <- feature_df[feature_df$Rank <= 10, ]
    
    plot_ly() %>%
      add_trace(
        data = gray_points,
        type = "scatter",
        mode = "markers",
        x = ~Rank,
        y = ~abs_coef,
        marker = list(color = "lightgray", size = 5),
        textposition = "outside",
        text = ~paste('Feature:', Features, '<br>Coefficient:', Coefficients, '<br>Rank:', Rank),
        hoverinfo = 'text',               # Show only the text info on hover
        name = "Lowly correlated Features"
      ) %>%
      add_trace(
        data = red_point,
        type = "scatter",
        mode = "markers",
        x = ~Rank,
        y = ~abs_coef,
        marker = list(color = "darkred", size = 10),
        textposition = "outside",
        text = ~paste('Feature:', Features, '<br>Coefficient:', Coefficients, '<br>Rank:', Rank),
        hoverinfo = 'text',               # Show only the text info on hover
        name = "Top correlated Features"
      ) %>%
      layout(
        title = paste0(input$selected_gene_feature, " Chronos Score Predictors"),
        xaxis = list(title = "Feature Rank"),
        yaxis = list(title = "Relative Feature Importance"),
        showlegend = T
      )
  })
  
  
  
  ################## Additional Tab for r
  output$scatter_plot <- renderPlotly({
    models_df_sorted <- models_df[order(-models_df$Correl_Coeffs_5CV), ]
    models_df_sorted$Rank <- 1:nrow(models_df_sorted)
    
    # Separate data for gray and red points
    gray_points <- models_df_sorted[models_df_sorted$Gene_Name != input$selected_gene_scatter, ]
    red_point <- models_df_sorted[models_df_sorted$Gene_Name == input$selected_gene_scatter, ]
    
    plot_ly() %>%
      add_trace(
        data = gray_points,
        type = "scatter",
        mode = "markers",
        x = ~Rank,
        y = ~Correl_Coeffs_5CV,
        marker = list(color = "lightgray", size = 5),
        textposition = "outside",
        text = ~paste('Gene:', Gene_Name, '<br>r:', Correl_Coeffs_5CV, '<br>Rank:', Rank),
        hoverinfo = 'text',               # Show only the text info on hover
        name = "Unselected Genes"
      ) %>%
      add_trace(
        data = red_point,
        type = "scatter",
        mode = "markers",
        x = ~Rank,
        y = ~Correl_Coeffs_5CV,
        marker = list(color = "darkred", size = 15),
        textposition = "outside",
        text = ~paste('Gene:', Gene_Name, '<br>r:', Correl_Coeffs_5CV, '<br>Rank:', Rank),
        hoverinfo = 'text',               # Show only the text info on hover
        name = input$selected_gene_scatter
      ) %>%
      layout(
        title = "Scatterplot of Averaged 5-Fold Cross-Validation Correlation Coefficients",
        xaxis = list(title = "Rank"),
        yaxis = list(title = "Correlation of Predicted vs. Observed Chronos Score in Test Data (R)"),
        showlegend = T
      )
  })
  
  output$dataSources_table <- renderDT({
    datatable(dataSources_df, escape = FALSE, options = list(pageLength=5))
  })
  
}


# Run the Shiny app
shinyApp(ui = ui, server = server)


