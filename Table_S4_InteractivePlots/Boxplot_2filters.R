set.seed(1234)


# R shiny for plots
if (!requireNamespace('shiny', quietly = TRUE)) {
  install.packages('shiny')
}
if (!requireNamespace('plotly', quietly = TRUE)) {
  install.packages('plotly')
}
library(shiny)
library(plotly)

# load data
plot_data <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657dfs_merged_across_23TumorTypes.rds")
genes <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657gene_names.rds")

ref_data <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657dfs_DepMap.rds")



# Define UI for the dashboard
ui <- fluidPage(
  headerPanel('Interactive Boxplots for Table S4 (select genes and tumor types)'),
  
  sidebarLayout(
    sidebarPanel(
      # Dropdown for gene selection
      selectInput("selected_gene", "Select a gene:", choices = genes, selected = genes[1]),
      
      # Dropdown for tumor type selection (initially populated based on the first gene)
      selectInput("selected_tumor_type", "Select a tumor type:", 
                  choices = unique(plot_data[[genes[1]]]$Tumor_type), 
                  selected = unique(plot_data[[genes[1]]]$Tumor_type)[1])
    ),
    
    mainPanel(
      plotlyOutput("boxplot_plot", height = "1200px", width = "900px")
    )
  )
)


# Define server logic for generating the boxplot
server <- function(input, output, session) {
  
  # Update the tumor type dropdown dynamically based on the selected gene
  observe({
    updateSelectInput(session, "selected_tumor_type",
                      choices = unique(plot_data[[input$selected_gene]]$Tumor_type),
                      selected = unique(plot_data[[input$selected_gene]]$Tumor_type)[1])
  })
  
  output$boxplot_plot <- renderPlotly({
    
    # Get the data for the selected gene
    gene_plot_data <- plot_data[[input$selected_gene]]
    gene_plot_data$Tumor_type <- as.factor(gene_plot_data$Tumor_type)
      
    ref_plot_data <- ref_data[[input$selected_gene]]
      
    # Create the first boxplot for gene_data
    gene_boxplot <- plot_ly(
      data = gene_plot_data[gene_plot_data$Tumor_type == input$selected_tumor_type, ],
      y = ~Chromos_score,               # Chromos scores on the y-axis
      type = 'box',                     # Boxplot type
      boxpoints = 'all',                # Show all data points on the same box
      jitter = 0.3,                     # Add some jitter to prevent points from overlapping
      pointpos = 0,                  # Position the points inside the box
      marker = list(size = 7, opacity = 0.5),  # Color points by Tumor_type
      text = ~paste('Tumor type:', Tumor_type, '<br>Tumor ID:', Tumor_ID, '<br>Chromos score:', Chromos_score),  # Tooltip info
      hoverinfo = 'text',               # Show only the text info on hover
      name = ~paste0('All Tumors in ', input$selected_tumor_type)
    )
      
      # Create the second boxplot for ref_plot_data
    ref_boxplot <- plot_ly(
      data = ref_plot_data,
      y = ~Chromos_score,               # Chromos scores on the y-axis for the reference data
      type = 'box',                     # Boxplot type
      boxpoints = 'all',                # Show all data points on the same box
      jitter = 0.3,                     # Add some jitter to prevent points from overlapping
      pointpos = 0,                  # Position the points inside the box
      marker = list(size = 7, opacity = 0.5),  # Uniform color for the reference plot
      text = ~paste('Model ID:', Model_ID, '<br>Chromos score:', Chromos_score),  # Tooltip info
      hoverinfo = 'text',               # Show only the text info on hover
      name = 'DepMap_Public_23Q2'
    )
      
    # Combine the two boxplots side by side using subplot
    subplot(gene_boxplot, ref_boxplot, nrows = 1, shareY = TRUE) %>%
      layout(
        title = paste0("Boxplots of Chromos_scores for ", input$selected_gene),
        yaxis = list(title = "Chromos Scores"),
        showlegend = TRUE,              # Show legend for Tumor_type colors
        margin = list(t = 60, b = 60)  # Adjust margins for better spacing
      )
  })
}


# Run the Shiny app
shinyApp(ui = ui, server = server)
















