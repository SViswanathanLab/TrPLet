# R shiny for plots
if (!requireNamespace('shiny', quietly = TRUE)) {
  install.packages('shiny')
}
if (!requireNamespace('plotly', quietly = TRUE)) {
  install.packages('plotly')
}
if (!requireNamespace('randomcoloR', quietly = TRUE)) {
  install.packages('randomcoloR')
}
library(shiny)
library(plotly)
library("randomcoloR") 

set.seed(1234)
color_panels <- distinctColorPalette(23)


# load data
plot_data <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657dfs_merged_across_23TumorTypes.rds")
genes <- readRDS("/Volumes/sviswanathan/users/ycui/TrPLet_Table_S4/0915_Boxplots/657gene_names.rds")


# Define UI for the dashboard
ui <- fluidPage(
  headerPanel('Interactive Boxplots for Table S4 (select genes)'),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_gene", "Select a gene:", choices = genes, selected = genes[1])
    ),
    
    mainPanel(
      plotlyOutput("boxplot_plot", height = "1200px", width = "3000px")
    )
  )
)


# Define server logic for generating the single boxplot
server <- function(input, output) {
  
  
  output$boxplot_plot <- renderPlotly({
    

    # Get the data for the selected gene
    gene_plot_data <- plot_data[[input$selected_gene]]
    gene_plot_data$Tumor_type <- as.factor(gene_plot_data$Tumor_type)
      
      
    # Create the interactive boxplot using plotly
    plot_ly(
      data = gene_plot_data,
      x = ~Tumor_type, 
      y = ~Chromos_score,               # Chromos scores on the y-axis
      type = 'box',                     # Boxplot type
      boxpoints = 'all',                # Show all data points on the same box
      jitter = 0.3,                     # Add some jitter to prevent points from overlapping
      pointpos = 0,                     # Position the points inside the box
      marker = list(size = 8, opacity = 0.8),  # Color points by Tumor_type
      color = ~Tumor_type,
      colors = color_panels,
      text = ~paste('Tumor type:', Tumor_type, '<br>Tumor ID:', Tumor_ID, '<br>Chromos score:', Chromos_score),  # Tooltip info
      hoverinfo = 'text',         # Show only the text info on hover
      name = ~Tumor_type
    ) %>%
      layout(
        title = paste0("Boxplot of Chromos_scores for ", input$selected_gene),
        yaxis = list(title = "Chromos Scores"),
        showlegend = TRUE,              # Show legend for Tumor_type colors
        margin = list(t = 60, b = 60)   # Adjust margins for better spacing
      )
  })
}



# Run the Shiny app
shinyApp(ui = ui, server = server)

