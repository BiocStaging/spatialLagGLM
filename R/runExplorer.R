#' Launch Interactive Spatial Data Explorer
#'
#' @param se An analyzed SpatialExperiment object.
#' @return Opens an interactive Shiny dashboard page in the web browser. None returned.
#' @export
#' @examples
#' \dontrun{
#' se <- make_mock_spatial()
#' analyzed <- spatialDESeq(se, design = ~ condition, radius = 2.5)
#' runExplorer(analyzed)
#' }
runExplorer <- function(se) {
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Please install the 'shiny' package to use this explorer.")
    }
    gene_list <- rownames(se)
    
    ui <- shiny::fluidPage(
        shiny::titlePanel("spatialDESeq Interactive Board"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::selectInput("gene", "Select Target Gene:", choices = gene_list),
                shiny::selectInput("metric", "Coordinate Visual Metric:", choices = c("Raw Expression" = "expression", "Spatial Effect" = "spatial_effect"))
            ),
            shiny::mainPanel(shiny::plotOutput("spatialPlot"), shiny::plotOutput("volcanoPlot"))
        )
    )
    
    server <- function(input, output, session) {
        output$spatialPlot <- shiny::renderPlot({ plotSpatialDESeq(se, gene_id = input$gene, metric = input$metric) })
        output$volcanoPlot <- shiny::renderPlot({ plotVolcano(se) })
    }
    shiny::runApp(list(ui = ui, server = server))
}
