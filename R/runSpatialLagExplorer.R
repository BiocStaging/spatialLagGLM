## Provenance: Gemini (original), Claude (revision). See inst/PROVENANCE.md.
#' Launch Interactive Spatial Lag Explorer
#'
#' @param se A \code{SpatialExperiment} analysed with
#'     \code{\link{spatialLagGLM}}.
#' @return Called for its side effect of starting a Shiny app; returns
#'     nothing useful.
#' @export
#' @examples
#' if (interactive()) {
#'     se <- makeMockSpatial()
#'     analyzed <- spatialLagGLM(se, design = ~condition, radius = 2.5)
#'     runSpatialLagExplorer(analyzed)
#' }
runSpatialLagExplorer <- function(se) {
    .checkSpatialExperiment(se)
    if (!requireNamespace("shiny", quietly = TRUE)) {
        stop("Please install the 'shiny' package to use this explorer.",
            call. = FALSE
        )
    }

    ui <- shiny::fluidPage(
        shiny::titlePanel("spatialLagGLM Interactive Board"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::selectInput("gene", "Select gene:",
                    choices = rownames(se)
                ),
                shiny::selectInput("metric", "Coordinate metric:",
                    choices = c(
                        "Raw expression" = "expression",
                        "Spatial lag" = "spatial_lag"
                    )
                )
            ),
            shiny::mainPanel(
                shiny::plotOutput("spatialPlot"),
                shiny::plotOutput("volcanoPlot")
            )
        )
    )

    server <- function(input, output, session) {
        output$spatialPlot <- shiny::renderPlot({
            plotSpatialLag(se, gene_id = input$gene, metric = input$metric)
        })
        output$volcanoPlot <- shiny::renderPlot({
            plotLagVolcano(se)
        })
    }
    shiny::runApp(shiny::shinyApp(ui, server))
}
