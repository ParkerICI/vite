

render_scaffold_ui <- function(working.directory, ...){renderUI({
    fluidPage(
        fluidRow(
            column(6,
                selectInput("scaffoldui_reference", "Choose a reference dataset:", choices = c("", list.files(path = working.directory, pattern = "*.clustered.txt$")), width = "100%"),
                selectInput("scaffoldui_markers", "Choose the markers for SCAFFoLD", choices = c(""), multiple = T, width = "100%"),
                selectInput("scaffoldui_ew_influence_type", "Edge weight influence", choices = c("Proportional", "Fixed"), width = "100%"),
                conditionalPanel(
                    condition = "input.scaffoldui_ew_influence_type == 'Fixed'",
                    numericInput("scaffoldui_ew_influence", "Specifiy Edge weight value", 12)
                ),
                checkboxInput("scaffoldui_inter_cluster_connections", "Add inter-cluster connections", value = TRUE),
                conditionalPanel(
                    condition = "input.scaffoldui_inter_cluster_connections == true",
                    selectInput("scaffoldui_markers_inter_cluster", "Markers for inter-cluster connections (if different)", choices = c(""), multiple = T, width = "100%"), 
                    numericInput("scaffoldui_inter_cluster_weight", "Weight factor for inter-cluster connections", 0.7, min = 0, max = 10, step = 0.1)
                ),
                numericInput("scaffoldui_asinh_cofactor", "asinh cofactor", 5),
                actionButton("scaffoldui_start_analysis", "Start analysis")
            )
        )
    )
})}

shinyServer(function(input, output, session) {
    working.directory <- dirname(file.choose())
    output$scaffoldUI <- render_scaffold_ui(working.directory, input, output, session)
    
    observe({
        if(!is.null(input$scaffoldui_reference) && input$scaffoldui_reference != "") {
                tab <- read.table(paste(working.directory, input$scaffoldui_reference, sep = "/"), header = T, sep = "\t", check.names = F)
                updateSelectInput(session, "scaffoldui_markers", choices = names(tab))
                updateSelectInput(session, "scaffoldui_markers_inter_cluster", choices = names(tab))
            }
    })


    observeEvent("scaffoldui_start_analysis", {



    })

    output$scaffoldui_empty <- renderText({
        if(!is.null(input$scaffoldui_start) && input$scaffoldui_start != 0)
            isolate({
                    if(!is.null(input$scaffoldui_reference) && input$scaffoldui_reference != "" &&
                        !is.null(input$scaffoldui_markers) && length(input$scaffoldui_markers) > 0)
                    {
                        files.analyzed <- NULL
                        ew_influence <- NULL
                        if(!is.null(input$scaffoldui_ew_influence_type)
                                && input$scaffoldui_ew_influence_type == 'Fixed')
                        {
                            if(!is.null(input$scaffoldui_ew_influence))
                                ew_influence <- input$scaffoldui_ew_influence
                        }
                        
                        
                        if(input$scaffoldui_mode == "Gated")
                        {
                            files.analyzed <- scaffold:::run_scaffold_gated(working.directory, input$scaffoldui_reference,
                                input$scaffoldui_markers, inter.cluster.connections = input$scaffoldui_inter_cluster_connections, col.names.inter_cluster = input$scaffoldui_markers_inter_cluster,
                                asinh.cofactor = input$scaffoldui_asinh_cofactor, ew_influence = ew_influence, inter_cluster.weight_factor = input$scaffoldui_inter_cluster_weight, overlap_method = "repel")
                        }
                        else if(input$scaffoldui_mode == "Existing")
                        {
                            files.analyzed <- scaffold:::run_scaffold_existing(working.directory, input$scaffoldui_reference,
                                input$scaffoldui_markers, inter.cluster.connections = input$scaffoldui_inter_cluster_connections, ew_influence = ew_influence)
                        }
                        if(input$scaffoldui_mode == "Unsupervised")
                        {
                            files.analyzed <- scaffold:::run_scaffold_unsupervised(working.directory, input$scaffoldui_reference,
                                input$scaffoldui_markers, inter.cluster.connections = input$scaffoldui_inter_cluster_connections, ew_influence = ew_influence)
                        }
                    }
                    updateSelectInput(session, "graphui_dataset", choices = c("", list.files(path = working.directory, pattern = "*.scaffold$")))
                    ret <- sprintf("Analysis completed with markers %s\n", paste(input$scaffoldui_markers, collapse = " "))
                    ret <- paste(ret, sprintf("Files analyzed:\n%s", paste(files.analyzed, collapse = "\n")), sep = "")
                    return(ret)
            })
            
            
    })
})