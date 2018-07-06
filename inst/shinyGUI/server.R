

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
                )
            ),
            column(6,
                p("Select landmarks data (select any file in the directory, and all the other FCS files will be loaded as well)"),
                actionButton("scaffoldui_select_landmarks_dir", "Select Landmarks directory"),
                verbatimTextOutput("scaffoldui_landmarks_dir"),
                fluidRow(
                    column(6,
                        checkboxInput("scaffoldui_transform_landmarks_data", "Transform landmarks data", value = TRUE)
                    ),
                    column(6, 
                        conditionalPanel(condition = "input.scaffoldui_transform_landmarks_data",
                            numericInput("scaffoldui_asinh_cofactor", "Cofactor for asinh transformation", 5)
                        )
                    )
                )
            )
        ),
        fluidRow(
            column(6,
                checkboxInput("scaffoldui_inter_cluster_connections", "Add inter-cluster connections", value = TRUE),
                conditionalPanel(
                    condition = "input.scaffoldui_inter_cluster_connections == true",
                    selectInput("scaffoldui_markers_inter_cluster", "Markers for inter-cluster connections (if different)", choices = c(""), multiple = T, width = "100%"), 
                    numericInput("scaffoldui_inter_cluster_weight", "Weight factor for inter-cluster connections", 0.7, min = 0, max = 10, step = 0.1)
                ),
                selectInput("scaffoldui_overlap_method", "Overlap resolution method", choices = c("Repel", "Expand")),
                actionButton("scaffoldui_start_analysis", "Start analysis")
            )
        )
    )
})}

shinyServer(function(input, output, session) {
    working.directory <- dirname(file.choose())
    output$scaffoldUI <- render_scaffold_ui(working.directory, input, output, session)
    
    scaffoldui.reactive.values <- reactiveValues(landmarks.dir = NULL)

    observe({
        if(!is.null(input$scaffoldui_reference) && input$scaffoldui_reference != "") {
                tab <- read.table(paste(working.directory, input$scaffoldui_reference, sep = "/"), header = T, sep = "\t", check.names = F)
                updateSelectInput(session, "scaffoldui_markers", choices = names(tab))
                updateSelectInput(session, "scaffoldui_markers_inter_cluster", choices = names(tab))
            }
    })


    observeEvent(input$scaffoldui_select_landmarks_dir, {
        scaffoldui.reactive.values$landmarks.dir <- dirname(file.choose())
    })

    output$scaffoldui_landmarks_dir <- renderText({
        if(is.null(scaffoldui.reactive.values$landmarks.dir))
            " "
        else
            scaffoldui.reactive.values$landmarks.dir
    })

    observeEvent(input$scaffoldui_start_analysis, {
        isolate({

            if(is.null(scaffoldui.reactive.values$landmarks.dir)) {
                showModal(modalDialog(
                    "Please select a directory for landmarks data"
                ))
                return(NULL)
            }

            showModal(modalDialog(
                "Analysis started, please wait..."
            ))

            files.list <- list.files(working.directory, pattern = "*.clustered.txt$", full.names = TRUE)
            landmarks.data <- load_landmarks_from_dir(scaffoldui.reactive.values$landmarks.dir, 
                input$scaffoldui_asinh_cofactor, input$scaffoldui_transform_landmarks_data)

            inter.cluster.col.names <- NULL
            
            if(input$scaffoldui_inter_cluster_connections) {
                if(length(input$scaffoldui_markers_inter_cluster) > 0)
                    inter.cluster.col.names <- input$scaffoldui_markers_inter_cluster
                else
                    inter.cluster.col.names <- input$scaffoldui_markers
            }


            out.dir <- file.path(working.directory, "scaffold_result")
            args.list <- list(
                    files.list = files.list,
                    ref.file = input$scaffoldui_reference,
                    landmarks.data = landmarks.data,
                    col.names = input$scaffoldui_markers,
                    inter.cluster.weight.factor = input$scaffoldui_inter_cluster_weight,
                    inter.cluster.col.names = inter.cluster.col.names,
                    overlap.method = tolower(input$scaffoldui_overlap_method),
                    out.dir = out.dir
            )


            if(input$scaffoldui_ew_influence_type == "Fixed")            
                args.list <- c(args.list, list(ew.influence = input$scaffoldui_ew_influence))
            
            do.call(scgraphs::run_scaffold_analysis, args.list)

            showModal(modalDialog(
                "Analysis Finished", br(),
                sprintf("Output files are located in %s", out.dir)
            ))
        })
    })
})