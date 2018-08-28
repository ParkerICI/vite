
render_unsupervised_ui <- function(working.directory, ...) {renderUI({
    files.list <- list.files(path = working.directory, pattern = "*.clustered.txt$")

    fluidPage(
        fluidRow(
            column(6,
                selectInput("unsupervisedui_markers_file", "Load marker names from file", choices = c("", list.files(path = working.directory, pattern = "*.clustered.txt$")), width = "100%"),
                selectInput("unsupervisedui_markers", "Choose the markers to use for the anlaysis", choices = c(""), multiple = T, width = "100%"),
                selectInput("unsupervisedui_files_list", label = "Files to include", choices = files.list, selected = files.list,
                               multiple = T, width = "100%")
            )
        ),
        fluidRow(
            column(6,
                checkboxInput("unsupervisedui_include_metadata", "Include metadata file"),
                conditionalPanel(condition = "input.unsupervisedui_include_metadata",
                    actionButton("unsupervisedui_select_metadata", "Select metadata file"),
                    verbatimTextOutput("unsupervisedui_metadata_file", placeholder = TRUE)
                )
            )
        ),
        fluidRow(
            column(6,
                numericInput("unsupervisedui_filtering_threshold", "Edge filtering threshold", 15, min = 1),
                textInput("unsupervisedui_out_name", "Output file name (the extension .graphml will be added automatically)", width = "100%"),
                actionButton("unsupervisedui_start_analysis", "Start analysis")
            )
        )
    )
})}


render_scaffold_ui <- function(working.directory, ...) {renderUI({
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
                verbatimTextOutput("scaffoldui_landmarks_dir", placeholder = TRUE),
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
    output$unsupervisedUI <- render_unsupervised_ui(working.directory, input, output, session)


    # Unsupervised UI functions
    unsupervisedui.reactive.values <- reactiveValues(metadata.file = NULL)

    output$unsupervisedui_metadata_file <- renderText({
        unsupervisedui.reactive.values$metadata.file
    })

    observeEvent(input$unsupervisedui_select_metadata, {
        unsupervisedui.reactive.values$metadata.file <- file.choose()
    })

    observe({
        if(!is.null(input$unsupervisedui_markers_file) && input$unsupervisedui_markers_file != "") {
                tab <- read.table(file.path(working.directory, input$unsupervisedui_markers_file), header = TRUE, sep = "\t", check.names = FALSE)
                updateSelectInput(session, "unsupervisedui_markers", choices = names(tab))
            }
    })

    observeEvent(input$unsupervisedui_start_analysis, {
        isolate({
            if(is.null(input$unsupervisedui_out_name) || input$unsupervisedui_out_name == "") {
                showModal(modalDialog(
                    "Please enter a name for the output"
                ))
                return(NULL)
            }

            showModal(modalDialog(
                "Analysis started, please wait..."
            ))


            files.list <- file.path(working.directory, input$unsupervisedui_files_list)
            metadata.tab <- NULL

            if(!is.null(unsupervisedui.reactive.values$metadata.file) && unsupervisedui.reactive.values$metadata.file != "")
                metadata.tab <- read.table(unsupervisedui.reactive.values$metadata.file,
                                            header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)


            G <- vite::get_unsupervised_graph_from_files(
                files.list = files.list,
                col.names = input$unsupervisedui_markers,
                filtering.threshold = input$unsupervisedui_filtering_threshold,
                metadata.tab = metadata.tab,
                metadata.filename.col = "filename",
                clusters.data.out.dir = working.directory
            )

            out.name <- file.path(working.directory, sprintf("%s.graphml", input$unsupervisedui_out_name))
            write_graph(G, out.name)

            showModal(modalDialog(
                "Analysis Finished", br(),
                sprintf("The ouptut file is located at are located at %s", out.name)
            ))

        })


    })




    # Scaffold UI functions
    scaffoldui.reactive.values <- reactiveValues(landmarks.dir = NULL)

    observe({
        if(!is.null(input$scaffoldui_reference) && input$scaffoldui_reference != "") {
                tab <- read.table(file.path(working.directory, input$scaffoldui_reference), header = T, sep = "\t", check.names = F)
                updateSelectInput(session, "scaffoldui_markers", choices = names(tab))
                updateSelectInput(session, "scaffoldui_markers_inter_cluster", choices = names(tab))
            }
    })


    observeEvent(input$scaffoldui_select_landmarks_dir, {
        scaffoldui.reactive.values$landmarks.dir <- dirname(file.choose())
    })

    output$scaffoldui_landmarks_dir <- renderText({
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
                    ref.file = file.path(working.directory, input$scaffoldui_reference),
                    landmarks.data = landmarks.data,
                    col.names = input$scaffoldui_markers,
                    inter.cluster.weight.factor = input$scaffoldui_inter_cluster_weight,
                    inter.cluster.col.names = inter.cluster.col.names,
                    overlap.method = tolower(input$scaffoldui_overlap_method),
                    out.dir = out.dir
            )


            if(input$scaffoldui_ew_influence_type == "Fixed")
                args.list <- c(args.list, list(ew.influence = input$scaffoldui_ew_influence))

            do.call(vite::run_scaffold_analysis, args.list)

            showModal(modalDialog(
                "Analysis Finished", br(),
                sprintf("Output files are located in %s", out.dir)
            ))
        })
    })
})
