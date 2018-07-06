shinyUI(
 navbarPage("scgraphs",
        tabPanel("Generate Scaffold map",
            uiOutput("scaffoldUI")
        ),
        tabPanel("Generate unsupervised graph",
            uiOutput("unsupervisedUI")
        )
    )

)