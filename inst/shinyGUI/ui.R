shinyUI(
 navbarPage("vite",
        tabPanel("Generate Scaffold map",
            uiOutput("scaffoldUI")
        ),
        tabPanel("Generate unsupervised graph",
            uiOutput("unsupervisedUI")
        )
    )

)