shinyUI(fluidPage(
  titlePanel(title = "Visualization of motifs"),
  hr(),
  br(),
  h3("Plot the intensity of values"),
  br(),
  fluidRow(
    column(12,
           sliderInput("rank", "Choose the rank:",min = 1, max = ifelse(length(listMotifShiny)<1000, length(listMotifShiny), 1000), value = 1, step = 1, width = "auto")
    )
  ),
  fluidRow(
    column(12,
           plotOutput('plot1',width = "auto", height = "800px")
    )
  ),
  h3("Plot the spatial-time series"),
  br(),
  fluidRow(
    column(12,
           sliderInput("rank2", "Choose the rank:",min = 1, max = ifelse(length(listMotifShiny)<1000, length(listMotifShiny), 1000), value = 1, step = 1, width = "auto")
    ),
    column(12,
           sliderInput("range", label = "Choose a range of columns", min = 1,  max = ncol(datasetShiny), value = c(1, 10), step = 1, width = "auto")
    )
  ),
  fluidRow(
    column(12,
           plotOutput('plot2',width = "auto", height = "1000px")
    )
  )
)
)
