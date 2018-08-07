function(input, output, session) {


  output$plot1 <- renderPlot({
    STMotif::intensityDataset(datasetShiny,listMotifShiny,input$rank,alphaShiny)
  })

  output$plot2 <- renderPlot({
    STMotif::displayPlotSeries(datasetShiny, listMotifShiny, input$rank2, input$range[1]:input$range[2])
  })

}
