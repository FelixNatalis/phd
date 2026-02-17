## Only run examples in interactive R sessions
if (interactive()) {
 
  
  # demoing group support in the `choices` arg
  shinyApp(
    ui = fluidPage(
      selectInput("state", "Choose a state:",
                  list(`East Coast` = list("NY", "NJ", "CT"),
                       `West Coast` = list("WA", "OR", "CA"),
                       `Midwest` = list("MN", "WI", "IA"))
      ),
      textOutput("result")
    ),
    server = function(input, output) {
      output$result <- renderText({
        paste("You chose", input$state)
      })
    }
  )
}
