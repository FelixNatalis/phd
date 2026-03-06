library(shiny)

ui <- fluidPage(
  titlePanel("Show UI on Button Click"),
  
  # Action button to trigger UI appearance
  actionButton("show_btn", "Show Extra Inputs"),
  
  # Placeholder for dynamic UI
  uiOutput("dynamic_ui")
)

server <- function(input, output, session) {
  
  output$dynamic_ui <- renderUI({
    # Only show UI if button has been clicked at least once
    if (input$show_btn > 0) {
      tagList(
        textInput("name", "Enter your name:"),
        numericInput("age", "Enter your age:", value = 30, min = 1, max = 120),
        actionButton("submit", "Submit")
      )
    }
  })
}

shinyApp(ui, server)


# role <- reactive({
#   input$kernel_label
# })
# 
# observeEvent(role(), {
#   if (role() %in% c("Periodic")) {
#     nav_show("nav", "C")
#   } else if (role() == "Linear") {
#     nav_hide("nav", "C")
#   } else {
#     stop(sprintf("user has unexpected role %s", role()))
#   }
# })
