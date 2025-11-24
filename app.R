# app.R
library(shiny)
library(ggplot2)

# squared exponential kernel
k_se <- function(x1, x2, ell, sigma2) {
  outer(x1, x2, function(a, b)
    sigma2 * exp(-(a - b)^2 / (2 * ell^2))
  )
}

# GP prior 
simulate_gp <- function(x, ell, sigma2) {
  K <- k_se(x, x, ell, sigma2)
  L <- chol(K + 1e-6 * diag(length(x)))
  drop(t(L) %*% rnorm(length(x)))
}

ui <- fluidPage(
  
  titlePanel("Gaussian Process Prior Simulator (Squared Exponential Kernel)"),
  
  sidebarLayout(
    
    sidebarPanel(
      sliderInput("ell", "Length Scale (ℓ):",
                  min = 0.1, max = 5, value = 1, step = 0.1),
      
      sliderInput("sigma2", "Variance (σ²):",
                  min = 0.1, max = 5, value = 1, step = 0.1),
      
      sliderInput("nfunc", "Number of Functions to Draw:",
                  min = 1, max = 20, value = 5, step = 1),
      
      sliderInput("npoints", "Number of X Points:",
                  min = 20, max = 400, value = 200, step = 10),
      
      sliderInput("xmax", "X Range (0 to Xmax):",
                  min = 2, max = 20, value = 10, step = 1),
      
      actionButton("draw", "Draw New Samples", class = "btn-primary")
    ),
    
    mainPanel(
      plotOutput("gpPlot", height = "600px")
    )
  )
)

server <- function(input, output) {
  
  gp_data <- eventReactive(input$draw, {
    
    x <- seq(0, input$xmax, length.out = input$npoints)
    
    funcs <- replicate(
      input$nfunc,
      simulate_gp(x, input$ell, input$sigma2)
    )
    
    df <- data.frame(
      x = rep(x, input$nfunc),
      f = as.vector(funcs),
      func = rep(1:input$nfunc, each = length(x))
    )
    
    df
  })
  
  output$gpPlot <- renderPlot({
    req(gp_data())
    
    ggplot(gp_data(), aes(x = x, y = f, group = func, color = factor(func))) +
      geom_line(alpha = 0.9, linewidth = 1) +
      scale_color_discrete(guide = "none") +
      labs(
        title = "Gaussian Process Prior Samples",
        subtitle = "Squared Exponential Kernel",
        x = "x", y = "f(x)"
      ) +
      theme_minimal(base_size = 16)
  })
}

shinyApp(ui, server)
