# app.R
library(shiny)
library(ggplot2)
library(digest)
library(VaRES)

set.seed(144)

#kernel
k_se <- function(x1, x2, lambda, sigma) {
  outer(x1, x2, function(a, b)
    (sigma^2) * exp(-(a - b)^2 / (2 * lambda^2))
  )
}

#GP prior 
simulate_gp <- function(x, lambda, sigma) {
  K <- k_se(x, x, lambda, sigma)
  L <- chol(K + 1e-6 * diag(length(x)))
  drop(t(L) %*% rnorm(length(x)))
}

ui <- fluidPage(
  
  titlePanel("Gaussian Process Prior Simulator (Squared Exponential Kernel)"),
  
  sidebarLayout(
    
    sidebarPanel(
      sliderInput("ig_alpha", "Inverse-Gamma alpha:", 
                  min = 1, max = 15, value = 1, step = 1),
      
      sliderInput("ig_beta", "Inverse-Gamma beta:", 
                  min = 1, max = 15, value = 1, step = 1),
      
      sliderInput("ht_mu", "Half-t mean:", 
                  min = 0, max = 15, value = 0, step = 1),
      
      sliderInput("ht_df", "Half-t degrees of freedom:", 
                  min = 1, max = 5, value = 4, step = 1),
      
      sliderInput("ht_scale", "Half-t scale:", 
                  min = 1, max = 15, value = 1, step = 1),
      
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
    
    seed_val <- digest(
      list(input$ig_alpha, input$ig_beta,
           input$ht_mu, input$ht_df, input$ht_scale, input$nfunc),
           #, input$npoints, input$xmax),
      algo = "xxhash32",
      serialize = TRUE
    ) |> substr(1, 7) |> strtoi(base = 16)
    
    set.seed(seed_val)
    
    x <- seq(0, input$xmax, length.out = input$npoints)
    
    lambda_draw   <- 1 / rgamma(1, shape = input$ig_alpha, rate = input$ig_beta) 
    sigma_draw <- input$ht_mu + input$ht_scale * varhalfT(runif(1), n = input$ht_df)
    
    funcs <- replicate(
      input$nfunc,
      simulate_gp(x, lambda_draw, sigma_draw)  
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
