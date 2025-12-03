# app.R
library(shiny)
library(ggplot2)
library(digest)
library(VaRES)
library(statmod)

set.seed(144)

#kernel
k_se <- function(x1, x2, lambda, sigma) {
  outer(x1, x2, function(a, b)
    (sigma^2) * exp(-(a - b)^2 / (2 * lambda^2))
  )
}

simulate_gp <- function(x, lambda, sigma, mean_fun = function(x) 0) {
  K <- k_se(x, x, lambda, sigma)
  L <- chol(K + 1e-6 * diag(length(x)))
  m <- mean_fun(x)
  f <- m + t(L) %*% rnorm(length(x))
  drop(f)
}


ui <- fluidPage(
  
  titlePanel("Gaussian Process Prior Simulator (Squared Exponential Kernel)"),
  
  sidebarLayout(
    
    sidebarPanel(
      #sliderInput("lambda", "lambda:", 
      #            min = 0, max = 15, value = 1, step = 1),
      
      #sliderInput("sigma", "sigma:", 
      #            min = 0, max = 50, value = 1, step = 5),

      
      #sliderInput("scale", "scale:", 
       #           min = 0, max = 10, value = 1, step = 1),
      
    #  sliderInput("df", "df:", 
     #             min = 0, max = 10, value = 1, step = 1),     
      
      sliderInput("mean", "mean:", 
                  min = 0, max = 10, value = 1, step = 1),
      
      sliderInput("shape", "shape:", 
                  min = 0, max = 10, value = 1, step = 1),
      
    ),
    
    mainPanel(
      #plotOutput("gpPlot", height = "600px"),
      #plotOutput("plot_ht", height = "600px"),
      plotOutput("plot_invgauss", height = "600px")
    )
  )
)

server <- function(input, output) {

  output$plot_ht <- renderPlot({
    
    x <- seq(0, 15, length.out = 400)
    df <- input$df
    mu <- input$mean
    sc <- input$scale
    
    dens <- 2 * dt((x - mu)/sc, df = df) / sc
    dens[x < mu] <- 0
    d <- data.frame(x=x, y=dens)
    
    ggplot(d, aes(x,y)) +
      geom_line(color="darkgreen", linewidth=1) +
      labs(title="Half-t prior for σ^2",
           y="density", x="σ^2") +
      theme_minimal(base_size=14)
  })
  
  output$plot_invgauss <- renderPlot({
    
    # parameters for inv-gaussian
    mean<-input$mean
    shape<-input$shape 
    
    # grid
    x <- seq(0, 10, length.out = 500)
    
    y <- dinvgauss(x, mean=mean, shape=shape)
    
    plot(x, y, type = "l", lwd = 2,
         main = paste0("Inv-Gaussian (mean = ", mean, ", shape = ", shape, ")"),
         xlab = "x", ylab = "density")
  })
    
  
  gp_data <- eventReactive(c(input$sigma, input$lambda), {
    
    x <- seq(0, 20, length.out = 200)
    
    funcs <- replicate(
      2,
      simulate_gp(x, input$lambda, input$sigma)#, mean_fun = function(x) 10 + 5 * x)  
    )
    
    df <- data.frame(
      x = rep(x, 2),
      f = as.vector(funcs),
      func = rep(1:2,each = length(x))
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
