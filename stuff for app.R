# imports
library(shiny)
library(ggplot2)
library(VaRES)
library(digest)
library(statmod)
library(hash)

# defaults parameters for GP prior predictive draws
n_functions = 10
x_min = -10
x_max = 10
n_points = 200
epsilon = 1e-6

# SE kernel function
squared_exponential_kernel <- function(x1, x2, lambda, sigma_2) {
  outer(x1, x2, function(a, b)
    (sigma_2 * exp(-(a - b)^2 / (2 * lambda^2)))
  )
}

# Linear kernel function
linear_kernel <- function(x1, x2, lambda, sigma_2) {
  outer(x1, x2, function(a, b)
    (sigma_2 * a * b)
  )
}

kernels <- hash(
  "Squared Exponential" = squared_exponential_kernel
  #,"Matérn" = 2, 
  ,"Linear" = linear_kernel#linear_kernel
  #,"Periodic" = 4
) 


## UI

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      tags$h4("Kernel"),
      selectInput("kernel_label", "Choose a kernel:",
                  list(`Simple kernels` = keys(kernels)
                       #, `Kernel combinations` = keys(kernel_combinations)
                  )
      ),
      tags$hr(),
      
      tags$h4("Hyperparameters for λ"),
      sliderInput("ig_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
      sliderInput("ig_beta", "Inverse-Gamma beta:", 1, 15, 1),
      #sliderInput("ig_mu", "Inverse-Gaussian mean:", 1, 15, 1),
      #sliderInput("ig_lambda", "Inverse-Gaussian shape:", 1, 15, 1),
      checkboxInput("lambda_mle", "Use MLE", FALSE),
      actionButton("draw_lambda", "Draw New Lambda"),
      tags$hr(),
      
      tags$h4("Hyperparameters for σ^2"),
      sliderInput("ht_mu", "Half-t mean:", 0, 15, 0),
      sliderInput("ht_df", "Half-t degrees of freedom:", 1, 5, 4),
      sliderInput("ht_scale", "Half-t scale:", 1, 15, 1),
      checkboxInput("sigma_2_mle", "Use MLE", FALSE),
      actionButton("draw_sigma_2", "Draw New Sigma^2"),
      tags$hr(),
    
    ),
    
    mainPanel(
      #actionButton("restart", "Restart Session"),
      fluidRow(
        column(6, plotOutput("plot_ig", height = "250px")),
        #column(6, plotOutput("plot_ig_gauss", height = "250px"))
        column(6, plotOutput("plot_ht", height = "250px"))
      ),
      tags$br(),
      plotOutput("kernelPlot", height = "250px")
      
    )
  )
)

server <- function(input, output) {
  
  # Kernel choice
  kernel_choice <- eventReactive(input$kernel_label, {
    return(kernels[[input$kernel_label]])
  })
  
  
  #Inv-Gamma
  lambda_draw <- eventReactive(input$draw_lambda, { 
    seed_val <- digest(list(input$ig_alpha, input$ig_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    if (isTRUE(input$lambda_mle)) {
      return(input$ig_beta / (input$ig_alpha + 1))   # MLE of InvGamma(α,β)
    } else {
      return(1 / rgamma(1, shape = input$ig_alpha, rate = input$ig_beta))
    }
  })
  
  sigma_2_draw <- eventReactive(input$draw_sigma_2, {
    seed_val <- digest(list(input$ht_mu, input$ht_df, input$ht_scale), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    if (isTRUE(input$sigma_2_mle)) {
      return(input$ht_mu)    # half-t MLE occurs at lower bound = μ
    } else {
      return(input$ht_mu + input$ht_scale * varhalfT(runif(1), n=input$ht_df))
    }
  })
  

  ### plots
  
  # Inverse-Gamma prior for lambda plot
  output$plot_ig <- renderPlot({
    req(lambda_draw())
    
    x <- seq(1e-6, 15, length.out = 400)
    alpha <- input$ig_alpha
    beta  <- input$ig_beta
    
    dens <- (beta^alpha / gamma(alpha)) * x^(-alpha-1) * exp(-beta / x)
    d <- data.frame(x=x, y=dens)
    
    ggplot(d, aes(x,y)) +
      geom_line(color="steelblue", linewidth=1) +
      geom_point(aes(x=lambda_draw(), y=0), color="red", size=3) +
      annotate("text",
               x = max(d$x) * 0.8,
               y = max(d$y) * 0.9,
               label = paste0("λ = ", signif(lambda_draw(),3)),
               color="red",
               size = 5,
               hjust = 0) +
      labs(title="Inverse-Gamma prior for λ",
           y="density", x="λ") +
      theme_minimal(base_size=14)
  })
  
  # Inverse-Gaussian prior for lambda plot
  #output$plot_ig_gauss <- renderPlot({
  #  req(lambda_draw_gauss())
  #  
  #  x <- seq(1e-6, 50, length.out = 400)
  #  
  #  dens <- dinvgauss(x, mean=input$ig_mu, shape=input$ig_lambda)
  #  d <- data.frame(x=x, y=dens)
  #  
  #  ggplot(d, aes(x,y)) +
  #    geom_line(color="steelblue", linewidth=1) +
  #    geom_point(aes(x=lambda_draw_gauss(), y=0), color="red", size=3) +
  #    annotate("text",
  #             x = max(d$x) * 0.8,
  #             y = max(d$y) * 0.9,
  #             label = paste0("λ = ", signif(lambda_draw_gauss(),3)),
  #             color="red",
  #             size = 5,
  #             hjust = 0) +
  #    labs(title="Inverse-Gaussian prior for λ",
  #         y="density", x="λ") +
  #    theme_minimal(base_size=14)
  #})
  
  # Half-t prior for sigma plot
  output$plot_ht <- renderPlot({
    req(sigma_2_draw())
    
    x <- seq(0, 15, length.out = 400)
    df <- input$ht_df
    mu <- input$ht_mu
    sc <- input$ht_scale
    
    dens <- 2 * dt((x - mu)/sc, df = df) / sc
    dens[x < mu] <- 0
    d <- data.frame(x=x, y=dens)
    
    ggplot(d, aes(x,y)) +
      geom_line(color="darkgreen", linewidth=1) +
      geom_point(aes(x=sigma_2_draw(), y=0), color="red", size=3) +
      annotate("text",
               x = max(d$x) * 0.8,
               y = max(d$y) * 0.9,
               label = paste0("σ^2 = ", signif(sigma_2_draw(),3)),
               color="red",
               size = 5,
               hjust = 0) +
      labs(title="Half-t prior for σ^2",
           y="density", x="σ^2") +
      theme_minimal(base_size=14)
  })
  
  # Kernel based on distance plot
  output$kernelPlot <- renderPlot({
    req(kernel_choice(), lambda_draw(), sigma_2_draw())
    
    dist <- seq(-3, 3, length.out = 300)
    x_o <- rep(0, length(dist))
    lambda <- lambda_draw()
    sigma_2  <- sigma_2_draw()
    
    
    k <- kernel_choice()(dist, x_o, lambda, sigma_2)
    
    ggplot(data.frame(dist=dist, k=k), aes(dist,k)) +
      geom_line(color="purple", linewidth=1.2) +
      labs(title=input$kernel_label, x="Distance", y="k(x)") +
      theme_minimal(base_size=14)
  })
  
}

shinyApp(ui, server)
