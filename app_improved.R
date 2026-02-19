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

# Kernels

# SE kernel function
squared_exponential_kernel <- function(x1, x2, length_scale = 1, variance = 1, roughness = 2.5, period = 1) {
  outer(x1, x2, function(a, b)
    (variance * exp(-(a - b)^2 / (2 * length_scale^2)))
  )
}

# Linear kernel function
linear_kernel <- function(x1, x2, length_scale = 1, variance = 1, roughness = 2.5, period = 1) {
  outer(x1, x2, function(a, b)
    (variance * a * b)
  )
}

# Matérn kernel function
matern_kernel <- function(x1, x2, length_scale = 1, variance = 1, roughness = 1.5, period = 1) {
  outer(x1, x2, function(a, b){
    length_scale = 1
    variance = 1
    roughness = 1.5
    distance <- (a - b)^2

    term <- sqrt(2 * roughness) * ((a - b)^2 + 1e-5) / length_scale^2
    matrix<- variance * 2^(1 - roughness) / gamma(roughness) * term^roughness * besselK(term, roughness)
    matrix[distance == 0] <- 1
    return (matrix)
  }
  )
}


kernels <- hash(
  "Squared Exponential" = squared_exponential_kernel
  ,"Matérn" = matern_kernel
  ,"Linear" = linear_kernel
  #,"Periodic" = 4
  ) 

# kernel_combinations = hash(
#           "Changepoint" = 5,
#           "Lin + SE" = 6,
#           "Lin * SE" = 7,
#           "Per + SE" = 8,
#           "Per * SE" = 9
# )

# Drawing GP
simulate_gp <- function(x, kernel, length_scale, variance, sigma_noise = 1e-3, mean_fun = function(x) 0) {
  
  K <- kernel(x, x, length_scale, variance)
  L <- chol(K + epsilon * diag(length(x)))
  
  m <- mean_fun(x)
  
  f <- m + t(L) %*% rnorm(length(x))
  
  # noise
  eps <- sigma_noise * rnorm(length(x))
  
  drop(f + eps)
}

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
      
      tags$h4("Hyperparameters for length scale"),
      sliderInput("ig_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
      sliderInput("ig_beta", "Inverse-Gamma beta:", 1, 15, 1),
      #sliderInput("ig_mu", "Inverse-Gaussian mean:", 1, 15, 1),
      #sliderInput("ig_length_scale", "Inverse-Gaussian shape:", 1, 15, 1),
      checkboxInput("length_scale_mle", "Use MLE", FALSE),
      actionButton("draw_length_scale", "Draw New Length Scale"),
      tags$hr(),
      
      tags$h4("Hyperparameters for variance"),
      sliderInput("ht_mu", "Half-t mean:", 0, 15, 0),
      sliderInput("ht_df", "Half-t degrees of freedom:", 1, 5, 4),
      sliderInput("ht_scale", "Half-t scale:", 1, 15, 1),
      checkboxInput("variance_mle", "Use MLE", FALSE),
      actionButton("draw_variance", "Draw New Variance"),
      tags$hr(),

      #tags$h4("Other parameters"),
      #sliderInput("sigma_n", "Noise amplitude:", 0, 15, 1),
      sliderInput("nfunc", "Number of Functions:", 1, n_functions, 3),
      #sliderInput("n_points", "Number of X Points:", 20, 400, 200),
      #sliderInput("x_max", "X Range:", 2, 20, 10),
      actionButton("draw_gp", "Draw GP")
    ),
    
    mainPanel(
      #actionButton("restart", "Restart Session"),
      fluidRow(
        column(6, plotOutput("plot_ig", height = "250px")),
        #column(6, plotOutput("plot_ig_gauss", height = "250px"))
        column(6, plotOutput("plot_ht", height = "250px"))
      ),
      tags$br(),
      plotOutput("kernelPlot", height = "250px"),
      tags$br(),
      plotOutput("gpPlot", height = "600px")
      
    )
  )
)

server <- function(input, output) {
  
  # Kernel choice
  kernel_choice <- eventReactive(input$kernel_label, {
    return(kernels[[input$kernel_label]])
  })
  

  #Inv-Gamma
  length_scale_draw <- eventReactive(input$draw_length_scale, { 
    seed_val <- digest(list(input$ig_alpha, input$ig_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    if (isTRUE(input$length_scale_mle)) {
      return(input$ig_beta / (input$ig_alpha + 1))   # MLE of InvGamma(α,β)
    } else {
      return(1 / rgamma(1, shape = input$ig_alpha, rate = input$ig_beta))
    }
  })
  
  # Inv-Gaussian
  #length_scale_draw_gauss <- eventReactive(input$draw_length_scale, { 
  #  seed_val <- digest(list(input$ig_mu, input$ig_length_scale), algo="xxhash32", serialize=TRUE) |> 
  #    substr(1,7) |> strtoi(base=16)
  #  set.seed(seed_val)
  #  
  #  if (isTRUE(input$length_scale_mle)) {
  #    return(input$ig_mu*(sqrt(1 + (9*input$ig_mu)/(4*input$ig_length_scale))-(3*input$ig_mu)/(2*input$ig_length_scale)))   # MLE of InvGaussian?
  #  } else {
  #    return(rinvgauss(1, mean=input$ig_mu, shape=input$ig_length_scale))
  #  }
  #})
  
  
  variance_draw <- eventReactive(input$draw_variance, {
    seed_val <- digest(list(input$ht_mu, input$ht_df, input$ht_scale), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    if (isTRUE(input$variance_mle)) {
      return(input$ht_mu)    # half-t MLE occurs at lower bound = μ
    } else {
      return(input$ht_mu + input$ht_scale * varhalfT(runif(1), n=input$ht_df))
    }
  })
  
  
  gp_funcs <- reactive({
    req(kernel_choice(), length_scale_draw(), variance_draw())
    
    x_orig <- seq(x_min, x_max, length.out = n_points)
    funcs <- replicate(input$nfunc,
                       simulate_gp(x_orig, kernel_choice(), length_scale_draw(), variance_draw()))#, mean_fun = function(x) 10 + 5 * x_orig) ))
    
    list(x_orig = x_orig, funcs = funcs)
  })

  
  last_params <- reactiveVal(NULL)
  last_pool   <- reactiveVal(NULL)
  
  gp_pool <- eventReactive(input$draw_gp, {
    req(kernel_choice(), length_scale_draw(), variance_draw())
    
    kernel <- kernel_choice()
    len <- length_scale_draw()
    var <- variance_draw()
    old_params <- last_params()
    old_pool   <- last_pool()
    
    # reuse pool when parameters unchanged
    if (!is.null(old_params) &&
        is.list(old_params) &&
        identical(old_params$kernel_prev, kernel) &&
        identical(signif(old_params$length_scale,10), signif(len,10)) &&
        identical(signif(old_params$variance,10),  signif(var,10)) &&
        !is.null(old_pool) &&
        is.list(old_pool)) {
      
      return(old_pool)
    }
    
    x_orig <- seq(x_min, x_max, length.out = n_points)
    funcs <- replicate(
      n_functions,
      simulate_gp(x_orig, kernel, len, var)#, mean_fun = function(x) 10 + 5 * x_orig)

    )
    
    new_pool <- list(x_orig = x_orig, funcs = funcs)
    
    last_params(list(kernel_prev = kernel, length_scale = len, variance = var))
    last_pool(new_pool)
    
    new_pool
  })
  
  
  gp_data <- reactive({
    req(gp_pool())
    
    idx <- 1:input$nfunc
    idx <- idx[idx <= 100]   # safety
    
    x_new <- seq(x_min, x_max, length.out = n_points)#input$x_max, length.out = input$n_points
    
    funcs_interp <- apply(gp_pool()$funcs[, idx, drop=FALSE], 2, function(f) {
      approx(gp_pool()$x_orig, f, xout = x_new)$y
    })
    
    data.frame(
      x = rep(x_new, length(idx)),
      f = as.vector(funcs_interp),
      func = rep(idx, each = length(x_new))
    )
  })
  
  
  
  ### plots
  
  # Inverse-Gamma prior for length_scale plot
  output$plot_ig <- renderPlot({
    req(length_scale_draw())
    
    x <- seq(1e-6, 15, length.out = 400)
    alpha <- input$ig_alpha
    beta  <- input$ig_beta
    
    dens <- (beta^alpha / gamma(alpha)) * x^(-alpha-1) * exp(-beta / x)
    d <- data.frame(x=x, y=dens)
    
    ggplot(d, aes(x,y)) +
      geom_line(color="steelblue", linewidth=1) +
      geom_point(aes(x=length_scale_draw(), y=0), color="red", size=3) +
      annotate("text",
               x = max(d$x) * 0.8,
               y = max(d$y) * 0.9,
               label = paste0("λ = ", signif(length_scale_draw(),3)),
               color="red",
               size = 5,
               hjust = 0) +
      labs(title="Inverse-Gamma prior for λ",
           y="density", x="λ") +
      theme_minimal(base_size=14)
  })
  
  # Inverse-Gaussian prior for length_scale plot
  #output$plot_ig_gauss <- renderPlot({
  #  req(length_scale_draw_gauss())
  #  
  #  x <- seq(1e-6, 50, length.out = 400)
  #  
  #  dens <- dinvgauss(x, mean=input$ig_mu, shape=input$ig_length_scale)
  #  d <- data.frame(x=x, y=dens)
  #  
  #  ggplot(d, aes(x,y)) +
  #    geom_line(color="steelblue", linewidth=1) +
  #    geom_point(aes(x=length_scale_draw_gauss(), y=0), color="red", size=3) +
  #    annotate("text",
  #             x = max(d$x) * 0.8,
  #             y = max(d$y) * 0.9,
  #             label = paste0("λ = ", signif(length_scale_draw_gauss(),3)),
  #             color="red",
  #             size = 5,
  #             hjust = 0) +
  #    labs(title="Inverse-Gaussian prior for λ",
  #         y="density", x="λ") +
  #    theme_minimal(base_size=14)
  #})
  
  # Half-t prior for variance plot
  output$plot_ht <- renderPlot({
    req(variance_draw())
    
    x <- seq(0, 15, length.out = 400)
    df <- input$ht_df
    mu <- input$ht_mu
    sc <- input$ht_scale
    
    dens <- 2 * dt((x - mu)/sc, df = df) / sc
    dens[x < mu] <- 0
    d <- data.frame(x=x, y=dens)
    
    ggplot(d, aes(x,y)) +
      geom_line(color="darkgreen", linewidth=1) +
      geom_point(aes(x=variance_draw(), y=0), color="red", size=3) +
      annotate("text",
               x = max(d$x) * 0.8,
               y = max(d$y) * 0.9,
               label = paste0("σ^2 = ", signif(variance_draw(),3)),
               color="red",
               size = 5,
               hjust = 0) +
      labs(title="Half-t prior for σ^2",
           y="density", x="σ^2") +
      theme_minimal(base_size=14)
  })
  
 # Kernel based on distance plot
  output$kernelPlot <- renderPlot({
    req(kernel_choice(), length_scale_draw(), variance_draw())
    
    dist <- seq(-3, 3, length.out = 300)
    x_o <- rep(0, length(dist))
    length_scale <- length_scale_draw()
    variance  <- variance_draw()
    
    k <- kernel_choice()(dist, x_o, length_scale, variance)
    
    ggplot(data.frame(dist=dist, k=k), aes(dist,k)) +
      geom_line(color="purple", linewidth=1.2) +
      labs(title=input$kernel_label, x="Distance", y="k(x)") +
      theme_minimal(base_size=14)
  })
  
  
  # GP prior draws plot
  output$gpPlot <- renderPlot({
    req(gp_data())
    
    ggplot(gp_data(), aes(x=x, y=f, group=func, color=factor(func))) +
      geom_line(alpha=0.9, linewidth=1) +
      scale_color_discrete(guide="none") +
      labs(title="Gaussian Process Prior Samples",
           subtitle="Squared Exponential Kernel",
           x="x", y="f(x)") +
      theme_minimal(base_size=16)
  })
  
}

shinyApp(ui, server)
