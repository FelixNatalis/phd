library(shiny)
library(ggplot2)
library(VaRES)
library(digest)
library(statmod)

n_functions = 10
xmax = 10
npoints = 200

# Kernel
k_se <- function(x1, x2, lambda, sigma_2) {
  outer(x1, x2, function(a, b)
    sigma_2 * exp(-(a - b)^2 / (2 * lambda^2))
  )
}

simulate_gp <- function(x, lambda, sigma_k, sigma_n = 1e-3, mean_fun = function(x) 0) {
  
  K <- k_se(x, x, lambda, sigma_k)
  L <- chol(K + 1e-6 * diag(length(x)))
  
  m <- mean_fun(x)
  
  f <- m + t(L) %*% rnorm(length(x))
  
  # noise
  eps <- sigma_n * rnorm(length(x))
  
  drop(f + eps)
}


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
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

      
      #tags$h4("Other parameters"),
      #sliderInput("sigma_n", "Noise amplitude:", 0, 15, 1),
      sliderInput("nfunc", "Number of Functions:", 1, n_functions, 3),
      #sliderInput("npoints", "Number of X Points:", 20, 400, 200),
      #sliderInput("xmax", "X Range:", 2, 20, 10),
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
  
  # Inv-Gaussian
  #lambda_draw_gauss <- eventReactive(input$draw_lambda, { 
  #  seed_val <- digest(list(input$ig_mu, input$ig_lambda), algo="xxhash32", serialize=TRUE) |> 
  #    substr(1,7) |> strtoi(base=16)
  #  set.seed(seed_val)
  #  
  #  if (isTRUE(input$lambda_mle)) {
  #    return(input$ig_mu*(sqrt(1 + (9*input$ig_mu)/(4*input$ig_lambda))-(3*input$ig_mu)/(2*input$ig_lambda)))   # MLE of InvGaussian?
  #  } else {
  #    return(rinvgauss(1, mean=input$ig_mu, shape=input$ig_lambda))
  #  }
  #})
  
  
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
  
  
  gp_funcs <- reactive({
    req(lambda_draw(), sigma_2_draw())
    
    x_orig <- seq(0, 10, length.out = 200)
    funcs <- replicate(input$nfunc,
                       simulate_gp(x_orig, lambda_draw(), sigma_2_draw()))#, mean_fun = function(x) 10 + 5 * x_orig) ))
    
    list(x_orig = x_orig, funcs = funcs)
  })

  
  last_params <- reactiveVal(NULL)
  last_pool   <- reactiveVal(NULL)
  
  gp_pool <- eventReactive(input$draw_gp, {
    req(lambda_draw(), sigma_2_draw())
    
    lam <- lambda_draw()
    sig <- sigma_2_draw()
    old_params <- last_params()
    old_pool   <- last_pool()
    
    # reuse pool when parameters unchanged
    if (!is.null(old_params) &&
        is.list(old_params) &&
        identical(signif(old_params$lambda,10), signif(lam,10)) &&
        identical(signif(old_params$sigma_2,10),  signif(sig,10)) &&
        !is.null(old_pool) &&
        is.list(old_pool)) {
      
      return(old_pool)
    }
    
    x_orig <- seq(0, 10, length.out = 200)
    funcs <- replicate(
      n_functions,
      simulate_gp(x_orig, lam, sig)#, mean_fun = function(x) 10 + 5 * x_orig)

    )
    
    new_pool <- list(x_orig = x_orig, funcs = funcs)
    
    last_params(list(lambda = lam, sigma_2 = sig))
    last_pool(new_pool)
    
    new_pool
  })
  
  
  gp_data <- reactive({
    req(gp_pool())
    
    idx <- 1:input$nfunc
    idx <- idx[idx <= 100]   # safety
    
    x_new <- seq(0, xmax, length.out = npoints)#input$xmax, length.out = input$npoints
    
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
  
  # Inverse-Gamma prior for lambda
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
  
  # Inverse-Gaussian prior for lambda
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
  
  # Half-t prior for sigma
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
  
  
  output$kernelPlot <- renderPlot({
    req(lambda_draw(), sigma_2_draw())
    
    x <- seq(-3, 3, length.out = 300)
    lambda <- lambda_draw()
    sigma_2  <- sigma_2_draw()
    
    k <- sigma_2 * exp(-x^2 / (2 * lambda^2))
    
    ggplot(data.frame(x=x, k=k), aes(x,k)) +
      geom_line(color="purple", linewidth=1.2) +
      labs(title="Squared-Exponential Kernel", x="Distance", y="k(x)") +
      theme_minimal(base_size=14)
  })
  
  
  
  # gp
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
