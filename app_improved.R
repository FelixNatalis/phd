library(shiny)
library(ggplot2)
library(VaRES)
library(digest)

# Kernel
k_se <- function(x1, x2, lambda, sigma) {
  outer(x1, x2, function(a, b)
    sigma^2 * exp(-(a - b)^2 / (2 * lambda^2))
  )
}

simulate_gp <- function(x, lambda, sigma) {
  K <- k_se(x, x, lambda, sigma)
  L <- chol(K + 1e-6 * diag(length(x)))
  drop(t(L) %*% rnorm(length(x)))
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      tags$h4("Hyperparameters for λ"),
      sliderInput("ig_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
      sliderInput("ig_beta", "Inverse-Gamma beta:", 1, 15, 1),
      actionButton("draw_lambda", "Draw New Lambda"),
      tags$hr(),
      
      tags$h4("Hyperparameters for σ"),
      sliderInput("ht_mu", "Half-t mean:", 0, 15, 0),
      sliderInput("ht_df", "Half-t degrees of freedom:", 1, 5, 4),
      sliderInput("ht_scale", "Half-t scale:", 1, 15, 1),
      actionButton("draw_sigma", "Draw New Sigma"),
      tags$hr(),

      
      tags$h4("Other parameters"),
      sliderInput("nfunc", "Number of Functions:", 1, 10, 3),
      sliderInput("npoints", "Number of X Points:", 20, 400, 200),
      sliderInput("xmax", "X Range:", 2, 20, 10),
      actionButton("draw_gp", "Draw GP")
    ),
    
    mainPanel(
      fluidRow(
        column(6, plotOutput("plot_ig", height = "250px")),
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
  
  # Lambda reactive (depends only on IG hyperparameters)
  lambda_draw <- eventReactive(input$draw_lambda, {
    seed_val <- digest(list(input$ig_alpha, input$ig_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    # draw 1/lambda ~ Gamma(...) ⇒ lambda ~ InvGamma(...)
    1 / rgamma(1, shape = input$ig_alpha, rate = input$ig_beta)
  })
  
  # Sigma reactive (depends only on Half-t hyperparameters)
  sigma_draw <- eventReactive(input$draw_sigma, {
    seed_val <- digest(list(input$ht_mu, input$ht_df, input$ht_scale), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    input$ht_mu + input$ht_scale * varhalfT(runif(1), n=input$ht_df)
  })
  
  # GP functions
  gp_funcs <- reactive({
    req(lambda_draw(), sigma_draw())
    
    x_orig <- seq(0, 10, length.out = 200)
    funcs <- replicate(input$nfunc,
                       simulate_gp(x_orig, lambda_draw(), sigma_draw()))
    
    list(x_orig = x_orig, funcs = funcs)
  })

  
  last_params <- reactiveVal(NULL)
  last_pool   <- reactiveVal(NULL)
  
  gp_pool <- eventReactive(input$draw_gp, {
    req(lambda_draw(), sigma_draw())
    
    lam <- lambda_draw()
    sig <- sigma_draw()
    old_params <- last_params()
    old_pool   <- last_pool()
    
    # reuse pool when parameters unchanged
    if (!is.null(old_params) &&
        is.list(old_params) &&
        identical(signif(old_params$lambda,10), signif(lam,10)) &&
        identical(signif(old_params$sigma,10),  signif(sig,10)) &&
        !is.null(old_pool) &&
        is.list(old_pool)) {
      
      return(old_pool)
    }
    
    x_orig <- seq(0, 10, length.out = 200)
    funcs <- replicate(
      100,
      simulate_gp(x_orig, lam, sig)
    )
    
    new_pool <- list(x_orig = x_orig, funcs = funcs)
    
    last_params(list(lambda = lam, sigma = sig))
    last_pool(new_pool)
    
    new_pool
  })
  
  
  gp_data <- reactive({
    req(gp_pool())
    
    idx <- 1:input$nfunc
    idx <- idx[idx <= 100]   # safety
    
    x_new <- seq(0, input$xmax, length.out = input$npoints)
    
    funcs_interp <- apply(gp_pool()$funcs[, idx, drop=FALSE], 2, function(f) {
      approx(gp_pool()$x_orig, f, xout = x_new)$y
    })
    
    data.frame(
      x = rep(x_new, length(idx)),
      f = as.vector(funcs_interp),
      func = rep(idx, each = length(x_new))
    )
  })
  
  
  
  ### --- PRIOR PLOTS ---------------------------------------------------------
  
  # Inverse-Gamma prior for lambda
  output$plot_ig <- renderPlot({
    req(lambda_draw())
    
    x <- seq(1e-6, 5, length.out = 400)
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
  
  # Half-t prior for sigma
  output$plot_ht <- renderPlot({
    req(sigma_draw())
    
    x <- seq(0, 15, length.out = 400)
    df <- input$ht_df
    mu <- input$ht_mu
    sc <- input$ht_scale
    
    dens <- 2 * dt((x - mu)/sc, df = df) / sc
    dens[x < mu] <- 0
    d <- data.frame(x=x, y=dens)
    
    ggplot(d, aes(x,y)) +
      geom_line(color="darkgreen", linewidth=1) +
      geom_point(aes(x=sigma_draw(), y=0), color="red", size=3) +
      annotate("text",
               x = max(d$x) * 0.8,
               y = max(d$y) * 0.9,
               label = paste0("σ = ", signif(sigma_draw(),3)),
               color="red",
               size = 5,
               hjust = 0) +
      labs(title="Half-t prior for σ",
           y="density", x="σ") +
      theme_minimal(base_size=14)
  })
  
  
  output$kernelPlot <- renderPlot({
    req(lambda_draw(), sigma_draw())
    
    x <- seq(-3, 3, length.out = 300)
    lambda <- lambda_draw()
    sigma  <- sigma_draw()
    
    k <- sigma^2 * exp(-x^2 / (2 * lambda^2))
    
    ggplot(data.frame(x=x, k=k), aes(x,k)) +
      geom_line(color="purple", linewidth=1.2) +
      labs(title="Squared-Exponential Kernel", x="Distance", y="k(x)") +
      theme_minimal(base_size=14)
  })
  
  
  
  ### --- GP PLOT -------------------------------------------------------------
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
