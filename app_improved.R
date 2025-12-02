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
      sliderInput("ig_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
      sliderInput("ig_beta", "Inverse-Gamma beta:", 1, 15, 1),
      sliderInput("ht_mu", "Half-t mean:", 0, 15, 0),
      sliderInput("ht_df", "Half-t degrees of freedom:", 1, 5, 4),
      sliderInput("ht_scale", "Half-t scale:", 1, 15, 1),
      sliderInput("nfunc", "Number of Functions:", 1, 10, 5),
      sliderInput("npoints", "Number of X Points:", 20, 400, 200),
      sliderInput("xmax", "X Range:", 2, 20, 10),
      actionButton("draw_sigma", "Draw New Sigma"),
      actionButton("draw_lambda", "Draw New Lambda")
    ),
    mainPanel(plotOutput("gpPlot", height = "600px"))
  )
)

server <- function(input, output) {
  
  # Lambda reactive (depends only on IG hyperparameters)
  lambda_draw <- eventReactive(input$draw_lambda, {
    seed_val <- digest(list(input$ig_alpha, input$ig_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    1 / rgamma(1, shape = input$ig_alpha, rate = input$ig_beta)
  })
  
  # Sigma reactive (depends only on Half-t hyperparameters)
  sigma_draw <- eventReactive(input$draw_sigma, {
    seed_val <- digest(list(input$ht_mu, input$ht_df, input$ht_scale), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    input$ht_mu + input$ht_scale * varhalfT(runif(1), n=input$ht_df)
  })
  
  # GP functions (depends on lambda, sigma)
  gp_funcs <- reactive({
    req(lambda_draw(), sigma_draw())
    
    # fixed grid
    x_orig <- seq(0, 10, length.out = 200)
    
    funcs <- replicate(input$nfunc,
                       simulate_gp(x_orig, lambda_draw(), sigma_draw()))
    
    list(x_orig = x_orig, funcs = funcs)
  })
  
  # Interpolate to current x-grid
  gp_data <- reactive({
    req(gp_funcs())
    
    x_new <- seq(0, input$xmax, length.out = input$npoints)
    
    funcs_interp <- apply(gp_funcs()$funcs, 2, function(f) {
      approx(gp_funcs()$x_orig, f, xout = x_new)$y
    })
    
    data.frame(
      x = rep(x_new, input$nfunc),
      f = as.vector(funcs_interp),
      func = rep(1:input$nfunc, each = length(x_new))
    )
  })
  
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
