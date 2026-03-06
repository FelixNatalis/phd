# imports
library(shiny)
library(ggplot2)
library(VaRES)
library(digest)
library(statmod)
library(hash)
library(shiny)
library(bslib)
library(gtools)
library(shinyjs)

#-------------------------------------------------------------------------------
################################################################################
# GLOBAL CONSTANTS
################################################################################
#-------------------------------------------------------------------------------

# defaults parameters for GP prior predictive draws
n_functions = 10
x_min = -10
x_max = 10
n_points = 200
epsilon = 1e-6

#-------------------------------------------------------------------------------
################################################################################
# LIBRARY
################################################################################
#-------------------------------------------------------------------------------

## Kernels

# SE kernel function
squared_exponential_kernel <- function(x1, x2, variance, length_scale) {
  outer(x1, x2, function(a, b)
    (variance^2 * exp(-(a - b)^2 / (2 * length_scale^2)))
  )
}

# Linear kernel function
linear_kernel <- function(x1, x2, variance) {
  outer(x1, x2, function(a, b)
    (variance^2 * a * b)
  )
}

# Matérn kernel function
matern_kernel <- function(x1, x2, variance, length_scale, roughness) {
  distance <- outer(x1, x2, function(a, b) abs(a - b))
  term <- sqrt(2 * roughness) * distance / length_scale
  K <- variance^2 * (2^(1 - roughness) / gamma(roughness)) * (term^roughness) * besselK(term, roughness)
  K[distance == 0] <- variance^2   # Handle the distance = 0 case explicitly
  K
}

# Periodic kernel function
periodic_kernel <- function(x1, x2, variance, length_scale, period) {
  outer(x1, x2, function(a, b)
    (variance^2 * exp(-2 * sin(pi * abs(a - b) / period)^2 / length_scale^2))
  )
}

# Periodic * squared exponential kernel function
per_and_se_kernel <- function(x1, x2, variance, length_scale, period) {
  outer(x1, x2, function(a, b)
    (variance^2 * exp(-2 * sin(pi * abs(a - b) / period)^2 / length_scale^2)) * 
      (variance^2 * exp(-(a - b)^2 / (2 * length_scale^2)))
  )
}

# Periodic + squared exponential kernel function
per_or_se_kernel <- function(x1, x2, variance, length_scale, period) {
  outer(x1, x2, function(a, b)
    (variance^2 * exp(-2 * sin(pi * abs(a - b) / period)^2 / length_scale^2)) + 
      (variance^2 * exp(-(a - b)^2 / (2 * length_scale^2)))
  )
}

kernels <- hash(
  "Squared Exponential" = squared_exponential_kernel
  ,"Matérn" = matern_kernel
  ,"Linear" = linear_kernel
  ,"Periodic" = periodic_kernel
  )

kernel_combinations = hash(
          #"Changepoint" = 5,
          #"Lin + SE" = 6,
          #"Lin * SE" = 7,
          "Per + SE" = per_or_se_kernel,
          "Per * SE" = per_and_se_kernel
)

kernel_wrapper <- function(kernel_label, x1, x2, params){ 
  variance = params[["variance"]]
  length_scale = params[["length_scale"]]
  roughness = params[["roughness"]]
  period = params[["period"]]
  
  if(kernel_label == "Linear"){
    if(!invalid(x1) & !invalid(x2) & !invalid(variance)){
      return(linear_kernel(x1 = x1, x2 = x2, variance = variance))
    }
  }
  
  if(kernel_label == "Squared Exponential"){
    if(!invalid(x1) & !invalid(x2) & !invalid(variance) & !invalid(length_scale)){
      return(squared_exponential_kernel(x1 = x1, x2 = x2, variance = variance, length_scale = length_scale))
    }
  }
  
  if(kernel_label == "Matérn"){
    if(!invalid(x1) & !invalid(x2) & !invalid(variance) & !invalid(length_scale) & !invalid(roughness)){
      return(matern_kernel(x1 = x1, x2 = x2, variance = variance, length_scale = length_scale, roughness = roughness))
    }
  }
  
  if(kernel_label == "Periodic"){
    if(!invalid(x1) & !invalid(x2) & !invalid(variance) & !invalid(length_scale) & !invalid(period)){
      return(periodic_kernel(x1 = x1, x2 = x2, variance = variance, length_scale = length_scale, period = period))
    }
  }
  
  if(kernel_label == "Per * SE"){
    if(!invalid(x1) & !invalid(x2) & !invalid(variance) & !invalid(length_scale) & !invalid(period)){
      return(per_and_se_kernel(x1 = x1, x2 = x2, variance = variance, length_scale = length_scale, period = period))
    }
  }
  
  if(kernel_label == "Per + SE"){
    if(!invalid(x1) & !invalid(x2) & !invalid(variance) & !invalid(length_scale) & !invalid(period)){
      return(per_or_se_kernel(x1 = x1, x2 = x2, variance = variance, length_scale = length_scale, period = period))
    }
  }
}

#-------------------------------------------------------------------------------
# Drawing GP
simulate_gp <- function(x, kernel_label, kernel_params, sigma_noise = 1e-3, mean_fun = function(x) 0) {
  K <- kernel_wrapper(kernel_label, x, x, params = kernel_params) 
  L <- chol(K + epsilon * diag(length(x)))
  m <- mean_fun(x)
  f <- m + t(L) %*% rnorm(length(x))
  
  # noise
  eps <- sigma_noise * rnorm(length(x))
  
  drop(f + eps)
}
#-------------------------------------------------------------------------------
## Distributions

# inverse gamma
inverse_gamma <- function(alpha, beta, x = NULL){
  if(invalid(x)){
    1 / rgamma(1, shape = alpha, rate = beta)
  }else{
    dgamma(1/x, shape = alpha, rate = beta) * 1/x^2
  }
}

# half-t
half_t <- function(df, scale, x = NULL){
  if(invalid(x)){
    scale * abs(rt(1, df = df))
  }else{
    dens <- 2 * dt((x)/scale, df = df) / scale
    dens[x < 0] <- 0
    dens
  }
}

#-------------------------------------------------------------------------------
################################################################################
# UI
################################################################################
#-------------------------------------------------------------------------------

ui <- page_fillable(
  tags$style(type='text/css', 
  ".selectize-input { 
    font-size:80%; 
    line-height: 16px;
  } 
  .selectize-dropdown { 
    font-size:80%; 
    line-height: 16px; 
  } 
  .control-label { 
    font-size:80%;
  }
  .card-header {
    padding: 4px;
    font-size:80%;
  }
  .btn {
    padding: 2px 6px;
    font-size: 80%;
  }
  .card-body {
    display: inline-block;
  }"),
  
  useShinyjs(),
  layout_columns(
    card( "Parameters",
            card( 
              card_header("Kernel", style='padding:4px; font-size:80%'),
              selectInput("kernel_label", "Choose a kernel:",
                          list(`Simple kernels` = keys(kernels), `Kernel combinations` = keys(kernel_combinations)
                          )),
              actionButton("draw_kernel", "Draw kernel")
            ), 
            
          card(
            card_header("Hyperparameters for variance"),
            # sliderInput("ht_mu", "Half-t mean:", 0, 15, 0),
            sliderInput("ht_df", "Half-t degrees of freedom:", 1, 5, 4),
            sliderInput("ht_scale", "Half-t scale:", 1, 15, 1),
           # checkboxInput("variance_mle", "Use MLE", FALSE),
            actionButton("draw_variance", "Draw New Variance"),

          ),
          uiOutput("dynamic_ls_choice"),
          uiOutput("dynamic_nu_choice"),
          uiOutput("dynamic_per_choice"),
          card(
            card_header("Other parameters"),
            #sliderInput("sigma_n", "Noise amplitude:", 0, 15, 1),
            sliderInput("nfunc", "Number of Functions:", 1, n_functions, 3),
            #sliderInput("n_points", "Number of X Points:", 20, 400, 200),
            #sliderInput("x_max", "X Range:", 2, 20, 10),
            disabled(actionButton("draw_gp", "Draw GP"))
          )
    ),
    card(
      layout_columns(
        card(plotOutput("plot_ht", height = "250px")),
        uiOutput("dynamic_ls_plot"),
        uiOutput("dynamic_per_plot"),
        
        card(plotOutput("kernelPlot", height = "250px"))
      ),
      card(plotOutput("gpPlot", height = "600px"))
      
    ),
 col_widths = c(2, 10) 
  )
)

#-------------------------------------------------------------------------------
################################################################################
# SERVER
################################################################################
#-------------------------------------------------------------------------------

server <- function(input, output) {
  
  ## Observable variables
  last_params   <- reactiveVal(NULL)
  last_pool     <- reactiveVal(NULL)
  variance      <- reactiveVal(NULL)
  length_scale  <- reactiveVal(NULL)
  roughness     <- reactiveVal(NULL)
  period        <- reactiveVal(NULL)

  rv <- reactiveValues(
    var = FALSE,
    ls  = FALSE,
    per = FALSE
  )
  
  #-------------------------------------------------------------------------------
  ## Dynamic ui elements 
  
  shinyjs::disable("draw_kernel")
  
  observe({
    if (input$kernel_label == "Linear" && rv$var) {
      shinyjs::enable("draw_kernel")
      
    } else if ((input$kernel_label == "Periodic" | input$kernel_label == "Per * SE"| input$kernel_label == "Per + SE") &&
               rv$var && rv$ls && rv$per) {
      shinyjs::enable("draw_kernel")
      
    } else if (input$kernel_label == "Squared Exponential" &&
               rv$var && rv$ls) {
      shinyjs::enable("draw_kernel")
      
    } else if (input$kernel_label == "Matérn" &&
               rv$var && rv$ls) {
      shinyjs::enable("draw_kernel")
      
    } else {
      shinyjs::disable("draw_kernel")
    }
  })
  
  
  output$dynamic_ls_choice <- renderUI({
    if (input$kernel_label != "Linear") {
      card(
        card_header("Hyperparameters for length scale"),
        sliderInput("ls_ig_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
        sliderInput("ls_ig_beta", "Inverse-Gamma beta:", 1, 15, 1),
        actionButton("draw_length_scale", "Draw New Length Scale"),
      )
    }
  })
  
  output$dynamic_per_choice <- renderUI({
    if (input$kernel_label == "Periodic" | input$kernel_label == "Per * SE"| input$kernel_label == "Per + SE") {
      card(
        card_header("Hyperparameters for period"),
        sliderInput("per_ig_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
        sliderInput("per_ig_beta", "Inverse-Gamma beta:", 1, 15, 1),
        actionButton("draw_period", "Draw New Period")
      )
    }
  })
  
  output$dynamic_nu_choice <- renderUI({
    if (input$kernel_label == "Matérn") {
      card(
        card_header("Hyperparameters for roughness"),
        sliderInput("nu", "Value of nu:", 0.5, 5.5, 0.5)
      )
    }
  })
  
  output$dynamic_ls_plot <- renderUI({
    if (input$kernel_label != "Linear") {
      card(plotOutput("plot_ig", height = "250px"))
    }
  })
  
  
  output$dynamic_per_plot <- renderUI({
    if (input$kernel_label == "Periodic" | input$kernel_label == "Per * SE"| input$kernel_label == "Per + SE") {
      card(plotOutput("plot_ig_per", height = "250px"))
    }
  })
  
  #-------------------------------------------------------------------------------
  ## Reactive events
  
  # drawing length scale on button click
  length_scale_draw <- eventReactive(input$draw_length_scale, { 
    seed_val <- digest(list(input$ls_ig_alpha, input$ls_ig_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$ls <- TRUE
    ls <- inverse_gamma(alpha = input$ls_ig_alpha, beta = input$ls_ig_beta)
    length_scale(ls)
  })
  
  # drawing period on button click
  period_draw <- eventReactive(input$draw_period, { 
    seed_val <- digest(list(input$per_ig_alpha, input$per_ig_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$per <- TRUE
    per <- inverse_gamma(alpha = input$per_ig_alpha, beta = input$per_ig_beta)
    period(per)
  })
  
  # drawing variance on button click
  variance_draw <- eventReactive(input$draw_variance, {
    seed_val <- digest(list(input$ht_df, input$ht_scale), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$var <- TRUE
    var <- half_t(df = input$ht_df, scale = input$ht_scale)
    variance(var) 
  })
  
  #-------------------------------------------------------------------------------
  ## GP redraw logic
  gp_pool <- eventReactive(input$draw_gp, {
    
    kernel <- input$kernel_label
    len <- length_scale()
    var <- variance()
    per <- period()
    ro <- input$nu
    old_params <- last_params()
    old_pool   <- last_pool()
    
    # reuse pool when parameters unchanged
    if (!is.null(old_params) &&
        is.list(old_params) &&
        identical(old_params$kernel_prev, kernel) &&
        identical(signif(old_params$length_scale,10), signif(len,10)) &&
        identical(signif(old_params$variance,10),  signif(var,10)) &&
        identical(signif(old_params$period,10), signif(per,10)) &&
        identical(signif(old_params$roughness,10), signif(ro,10)) &&
        !is.null(old_pool) &&
        is.list(old_pool)) {
      
      return(old_pool)
    }
    
    x_orig <- seq(x_min, x_max, length.out = n_points)
    
    kernel_params <- hash(
      "variance" = var
      ,"length_scale" = len
      ,"period" = per
      ,"roughness" = ro
    )
    
    funcs <- replicate(
      n_functions,
      simulate_gp(x_orig, input$kernel_label, kernel_params)
    )
    
    new_pool <- list(x_orig = x_orig, funcs = funcs)
    
    last_params(list(kernel_prev = kernel, length_scale = len, variance = var, period = per, roughness = ro))
    last_pool(new_pool)
    
    new_pool
  })
  
  # For storing and reproducing draws that were already computed
  gp_data <- reactive({
    req(gp_pool())
    
    idx <- 1:input$nfunc
    idx <- idx[idx <= 100]   # safety
    
    x_new <- seq(x_min, x_max, length.out = n_points)
    
    funcs_interp <- apply(gp_pool()$funcs[, idx, drop=FALSE], 2, function(f) {
      approx(gp_pool()$x_orig, f, xout = x_new)$y
    })
    
    data.frame(
      x = rep(x_new, length(idx)),
      f = as.vector(funcs_interp),
      func = rep(idx, each = length(x_new))
    )
  })
  
  observeEvent(input$draw_kernel, {
    enable("draw_gp")
  })
  
  #-------------------------------------------------------------------------------
  ## plots
  
  # Inverse-Gamma prior for length_scale plot
  output$plot_ig <- renderPlot({
    observe(length_scale_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$ls_ig_alpha
    beta  <- input$ls_ig_beta
    
    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq) 
    d <- data.frame(x=x_seq, y=dens)
    
    plot <- ggplot(d, aes(x,y)) +
      geom_line(color="steelblue", linewidth=1) +
      labs(title="Inverse-Gamma prior for length scale",
           y="density", x="length scale") +
      theme_minimal(base_size=14)
    
    if(!invalid(length_scale())){
      plot <- plot + annotate("point", x = length_scale(), y = 0, colour = "red", size = 3) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 hjust = 1.1,  
                 vjust = 1.5,
                 label = paste0("draw = ", signif(length_scale(),3)),
                 color="red",
                 size = 5)
    }
    
    plot
  })
  
  # Half-t prior for variance plot
  output$plot_ht <- renderPlot({
    observe(variance_draw())
    
    x_seq <- seq(0, 15, length.out = 400)
    df <- input$ht_df
    sc <- input$ht_scale
    
    dens <- half_t(df = df, scale = sc, x = x_seq)

    d <- data.frame(x=x_seq, y=dens)
    
    plot <- ggplot(d, aes(x,y)) +
      geom_line(color="darkgreen", linewidth=1) +
      labs(title="Half-t prior for variance",
           y="density", x="variance") +
      theme_minimal(base_size=14)
    
    if(!invalid(variance())){
      plot <- plot + 
        annotate("point", x = variance(), y = 0, colour = "red", size = 3) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 hjust = 1.1,  
                 vjust = 1.5,
                 label = paste0("draw = ", signif(variance(),3)),
                 color="red",
                 size = 5)
    }
    
    plot
  })
  
  # Inverse-Gamma prior for period plot
  output$plot_ig_per <- renderPlot({
    observe(period_draw())

    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$per_ig_alpha
    beta  <- input$per_ig_beta

    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq) 
    d <- data.frame(x=x_seq, y=dens)

    plot <- ggplot(d, aes(x,y)) +
      geom_line(color="steelblue", linewidth=1) +
      labs(title="Inverse-Gamma prior for period",
           y="density", x="period") +
      theme_minimal(base_size=14)
    
    if(!invalid(period())){
      plot <- plot + 
        annotate("point", x = period(), y = 0, colour = "red", size = 3) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 hjust = 1.1,  
                 vjust = 1.5,
                 label = paste0("draw = ", signif(period(),3)),
                 color="red",
                 size = 5)
    }
    
    plot
 })
  
 # Kernel based on distance plot
  output$kernelPlot <- renderPlot({
    req(input$draw_kernel)
    
    dist <- seq(-3, 3, length.out = 300)
    x_o <- rep(0, length(dist))

    kernel_params <- hash(
      "variance" = variance()
      ,"length_scale" = length_scale()
      ,"period" = period()
      ,"roughness" = input$nu
    )
    
    k <- kernel_wrapper(input$kernel_label, dist, x_o, kernel_params)[,1]
    
    ggplot(data.frame(dist=dist, k=k), aes(dist,k)) +
      geom_line(color="purple", linewidth=1.2) +
      labs(title=input$kernel_label, x="Distance", y="k(x)") +
      theme_minimal(base_size=14)
  })
  
  
  # GP prior draws plot
  output$gpPlot <- renderPlot({
    req(input$draw_gp)
    observe(gp_data())#, period_draw(), length_scale_draw(), variance_draw())
    ggplot(gp_data(), aes(x=x, y=f, group=func, color=factor(func))) +
      geom_line(alpha=0.9, linewidth=1) +
      scale_color_discrete(guide="none") +
      labs(title="Gaussian Process Prior Samples",
           subtitle= paste(input$kernel_label, "Kernel"),
           x="x", y="f(x)") +
      theme_minimal(base_size=16)
  })
  
}

shinyApp(ui, server)
