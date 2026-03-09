# imports
library(shiny)
library(ggplot2)
library(VaRES)
library(digest)
library(statmod)
library(hash)
library(shinyWidgets)
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

# # Periodic * squared exponential kernel function
# per_and_se_kernel <- function(x1, x2, variance_1, length_scale_1, period, variance_2, length_scale_2) {
#   return(periodic_kernel(x1, x2, variance = variance_1, length_scale = length_scale_1, period = period)*
#     squared_exponential_kernel(x1, x2, variance = variance_2, length_scale = length_scale_2))
#   # outer(x1, x2, function(a, b)
#   #   (variance_1^2 * exp(-2 * sin(pi * abs(a - b) / period)^2 / length_scale_1^2)) * 
#   #     (variance_2^2 * exp(-(a - b)^2 / (2 * length_scale_2^2)))
#   # )
# }
# 
# # Periodic + squared exponential kernel function
# per_or_se_kernel <- function(x1, x2, variance_1, length_scale_1, period, variance_2, length_scale_2) {
#   return(periodic_kernel(x1, x2, variance = variance_1, length_scale = length_scale_1, period = period)+
#            squared_exponential_kernel(x1, x2, variance = variance_2, length_scale = length_scale_2))
# }

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
          "Per + SE" = 3,
          "Per * SE" = 3
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
    if(!invalid(x1) & !invalid(x2) 
       & !invalid(variance[1]) & !invalid(length_scale[1]) & !invalid(period) 
       & !invalid(variance[2]) & !invalid(length_scale[2])){
      # return(per_and_se_kernel(x1 = x1, x2 = x2, variance_1 = variance[1], length_scale_1 = length_scale[1], period = period, 
      #                           variance_2 = variance[2], length_scale_2 = length_scale[2]))
      return(periodic_kernel(x1 = x1, x2 = x2, variance = variance[1], length_scale = length_scale[1], period = period)*
        squared_exponential_kernel(x1 = x1, x2 = x2, variance = variance[2], length_scale = length_scale[2]))
    }
  }
  
  if(kernel_label == "Per + SE"){
    if(!invalid(x1) & !invalid(x2) 
       & !invalid(variance[1]) & !invalid(length_scale[1]) & !invalid(period) 
       & !invalid(variance[2]) & !invalid(length_scale[2])){
      return(periodic_kernel(x1 = x1, x2 = x2, variance = variance[1], length_scale = length_scale[1], period = period)+
               squared_exponential_kernel(x1 = x1, x2 = x2, variance = variance[2], length_scale = length_scale[2]))
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
  navset_tab( 
    nav_panel("Kernel", 
              card( "Parameters",
                    card( 
                      card_header("Kernel", style='padding:4px; font-size:80%'),
                      layout_columns(
                        column(width = 8,
                          selectInput("kernel_label", "Choose a kernel:",
                                      list(`Simple kernels` = keys(kernels), 
                                           `Kernel combinations` = keys(kernel_combinations))),
                          actionButton("draw_kernel", "Draw kernel")),
                        plotOutput("kernelPlot", height = "250px")
                      )
                    ), 
                    
                    card(
                      card_header("Hyperparameters for variance"),
                      layout_columns(
                        column(width = 8, 
                          sliderInput("variance_df", "Half-t degrees of freedom:", 1, 5, 4),
                          sliderInput("variance_scale", "Half-t scale:", 1, 15, 1),
                          actionButton("draw_variance", "Draw New Variance"))
                      ,plotOutput("plot_variance", height = "250px")
                      )
                    ),
                    uiOutput("dynamic_length_scale_choice"),
                    uiOutput("dynamic_roughness_choice"),
                    uiOutput("dynamic_period_choice"),
                    uiOutput("dynamic_variance_2_choice"),#!
                    uiOutput("dynamic_length_scale_2_choice"),#!
                    uiOutput("dynamic_roughness_2_choice"),#!
                    uiOutput("dynamic_period_2_choice"),#!
                    uiOutput("dynamic_changepoint_location_choice"),#!
                    uiOutput("dynamic_steepness_choice"),#!
              )), 
    nav_panel("GP", 
              
              card(

                card(
                  card_header("Other parameters"),
                  #sliderInput("sigma_n", "Noise amplitude:", 0, 15, 1),
                  sliderInput("nfunc", "Number of Functions:", 1, n_functions, 3),
                  sliderInput("n_points", "Number of X Points:", value = n_points, min = 2, max = 400),
                  sliderInput("x_range", "x range", min = -100, max = 100, value = c(x_min, x_max)), 
                  #sliderInput("x_range", "X Range:", 2, 20, 10),
                  #sliderInput("x_max", "X Range:", 2, 20, 10),
                  disabled(actionButton("draw_gp", "Draw GP"))
                ),
                card(plotOutput("gpPlot", height = "600px"))
                
              ),
              ), 
    id = "nav"),
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
  variance_2      <- reactiveVal(NULL)
  length_scale_2  <- reactiveVal(NULL)
  roughness_2     <- reactiveVal(NULL)
  period_2        <- reactiveVal(NULL)

  rv <- reactiveValues(
    var = FALSE,
    ls  = FALSE,
    per = FALSE,
    var_2 = FALSE,
    ls_2  = FALSE,
    per_2 = FALSE
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
  
  
  output$dynamic_length_scale_choice <- renderUI({
    if (input$kernel_label != "Linear") {
      card(
        card_header("Hyperparameters for length scale"),
        layout_columns(
          column(width = 8,
        sliderInput("length_scale_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
        sliderInput("length_scale_beta", "Inverse-Gamma beta:", 1, 15, 1),
        actionButton("draw_length_scale", "Draw New Length Scale")),
        plotOutput("plot_length_scale", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_period_choice <- renderUI({
    if (input$kernel_label == "Periodic" || input$kernel_label == "Per + SE" || input$kernel_label == "Per * SE" ) {
      card(
        card_header("Hyperparameters for period"),
        layout_columns(
          column(width = 8,
        sliderInput("period_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
        sliderInput("period_beta", "Inverse-Gamma beta:", 1, 15, 1),
        actionButton("draw_period", "Draw New Period")),
        plotOutput("plot_period", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_roughness_choice <- renderUI({
    if (input$kernel_label == "Matérn") {
      card(
        card_header("Hyperparameters for roughness"),
        sliderTextInput("nu", "Value of roughness:", choices = c(0.5, 1.5, 2.5, 3.5),
                        grid = TRUE, selected = 1.5)
      )
    }
  })
  #-------------------------------------------------------------------------------
  
  # Parameters for the second kernel
  
  output$dynamic_variance_2_choice <- renderUI({
    
    if (has.key(input$kernel_label , kernel_combinations)) {
      card(
        card_header("Hyperparameters for variance of the second kernel"),
        layout_columns(
          column(width = 8, 
                 sliderInput("variance_2_df", "Half-t degrees of freedom:", 1, 5, 4),
                 sliderInput("variance_2_scale", "Half-t scale:", 1, 15, 1),
                 actionButton("draw_variance_2", "Draw New Variance"))
          ,plotOutput("plot_variance_2", height = "250px")
        )
      )
    }
  })  
  
  output$dynamic_length_scale_2_choice <- renderUI({
    if (has.key(input$kernel_label , kernel_combinations)) {#! different condition
      card(
        card_header("Hyperparameters for length scale of the second kernel"),
        layout_columns(
          column(width = 8,
                 sliderInput("length_scale_2_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
                 sliderInput("length_scale_2_beta", "Inverse-Gamma beta:", 1, 15, 1),
                 actionButton("draw_length_scale_2", "Draw New Length Scale")),
          plotOutput("plot_length_scale_2", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_period_2_choice <- renderUI({
    # !
    # if (input$kernel_label == "Periodic" | input$kernel_label == "Per * SE"| input$kernel_label == "Per + SE") {
    #   card(
    #     card_header("Hyperparameters for period of the second kernel"),
    #     layout_columns(
    #       column(width = 8,
    #              sliderInput("period_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
    #              sliderInput("period_beta", "Inverse-Gamma beta:", 1, 15, 1),
    #              actionButton("draw_period", "Draw New Period")),
    #       plotOutput("plot_period", height = "250px")
    #     )
    #   )
    # }
  })
  
  output$dynamic_roughness_2_choice <- renderUI({
    # !
    # if (input$kernel_label == "Matérn") {
    #   card(
    #     card_header("Hyperparameters for roughness of the second kernel"),
    #     sliderTextInput("nu", "Value of roughness:", choices = c(0.5, 1.5, 2.5, 3.5),
    #                     grid = TRUE, selected = 1.5)
    #   )
    # }
  })
  
  output$dynamic_changepoint_location_choice <- renderUI({
    # !
    # if (input$kernel_label == "Matérn") {
    #   card(
    #     card_header("Hyperparameters for roughness of the second kernel"),
    #     sliderTextInput("nu", "Value of roughness:", choices = c(0.5, 1.5, 2.5, 3.5),
    #                     grid = TRUE, selected = 1.5)
    #   )
    # }
  })
  
  output$dynamic_steepness_choice <- renderUI({
    #!
    # if (input$kernel_label == "Matérn") {
    #   card(
    #     card_header("Hyperparameters for roughness of the second kernel"),
    #     sliderTextInput("nu", "Value of roughness:", choices = c(0.5, 1.5, 2.5, 3.5),
    #                     grid = TRUE, selected = 1.5)
    #   )
    # }
  })
  #-------------------------------------------------------------------------------
  ## Reactive events
  
  # drawing length scale on button click
  length_scale_draw <- eventReactive(input$draw_length_scale, { 
    seed_val <- digest(list(input$length_scale_alpha, input$length_scale_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$ls <- TRUE
    ls <- inverse_gamma(alpha = input$length_scale_alpha, beta = input$length_scale_beta)
    length_scale(ls)
  })
  
  # drawing period on button click
  period_draw <- eventReactive(input$draw_period, { 
    seed_val <- digest(list(input$period_alpha, input$period_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$per <- TRUE
    per <- inverse_gamma(alpha = input$period_alpha, beta = input$period_beta)
    period(per)
  })
  
  # drawing variance on button click
  variance_draw <- eventReactive(input$draw_variance, {
    seed_val <- digest(list(input$variance_df, input$variance_scale), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$var <- TRUE
    var <- half_t(df = input$variance_df, scale = input$variance_scale)
    variance(var) 
  })
  
  #-------------------------------------------------------------------------------
  # second kernel parameters
  # drawing length scale on button click
  length_scale_2_draw <- eventReactive(input$draw_length_scale_2, { 
    seed_val <- digest(list(input$length_scale_2_alpha, input$length_scale_2_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$ls_2 <- TRUE
    ls_2 <- inverse_gamma(alpha = input$length_scale_2_alpha, beta = input$length_scale_2_beta)
    length_scale_2(ls_2)
  })
  
  # drawing period on button click
  period_2_draw <- eventReactive(input$draw_period_2, { 
    seed_val <- digest(list(input$period_2_alpha, input$period_2_beta), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$per_2 <- TRUE
    per_2 <- inverse_gamma(alpha = input$period_2_alpha, beta = input$period_2_beta)
    period_2(per_2)
  })
  
  # drawing variance on button click
  variance_2_draw <- eventReactive(input$draw_variance_2, {
    seed_val <- digest(list(input$variance_2_df, input$variance_2_scale), algo="xxhash32", serialize=TRUE) |> 
      substr(1,7) |> strtoi(base=16)
    set.seed(seed_val)
    
    rv$var_2 <- TRUE
    var_2 <- half_t(df = input$variance_2_df, scale = input$variance_2_scale)
    variance_2(var_2) 
  })
  #-------------------------------------------------------------------------------
  ## GP redraw logic
  gp_pool <- eventReactive(input$draw_gp, {
    
    kernel <- input$kernel_label
    len <- length_scale()
    var <- variance()
    per <- period()
    ro <- as.numeric(input$nu)
    len_2 <- length_scale_2()
    var_2 <- variance_2()
    per_2 <- period_2()
    ro_2 <- as.numeric(input$nu_2)
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
        identical(signif(old_params$length_scale_2,10), signif(len_2,10)) &&
        identical(signif(old_params$variance_2,10),  signif(var_2,10)) &&
        identical(signif(old_params$period_2,10), signif(per_2,10)) &&
        identical(signif(old_params$roughness_2,10), signif(ro_2,10)) &&
        !is.null(old_pool) &&
        is.list(old_pool)) {
      
      return(old_pool)
    }
    
    x_orig <- seq(input$x_range[1], input$x_range[2], length.out = input$n_points)
    
    kernel_params <- hash(
      "variance" = c(var, var_2)
      ,"length_scale" =c(len, len_2) 
      ,"period" = c(per, per_2)
      ,"roughness" = c(ro, ro_2)
    )
    
    funcs <- replicate(
      n_functions,
      simulate_gp(x_orig, input$kernel_label, kernel_params)
    )
    
    new_pool <- list(x_orig = x_orig, funcs = funcs)
    
    last_params(list(kernel_prev = kernel, length_scale = c(len, len_2), variance = c(var, var_2), period = c(per, per_2), roughness = c(ro, ro_2)))
    last_pool(new_pool)
    
    new_pool
  })
  
  # For storing and reproducing draws that were already computed
  gp_data <- reactive({
    req(gp_pool())
    
    idx <- 1:input$nfunc
    idx <- idx[idx <= 100]   # safety
    
    x_new <- seq(input$x_range[1], input$x_range[2], length.out = input$n_points)
    
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
  output$plot_length_scale <- renderPlot({
    observe(length_scale_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$length_scale_alpha
    beta  <- input$length_scale_beta
    
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
  output$plot_variance <- renderPlot({
    observe(variance_draw())
    
    x_seq <- seq(0, 15, length.out = 400)
    df <- input$variance_df
    sc <- input$variance_scale
    
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
  output$plot_period <- renderPlot({
    observe(period_draw())

    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$period_alpha
    beta  <- input$period_beta

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
  #-------------------------------------------------------------------------------
  ## plots
  
  # Inverse-Gamma prior for length_scale plot
  output$plot_length_scale_2 <- renderPlot({
    observe(length_scale_2_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$length_scale_2_alpha
    beta  <- input$length_scale_2_beta
    
    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq) 
    d <- data.frame(x=x_seq, y=dens)
    
    plot <- ggplot(d, aes(x,y)) +
      geom_line(color="steelblue", linewidth=1) +
      labs(title="Inverse-Gamma prior for length scale of the second kernel",
           y="density", x="length scale") +
      theme_minimal(base_size=14)
    
    if(!invalid(length_scale_2())){
      plot <- plot + annotate("point", x = length_scale_2(), y = 0, colour = "red", size = 3) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 hjust = 1.1,  
                 vjust = 1.5,
                 label = paste0("draw = ", signif(length_scale_2(),3)),
                 color="red",
                 size = 5)
    }
    
    plot
  })
  
  # Half-t prior for variance plot
  output$plot_variance_2 <- renderPlot({
    observe(variance_2_draw())
    
    x_seq <- seq(0, 15, length.out = 400)
    df <- input$variance_2_df
    sc <- input$variance_2_scale
    
    dens <- half_t(df = df, scale = sc, x = x_seq)
    
    d <- data.frame(x=x_seq, y=dens)
    
    plot <- ggplot(d, aes(x,y)) +
      geom_line(color="darkgreen", linewidth=1) +
      labs(title="Half-t prior for variance of the second kernel",
           y="density", x="variance") +
      theme_minimal(base_size=14)
    
    if(!invalid(variance_2())){
      plot <- plot + 
        annotate("point", x = variance_2(), y = 0, colour = "red", size = 3) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 hjust = 1.1,  
                 vjust = 1.5,
                 label = paste0("draw = ", signif(variance_2(),3)),
                 color="red",
                 size = 5)
    }
    
    plot
  })
  
  # Inverse-Gamma prior for period plot
  output$plot_period_2 <- renderPlot({
    observe(period_2_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$period_2_alpha
    beta  <- input$period_2_beta
    
    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq) 
    d <- data.frame(x=x_seq, y=dens)
    
    plot <- ggplot(d, aes(x,y)) +
      geom_line(color="steelblue", linewidth=1) +
      labs(title="Inverse-Gamma prior for period of the second kernel",
           y="density", x="period") +
      theme_minimal(base_size=14)
    
    if(!invalid(period_2())){
      plot <- plot + 
        annotate("point", x = period_2(), y = 0, colour = "red", size = 3) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 hjust = 1.1,  
                 vjust = 1.5,
                 label = paste0("draw = ", signif(period_2(),3)),
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
      "variance" = c(variance(), variance_2())
      ,"length_scale" = c(length_scale(), length_scale_2())
      ,"period" = c(period(), period_2())
      ,"roughness" = as.numeric(input$nu)
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
