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

kernel_labels <- c("Squared Exponential", "Matérn", "Linear", "Periodic")

kernel_operation_labels <- c("add", "multiply", "changepoint")

combine_kernels<- function(label_1, label_2, params_1, params_2, operation, x1, x2, params){
  k_1 <- simple_kernel_wrapper(label_1, x1, x2, params = params_1) 
  k_2 <- simple_kernel_wrapper(label_2, x1, x2, params = params_2)
  if(operation == "add"){
    return(k_1 + k_2)
  }else if (operation == "multiply"){
    return(k_1 * k_2)
  }
  else if (operation == "changepoint"){
    # TODO params 
    return(1)
  }
}

kernel_wrapper <- function(is_combination, kernel_label, x1, x2, params){
  if(is_combination){
    extra_params <- params[["extra"]]
    operation <- extra_params[["operation"]]
    additional_params <- extra_params[["additional"]]
    return(combine_kernels(kernel_label[1], kernel_label[2], params[["kernel_1"]], params[["kernel_2"]], operation, x1, x2, additional_params))
  }else{
    return(simple_kernel_wrapper(kernel_label, x1, x2, params[["kernel_1"]]))
  }
}

simple_kernel_wrapper <- function(kernel_label, x1, x2, params){ 
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
}

#-------------------------------------------------------------------------------
# Drawing GP
simulate_gp <- function(x, is_combination, kernel_label, kernel_params, sigma_noise = 1e-3, mean_fun = function(x) 0) {
  K <- kernel_wrapper(is_combination, kernel_label, x, x, params = kernel_params) 
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
  .checkbox { 
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
                               
                               checkboxInput("is_combination", "Use kernel combination?", value = FALSE, width = NULL),
                               uiOutput("kernel_label_block"),
                          
                          actionButton("draw_kernel", "Draw kernel")),
                        plotOutput("plot_kernel", height = "250px")
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
                    uiOutput("dynamic_changepoint_choice"),#!
              )), 
    nav_panel("GP", 
              card(
                card(
                  card_header("Other parameters"),
                  #sliderInput("sigma_n", "Noise amplitude:", 0, 15, 1),
                  sliderInput("nfunc", "Number of Functions:", 1, n_functions, 3),
                  sliderInput("n_points", "Number of X Points:", value = n_points, min = 2, max = 400),
                  sliderInput("x_range", "x range", min = -100, max = 100, value = c(x_min, x_max)), 
                  disabled(actionButton("draw_gp", "Draw GP"))
                ),
                card(plotOutput("plot_gp", height = "600px"))
              ),), 
    id = "nav"),
)

#-------------------------------------------------------------------------------
################################################################################
# SERVER
################################################################################
#-------------------------------------------------------------------------------

server <- function(input, output) {
  
  ## Observable variables
  last_params     <- reactiveVal(NULL)
  last_pool       <- reactiveVal(NULL)
  variance        <- reactiveVal(NULL)
  length_scale    <- reactiveVal(NULL)
  roughness       <- reactiveVal(NULL)
  period          <- reactiveVal(NULL)
  variance_2      <- reactiveVal(NULL)
  length_scale_2  <- reactiveVal(NULL)
  roughness_2     <- reactiveVal(NULL)
  period_2        <- reactiveVal(NULL)
  operation       <- reactiveVal(NULL)
  location        <- reactiveVal(NULL)
  steepness       <- reactiveVal(NULL)

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
  
  output$kernel_label_block <- renderUI({
    if(input$is_combination){
      layout_column_wrap(
        selectInput("kernel_label", "Choose the first kernel:",
                    kernel_labels),
        selectInput("operation", "Choose a combining operation:",
                    kernel_operation_labels),
        selectInput("kernel_label_2", "Choose the second kernel:",
                    kernel_labels))
    }else{
      selectInput("kernel_label", "Choose a kernel:",
                  choices = kernel_labels)
    }
  })
  
  shinyjs::disable("draw_kernel")
  
  condition_parameters_enough<- function(){
    return((input$kernel_label == "Linear" 
            && rv$var 
            || input$kernel_label == "Periodic" 
            && rv$var && rv$ls && rv$per 
            || input$kernel_label == "Squared Exponential" 
            && rv$var && rv$ls 
            || input$kernel_label == "Matérn" 
            && rv$var && rv$ls) 
           && (!input$is_combination ||# conditions on kernel 1
                 input$is_combination && 
                 (input$kernel_label_2 == "Linear" 
                  && rv$var_2 
                  || input$kernel_label_2 == "Periodic" 
                  && rv$var_2 && rv$ls_2 && rv$per_2 
                  || input$kernel_label_2 == "Squared Exponential" 
                  && rv$var_2 && rv$ls_2 
                  || input$kernel_label_2 == "Matérn" 
                  && rv$var_2 && rv$ls_2)))
  }
  
  observe({                   
      if (condition_parameters_enough()) {
        shinyjs::enable("draw_kernel")
      } else {
        shinyjs::disable("draw_kernel")
      }
  })
  
  output$dynamic_length_scale_choice <- renderUI({
    if (!invalid(input$kernel_label) && (input$kernel_label != "Linear")) {
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
    if (!invalid(input$kernel_label) && (input$kernel_label == "Periodic")) {
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
    if (!invalid(input$kernel_label) && (input$kernel_label == "Matérn")) {
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
    if (input$is_combination && !invalid(input$kernel_label_2)) {
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
    if (!invalid(input$kernel_label_2) && (input$kernel_label_2 != "Linear")) {#! different condition
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
    if (!invalid(input$kernel_label_2) && (input$kernel_label_2 == "Periodic")) {
      card(
        card_header("Hyperparameters for period of the second kernel"),
        layout_columns(
          column(width = 8,
                 sliderInput("period_2_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
                 sliderInput("period_2_beta", "Inverse-Gamma beta:", 1, 15, 1),
                 actionButton("draw_period_2", "Draw New Period")),
          plotOutput("plot_period_2", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_roughness_2_choice <- renderUI({
    if (!invalid(input$kernel_label_2) && (input$kernel_label_2 == "Matérn")) {
      card(
        card_header("Hyperparameters for roughness of the second kernel"),
        sliderTextInput("nu_2", "Value of roughness:", choices = c(0.5, 1.5, 2.5, 3.5),
                        grid = TRUE, selected = 1.5)
      )
    }
  })
  
  output$dynamic_changepoint_choice <- renderUI({ 
    if (input$is_combination && !invalid( input$operation )&& input$operation == "changepoint") {
       card(
         card_header("Changepoint kernel parameters"),
         sliderInput("location", "Location of changepoint:", min = -100, max = 100, value = 0),
         sliderInput("steepness", "Steepness of changepoint:", 0.1, 10, 0.1)
       )
     }
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
    tryCatch({
    kernel <- c(input$kernel_label, input$kernel_label_2)
    len <- length_scale()
    var <- variance()
    per <- period()
    ro <- as.numeric(input$nu)
    len_2 <- length_scale_2()
    var_2 <- variance_2()
    per_2 <- period_2()
    ro_2 <- as.numeric(input$nu_2)
    loc <- location()
    ste <- steepness()
    ope <- operation()
    old_params <- last_params()
    old_pool   <- last_pool()
   
    # cat(paste("\n kernel_prev", old_params$kernel_prev, "\n")) 
    # cat(paste("\n length_scale", old_params$length_scale, "\n"))
    # cat(paste("\n variance", old_params$variance, "\n"))
    # cat(paste("\n period", old_params$period, "\n"))
    # cat(paste("\n roughness", old_params$roughness, "\n"))
    
    # reuse pool when parameters unchanged
    if (!is.null(old_params) &&
        is.list(old_params) &&
        identical(old_params$kernel_prev, kernel) &&
        identical(old_params$length_scale, len) &&
        identical(old_params$variance,  var) &&
        identical(old_params$period, per) &&
        identical(old_params$roughness, ro) &&
        identical(old_params$length_scale_2, len_2) &&
        identical(old_params$variance_2,  var_2) &&
        identical(old_params$period_2, per_2) &&
        identical(old_params$roughness_2, ro_2) &&
        identical(old_params$location, loc) &&
        identical(old_params$steepness, ste) &&
        identical(old_params$operation, ope) &&
        !is.null(old_pool) &&
        is.list(old_pool)) {
      
      return(old_pool)
    }
    
    x_orig <- seq(input$x_range[1], input$x_range[2], length.out = input$n_points)
    
    if(input$is_combination){
      kernel_label <- c(input$kernel_label, input$kernel_label_2)
    }else{
      kernel_label <- input$kernel_label
    }
    
    kernel_params <- hash(
      "kernel_1" = hash(
        "variance" = var
        ,"length_scale" = len
        ,"period" = per
        ,"roughness" = ro
      ),
      "kernel_2" = hash(
        "variance" = var_2
        ,"length_scale" = len_2
        ,"period" = per_2
        ,"roughness" = ro_2
      ),
      "extra" = hash(
        "operation" = input$operation,
        "additional" = hash(
          "location" = input$location,
          "steepness" = input$steepness
        )))
    
    funcs <- replicate(
      n_functions,
      simulate_gp(x_orig, input$is_combination,kernel_label, kernel_params)
    )
    
    new_pool <- list(x_orig = x_orig, funcs = funcs)
    
    last_params(list(kernel_prev = kernel, length_scale = len, variance = var, period = per, roughness = ro
                     ,length_scale_2 = len_2, variance_2 = var_2, period_2 = per_2, roughness_2 = ro_2
                     ,location = loc, steepness = ste, operation = ope
                     ))
    last_pool(new_pool)
    
    new_pool
    }, error=function(e) {
      cat(paste("\nerror in gp draw\n",e,"\n","\n"  #TODO
                ))
    }, warning=function(w) {
      cat(paste("\nwarning in gp draw\n"))
    })
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
  
  param_plot<- function(dataframe, title, y_label, x_label, observation_value){
    
    plot <- ggplot(dataframe, aes(x,y)) +
      geom_line(color="steelblue", linewidth=1) +
      labs(title=title,
           y=y_label, x=x_label) +
      theme_minimal(base_size=14)
    
    if(!invalid(observation_value)){
      plot <- plot + 
        annotate("point", x = observation_value, y = 0, colour = "red", size = 3) +
        annotate("text",
                 x = Inf,
                 y = Inf,
                 hjust = 1.1,  
                 vjust = 1.5,
                 label = paste0("draw = ", signif(observation_value,3)),
                 color="red",
                 size = 5)
    }
    
    plot
    
  }
  
  # Inverse-Gamma prior for length_scale plot
  output$plot_length_scale <- renderPlot({
    observe(length_scale_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$length_scale_alpha
    beta  <- input$length_scale_beta
    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq)
    d <- data.frame(x=x_seq, y=dens)
    
    param_plot(d, "Inverse-Gamma prior for length scale", "density", "length scale", length_scale())
  })
  
  # Half-t prior for variance plot
  output$plot_variance <- renderPlot({
    observe(variance_draw())
    
    x_seq <- seq(0, 15, length.out = 400)
    df <- input$variance_df
    sc <- input$variance_scale
    dens <- half_t(df = df, scale = sc, x = x_seq)
    d <- data.frame(x=x_seq, y=dens)
    
    param_plot(d, "Half-t prior for variance", "density", "variance", variance())
  })
  
  # Inverse-Gamma prior for period plot
  output$plot_period <- renderPlot({
    observe(period_draw())

    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$period_alpha
    beta  <- input$period_beta

    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq) 
    d <- data.frame(x=x_seq, y=dens)

    param_plot(d, "Inverse-Gamma prior for period", "density", "period", period())
 })
  #-------------------------------------------------------------------------------
  ## plots for second kernel parameters
  
  # Inverse-Gamma prior for length_scale plot
  output$plot_length_scale_2 <- renderPlot({
    observe(length_scale_2_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$length_scale_2_alpha
    beta  <- input$length_scale_2_beta
    
    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq) 
    d <- data.frame(x=x_seq, y=dens)
    
    param_plot(d, "Inverse-Gamma prior for length scale of the second kernel", "density", "length scale", length_scale_2())    
  })
  
  # Half-t prior for variance plot
  output$plot_variance_2 <- renderPlot({
    observe(variance_2_draw())
    
    x_seq <- seq(0, 15, length.out = 400)
    df <- input$variance_2_df
    sc <- input$variance_2_scale
    dens <- half_t(df = df, scale = sc, x = x_seq)
    d <- data.frame(x=x_seq, y=dens)
    
    param_plot(d, "Half-t prior for variance of the second kernel", "density", "variance", variance_2())
  })
  
  # Inverse-Gamma prior for period plot
  output$plot_period_2 <- renderPlot({
    observe(period_2_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$period_2_alpha
    beta  <- input$period_2_beta
    dens <- inverse_gamma(alpha = alpha, beta = beta, x = x_seq) 
    d <- data.frame(x=x_seq, y=dens)
    
    param_plot(d, "Inverse-Gamma prior for period of the second kernel", "density", "period", period_2())
  })
  
 # Kernel based on distance plot
  output$plot_kernel <- renderPlot({
    req(input$draw_kernel)
    tryCatch({
      if(condition_parameters_enough()){
    
    dist <- seq(-3, 3, length.out = 300)
    x_o <- rep(0, length(dist))

    if(input$is_combination){
      kernel_label <- c(input$kernel_label, input$kernel_label_2)
    }else{
      kernel_label <- input$kernel_label
    }
    
    
    kernel_params <- hash(
      "kernel_1" = hash(
        "variance" = variance()
        ,"length_scale" = length_scale()
        ,"period" = period()
        ,"roughness" = as.numeric(input$nu)
      ),
      "kernel_2" = hash(
        "variance" = variance_2()
        ,"length_scale" = length_scale_2()
        ,"period" = period_2()
        ,"roughness" = as.numeric(input$nu_2)
      ),
      "extra" = hash(
        "operation" = input$operation,
        "additional" = hash(
          "location" = input$location,
          "steepness" = input$steepness
        )))
    
    k <- kernel_wrapper(input$is_combination, kernel_label, dist, x_o, kernel_params)[,1]
    #cat(paste("\n",length(dist), " ", length(k), "\n"))
    kernel_title <- paste(input$kernel_label, input$operation, input$kernel_label_2)
    ggplot(data.frame(dist=dist, k=k), aes(dist,k)) +
      geom_line(color="purple", linewidth=1.2) +
      labs(title=paste(kernel_title, "Kernel"), x="Distance", y="k(x)") +
      theme_minimal(base_size=14)
      }
    }, error=function(e) {
      cat(paste("\nError in kernel plot\n", e, "\n"))
    }, warning=function(w) {
      cat(paste("\nWarning in kernel plot\n", w, "\n"))
    })
  })
  
  
  # GP prior draws plot
  output$plot_gp <- renderPlot({
    tryCatch({
    req(input$draw_gp)
    observe(gp_data())#, period_draw(), length_scale_draw(), variance_draw())
      kernel_title <- paste(input$kernel_label, input$operation, input$kernel_label_2)
    ggplot(gp_data(), aes(x=x, y=f, group=func, color=factor(func))) +
      geom_line(alpha=0.9, linewidth=1) +
      scale_color_discrete(guide="none") +
      labs(title="Gaussian Process Prior Samples",
           subtitle= paste(kernel_title, "Kernel"),
           x="x", y="f(x)") +
      theme_minimal(base_size=16)
    }, error=function(e) {
      cat(paste("\nerror in gp plot\n", e, "\n"))
    }, warning=function(w) {
      cat(paste("\nwarning in gp plot\n", w, "\n"))
    })
  })
  
}

shinyApp(ui, server)
