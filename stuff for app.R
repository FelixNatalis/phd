library(shiny)

ui <- fluidPage(
  titlePanel("Show UI on Button Click"),
  
  # Action button to trigger UI appearance
  actionButton("show_btn", "Show Extra Inputs"),
  
  # Placeholder for dynamic UI
  uiOutput("dynamic_ui")
)

server <- function(input, output, session) {
  
  output$dynamic_ui <- renderUI({
    # Only show UI if button has been clicked at least once
    if (input$show_btn > 0) {
      tagList(
        textInput("name", "Enter your name:"),
        numericInput("age", "Enter your age:", value = 30, min = 1, max = 120),
        actionButton("submit", "Submit")
      )
    }
  })
}

shinyApp(ui, server)


# role <- reactive({
#   input$kernel_label
# })
# 
# observeEvent(role(), {
#   if (role() %in% c("Periodic")) {
#     nav_show("nav", "C")
#   } else if (role() == "Linear") {
#     nav_hide("nav", "C")
#   } else {
#     stop(sprintf("user has unexpected role %s", role()))
#   }
# })



#________________________________________________________________________________________


## GP redraw logic
gp_pool <- eventReactive(input$draw_gp, {
  kernel    <- c(input$kernel_label, input$kernel_label_2)
  len       <- length_scale()   # scalar in single mode, array in multiple mode
  var       <- variance()       # same
  old_params <- last_params()
  old_pool   <- last_pool()
  multi      <- isTRUE(input$multiple_draws_switch)
  
  # reuse pool when parameters unchanged
  if (!is.null(old_params) &&
      is.list(old_params) &&
      identical(old_params$kernel_prev, kernel) &&
      identical(old_params$length_scale, len) &&
      identical(old_params$variance, var) &&
      identical(old_params$multi, multi) &&
      !is.null(old_pool) &&
      is.list(old_pool)) {
    return(old_pool)
  }
  
  x_orig <- seq(input$x_range[1], input$x_range[2], length.out = input$n_points)
  
  if (!multi) {
    # ── Original mode: one param set, n_functions random draws ──────────────
    kernel_params <- hash(
      "kernel_1" = hash(
        "variance"     = var,
        "length_scale" = len,
        ...
      )
    )
    
    funcs <- replicate(
      n_functions,
      simulate_gp(
        x_orig,
        input$is_combination,
        kernel_label,
        kernel_params
      )
    )
    
    new_pool <- list(
      mode   = "single",
      x_orig = x_orig,
      funcs  = funcs          # matrix: n_points × n_functions
    )
    
  } else {
    # ── New mode: n_to_draw param sets, one GP draw each ────────────────────
    n_draw <- min(input$n_to_draw, 100)   # guard against runaway draws
    
    # len / var are now vectors of length n_to_draw (produced by the
    # reactive sliders / samplers upstream)
    stopifnot(length(len) == n_draw, length(var) == n_draw)
    
    funcs <- vapply(seq_len(n_draw), function(i) {
      kernel_params_i <- hash(
        "kernel_1" = hash(
          "variance"     = var[i],
          "length_scale" = len[i],
          ...
        )
      )
      simulate_gp(
        x_orig,
        input$is_combination,
        kernel_label,
        kernel_params_i
      )
    }, numeric(length(x_orig)))
    # result is still an n_points × n_draw matrix, one column per param set
    
    new_pool <- list(
      mode   = "multi",
      x_orig = x_orig,
      funcs  = funcs          # matrix: n_points × n_draw
    )
  }
  
  last_params(list(
    kernel_prev  = kernel,
    length_scale = len,
    variance     = var,
    multi        = multi      # include mode in cache key
  ))
  last_pool(new_pool)
  
  new_pool
})

# For storing and reproducing draws that were already computed
gp_data <- reactive({
  idx <- seq_len(input$nfunc)
  idx <- idx[idx <= 100]
  req(gp_pool())
  pool  <- gp_pool()
  multi <- isTRUE(input$multiple_draws_switch)
  
  x_new <- seq(input$x_range[1], input$x_range[2], length.out = input$n_points)
  
  if (!multi) {
    # ── Original mode: user-selected slice of the pre-drawn pool ────────────

    
    funcs_interp <- apply(pool$funcs[, idx, drop = FALSE], 2, function(f) {
      approx(pool$x_orig, f, xout = x_new)$y
    })
    
    data.frame(
      x    = rep(x_new, length(idx)),
      f    = as.vector(funcs_interp),
      func = rep(idx, each = length(x_new))
    )
    
  } else {
    # ── New mode: every column is one param-set draw; show all of them ──────
    n_draw <- ncol(pool$funcs)
    
    funcs_interp <- apply(pool$funcs, 2, function(f) {
      approx(pool$x_orig, f, xout = x_new)$y
    })
    
    data.frame(
      x    = rep(x_new, n_draw),
      f    = as.vector(funcs_interp),
      func = rep(seq_len(n_draw), each = length(x_new))
    )
  }
})