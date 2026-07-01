
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
library(lineqGPR)
library(shinysurveys)

source("GreatPlains.R")

#-------------------------------------------------------------------------------
################################################################################
# GLOBAL CONSTANTS
################################################################################
#-------------------------------------------------------------------------------

# defaults parameters for GP prior predictive draws
n_functions = 10
n_func = 1 # how many GP functions are drawn per parameter combinations
x_min = -10
x_max = 10
n_points = 200
epsilon = 1e-6
n_draws = 6
colors = c( "indianred1", "goldenrod2", "seagreen2", "turquoise2","royalblue1", "mediumorchid1")

#-------------------------------------------------------------------------------
################################################################################
# UI
################################################################################
#-------------------------------------------------------------------------------

ui <- page_fillable(
  tags$style(
    type = 'text/css',
    ".selectize-input {
    font-size:80%;
    line-height: 16px;
  }
  .checkbox {
      font-size:85%;
      line-height: 16px;
      padding: 4px;
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
    font-size:85%;
  }
  .btn {
    padding: 2px 6px;
    font-size: 80%;
  }
  .card-body {
    display: inline-block;
  }
"
  ),
  
  useShinyjs(),
  navset_tab(
    nav_panel(
      "0. Priors",
        card(card_header("1. General"),
          layout_columns(
            col_widths = c(6, 6),
            gap = "0.1rem",
            column(
              width = 10,
          checkboxInput("is_formula", "Display hints and kernel formula?", value = FALSE),
          sliderInput("n_to_draw", "Number of parameter sets to draw:", 1, n_draws, 1),
          numberInput("seed_value", "Choose seed value:", width = 120), 
          actionButton("set_seed", "Set seed")),
          column(width = 10,
                 div(
                  style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
                   uiOutput("hint_general")
                 ) 
)
        )),
       card(
         card_header("2. Kernel"),
         layout_columns(
           col_widths = c(3, 3, 6),
           gap = "0.1rem",
           column(
             width = 10,
             checkboxInput("is_combination", "Use kernel combination?", value = FALSE),
             uiOutput("kernel_label_block"),
             actionButton("draw_kernel", "4. Draw kernel")
           ),
           column(
             width = 10,
             div(
               style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
               uiOutput("kernel_formula_block")
             )
           ),
           plotOutput("plot_kernel", height = "250px")
         )
       ),
card(card_header("3. Kernel parameters",
        card(
          card_header("Magnitude (σ)", style = 'padding:4px; font-size:90%'),
          layout_columns(
            col_widths = c(3, 3, 6),
            gap = "0.1rem",
            column(
              width = 10,
              sliderInput("magnitude_df", "Half-t degrees of freedom:", 1, 5, 4),
              sliderInput("magnitude_scale", "Half-t scale:", 1, 15, 1),
              actionButton("draw_magnitude", "Draw magnitude"),
              checkboxInput(inputId = "is_fix_magnitude", "Fix magnitude at a constant value?", value = FALSE),
              sliderInput(inputId = "magnitude_fixed_value", label = NULL, value = 0, min = 0, max = 5, step = 0.1)
            ),
            div(
              style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
              uiOutput("hint_magnitude")
            )
            ,
            plotOutput("plot_magnitude", height = "250px")
          )
        ),
        uiOutput("dynamic_length_scale_choice"),
        uiOutput("dynamic_roughness_choice"),
        uiOutput("dynamic_period_choice"),
        uiOutput("dynamic_magnitude_2_choice"),
        uiOutput("dynamic_length_scale_2_choice"),
        uiOutput("dynamic_roughness_2_choice"),
        uiOutput("dynamic_period_2_choice"),
        uiOutput("dynamic_changepoint_choice"),
      )
    )),
    nav_panel(
      "5. Constraints",
      card(
        card_header("Function characteristics"),
        div(
          style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
          uiOutput("hint_constraints")
        ),
        layout_columns(
          column(
            width = 8,
            checkboxInput(
              "is_upper_bound",
              "Upper Y-scale bound",
              value = FALSE,
              width = NULL
            ),
            div(
              style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
              uiOutput("hint_upper_bound")
            ),
            uiOutput("dynamic_upper_bound_choice")
          ),
          column(
            width = 8,
            checkboxInput(
              "is_lower_bound",
              "Lower Y-scale bound",
              value = FALSE,
              width = NULL
            ),
            div(
              style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
              uiOutput("hint_lower_bound")
            ),
            uiOutput("dynamic_lower_bound_choice")
          ),
          column(
            width = 8,
            checkboxInput(
              "is_monotonicity",
              "Non-decreasing monotonicity",
              value = FALSE,
              width = NULL
            ),
            div(
              style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
              uiOutput("hint_monotonicity")
            )#,
        #    uiOutput("dynamic_monotonicity_choice")
          ),
          column(
            width = 8,
            checkboxInput(
              "is_convexity",
              "Convexity",
              value = FALSE,
              width = NULL
            ),
        
          div(
            style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
            uiOutput("hint_convexity")
          ))
        )
      ),
      card(
        card_header("Specific function values"),
        div(
          style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
          uiOutput("hint_values")
        ),
        layout_columns(column(width = 2,numberInput("x_1", "x₁")),column(width = 2, numberInput("y_1", "y₁"))),
        layout_columns(column(width = 2,numberInput("x_2", "x₂")), column(width = 2,numberInput("y_2", "y₂"))),
        layout_columns(column(width = 2,numberInput("x_3", "x₃")),column(width = 2, numberInput("y_3", "y₃")))
      )
    ),
    
    nav_panel("6. GP", card(
      card(
        card_header("Other parameters"),
        #sliderInput("sigma_n", "Noise amplitude:", 0, 15, 1),
        sliderInput(
          "n_points",
          "Number of grid points on the X-scale:",
          value = n_points,
          min = 2,
          max = 400
        ),
        sliderInput(
          "x_range",
          "Range of the X-scale:",
          min = -100,
          max = 100,
          value = c(x_min, x_max)
        ),
        disabled(actionButton("draw_gp", "Draw GP"))
      ),
      card(shinycssloaders::withSpinner(
        plotOutput("plot_gp", height = "600px"), type = 5
      ))
    ), ),
    id = "nav"
  ),
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
  magnitude        <- reactiveVal(NULL)
  length_scale    <- reactiveVal(NULL)
  roughness       <- reactiveVal(NULL)
  period          <- reactiveVal(NULL)
  magnitude_2      <- reactiveVal(NULL)
  length_scale_2  <- reactiveVal(NULL)
  roughness_2     <- reactiveVal(NULL)
  period_2        <- reactiveVal(NULL)
  operation       <- reactiveVal(NULL)
  location        <- reactiveVal(NULL)
  steepness       <- reactiveVal(NULL)
  x_fixed         <- reactiveVal(NULL)
  y_fixed         <- reactiveVal(NULL)
  draws           <- reactiveVal(NULL)
  
  rv <- reactiveValues(
    mag = FALSE,
    ls  = FALSE,
    per = FALSE,
    mag_2 = FALSE,
    ls_2  = FALSE,
    per_2 = FALSE
  )
  
  condition_parameters_check <- function() {
    return((
      input$kernel_label == "Linear"
      && rv$mag
      || input$kernel_label == "Periodic"
      && rv$mag && rv$ls && rv$per
      || input$kernel_label == "Squared Exponential"
      && rv$mag && rv$ls
      || input$kernel_label == "Matérn"
      && rv$mag && rv$ls
      )
    && (
      !input$is_combination || # conditions on kernel 1
        input$is_combination &&
        (
          input$kernel_label_2 == "Linear"
          && rv$mag_2
          || input$kernel_label_2 == "Periodic"
          && rv$mag_2 && rv$ls_2 && rv$per_2
          || input$kernel_label_2 == "Squared Exponential"
          && rv$mag_2 && rv$ls_2
          || input$kernel_label_2 == "Matérn"
          && rv$mag_2 && rv$ls_2
        )
    )
    )
  }
  
  constrained_check <- function() {
    return(
      input$is_upper_bound ||
        input$is_lower_bound ||
        input$is_monotonicity || 
        input$is_convexity
   #   || (input$x_1 && input$y_1) 
    #  || (input$x_2 && input$y_2) 
     # || (input$x_3 && input$y_3) 
    )
  }
  
  condition_constrained_parameters_check <- function() {
    # TODO check that everything is filled
    return(
      input$is_upper_bound ||
        input$is_lower_bound ||
        input$is_monotonicity || input$is_convexity
    )
  }
  
  observeEvent(input$set_seed, {
    set.seed(input$seed_value)
  })
  
#-------------------------------------------------------------------------------
## Dynamic ui elements
  
# Hints and text
  
  output$kernel_label_block <- renderUI({
    if (input$is_combination) {
      layout_column_wrap(
        selectInput("kernel_label", "Choose the first kernel:", kernel_labels),
        selectInput(
          "operation",
          "Choose a combining operation:",
          kernel_operation_labels
        ),
        selectInput(
          "kernel_label_2",
          "Choose the second kernel:",
          kernel_labels
        )
      )
    } else{
      selectInput("kernel_label", "Choose a kernel:", choices = kernel_labels)
    }
  })
  
  output$hint_general <- renderUI({
    if (input$is_formula) {
      general <- "This app helps you explore the prior distribution formulation for Gaussian processes. Follow these steps:"
      step_1<-"1. Set the general parameters in this block"
      step_2<-"2. Choose the kernel structure in the \"Kernel\" block"
      step_3<-"3. In the \"Kernel parameters\" block, fill in the parameters for the chosen kernel, selecting prior hyperparameters and pressing the \"Draw parameter\" buttons or setting them as fixed values via checking the corresponding check boxes"
      step_4<-"4. Examine the kernel by pressing the \"Draw kernel\" button"
      step_5<-"5. Navigate to the \"Constraints\" page and set constraints (optional)"
      step_6<-"6. Navigate to the \"GP\" page, adjust GP drawing settings and draw GPs by pressing the \"Draw GP\" button"
      
      withMathJax(
        p(general),
        p(step_1),
        p(step_2),
        p(step_3),
        p(step_4),
        p(step_5),
        p(step_6)
      )
      }
  })
  
  output$hint_magnitude <- renderUI({
    if (input$is_formula) {
        hint <- "Magnitude refers to the amplitude in the \\(y\\) values. The prior is set on the \\(standard\\) \\(deviation\\), which is squared in most kernels, making prior values close to \\(0\\) shift towards \\(0\\) even more."
        
        withMathJax(
          p(hint)
        )
    }
  })
  
  output$hint_length_scale <- renderUI({
    if (input$is_formula) {
      hint <- "Length scale refers to the scale of the \\(x\\) values after which the function fluctuates. Shorter length scales correspond to more wiggly functions, and longer - to flatter functions. The length scale is always relative to the scale of \\(x\\)."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_period <- renderUI({
    if (input$is_formula) {
      hint <- "Period determines the scale of the \\(x\\) values after which the periodic function repeats its pattern. Shorter periods correspond to more frequent oscillations, and longer - to slower oscillations. The period is always relative to the scale of \\(x\\)."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_roughness <- renderUI({
    if (input$is_formula) {
      hint <- "Roughness of the Matérn kernel represents how smooth the function is. Higher values correspond to higher differentiability and smoothness, and lower produce more jagged functions. As \\(\\nu\\) approaches infinity, Matérn kernel converges to the Squared Exponential kernel."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_roughness_2 <- renderUI({
    if (input$is_formula) {
      hint <- "Roughness of the Matérn kernel represents how smooth the function is. Higher values correspond to higher differentiability and smoothness, and lower produce more jagged functions. As \\(\\nu\\) approaches infinity, Matérn kernel converges to the Squared Exponential kernel."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_period_2 <- renderUI({
    if (input$is_formula) {
      hint <- "Period determines the scale of the \\(x\\) values after which the periodic function repeats its pattern. Shorter periods correspond to more frequent oscillations, and longer - to slower oscillations. The period is always relative to the scale of \\(x\\)."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_length_scale_2 <- renderUI({
    if (input$is_formula) {
      hint <- "Length scale refers to the scale of the \\(x\\) values after which the function fluctuates. Shorter length scales correspond to more wiggly functions, and longer - to flatter functions. The length scale is always relative to the scale of \\(x\\)."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_magnitude_2 <- renderUI({
    if (input$is_formula) {
      hint <- "Magnitude refers to the amplitude in the \\(y\\) values. The prior is set on the \\(standard\\) \\(deviation\\), which is squared in most kernels, making prior values close to \\(0\\) shift towards \\(0\\) even more."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_changepoint <- renderUI({
    if (input$is_formula) {
      hint <- "Parameters of the changepoint kernel govern how the two base kernels replace each other. Location of the changepoint determines the \\(x\\) value where the change in base kernels takes place. Steepness regulates how rapidly the fade-out of the first kernel and the fade-in of the second kernel happen."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_lower_bound <- renderUI({
    if (input$is_formula) {
      hint <- "Lower bound constraint sets the minimum value that the \\(y\\) is allowed to take."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_upper_bound <- renderUI({
    if (input$is_formula) {
      hint <- "Upper bound constraint sets the maximum value that the \\(y\\) is allowed to take."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_monotonicity <- renderUI({
    if (input$is_formula) {
      hint <- "Non-decreasing monotonicity constraint forces the function to only increase its values or keep them at the same level. This also means that the \\(first\\) \\(derivative\\) of the function is \\(equal\\) \\(to\\) \\(or\\) \\(greater\\) \\(than\\) \\(0\\)"
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_convexity <- renderUI({
    if (input$is_formula) {
      hint <- "Convexity forces the \\(second\\) \\(derivative\\) of the function to be \\(equal\\) \\(to\\) \\(or\\) \\(greater\\) \\(than\\) \\(0\\)."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_constraints <- renderUI({
    if (input$is_formula) {
      hint <- "Constraints allow to fix the function specifically in qualitative ways informed by the domain knowledge. Due to the limitations of the package used for this, setting constraints is only possible for the \\(Squared\\) \\(exponential\\) \\(kernel\\), it only allows the \\(x\\) range of \\([0;1]\\), and only makes one GP draw."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$hint_values <- renderUI({
    if (input$is_formula) {
      hint <- "It is possible to set several specific points as \\((x, y)\\) pairs."
      
      withMathJax(
        p(hint)
      )
    }
  })
  
  output$kernel_formula_block <- renderUI({
    if (input$is_formula) {
      if (input$is_combination) {
        kernel_1_formula <- paste("\\(k_1(x, x') = ", kernel_formulae[[input$kernel_label]])
        kernel_2_formula <- paste("\\(k_2(x, x') = ", kernel_formulae[[input$kernel_label_2]])
        operation_formula <- kernel_operation_formulae[[input$operation]]
        
        withMathJax(
          p(operation_formula),
          p(kernel_1_formula),
          p(kernel_2_formula)
        )
      } else {
        kernel_1_formula <- paste("\\(k(x, x') = ", kernel_formulae[[input$kernel_label]])
        
        withMathJax(
          p(kernel_1_formula)
        )
      }
    }
  })
  
#-------------------------------------------------------------------------------
# Enabling and disabling elements

  shinyjs::disable("draw_kernel")
  
  observe({
    if (condition_parameters_check()) {
      shinyjs::enable("draw_kernel")
    } else {
      shinyjs::disable("draw_kernel")
    }
  })
  
  shinyjs::disable("magnitude_fixed_value")
  shinyjs::disable("length_scale_fixed_value")
  shinyjs::disable("period_fixed_value")
  shinyjs::disable("magnitude_2_fixed_value")   
  shinyjs::disable("period_2_fixed_value")
  shinyjs::disable("length_scale_2_fixed_value")
  
  observe({
    if (input$is_fix_magnitude) {
      shinyjs::disable("draw_magnitude")
      shinyjs::disable("magnitude_df")
      shinyjs::disable("magnitude_scale")
      shinyjs::enable("magnitude_fixed_value")
    } else {
      shinyjs::enable("draw_magnitude")
      shinyjs::enable("magnitude_df")
      shinyjs::enable("magnitude_scale")
      shinyjs::disable("magnitude_fixed_value")
    }
  })
  
  observe({
    req(!is.null(input$is_fix_length_scale))
    if (input$is_fix_length_scale) {
      shinyjs::disable("draw_length_scale")
      shinyjs::disable("length_scale_alpha")
      shinyjs::disable("length_scale_beta")
      shinyjs::enable("length_scale_fixed_value")
    } else {
      shinyjs::enable("draw_length_scale")
      shinyjs::enable("length_scale_alpha")
      shinyjs::enable("length_scale_beta")
      shinyjs::disable("length_scale_fixed_value")
    }
  })
  
  observe({
    req(!is.null(input$is_fix_period))
    if (input$is_fix_period) {
      shinyjs::disable("draw_period")
      shinyjs::disable("period_alpha")
      shinyjs::disable("period_beta")
      shinyjs::enable("period_fixed_value")
    } else {
      shinyjs::enable("draw_period")
      shinyjs::enable("period_alpha")
      shinyjs::enable("period_beta")
      shinyjs::disable("period_fixed_value")
    }
  })
  
  observe({
    req(!is.null(input$is_fix_magnitude_2))
    if (input$is_fix_magnitude_2) {
      shinyjs::disable("draw_magnitude_2")
      shinyjs::disable("magnitude_2_df")
      shinyjs::disable("magnitude_2_scale")
      shinyjs::enable("magnitude_2_fixed_value")
    } else {
      shinyjs::enable("draw_magnitude_2")
      shinyjs::enable("magnitude_2_df")
      shinyjs::enable("magnitude_2_scale")
      shinyjs::disable("magnitude_2_fixed_value")
    }
  })
  
  observe({
    req(!is.null(input$is_fix_length_scale_2))
    if (input$is_fix_length_scale_2) {
      shinyjs::disable("draw_length_scale_2")
      shinyjs::disable("length_scale_2_alpha")
      shinyjs::disable("length_scale_2_beta")
      shinyjs::enable("length_scale_2_fixed_value")
    } else {
      shinyjs::enable("draw_length_scale_2")
      shinyjs::enable("length_scale_2_alpha")
      shinyjs::enable("length_scale_2_beta")
      shinyjs::disable("length_scale_2_fixed_value")
    }
  })
  
  observe({
    req(!is.null(input$is_fix_period_2))
    if (input$is_fix_period_2) {
      shinyjs::disable("draw_period_2")
      shinyjs::disable("period_2_alpha")
      shinyjs::disable("period_2_beta")
      shinyjs::enable("period_2_fixed_value")
    } else {
      shinyjs::enable("draw_period_2")
      shinyjs::enable("period_2_alpha")
      shinyjs::enable("period_2_beta")
      shinyjs::disable("period_2_fixed_value")
    }
  })
  
  observe({
    if (constrained_check()) {
      shinyjs::disable("x_range")
    }
  })
  
#-------------------------------------------------------------------------------
# Parameters for the first kernel
  
  output$dynamic_length_scale_choice <- renderUI({
    if (!invalid(input$kernel_label) &&
        (input$kernel_label != "Linear")) {
      card(
        card_header("Length scale (λ)", style = 'padding:4px; font-size:90%'),
        layout_columns(
          col_widths = c(3, 3, 6),
          gap = "0.1rem",
          column(
            width = 10,
            sliderInput("length_scale_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
            sliderInput("length_scale_beta", "Inverse-Gamma beta:", 1, 15, 1),
            actionButton("draw_length_scale", "Draw length scale"),
            checkboxInput(inputId = "is_fix_length_scale", "Fix length scale at a constant value?", value = FALSE),
            sliderInput(inputId = "length_scale_fixed_value", label = NULL, value = 0, min = 0, max = 5, step = 0.1)
          ),
          div(
            style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
            uiOutput("hint_length_scale")
          ),
          plotOutput("plot_length_scale", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_period_choice <- renderUI({
    if (!invalid(input$kernel_label) &&
        (input$kernel_label == "Periodic")) {
      card(
        card_header("Period (p)", style = 'padding:4px; font-size:90%'),
        layout_columns(
          col_widths = c(3, 3, 6),
          gap = "0.1rem",
          column(
            width = 10,
            sliderInput("period_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
            sliderInput("period_beta", "Inverse-Gamma beta:", 1, 15, 1),
            actionButton("draw_period", "Draw period"),
            checkboxInput(inputId = "is_fix_period", "Fix period at a constant value?", value = FALSE),
            sliderInput(inputId = "period_fixed_value", label = NULL, value = 0, min = 0, max = 5, step = 0.1)
          ),
          div(
            style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
            uiOutput("hint_period")
          ),
          plotOutput("plot_period", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_roughness_choice <- renderUI({
    if (!invalid(input$kernel_label) &&
        (input$kernel_label == "Matérn")) {
      card(
        card_header("Roughness (ν)", style = 'padding:4px; font-size:90%'),
        layout_columns(
          col_widths = c(3, 3),
          gap = "0.1rem",
          column(
            width = 10,
        sliderTextInput(
          "nu",
          "Value of roughness:",
          choices = c(0.5, 1.5, 2.5, 3.5),
          grid = TRUE,
          selected = 1.5
        )
          ),
        div(
          style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
          uiOutput("hint_roughness")
        )
        )
      )
    }
  })

#-------------------------------------------------------------------------------
# Parameters for the second kernel
  
  output$dynamic_magnitude_2_choice <- renderUI({
    if (input$is_combination && !invalid(input$kernel_label_2)) {
      card(
        card_header("Magnitude (σ) of the second kernel", style = 'padding:4px; font-size:90%'),
        layout_columns(
          col_widths = c(3, 3, 6),
          gap = "0.1rem",
          column(
            width = 10,
            sliderInput("magnitude_2_df", "Half-t degrees of freedom:", 1, 5, 4),
            sliderInput("magnitude_2_scale", "Half-t scale:", 1, 15, 1),
            actionButton("draw_magnitude_2", "Draw magnitude"),
            checkboxInput(inputId = "is_fix_magnitude_2", "Fix magnitude at a constant value?", value = FALSE),
            sliderInput(inputId = "magnitude_2_fixed_value", label = NULL, value = 0, min = 0, max = 5, step = 0.1)
          ),
          div(
            style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
            uiOutput("hint_magnitude_2")
          )
          ,
          plotOutput("plot_magnitude_2", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_length_scale_2_choice <- renderUI({
    if (!invalid(input$kernel_label_2) &&
        (input$kernel_label_2 != "Linear")) {
      card(
        card_header("Length scale (λ) of the second kernel", style = 'padding:4px; font-size:90%'),
        layout_columns(
          col_widths = c(3, 3, 6),
          gap = "0.1rem",
          column(
            width = 10,
            sliderInput("length_scale_2_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
            sliderInput("length_scale_2_beta", "Inverse-Gamma beta:", 1, 15, 1),
            actionButton("draw_length_scale_2", "Draw length scale"),
            checkboxInput(inputId = "is_fix_length_scale_2", "Fix length scale at a constant value?", value = FALSE),
            sliderInput(inputId = "length_scale_2_fixed_value", label = NULL, value = 0, min = 0, max = 5, step = 0.1)
          ),
          div(
            style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
            uiOutput("hint_length_scale_2")
          ),
          plotOutput("plot_length_scale_2", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_period_2_choice <- renderUI({
    if (!invalid(input$kernel_label_2) &&
        (input$kernel_label_2 == "Periodic")) {
      card(
        card_header("Period (p) of the second kernel", style = 'padding:4px; font-size:90%'),
        layout_columns(
          col_widths = c(3, 3, 6),
          gap = "0.1rem",
          column(
            width = 10,
            sliderInput("period_2_alpha", "Inverse-Gamma alpha:", 1, 15, 1),
            sliderInput("period_2_beta", "Inverse-Gamma beta:", 1, 15, 1),
            actionButton("draw_period_2", "Draw period"),
            checkboxInput(inputId = "is_fix_period_2", "Fix period at a constant value?", value = FALSE),
            sliderInput(inputId = "period_2_fixed_value", label = NULL, value = 0, min = 0, max = 5, step = 0.1)
          ),
          div(
            style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
            uiOutput("hint_period_2")
          ),
          plotOutput("plot_period_2", height = "250px")
        )
      )
    }
  })
  
  output$dynamic_roughness_2_choice <- renderUI({
    if (!invalid(input$kernel_label_2) &&
        (input$kernel_label_2 == "Matérn")) {
      card(
        card_header("Roughness (ν) of the second kernel", style = 'padding:4px; font-size:90%'),
        layout_columns(
          col_widths = c(3, 3),
          gap = "0.1rem",
          column(
            width = 10,
        sliderTextInput(
          "nu_2",
          "Value of roughness:",
          choices = c(0.5, 1.5, 2.5, 3.5),
          grid = TRUE,
          selected = 1.5
        )
        ),
        div(
          style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
          uiOutput("hint_roughness_2")
        )
        )
      )
    }
  })
  
  output$dynamic_changepoint_choice <- renderUI({
    if (input$is_combination &&
        !invalid(input$operation) &&
        input$operation == "changepoint") {
      card(
        card_header("Changepoint parameters"),
        layout_columns(
          col_widths = c(3, 3),
          gap = "0.1rem",
          column(
            width = 10,
        sliderInput(
          "location",
          "Location of changepoint (x₀):",
          min = -100,
          max = 100,
          value = 0
        ),
        sliderInput("steepness", "Steepness of changepoint (s):", 0.1, 10, 0.1)
      ),
      div(
        style = "font-size: 70%; overflow: hidden; overflow-x: auto; word-break: break-word; max-width: 100%;",
        uiOutput("hint_changepoint")
      )
      )
      )
    }
  })
  
#-------------------------------------------------------------------------------
# Constraint parameters
  
  output$dynamic_upper_bound_choice <- renderUI({
    if (input$is_upper_bound) {
      numericInput(
        "upper_bound",
        "",
        value = 1,
        min = -10,
        max = 10
      )
    }
  })
  
  output$dynamic_lower_bound_choice <- renderUI({
    if (input$is_lower_bound) {
      numericInput(
        "lower_bound",
        "",
        value = 0,
        min = -10,
        max = 10
      )
    }
  })
  
  output$dynamic_monotonicity_choice <- renderUI({
    if (input$is_monotonicity) {
      selectInput(
        "monotonicity",
        label = "Monotonicity type",
        choices = c("nondecreasing", "nonincreasing")
      )
    }
  })
  
#-------------------------------------------------------------------------------
## Reactive events
  
  # drawing length scale on button click
  length_scale_draw <- eventReactive(input$draw_length_scale, {
    draws <- c()
    for (i in 1:input$n_to_draw) {
      draw <- inverse_gamma(alpha = input$length_scale_alpha,
                            beta = input$length_scale_beta)
      draws <- append(draws, draw)
    }
    rv$ls <- TRUE
    length_scale(draws)
  })
  
  # fixing length scale in the constant regime
  length_scale_fixed <- observeEvent(input$length_scale_fixed_value, {
    req(input$is_fix_length_scale)
    draws <- rep(input$length_scale_fixed_value, input$n_to_draw)
    rv$ls <- TRUE
    length_scale(draws)
  })
  
  # drawing period on button click
  period_draw <- eventReactive(input$draw_period, {
    draws <- c()
    for (i in 1:input$n_to_draw) {
      draw <- inverse_gamma(alpha = input$period_alpha,
                            beta = input$period_beta)
      draws <- append(draws, draw)
    }
    rv$per <- TRUE
    period(draws)
  })
  
  # fixing period in the constant regime
  period_fixed <- observeEvent(input$period_fixed_value, {
    req(input$is_fix_period)
    draws <- rep(input$period_fixed_value, input$n_to_draw)
    rv$per <- TRUE
    period(draws)
  })
  
  # drawing magnitude on button click
  magnitude_draw <- eventReactive(input$draw_magnitude, {
      draws <- c()
      for (i in 1:input$n_to_draw) {
        draw <- half_t(df = input$magnitude_df,
                       scale = input$magnitude_scale)
        draws <- append(draws, draw)
      }
    
    rv$mag <- TRUE
    magnitude(draws)
  })
  
  # fixing magnitude in the constant regime
  magnitude_fixed <- observeEvent(input$magnitude_fixed_value, {
    req(input$is_fix_magnitude)
    draws <- rep(input$magnitude_fixed_value, input$n_to_draw)
    rv$mag <- TRUE
    magnitude(draws)
  })
  
#-------------------------------------------------------------------------------
# Second kernel parameters
  
  # drawing length scale on button click
  length_scale_2_draw <- eventReactive(input$draw_length_scale_2, {
    draws <- c()
    for (i in 1:input$n_to_draw) {
      draw <- inverse_gamma(
        alpha = input$length_scale_2_alpha,
        beta = input$length_scale_2_beta
      )
      draws <- append(draws, draw)
    }
    rv$ls_2 <- TRUE
    length_scale_2(draws)
  })
  
  # fixing length scale 2 in the constant regime
  length_scale_2_fixed <- observeEvent(input$length_scale_2_fixed_value, {
    req(input$is_fix_length_scale_2)
    draws <- rep(input$length_scale_2_fixed_value, input$n_to_draw)
    rv$ls_2 <- TRUE
    length_scale_2(draws)
  })
  
  
  # drawing period on button click
  period_2_draw <- eventReactive(input$draw_period_2, {
    draws <- c()
    for (i in 1:input$n_to_draw) {
      draw <- inverse_gamma(alpha = input$period_2_alpha,
                            beta = input$period_2_beta)
      draws <- append(draws, draw)
    }
    rv$per_2 <- TRUE
    period_2(draws)
  })
  
  # fixing period 2 in the constant regime
  period_2_fixed <- observeEvent(input$period_2_fixed_value, {
    req(input$is_fix_period_2)
    draws <- rep(input$period_2_fixed_value, input$n_to_draw)
    rv$per_2 <- TRUE
    period_2(draws)
  })
  
  # drawing magnitude on button click
  magnitude_2_draw <- eventReactive(input$draw_magnitude_2, {
    draws <- c()
    for (i in 1:input$n_to_draw) {
      draw <- half_t(df = input$magnitude_2_df,
                     scale = input$magnitude_2_scale)
      draws <- append(draws, draw)
    }
    rv$mag_2 <- TRUE
    magnitude_2(draws)
  })
  
  # fixing magnitude 2 in the constant regime 
  magnitude_2_fixed <- observeEvent(input$magnitude_2_fixed_value, {
    req(input$is_fix_magnitude_2)
    draws <- rep(input$magnitude_2_fixed_value, input$n_to_draw)
    rv$mag_2 <- TRUE
    magnitude_2(draws)
  })
  
#-------------------------------------------------------------------------------
## GP redraw logic
  
  gp_pool <- eventReactive(input$draw_gp, {
    if (!constrained_check()) {
      tryCatch({
        kernel <- c(input$kernel_label, input$kernel_label_2)
        len <- length_scale()
        mag <- magnitude()
        per <- period()
        ro <- as.numeric(input$nu)
        len_2 <- length_scale_2()
        mag_2 <- magnitude_2()
        per_2 <- period_2()
        ro_2 <- as.numeric(input$nu_2)
        loc <- location()
        ste <- steepness()
        ope <- operation()
        multi      <- isTRUE(TRUE)#input$multiple_draws_switch)
        old_params <- last_params()
        old_pool   <- last_pool()
        #TODO add bound parameters
        
        # reuse pool when parameters unchanged
        if (!is.null(old_params) &&
            is.list(old_params) &&
            identical(old_params$kernel_prev, kernel) &&
            identical(old_params$length_scale, len) &&
            identical(old_params$magnitude, mag) &&
            identical(old_params$period, per) &&
            identical(old_params$roughness, ro) &&
            identical(old_params$length_scale_2, len_2) &&
            identical(old_params$magnitude_2, mag_2) &&
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
        
        if (input$is_combination) {
          kernel_label <- c(input$kernel_label, input$kernel_label_2)
        } else{
          kernel_label <- input$kernel_label
        }
        
          n_draw <- input$n_to_draw
          
          funcs <- vapply(seq_len(n_draw), function(i) {
            kernel_params_i <- hash(
              "kernel_1" = hash(
                "magnitude" = mag[i]
                ,
                "length_scale" = len[i]
                ,
                "period" = per[i]
                ,
                "roughness" = ro
              ),
              "kernel_2" = hash(
                "magnitude" = mag_2[i]
                ,
                "length_scale" = len_2[i]
                ,
                "period" = per_2[i]
                ,
                "roughness" = ro_2
              ),
              "extra" = hash(
                "operation" = input$operation,
                "additional" = hash(
                  "location" = input$location,
                  "steepness" = input$steepness
                )
              )
            )
            
            simulate_gp(x_orig,
                        input$is_combination,
                        kernel_label,
                        kernel_params_i)
          }, numeric(length(x_orig)))
          
          new_pool <- list(
            mode   = "multi",
            x_orig = x_orig,
            funcs  = funcs         
          )

        last_params(
          list(
            kernel_prev = kernel,
            length_scale = len,
            magnitude = mag,
            period = per,
            roughness = ro ,
            length_scale_2 = len_2,
            magnitude_2 = mag_2,
            period_2 = per_2,
            roughness_2 = ro_2,
            location = loc,
            steepness = ste,
            operation = ope,
            multi        = multi
          )
        )
        last_pool(new_pool)
        
        new_pool
      }, error = function(e) {
        cat(paste("\nError in gp draw\n", e, "\n"))
      }, warning = function(w) {
        cat(paste("\nWarning in gp draw\n", w, "\n"))
      })
      
    } else{
      x <- c()
      y <- c()
      if (!invalid(input$x_1) && !invalid(input$y_1)) {
        x <- append(x, input$x_1)
        y <- append(y, input$y_1)
      }
      if (!invalid(input$x_2) && !invalid(input$y_2)) {
        x <- append(x, input$x_2)
        y <- append(y, input$y_2)
      }
      if (!invalid(input$x_3) && !invalid(input$y_3)) {
        x <- append(x, input$x_3)
        y <- append(y, input$y_3)
      }
      if (length(x) > 0) {
        x_train <- x
        y_train <- y
        x_draw <- seq(0, 1, 0.1)
        data_noise <- 1e-4
        x_fixed(x_train)
        y_fixed(y_train)
      } else{
        x_train <- seq(0, 1, 0.1)
        y_train <- x_train
        x_draw <- x_train
        data_noise <- 1.5
      }
      
      
      x_orig <-  seq(0, 1, 0.1) # TODO fix to normal
      y <- x_orig
      kernel_params <- hash(
        "kernel_1" = hash(
          "magnitude" = magnitude()
          ,
          "length_scale" = length_scale()
          ,
          "period" = period()
          ,
          "roughness" = as.numeric(input$nu)
        )
      )
      # TODO draw constrained
      constraints <- c()
      if (input$is_upper_bound || input$is_lower_bound) {
        constraints <- append(constraints, "boundedness")
      }
      if (input$is_monotonicity) {
        constraints <- append(constraints, "monotonicity")
      }
      if (input$is_convexity) {
        constraints <- append(constraints, "convexity")
      }
      
      constraint_params <- hash(
        "lower_bound" = input$lower_bound,
        "upper_bound" = input$upper_bound
        #, "monotonicity_type" = input$monotonicity # TODO is there a monotonically decreasing option?
      )
      
      funcs <- simulate_constrained_gp(
        x_train = x_train,
        y_train = y_train,
        kernel_label = input$kernel_label,
        kernel_params = kernel_params,
        constraints = constraints,
        constraint_params = constraint_params,
        n_functions = n_func,
        x_draw = x_draw,
        data_noise = data_noise
      )

      funcs
    }
  })
  
  # For storing and reproducing draws that were already computed
  gp_data <- reactive({
   
    req(gp_pool())
    pool  <- gp_pool()
    func_sequence <- seq_len(input$n_to_draw)
    
    if (!constrained_check()) {
      x_new <- seq(input$x_range[1], input$x_range[2], length.out = input$n_points)
      
        funcs_interp <- apply(gp_pool()$funcs, 2, function(f) {
          approx(gp_pool()$x_orig, f, xout = x_new)$y
        })
    }
    else{
      x_new <- seq(0, 1, 0.1)
      funcs_interp <- gp_pool()
    }
    
    data <- data.frame(
      x = rep(x_new, length(func_sequence)),
      f = as.vector(funcs_interp),
      func = rep(func_sequence, each = length(x_new))
    )
    
    data
  })
  
  observeEvent(input$draw_kernel, {
    enable("draw_gp")
  })
  
#-------------------------------------------------------------------------------
## Plots
  
  param_plot <- function(dataframe,
                         title,
                         y_label,
                         x_label,
                         observation_value) {
    plot <- ggplot(dataframe, aes(x, y)) +
      geom_line(color = "steelblue", linewidth = 1) +
      labs(title = title, y = y_label, x = x_label) +
      theme_minimal(base_size = 14)
    
    if (!invalid(observation_value)) {
      for (i in 1:length(observation_value)) {
        plot <- plot +
          annotate(
            "point",
            x = observation_value[i],
            y = 0,
            colour = colors[i],
            size = 2.5
          ) + annotate(
            "text",
            x = Inf,
            y = Inf,
            hjust = 1.1,
            vjust = 1.5 + i * 1.3,
            color = colors[i],
            size = 5,
            label = paste0("draw = ", signif(observation_value[i], 3))
          )
      }
    }
    plot
  }
  
  # Inverse-Gamma prior for length_scale plot
  output$plot_length_scale <- renderPlot({
    observe(length_scale_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$length_scale_alpha
    beta  <- input$length_scale_beta
    dens <- inverse_gamma(alpha = alpha,
                          beta = beta,
                          x = x_seq)
    d <- data.frame(x = x_seq, y = dens)
    
    param_plot(
      d,
      "Inverse-Gamma prior for length scale",
      "density",
      "length scale",
      length_scale()
    )
  })
  
  # Half-t prior for magnitude plot
  output$plot_magnitude <- renderPlot({
    observe(magnitude_draw())
    
    x_seq <- seq(0, 15, length.out = 400)
    df <- input$magnitude_df
    sc <- input$magnitude_scale
    dens <- half_t(df = df, scale = sc, x = x_seq)
    d <- data.frame(x = x_seq, y = dens)
    
    param_plot(d,
               "Half-t prior for magnitude",
               "density",
               "magnitude",
               magnitude())
  })
  
  # Inverse-Gamma prior for period plot
  output$plot_period <- renderPlot({
    observe(period_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$period_alpha
    beta  <- input$period_beta
    
    dens <- inverse_gamma(alpha = alpha,
                          beta = beta,
                          x = x_seq)
    d <- data.frame(x = x_seq, y = dens)
    
    param_plot(d,
               "Inverse-Gamma prior for period",
               "density",
               "period",
               period())
  })
  
#-------------------------------------------------------------------------------
## Plots for second kernel parameters
  
  # Inverse-Gamma prior for length_scale plot
  output$plot_length_scale_2 <- renderPlot({
    observe(length_scale_2_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$length_scale_2_alpha
    beta  <- input$length_scale_2_beta
    
    dens <- inverse_gamma(alpha = alpha,
                          beta = beta,
                          x = x_seq)
    d <- data.frame(x = x_seq, y = dens)
    
    param_plot(
      d,
      "Inverse-Gamma prior for length scale of the second kernel",
      "density",
      "length scale",
      length_scale_2()
    )
  })
  
  # Half-t prior for magnitude plot
  output$plot_magnitude_2 <- renderPlot({
    observe(magnitude_2_draw())
    
    x_seq <- seq(0, 15, length.out = 400)
    df <- input$magnitude_2_df
    sc <- input$magnitude_2_scale
    dens <- half_t(df = df, scale = sc, x = x_seq)
    d <- data.frame(x = x_seq, y = dens)
    
    param_plot(
      d,
      "Half-t prior for magnitude of the second kernel",
      "density",
      "magnitude",
      magnitude_2()
    )
  })
  
  # Inverse-Gamma prior for period plot
  output$plot_period_2 <- renderPlot({
    observe(period_2_draw())
    
    x_seq <- seq(1e-6, 15, length.out = 400)
    alpha <- input$period_2_alpha
    beta  <- input$period_2_beta
    dens <- inverse_gamma(alpha = alpha,
                          beta = beta,
                          x = x_seq)
    d <- data.frame(x = x_seq, y = dens)
    
    param_plot(
      d,
      "Inverse-Gamma prior for period of the second kernel",
      "density",
      "period",
      period_2()
    )
  })
  
  # Kernel based on distance plot
  output$plot_kernel <- renderPlot({
    req(input$draw_kernel)
    tryCatch({
      if (condition_parameters_check()) {
        dist <- seq(-3, 3, length.out = 300)
        x_o <- rep(0, length(dist))
        
        if (input$is_combination) {
          kernel_label <- c(input$kernel_label, input$kernel_label_2)
        } else{
          kernel_label <- input$kernel_label
        }
        kernel_title <- paste(input$kernel_label,
                              input$operation,
                              input$kernel_label_2)
        
        kernel_data <- data.frame(dist = dist)
        
        for (i in 1:input$n_to_draw) {
          kernel_params <- hash(
            "kernel_1" = hash(
              "magnitude" = magnitude()[i]
              ,
              "length_scale" = length_scale()[i]
              ,
              "period" = period()[i]
              ,
              "roughness" = as.numeric(input$nu)
            ),
            "kernel_2" = hash(
              "magnitude" = magnitude_2()[i]
              ,
              "length_scale" = length_scale_2()[i]
              ,
              "period" = period_2()[i]
              ,
              "roughness" = as.numeric(input$nu_2)
            ),
            "extra" = hash(
              "operation" = input$operation,
              "additional" = hash(
                "location" = input$location,
                "steepness" = input$steepness
              )
            )
          )
          
          k <- kernel_wrapper(input$is_combination,
                              kernel_label,
                              dist,
                              x_o,
                              kernel_params)[, 1]
          
          kernel_data[paste("k", sep = "", i)] <- k
        }
        
        plot <- ggplot()
        
        for (i in 1:input$n_to_draw) {
          col_name <- paste0("k", i)
          
          layer_data <- data.frame(x = kernel_data$dist, y = kernel_data[[col_name]])
          
          plot <- plot + geom_line(
            data = layer_data,
            aes(x = x, y = y),
            color = colors[i],
            linewidth = 1.2
          )
        }
        
        plot +  labs(
          title = paste(kernel_title, "Kernel"),
          x = "Distance",
          y = "k(x)"
        ) +
          theme_minimal(base_size = 14)
      }
    }, error = function(e) {
      cat(paste("\nError in kernel plot\n", e, "\n"))
    }, warning = function(w) {
      cat(paste("\nWarning in kernel plot\n", w, "\n"))
    })
  })
  
  # GP prior draws plot
  output$plot_gp <- renderPlot({
    tryCatch({
      req(input$draw_gp)
      if (1 ||
          condition_parameters_check() &&
          !constrained_check() 
          || condition_constrained_parameters_check()) {
        observe(gp_data())
        kernel_title <- paste(input$kernel_label,
                              input$operation,
                              input$kernel_label_2)
        
        lines <- geom_line(
          data = gp_data(),
          aes(
            x = x,
            y = f,
            group = func,
            color = factor(func)
          ),
          alpha     = 0.9,
          linewidth = 1
        )
        
        ribbon_data <- NULL
        CI <- NULL
        fixed_points <- NULL
        bounds <- NULL
        
        if (constrained_check()) {
          x <- as.vector(gp_data()$x)
          x_m <- matrix(x, length(x) / n_func, n_func)
          x_draw <- x_m[, 1]
          dat <- gp_data()$f
          funcs <- matrix(dat, length(dat) / n_func, n_func)
          qtls <- apply(funcs, 1, quantile, probs =  c(0.05, 0.95))
          # confidence interval
          ribbon_data <- data.frame(x    = x_draw,
                                    ymin = qtls[1, ],
                                    ymax = qtls[2, ])
          
          CI <- geom_ribbon(
            data = ribbon_data,
            aes(
              x = x,
              ymin = ymin,
              ymax = ymax
            ),
            fill  = "gray80",
            alpha = 0.6,
            inherit.aes = FALSE
          )
          
          if (!invalid(x_fixed()) && !invalid(y_fixed())) {
            fixed_points <- geom_point(
              data = data.frame(x = x_fixed(), y = y_fixed()),
              aes(x = x, y = y),
              color = "black",
              size = 2
            )
          }
          
          if (input$is_lower_bound || input$is_upper_bound) {
            bounds <- ylim(input$lower_bound, input$upper_bound)
          }
        }
        
        ggplot() +
          # CI first so lines are drawn on top
          CI +
          # The sampled function lines
          lines +
          # Training points
          fixed_points +
          bounds +
          labs(
            x     = "x",
            y     = "f(x)",
            title = "Gaussian Process Prior Samples",
            subtitle = paste(kernel_title, "Kernel"),
            color = "Function"
          ) +
          scale_color_manual(values = setNames(colors, levels(factor(
            gp_data()$func
          )))) +
          # Manual legend entries for ribbon + dashed line + points
          guides(color = "none") +          # drop per-function color legend if not needed
          annotate(
            "rect",
            # proxy for the ribbon in the legend
            xmin = -Inf,
            xmax = -Inf,
            ymin = -Inf,
            ymax = -Inf,
            fill = "gray80"
          ) +
          theme_minimal(base_size = 16)
      }
    }, error = function(e) {
      cat(paste("\nError in gp plot\n", e, "\n"))
    }, warning = function(w) {
      cat(paste("\nWarning in gp plot\n", w, "\n"))
    })
  })
}

shinyApp(ui, server)
