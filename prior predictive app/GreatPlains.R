library(gtools)

## Kernels

# SE kernel function
squared_exponential_kernel <- function(x1, x2, magnitude, length_scale) {
  outer(x1, x2, function(a, b)
    (magnitude^2 * exp(-(a - b)^2 / (2 * length_scale^2))))
}

# Linear kernel function
linear_kernel <- function(x1, x2, magnitude) {
  outer(x1, x2, function(a, b)
    (magnitude^2 * a * b))
}

# Matérn kernel function
matern_kernel <- function(x1, x2, magnitude, length_scale, roughness) {
  distance <- outer(x1, x2, function(a, b)
    abs(a - b))
  term <- sqrt(2 * roughness) * distance / length_scale
  K <- magnitude^2 * (2^(1 - roughness) / gamma(roughness)) * (term^roughness) * besselK(term, roughness)
  K[distance == 0] <- magnitude^2   # Handle the distance = 0 case explicitly
  K
}

# Periodic kernel function
periodic_kernel <- function(x1, x2, magnitude, length_scale, period) {
  outer(x1, x2, function(a, b)
    (magnitude^2 * exp(
      -2 * sin(pi * abs(a - b) / period)^2 / length_scale^2
    )))
}

kernel_labels <- c("Squared Exponential", "Matérn", "Linear", "Periodic")

kernel_operation_labels <- c("add", "multiply", "changepoint")

kernel_formulae <- hash(kernel_labels, c("\\sigma^2\\exp\\left(\\frac{(x - x')^2}{\\lambda^2}\\right)\\)",
                                         "\\sigma^2\\frac{2^{1 - \\nu}}{\\Gamma(\\nu)} \\left(\\frac{\\sqrt{2\\nu}|x - x'|}{\\lambda}\\right)^\\nu K_{\\nu}\\left(\\frac{\\sqrt{2\\nu}|x - x'|}{\\lambda}\\right)\\)
                                         where \\(\\Gamma\\) is the gamma function, \\(K_{\\nu}\\) is a modified Bessel function", 
                                         "\\sigma^2xx'\\)", 
                                         "\\sigma^2\\exp\\left(-\\frac{2\\sin^2{\\frac{\\pi|x - x'|}{p}}}{\\lambda^2}\\right) \\)"))

kernel_operation_formulae = hash(kernel_operation_labels, c("\\( k(x, x') = k_1(x, x') + k_2(x, x') \\)", 
"\\( k(x, x') = k_1(x, x')k_2(x, x') \\)",
"\\( k(x, x') = k_1(x, x')\\psi + k_2(x, x')\\hat\\psi\\)
\\(\\psi = \\psi(x)\\psi(x')\\)
\\(\\hat\\psi = (1 - \\psi(x))(1 - \\psi(x'))\\)
\\(\\psi(x) = \\frac{1}{1 + \\exp\\left(-s(x-x_{0})\\right)}\\)"))

psi <- function(x, steepness, location) {
  return(1 / (1 + exp(-steepness * (x - location))))
}

combine_kernels <- function(label_1,
                            label_2,
                            params_1,
                            params_2,
                            operation,
                            x1,
                            x2,
                            additional_params) {
  k_1 <- simple_kernel_wrapper(label_1, x1, x2, params = params_1)
  k_2 <- simple_kernel_wrapper(label_2, x1, x2, params = params_2)
  if (operation == "add") {
    return(k_1 + k_2)
  } else if (operation == "multiply") {
    return(k_1 * k_2)
  }
  else if (operation == "changepoint") {
    steepness <- additional_params[["steepness"]]
    location <- additional_params[["location"]]
    psi_total <- outer(psi(x1, steepness, location), psi(x2, steepness, location))
    psi_total_rev <- outer((1 - psi(x1, steepness, location)), (1 - psi(x2, steepness, location)))
    return(
      k_1 * psi_total+ k_2 *psi_total_rev 
    )
  }
}

kernel_wrapper <- function(is_combination,
                           kernel_label,
                           x1,
                           x2,
                           params) {
  
  if (is_combination) {
    extra_params <- params[["extra"]]
    operation <- extra_params[["operation"]]
    additional_params <- extra_params[["additional"]]
    return(
      combine_kernels(
        kernel_label[1],
        kernel_label[2],
        params[["kernel_1"]],
        params[["kernel_2"]],
        operation,
        x1,
        x2,
        additional_params
      )
    )
  } else{
    return(simple_kernel_wrapper(kernel_label, x1, x2, params[["kernel_1"]]))
  }
}

simple_kernel_wrapper <- function(kernel_label, x1, x2, params) {
  magnitude = params[["magnitude"]]
  length_scale = params[["length_scale"]]
  roughness = params[["roughness"]]
  period = params[["period"]]
  
  if (kernel_label == "Linear") {
    if (!invalid(x1) & !invalid(x2) & !invalid(magnitude)) {
      return(linear_kernel(
        x1 = x1,
        x2 = x2,
        magnitude = magnitude
      ))
    }
  }
  
  if (kernel_label == "Squared Exponential") {
    if (!invalid(x1) &
        !invalid(x2) &
        !invalid(magnitude) & !invalid(length_scale)) {
      return(
        squared_exponential_kernel(
          x1 = x1,
          x2 = x2,
          magnitude = magnitude,
          length_scale = length_scale
        )
      )
    }
  }
  
  if (kernel_label == "Matérn") {
    if (!invalid(x1) &
        !invalid(x2) &
        !invalid(magnitude) &
        !invalid(length_scale) & !invalid(roughness)) {
      return(
        matern_kernel(
          x1 = x1,
          x2 = x2,
          magnitude = magnitude,
          length_scale = length_scale,
          roughness = roughness
        )
      )
    }
  }
  
  if (kernel_label == "Periodic") {
    if (!invalid(x1) &
        !invalid(x2) &
        !invalid(magnitude) &
        !invalid(length_scale) & !invalid(period)) {
      return(
        periodic_kernel(
          x1 = x1,
          x2 = x2,
          magnitude = magnitude,
          length_scale = length_scale,
          period = period
        )
      )
    }
  }
  print("Some of the kernel parameters are invalid, correct them and try again.")
}

#-------------------------------------------------------------------------------
# Drawing GP
make_psd <- function(K, jitter = 1e-6) {
  K <- (K + t(K)) / 2          # force exact symmetry
  eig <- eigen(K, symmetric = TRUE)
  
  min_eig <- min(eig$values)
  if (min_eig < 0) {
    message(sprintf("Kernel has negative eigenvalue: %.4e — repairing", min_eig))
    # Clip negative eigenvalues to zero, add jitter
    eig$values <- pmax(eig$values, 0) + jitter
    K <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  } else {
    K <- K + jitter * diag(nrow(K))
  }
  K
}

simulate_gp <- function(x,
                        is_combination,
                        kernel_label,
                        kernel_params,
                        sigma_noise = 1e-3,
                        mean_fun = function(x)
                          0) {
  K <- kernel_wrapper(is_combination, kernel_label, x, x, params = kernel_params)
  K <- make_psd(K) 
  L <- chol(K)
  m <- mean_fun(x)
  f <- m + t(L) %*% rnorm(length(x))
  #cat(paste("\neeee\n"))
  #cat(paste(f))
  
  # noise
  eps <- sigma_noise * rnorm(length(x))
  
  drop(f + eps)
}

simulate_constrained_gp <- function(x_train,
                                    y_train,
                                    kernel_label,
                                    kernel_params,
                                    constraints,
                                    constraint_params,
                                    n_functions,
                                    x_draw,
                                    n_points = 50,
                                    data_noise = 1.5) {
  k_1 <- kernel_params[["kernel_1"]]$kernel_1
  magnitude <- k_1[["magnitude"]]
  length_scale <- k_1[["length_scale"]]
  roughness <- k_1[["roughness"]]
  
  if (kernel_label == "Squared Exponential") {
    kernel_type <- "gaussian"
  } else if (kernel_label == "Matérn") {
    if (roughness == 2.5) {
      kernel_type <- "matern52"
    } else if (roughness == 1.5) {
      kernel_type <- "matern32"
    }
  }
  
  model <- create(
    class = "lineqGP",
    x = x_train,
    y = y_train,
    constrType = constraints
  )
  model$localParam$m <- n_points
  if ("boundedness" %in% constraints) {
    model$bounds[1, ] <- c(constraint_params[["lower_bound"]], constraint_params[["upper_bound"]])
  }
  
  model$kernParam$type = kernel_type
  model$kernParam$par <- c(magnitude, length_scale)
  model$varnoise <- data_noise
  model <- augment(model)
  
  # sampling from the model
  sim.model <- simulate(model,
                        nsim = n_functions,
                        seed = 1,
                        xtest = x_draw)
  
  return(sim.model$ysim)
}
#-------------------------------------------------------------------------------
## Distributions

# inverse gamma
inverse_gamma <- function(alpha, beta, x = NULL) {
  if (invalid(x)) {
    1 / rgamma(1, shape = alpha, rate = beta)
  } else{
    dgamma(1 / x, shape = alpha, rate = beta) * 1 / x^2
  }
}

# half-t
half_t <- function(df, scale, x = NULL) {
  if (invalid(x)) {
    scale * abs(rt(1, df = df))
  } else{
    dens <- 2 * dt((x) / scale, df = df) / scale
    dens[x < 0] <- 0
    dens
  }
}