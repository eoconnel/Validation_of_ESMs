


gev_density <- function(params, x) {
  mu <- params[1]
  sigma <- params[2]
  gamma <- params[3]
  
  if (gamma == 0) {
    f <- (1 / sigma) * exp(-((x - mu) / sigma) - exp(-((x - mu) / sigma)))
  } else {
    f <- (1 / sigma) * (1 + gamma * ((x - mu) / sigma))^(-(1 / gamma + 1)) *
      exp(- (1 + gamma * ((x - mu) / sigma))^(-(1 / gamma)))
  }
  
  f[(1 + gamma * (x - mu) / sigma) <= 0] <- 0
  f
}


gp_density <- function(params, x) {
  sigma <- params[1]
  gamma <- params[2]
  
  if (gamma == 0) {
    f <- (1 / sigma)*exp(-x / sigma)
  } else {
    f <- (1 / sigma) * (1 + gamma * x / sigma) ^ (-(1/gamma + 1))
  }
  
  f[x <= 0] <- 0
  f[(1 + gamma * x / sigma) <= 0] <- 0
  f
}


gev_log_like <- function(params, data, transform = TRUE) {
  n <- length(data)
  mu <- params[1]
  sigma <- params[2]
  gamma <- params[3]
  
  if (transform) {
    sigma <- exp(sigma)
  }
  
  if (any((1 + gamma * (data - mu) / sigma) <= 0)) {
    ll <- - Inf
  } else if (gamma == 0) {
    ll <- -n*log(sigma) - (1/sigma)*sum(data) + n*mu/sigma - sum(exp(-((x-mu)/sigma)))
  } else {
    ll <- -n*log(sigma) - (1/gamma + 1)*sum(log(1 + gamma*(data-mu)/sigma)) - 
      sum((1 + gamma*(data-mu)/sigma)^(-(1/gamma)))
  }
  ll
}


gp_log_like <- function(params, data, transform = TRUE) {
  n <- length(data)
  sigma <- params[1]
  gamma <- params[2]
  
  if (transform) {
    sigma <- exp(sigma)
  }
  
  if (any((1 + gamma * data / sigma) <= 0)) {
    ll <- - Inf
  } else if (gamma == 0) {
    ll <- - n*log(sigma) - (1/sigma)*sum(data)
  } else {
    ll <- -n*log(sigma) - (1 + 1/gamma)*sum(log(1 + gamma*data/sigma))
  }
  ll
}


mle_transform <- function(x, init, log_like, transform, jacobian) {
  
  # Maximizing log-likelihood function
  mle <- optim(init, log_like, data = x, control = list(fnscale = -1), hessian = TRUE)
  
  # Transforming parameters
  mle$par <- do.call(transform, list(params = mle$par))
  
  # Extracting covariance matrix and defining jacobian
  cov_mat <- solve(-mle$hessian)
  J_inv <- solve(do.call(jacobian, list(params = mle$par)))
  
  # Calculating transformed covariance matrix
  mle$cov <- J_inv %*% cov_mat %*% t(J_inv)
  
  # Calculating parameter standard errors
  mle$se <- sqrt(diag(mle$cov))
  
  # Returning maximum likelihood estimation
  mle
}


gev_mle <- function(x, init) {
  
  transform_func <- function(params) {
    params[2] <- exp(params[2])
    params
  }
  
  jacobian_func <- function(params) {
    matrix(c(1, 0, 0,
             0, 1/params[2], 0,
             0, 0, 1), nrow = 3, byrow = TRUE)
  }
  
  mle_transform(
    x = x,
    init = init,
    log_like = gev_log_like,
    transform = transform_func,
    jacobian = jacobian_func
  )
}


gp_mle <- function(x, init) {
  
  transform_func <- function(params) {
    params[1] <- exp(params[1])
    params
  }
  
  jacobian_func <- function(params) {
    matrix(c(1/params[1], 0,
             0, 1), nrow = 2, byrow = TRUE)
  }
  
  mle_transform(
    x = x,
    init = init,
    log_like = gp_log_like,
    transform = transform_func,
    jacobian = jacobian_func
  )
}


plot_gp_hist <- function(x, mle, color = "blue", line_width = 0.5, bins = 20) {
  
  x_range <- seq(from = 0, to = max(x), length.out = 100)
  f <- gp_density(mle$par, x_range)
  
  hist_data <- tibble(x)
  dist_data <- tibble("x" = x_range, "f" = f)
  
  hist_data |>
    ggplot(aes(x)) +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = seq(0, max(x), length.out = bins),
                   fill = "lightgrey", alpha = 0.5, color = "black") +  
    geom_line(data = dist_data, aes(x = x, y = f), color = color, lwd = line_width) +
    labs(
      title = "Generalized Pareto Modeled Distribution Function",
      y = "Density"
    ) +
    theme_bw()
}


plot_gp_empirical <- function(x, mle) {
  
  x_range <- seq(from = 0, to = max(x), length.out = 100)
  f <- gp_density(mle$par, x_range)
  
  empirical <- logdensity(x, from = 0, to = max(x))
  empirical$y[is.na(empirical$y)] = 0
  
  tibble(
    "x" = c(x_range, empirical$x),
    "y" = c(f, empirical$y),
    "label" = c(rep("Modeled", length(x_range)), rep("Empirical", length(empirical$x)))
  ) |>
    ggplot(aes(x = x, y = y, col = label, lty = label)) +
    geom_line() +
    labs(
      title = "Generalized Pareto Distribution Functions",
      x = "",
      y = "Denisty",
      color = "Distribution", 
      lty = "Distribution") +
    scale_color_manual(
      values = c("black", "blue"),
      breaks = c("Empirical", "Modeled"),
      labels = c("Empirical", "Modeled")
    ) +
    theme_bw()
}


fit_thresholds <- function(x, range, nout = 10, alpha = NULL, init = c(0,0)) {
  
  # Defining thresholds
  u <- seq(from = range[1], to = range[2], length.out = nout)
  
  # Initializing parameters
  params <- matrix(nrow = nout, ncol = 2)
  se <- matrix(nrow = nout, ncol = 2)
  param_labels <- c("Transformed scale", "Shape")
  error_width <- ((range[2] - range[1]) / nout) / 10
  crit <- 0
  if (!is.null(alpha)) {
    crit <- qnorm(1 - alpha / 2)
  }
  
  # Looping through each threshold
  for (i in 1:nout) {
    
    # Getting exceedances
    Z <- x[x > u[i]] - u[i]
    
    # Fitting maximum likelihood
    mle <- gp_mle(Z, init)
    
    # Printing failed convergence code
    if (mle$convergence != 0) {
      print(paste("Convergence code (threshold = ", u[i], "): ", mle$convergence, sep = ""))
    }
    
    # Extracting parameters
    params[i,] = mle$par
    se[i,] = mle$se
    
    # Transforming scale
    new_var <- matrix(c(1,-u[i]), nrow = 1) %*% mle$cov %*% t(matrix(c(1,-u[i]), nrow = 1))
    params[i,1] <- params[i,1] - params[i,2]*u[i]
    se[i,1] <- sqrt(new_var)
    
  }
  
  # Plotting thresholds
  tibble(
    "threshold" = rep(u, 2),
    "mle" = as.vector((params)),
    "se" = as.vector((se)),
    "param" = factor(rep(param_labels, each = nout), levels = param_labels)
  ) |>
    ggplot(aes(x = threshold, y = mle)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mle - crit*se, ymax = mle + crit*se), width = error_width,
                  position = position_dodge(0.05)) +
    facet_wrap(~param, ncol = 1, scales = "free") +
    labs(
      title = "Threshold Selection Plot",
      y = "",
      x = "Threshold"
    )
}
