library(tidyverse)

# Function to extract data from a single site and calculate
# seasonal cycle + anomalies
# @param data: time series data with OBS and FVCOM columns
# @param temp_site: temperature site to subset from the data
# @return time series tibble with seasonal cycle + anomalies
calculate_anomalies <- function(data, temp_site) {
  
  fvcom <- data |>
    filter(site == temp_site) |>
    select(time, FVCOM) |>
    rename("temp" = FVCOM) |>
    mutate(seas = smooth.spline(temp, nknots = 15)$y) |>
    mutate(anomaly = (temp - seas)) |>
    mutate(type = "FVCOM")
  
  obs <- data |>
    filter(site == temp_site) |>
    select(time, OBS) |>
    rename("temp" = OBS) |>
    mutate(seas = smooth.spline(temp, nknots = 15)$y) |>
    mutate(anomaly = (temp - seas)) |>
    mutate(type = "OBS")
  
  bind_rows(obs, fvcom)
}


# Function to filter time series in the frequency domain
# @param x: time series to filter
# @param dt: sampling interval
# @param fc: cutoff frequency (scalar for low/high pass, vector for band pass)
# @param type: filter type ("lowpass", "highpass", or "bandpass")
# @return filtered time series
freq_filter <- function(x, fc, dt = 1, type = "lowpass") {
  
  # Padding to ensure even time series length
  odd_length <- (length(x) %% 2) == 1
  if (odd_length) {
    x[length(x) + 1] <- x[length(x)]
  }
  
  # Compute Fourier Transform of the input
  x_fft <- fft(x)
  nt <- length(x)
  nf <- length(x_fft) / 2
  
  # Normalizing fft and constructing frequency vector
  x_fft <- x_fft[1:(nf+1)] / (nt / 2)
  freq <- (1 / nt) * seq(0, nf)
  freq <- (1 / dt) * freq
  
  # Filtering FFT
  idx <- switch(
    type,
    "lowpass" = which(freq <= fc[1]),
    "highpass" = which(freq > fc[1]),
    "bandpass" = which(freq >= fc[1] & freq <= fc[2])
  )
  filtered_fft <- rep(0, nf + 1)
  filtered_fft[idx] <- x_fft[idx]
  
  # Reconstructing signal from filtered FFT
  filtered_fft_rev <- c(filtered_fft, rev(Conj(filtered_fft[2:nf])))
  x_filtered <- Re(0.5*fft(filtered_fft_rev, inverse=TRUE))
  
  # Removing padding for odd time series
  if (odd_length) {
    x_filtered <- x_filtered[-length(x_filtered)]
  }
  
  x_filtered
}


# Function to conduct matrix-based ordinary least squares regression
# @param x: time series vector of response variable
# @param Z: regression design matrix
# @return betahat: estimated regression parameters
# @return xhat: estimated response variables
# @return e: error (x - xhat)
# @return R2: coefficient of determination
matrixOLS <- function(x, Z) {
  
  # Compute Z'Z and its inverse
  ZtZ <- t(Z) %*% Z
  ZtZi <- solve(ZtZ)
  
  # Compute Z'x
  Zty <- t(Z) %*% x
  
  # Estimating coefficients and residuals
  betahat <- as.vector(ZtZi %*% Zty)
  xhat <- as.vector(Z %*% betahat)
  e <- x - xhat
  
  # Calculating R^2
  SSE=sum(e^2)
  SST=sum(x^2)
  R2 = 1 - (SSE/SST)
  
  # Returning results
  results <- list(betahat = betahat, xhat = xhat, e = e, R2 = R2)
}


# Function to conduct a harmonic regression analysis for a given set
# of periods. Note that the CCF is calculated as model x observations.
# @param time: common time vector
# @param x_model: model time series data
# @param x_obs: observation time series data
# @param dt: time difference between observations
# @param periods: scalar or vector of periods (in the units of dt)
#                 on which to regress
# @param lags: maximum lag to use for the cross-correlation function
harmonic_regression <- function(time, x_model, x_obs, dt, periods, lags = 250) {
  
  # Conduct harmonic regression
  t <- (0:(length(time) - 1)) * dt
  omegas <- 2 * pi / periods
  
  Z <- NULL
  for (omega in omegas) {
    Z <- cbind(Z, cos(omega * t), sin(omega * t))
  } 
  
  model_regression <- matrixOLS(x_model, Z)
  obs_regression <- matrixOLS(x_obs, Z)
  
  # Calculating metrics -----------------------------------------------------
  model_amplitude <- norm(model_regression$betahat, type = "2")
  obs_amplitude <- norm(obs_regression$betahat, type = "2")
  amplitudes <- list("Model" = model_amplitude, "Observation" = obs_amplitude)
  
  reg_ccf <- ccf(model_regression$xhat, obs_regression$xhat,
                 lag.max = lags, type = "correlation", plot = FALSE)
  
  # Creating dataframes for plotting
  regression <- tibble(
    "time" = rep(time, 2),
    "temp" = c(x_model, x_obs),
    "reg" = c(model_regression$xhat, obs_regression$xhat),
    "type" = rep(c("Model", "Observation"), each = length(x_model))
  )
  reg_ccf_df <- tibble(
    "lag" = reg_ccf$lag[, 1, 1],
    "ccf" = reg_ccf$acf[, 1, 1]
  )
  ccf_peak <- reg_ccf_df$lag[which.max(reg_ccf_df$ccf)]
  
  # Constructing plots
  reg_plot <- regression |>
    ggplot(aes(x = time, y = reg, col = type)) +
    geom_line() +
    labs(
      title = "Harmonic Regression Estimates",
      x = "Time",
      y = "Temp (C)"
    )
  
  ccf_plot <- reg_ccf_df |>
    ggplot(aes(x = lag, y = ccf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    theme_light() +
    labs(
      title = "Temperature Regression CCF",
      x = "Lag",
      y = "CCF"
    )
  
  # Returning results
  list(
    "amplitudes" = amplitudes,
    "ccf_peak" = ccf_peak,
    "regression_plot" = reg_plot,
    "ccf_plot" = ccf_plot
  )
}


# Function to conduct a frequency band filtered analysis for a given type
# of frequency filter. Note that the CCF is calculated as model x
# observations.
# @param time: common time vector
# @param x_model: model time series data
# @param x_obs: observation time series data
# @param dt: time step between observations
# @param fc: filter cutoff frequency (scalar for low- or high-pass filter,
#            vector for bandpass filter)
# @param type: type of filter ("lowpass", "highpass", or "bandpass")
# @param lags: maximum lag to use for the cross-correlation function
# @param spec_spans: degree of smoothing for the spectral density calculation
band_analysis <- function(time, x_model, x_obs, dt, fc, type,
                          lags = 250, spec_spans = NULL) {
  
  # Ensuring even number of observations
  n <- length(time)
  if (n %% 2 == 1) {
    time <- time[1:(n-1)]
    x_model <- x_model[1:(n-1)]
    x_obs <- x_obs[1:(n-1)]
  }
  
  # Filtering time series
  model <- tibble(
    "time" = time,
    "temp" = x_model,
    "filtered" = freq_filter(x_model, dt, fc, type = type),
    "type" = "Model"
  )
  obs <- tibble(
    "time" = time,
    "temp" = x_obs,
    "filtered" = freq_filter(x_obs, dt, fc, type = type),
    "type" = "Observations"
  )
  filtered <- bind_rows(model, obs)
  
  # Calculating CCF and spectrum
  band_ccf <- ccf(model$filtered, obs$filtered, lag.max = lags,
                  type = "correlation", plot = FALSE)
  band_spec <- spec.pgram(cbind(model$filtered, obs$filtered), 
                          spans = spec_spans, plot = FALSE)
  
  # Filtering spectrum
  idx <- switch(
    type,
    "lowpass" = which(band_spec$freq <= fc[1]),
    "highpass" = which(band_spec$freq > fc[1]),
    "bandpass" = which(band_spec$freq >= fc[1] & band_spec$freq <= fc[2])
  )
  
  # Fitting spectral slopes
  model_lm <- matrixOLS(
    x = log(band_spec$spec[idx, 1]),
    Z = cbind(1, log(band_spec$freq[idx]))
  )
  obs_lm <- matrixOLS(
    x = log(band_spec$spec[idx, 2]),
    Z = cbind(1, log(band_spec$freq[idx]))
  )
  
  # Constructing dataframes for plotting 
  band_ccf_df <- tibble(
    lag = band_ccf$lag[, 1, 1],
    ccf = band_ccf$acf[, 1, 1]
  )
  band_spec_df <- tibble(
    freq = rep(band_spec$freq[idx], 2),
    power = c(band_spec$spec[idx, 1], band_spec$spec[idx, 2]),
    lm = c(model_lm$xhat, obs_lm$xhat),
    ID = rep(c("Model", "Observations"), each = length(idx))
  )
  
  # Calculating metrics
  rmse <- sqrt(mean((model$filtered - obs$filtered)^2))
  sd_ratio <- sd(obs$filtered) / sd(model$filtered)
  ccf_peak <- band_ccf_df$lag[which.max(band_ccf_df$ccf)]
  slopes <- list("Model" = model_lm$betahat[2], 
                 "Observation" = obs_lm$betahat[2])
  
  # Constructing plot of filtered time series
  filtered_plot <- ggplot(filtered, aes(x = time, y = filtered, 
                                        col = type)) +
    geom_line() +
    labs(
      title = "Filtered Temperature Anomalies",
      x = "Time",
      y = "Temp (C)",
      col = "Data"
    )
  
  # Constructing CCF plot
  ccf_plot <- ggplot(band_ccf_df, aes(x = lag, y = ccf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    theme_light() +
    labs(
      title = "Filtered Temperature Band CCF",
      x = "Lag",
      y = "CCF"
    )
  
  # Constructing power spectrum plot
  spectrum_plot <- ggplot(band_spec_df, aes(x = freq, y = log(power), 
                                            col = ID, group = ID)) +
    geom_line() +
    labs(
      title = "Auto-Spectrum of Filtered Temperatures",
      x = "Frequency",
      y = "Log power",
      col = "Data"
    ) +
    theme_light()
  
  # Constructing spectral slope plot
  slope_plot <- ggplot(band_spec_df, aes(x = log(freq), y = log(power), 
                                         col = ID, group = ID)) +
    geom_line() +
    geom_line(aes(y = lm)) +
    labs(
      title = "Auto-Spectrum of Filtered Temperatures",
      x = "Log Frequency",
      y = "Log power",
      col = "Data"
    ) +
    theme_light()
  
  # Returning results
  list(
    "rmse" = rmse,
    "sd_ratio" = sd_ratio,
    "ccf_peak" = ccf_peak,
    "spectral_slopes" = slopes,
    "filtered_plot" = filtered_plot,
    "ccf_plot" = ccf_plot,
    "spectrum_plot" = spectrum_plot,
    "slope_plot" = slope_plot
  )
}


# Function to create a spectrogram for a given time series. Function returns
# the spectrogram plot
# @param time: time vector for time series
# @param x: data vector for time series
# @param dt: time step between observations
# @param width: number of observations in each segment
# @param overlap: overlap between consecutive segments
# @param spec_spans: spans argument used for the spec.pgram function
# @param spec_taper: taper argument used for the spec.pgram function
# @param noisefloor: smallest possible value for the periodograms
# @param fc: cutoff frequency for plotting the periodograms
# @param log_power: set to true to return the log power of the spectrogram
# @param normalize: set to true to normalize each segment to a standard 
#                   deviation of 1
evolutionary_spectrum <- function(time, x, dt, width, overlap, spec_spans,
                                  spec_taper, noisefloor, fc, 
                                  log_power = TRUE, normalize = FALSE) {
  
  n <- length(time)
  nseg <- floor(n / (width - overlap)) - 1
  
  # Extracting first segment
  initial_segment <- x[1:width] - smooth.spline(x[1:width], nknots = 9)$y
  
  if (normalize) {
    initial_segment <- initial_segment / sd(initial_segment)
  }
  
  # Calculating first spectrum
  initial_spec <- spec.pgram(initial_segment, plot = FALSE,
                             spans = spec_spans, taper = spec_taper)
  
  # Initializing result matrices
  spec_matrix <- matrix(nrow = length(initial_spec$spec), ncol = nseg)
  spec_time <- array(NA, dim = nseg)
  freq <- (1 / dt) * initial_spec$freq
  
  # Saving first spectrum
  spec_matrix[, 1] <- initial_spec$spec
  spec_time[1] <- median(time[1:width])
  
  # Looping over remaining segments
  for (i in 2:nseg) {
    start_ind <- (i - 1) * (width - overlap) + 1
    end_ind <- (i - 1) * (width - overlap) + width
    
    # Extracting and detrending segment
    segment <- x[start_ind:end_ind]
    segment_trend <- smooth.spline(segment, nknots = 9)$y
    segment_anomaly <- segment - segment_trend
    
    if (normalize) {
      segment_anomaly <- segment_anomaly / sd(segment_anomaly)
    }
    
    # Calculating spectrum
    segment_spec <- spec.pgram(segment_anomaly, plot = FALSE,
                               spans = spec_spans, taper = spec_taper)
    spec_matrix[, i] <- segment_spec$spec
    spec_time[i] <- median(time[start_ind:end_ind])
  }
  
  # Setting noisefloor and cutting off frequencies
  spec_matrix[spec_matrix < noisefloor] <- noisefloor
  spec_matrix <- spec_matrix[freq < 0.25,]
  freq <- freq[freq < 0.25]
  
  legend_name <- "Power"
  if (log_power) {
    spec_matrix <- log(spec_matrix)
    legend_name <- "Log power"
  }
  
  spec_matrix |>
    as_tibble(.name_repair = NULL) |>
    mutate("freq" = freq) |>
    pivot_longer(
      cols = starts_with("V"),
      names_to = "time",
      values_to = "power"
    ) |>
    mutate(
      "time" = as.numeric(gsub("V", "", time))
    ) |>
    mutate(
      "time" = as.POSIXct(spec_time[time], origin = "1970-01-01")
    ) |>
    ggplot(aes(x = time, y = freq, fill = power)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("darkblue", "blue", "lightblue1",
                                    "green", "yellow", "red", "darkred"),
                         name = legend_name)
}


evolutionary_spectrum2 <- function(time, x, dt = 1, width = 24*30, 
                                   overlap = 24*15, spec_spans = 5,
                                   spec_taper = 0.1, noisefloor = 1e-3, 
                                   normalize = FALSE) {
  
  n <- length(time)
  nseg <- floor(n / (width - overlap)) - 1
  
  # Extracting first segment
  initial_segment <- x[1:width]
  
  if (normalize) {
    initial_segment <- initial_segment / sd(initial_segment)
  }
  
  # Calculating first spectrum
  initial_spec <- spec.pgram(initial_segment, plot = FALSE,
                             spans = spec_spans, taper = spec_taper)
  spec_width <- length(initial_spec$spec)
  
  # Initializing result matrices
  spec_vec <- rep(NA, nseg*spec_width)
  spec_time <- rep(NA, nseg)
  freq <- (1 / dt) * initial_spec$freq
  
  # Saving first spectrum
  spec_vec[1:spec_width] <- initial_spec$spec
  spec_time[1] <- median(time[1:width])
  
  # Looping over remaining segments
  for (i in 2:nseg) {
    start_ind <- (i - 1) * (width - overlap) + 1
    end_ind <- (i - 1) * (width - overlap) + width
    
    # Extracting and detrending segment
    segment <- x[start_ind:end_ind]
    
    if (normalize) {
      segment <- segment / sd(segment)
    }
    
    # Calculating spectrum
    segment_spec <- spec.pgram(segment, plot = FALSE,
                               spans = spec_spans, taper = spec_taper)
    spec_vec[(1+360*(i-1)):(360*i)] <- segment_spec$spec
    spec_time[i] <- median(time[start_ind:end_ind])
  }
  
  # Converting times
  spec_time <- as.POSIXct(spec_time, origin = "1970-01-01", tz = "UTC")
  
  # Setting noisefloor
  spec_vec[spec_vec < noisefloor] <- noisefloor
  
  # Creating dataframe
  tibble(
    "t" = rep(spec_time, each = spec_width),
    "freq" = rep(freq, nseg),
    "power" = spec_vec
  )
}


# Function to calculate a metric for a given time series of model and 
# observations, by segmenting the time series.
# @param time: time vector for time series
# @param x: first data vector for time series
# @param y: second data vector for time series (can be NULL if metric is 
#           univariate)
# @param width: number of observations in each segment
# @param overlap: overlap between consecutive segments
# @param metric: metric to calculate between x and y (or just on x vector
#                if y is NULL)
# @param metric_params: any additional paramaters for calculating metrics
# @param detrend: set to TRUE to remove the trend in each segment (via 
#                 smoothing spline function)
evolutionary_metrics <- function(time, obs, pred, metric, metric_params = NULL,
                                 width = 24*30, overlap = 24*15) {
  
  # Initializing variables
  n <- length(time)
  nseg <- floor(n / (width - overlap)) - 1

  metrics_vec <- rep(NA, nseg)
  metrics_time <- rep(NA, nseg)
  
  # Looping through each segment
  for (i in 1:nseg) {
    
    # Extracting segment
    start_ind <- (i - 1) * (width - overlap) + 1
    end_ind <- (i - 1) * (width - overlap) + width
    obs_seg <- obs[start_ind:end_ind]
    pred_seg <- pred[start_ind:end_ind]
    
    # Constructing segment list
    segments <- list("obs" = obs_seg)
    if (!is.null(pred)) {
      segments <- c(segments, list("pred" = pred_seg))
    }
    
    # Calculating metrics
    metrics_vec[i] <- do.call(metric, c(segments, metric_params))
    metrics_time[i] <- median(time[start_ind:end_ind])
  }
  
  # Contructing tibble to return
  tibble(
    "t" = as.POSIXct(metrics_time, origin = "1970-01-01", tz = "UTC"),
    "metric" = metrics_vec
  )
}

aae <- function(obs, pred) {
  mean(abs(pred - obs))
}

mse <- function(obs, pred) {
  mean((pred - obs)^2)
}

rmse <- function(obs, pred) {
  sqrt(mean((pred - obs)^2))
}

bias <- function(obs, pred) {
  mean(pred - obs)
}

r <- function(obs, pred) {
  obs_bar <- mean(obs)
  pred_bar <- mean(pred)
  
  sum((obs - obs_bar) * (pred - pred_bar)) / 
    sqrt(sum((obs - obs_bar)^2) * sum((pred - pred_bar)^2))
}

willmott <- function(obs, pred) {
  1 - rmse(obs, pred)^2 / mean((abs(pred - mean(obs)) + abs(obs - mean(obs)))^2)
}

nse <- function(obs, pred) {
  1 - rmse(obs, pred)^2 / var(obs)
}

sig_star <- function(obs, pred) {
  sd(pred) / sd(obs)
}

yhour <- function(time) {
  24 * (yday(time) - 1) + hour(time) + 1
}

remove_missing_times <- function(data) {
  
  sites <- unique(data$site)
  filtered <- NULL
  
  for (i in 1:length(sites)) {
    
    site_data <- data |>
      filter(site == sites[i])
    
    site_obs <- site_data |>
      filter(type == "OBS")
    
    missing_times <- site_obs |>
      filter(is.na(temp)) |>
      pull(t)
    
    site_fvcom <- site_data |>
      filter(type == "FVCOM")
    site_fvcom2 <- site_data |>
      filter(type == "FVCOM2")
    site_glorys <- site_data |>
      filter(type == "GLORYS")
    site_cmip6 <- site_data |>
      filter(type == "CMIP6")
    
    site_fvcom$temp[site_fvcom$t %in% missing_times] <- NA
    site_fvcom2$temp[site_fvcom2$t %in% missing_times] <- NA
    site_glorys$temp[site_glorys$t %in% missing_times] <- NA
    site_cmip6$temp[site_cmip6$t %in% missing_times] <- NA
    
    filtered_site <- bind_rows(
      site_obs,
      site_fvcom,
      site_fvcom2,
      site_glorys,
      site_cmip6
    )
    filtered <- bind_rows(filtered, filtered_site)
  }
  
  filtered
}

exceedance_area <- function(clm) {
  
  intersections <- clm |>
    mutate("t" = as.numeric(as_datetime(t))) |>
    get_intersections("t", "temp", "thresh") |>
    mutate("t" = as_datetime(t))
  
  clm |>
    mutate(t = as_datetime(t)) |>
    bind_rows(intersections) |>
    arrange(t) |>
    mutate(
      "event_no" = pmax(event_no, lag(event_no), lead(event_no), na.rm = TRUE)
    ) |>
    drop_na(event_no) |>
    group_by(event_no) |>
    summarize("area" = trapezoid_method(t, temp - thresh)) |>
    pull(area) |>
    sum()
}

trapezoid_method <- function(t, data, unit = "days") {
  
  h <- as.numeric(difftime(t[-1], t[-length(t)], units = unit))
  b1 <- data[-1]
  b2 <- data[-length(data)]
  sum(pmax(h * (b1 + b2) / 2, 0))
}

apply_metrics <- function(obs, pred, pred_label, metric_funcs, metric_labels) {
  
  metrics <- rep(NA, length(metric_funcs))
  for (i in 1:length(metric_funcs)) {
    metrics[i] <- metric_funcs[[i]](obs, pred)
  }
  
  tibble(
    "metric" = metric_labels,
    !!pred_label := metrics
  )
}

get_slope <- function(x, y) {
  data <- data.frame(x, y)
  fit <- lm(y ~ x, data)
  fit$coefficients[2]
}
