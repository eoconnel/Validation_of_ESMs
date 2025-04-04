library(tidyverse)

event_funcs <- list(
  "mean" = ~ mean(.x, na.rm = TRUE),
  "max" = ~ max(.x, na.rm = TRUE),
  "cumulative" = ~ sum(.x, na.rm = TRUE)
)

mhw_climatology <- function(ts_data, 
                            hourly = FALSE, 
                            start_time = NULL,
                            end_time = NULL, 
                            filter_function = NULL,
                            filter_params = NULL,
                            clm_filter_width = 31, 
                            window_halfwidth = 5,
                            thresh_quantile = 0.9) {
  
  # Setting parameters for daily or hourly data
  if (hourly) {
    leap_start <- 1417
    leap_end <- 1440
    ind_unit <- 24
  } else {
    leap_start <- 60
    leap_end <- 60
    ind_unit <- 1
  }
  ind_max <- 365*ind_unit
  window_halfwidth <- window_halfwidth*ind_unit
  clm_filter_width <- clm_filter_width*ind_unit
  
  # Setting default filter function
  if (is.null(filter_function)) {
    filter_function <- stats::filter
    filter_params <- list(
      "filter" = rep(1/clm_filter_width, clm_filter_width),
      "method" = "convolution",
      "sides" = 2,
      "circular" = TRUE
    )
  }
  
  # Getting time frame for the climatology
  if (is.null(start_time)) {
    start_time <- min(ts_data$t)
  }
  if (is.null(end_time)) {
    end_time <- max(ts_data$t)
  }
  
  # Calculating daily or hourly index
  clm <- ts_data |>
    mutate(
      "included" = between(t, start_time, end_time),
      "hourly" = hourly,
      "year" = year(t),
      "ind" = if_else(hourly, yhour(t), yday(t))
    )
  
  # Extracting leap days
  leap_days <- clm |>
    filter((month(t) == 2) & day(t) == 29)
  
  # Filtering out leap days
  clm <- clm |>
    filter(!((month(t) == 2) & day(t) == 29)) |>
    mutate(
      "ind" = if_else((year %in% leap_days$year) & (ind > leap_end), 
                      ind - ind_unit, ind)
    )
  
  # Initializing seasonal and threshold vectors
  seas <- rep(NA, ind_max)
  thresh <- rep(NA, ind_max)
  
  # Looping over each day of year
  for (i in 1:ind_max) {
    
    # Determining 11 day window
    if (i <= window_halfwidth) {
      inds <- ((i - window_halfwidth - 1):(i + window_halfwidth)) %% (ind_max + 1)
      inds <- inds[-which(inds == 0)]
    } else if (i > (ind_max - window_halfwidth)) {
      inds <- ((i - window_halfwidth):(i + window_halfwidth + 1)) %% (ind_max + 1)
      inds <- inds[-which(inds == 0)]
    } else {
      inds <- ((i - window_halfwidth):(i + window_halfwidth)) %% (ind_max + 1)
    }
    
    # Getting temperatures
    temps <- clm$temp[clm$included & (clm$ind %in% inds)]
    
    # Calculating climatology and threshold
    seas[i] <- mean(temps, na.rm = TRUE)
    thresh[i] <- quantile(temps, thresh_quantile, na.rm = TRUE)
  }
  
  # Smoothing seasonal and threshold vectors
  smooth_seas <- as.vector(do.call(filter_function, c(list("x" = seas), filter_params)))
  smooth_thresh <- as.vector(do.call(filter_function, c(list("x" = thresh), filter_params)))
  
  # Interpolating over leap days
  smooth_seas[(leap_end + 1):(ind_max + ind_unit)] <- smooth_seas[leap_start:ind_max]
  smooth_thresh[(leap_end + 1):(ind_max + ind_unit)] <- smooth_thresh[leap_start:ind_max]
  
  for (i in leap_start:leap_end) {
    smooth_seas[i] <- (smooth_seas[i - ind_unit] + smooth_seas[i + ind_unit]) / 2
    smooth_thresh[i] <- (smooth_thresh[i - ind_unit] + smooth_thresh[i + ind_unit]) / 2
  }

  # Constructing and returning dataframe
  clm <- clm |>
    add_row(leap_days) |>
    arrange(t) |>
    mutate(
      "ind" = if_else(hourly, yhour(t), yday(t)),
      "modified_ind" = if_else(!(year %in% leap_days$year) & (ind >= leap_start), 
                               ind + ind_unit, ind),
      "seas" = smooth_seas[modified_ind],
      "thresh" = smooth_thresh[modified_ind],
    ) |>
    select(-year, -modified_ind, -included)
}

periodic_climatology <- function(ts_data, 
                                 filter_func,
                                 filter_params,
                                 extract = NULL,
                                 hourly = FALSE, 
                                 start_time = NULL,
                                 end_time = NULL,
                                 thresh_quantile = 0.9) {
  
  # Setting parameters for daily or hourly data
  if (hourly) {
    leap_start <- 1417
    leap_end <- 1440
    ind_unit <- 24
  } else {
    leap_start <- 60
    leap_end <- 60
    ind_unit <- 1
  }
  ind_max <- 365*ind_unit
  
  # Getting time frame for the climatology
  if (is.null(start_time)) {
    start_time <- min(ts_data$t)
  }
  if (is.null(end_time)) {
    end_time <- max(ts_data$t)
  }
  
  # Calculating daily or hourly index
  clm <- ts_data |>
    mutate(
      "included" = between(t, start_time, end_time),
      "hourly" = hourly,
      "year" = year(t),
      "ind" = if_else(hourly, yhour(t), yday(t))
    )
  
  # Extracting leap days
  leap_days <- clm |>
    filter((month(t) == 2) & day(t) == 29)
  
  # Filtering out leap days
  clm <- clm |>
    filter(!((month(t) == 2) & day(t) == 29)) |>
    mutate(
      "ind" = if_else((year %in% leap_days$year) & (ind > leap_end), 
                      ind - ind_unit, ind)
    )
  
  # Getting mean temperature for each index
  mean_temp <- clm |>
    filter(included) |>
    group_by(ind) |>
    summarize(
      "mean_temp" = mean(temp, na.rm = TRUE)
    ) |>
    pull(mean_temp)

  # Calculating seasonal cycle
  seas <- do.call(filter_func, c(list("x" = mean_temp), filter_params))
  
  # Extracting seasonal cycle if specified
  if (!is.null(extract)) {
    seas <- seas[[extract]]
  }
  
  # Ensuring seas is a vector
  seas <- as.vector(seas)
  
  # Interpolating over leap days
  seas[(leap_end + 1):(ind_max + ind_unit)] <- seas[leap_start:ind_max]
  for (i in leap_start:leap_end) {
    seas[i] <- (seas[i - ind_unit] + seas[i + ind_unit]) / 2
  }
  
  # Constructing and returning dataframe
  clm <- clm |>
    add_row(leap_days) |>
    arrange(t) |>
    mutate(
      "ind" = if_else(hourly, yhour(t), yday(t)),
      "modified_ind" = if_else(!(year %in% leap_days$year) & (ind >= leap_start), 
                               ind + ind_unit, ind),
      "seas" = seas[modified_ind],
      "anomaly" = temp - seas,
      "thresh" = seas + quantile(anomaly, probs = thresh_quantile, na.rm = TRUE)
    ) |>
    select(-year, -modified_ind, -included, -hourly, -anomaly)
}


periodic_climatology2 <- function(ts_data, 
                                  filter_func,
                                  filter_params,
                                  extract = NULL,
                                  hourly = FALSE, 
                                  start_time = NULL,
                                  end_time = NULL,
                                  padding = 0,
                                  thresh_quantile = 0.9) {
  
  
  # Setting parameters for daily or hourly data
  if (hourly) {
    leap_start <- 1417
    leap_end <- 1440
    ind_unit <- 24
  } else {
    leap_start <- 60
    leap_end <- 60
    ind_unit <- 1
  }
  ind_max <- 365*ind_unit
  
  # Getting time frame for the climatology
  if (is.null(start_time)) {
    start_time <- min(ts_data$t)
  }
  if (is.null(end_time)) {
    end_time <- max(ts_data$t)
  }
  
  # Calculating daily or hourly index
  clm <- ts_data |>
    mutate(
      "included" = between(t, start_time, end_time),
      "hourly" = hourly,
      "year" = year(t),
      "ind" = if_else(hourly, yhour(t), yday(t))
    )
  
  # Extracting leap days
  leap_days <- clm |>
    filter((month(t) == 2) & day(t) == 29)
  
  # Filtering out leap days
  clm <- clm |>
    filter(!((month(t) == 2) & day(t) == 29)) |>
    mutate(
      "ind" = if_else((year %in% leap_days$year) & (ind > leap_end), 
                      ind - ind_unit, ind)
    )
  
  # Getting mean temperature for each index
  mean_temp <- clm |>
    filter(included) |>
    group_by(ind) |>
    summarize(
      "mean_temp" = mean(temp, na.rm = TRUE)
    ) |>
    pull(mean_temp)
  
  # Adding window padding
  if (padding > 0) {
    mean_temp <- c(mean_temp[(ind_max-padding+1):ind_max], mean_temp, mean_temp[1:padding])
  }
  
  # Calculating seasonal cycle
  seas <- do.call(filter_func, c(list("x" = mean_temp), filter_params))
  
  # Extracting seasonal cycle if specified
  if (!is.null(extract)) {
    seas <- seas[[extract]]
  }
  
  # Ensuring seas is a vector
  seas <- as.vector(seas)
  
  # Removing window padding
  if (padding > 0) {
    seas <- seas[(padding+1):(ind_max+padding)]
  }
  
  # Interpolating over leap days
  seas[(leap_end + 1):(ind_max + ind_unit)] <- seas[leap_start:ind_max]
  for (i in leap_start:leap_end) {
    seas[i] <- (seas[i - ind_unit] + seas[i + ind_unit]) / 2
  }
  
  # Constructing and returning dataframe
  clm <- clm |>
    add_row(leap_days) |>
    arrange(t) |>
    mutate(
      "ind" = if_else(hourly, yhour(t), yday(t)),
      "modified_ind" = if_else(!(year %in% leap_days$year) & (ind >= leap_start), 
                               ind + ind_unit, ind),
      "seas" = seas[modified_ind],
      "anomaly" = temp - seas,
      "thresh" = seas + quantile(anomaly, probs = thresh_quantile, na.rm = TRUE)
    ) |>
    select(-year, -modified_ind, -included, -hourly, -anomaly)
}


fill_gaps <- function(x, max_gap) {
  
  # Initializing data
  n <- length(x)
  y <- rep(FALSE, n)
  gap_count <- Inf
  
  # Getting index of first true element
  st_ind <- which(x)[1]
  
  # Returning false vector if no true elements in x
  if(is.na(st_ind)) {
    return(y)
  }
  
  # Looping through remaining elements
  for (i in st_ind:n) {
    if (x[i]) {
      y[i] <- TRUE
      if (gap_count <= max_gap) {
        y[(i-gap_count):(i-1)] <- TRUE
      }
      gap_count <- 0
    } else {
      gap_count <- gap_count + 1
    }
  }
  
  # Filling final gap
  if (!x[i] & (gap_count <= max_gap)) {
    y[(i-gap_count+1):i] <- TRUE
  }
  
  # Returning output
  y
}

detect_events <- function(clm, min_duration = 5, max_gap = 2, exceedance = FALSE) {
  
  if (exceedance) {
    min_duration <- 1
    max_gap <- 0
  }
  
  clm |>
    mutate(
      "exceedanceCriterion" = exceedance,
      "threshCriterion" = if_else(!is.na(temp), temp > thresh, FALSE),
      "sums" = cumsum(!threshCriterion | 
                        (threshCriterion & !lag(threshCriterion, default = FALSE))),
      "sums" = if_else(exceedanceCriterion, 
                       cumsum(!threshCriterion | threshCriterion), sums)
    ) |>
    group_by(sums) |>
    mutate("event_length" = n()) |>
    ungroup() |>
    mutate(
      "event_length" = if_else(threshCriterion, event_length, 0),
      "durationCriterion" = event_length >= min_duration,
      "event" = fill_gaps(durationCriterion, max_gap = max_gap),
      "event_no" = if_else(exceedanceCriterion, 
                           cumsum(event), cumsum(event & !lag(event, default = FALSE))),
      "event_no" = if_else(event, event_no, NA),
    ) |>
    select(-exceedanceCriterion)
  
}

calculate_events <- function(clm, funcs, target_col = "intensity") {
  
  clm |>
    mutate(
      "intensity" = temp - seas,
      "exceedance" = temp - thresh
      ) |>
    group_by(event_no) |>
    summarize(
      "duration" = n(),
      "date_start" = min(t, na.rm = TRUE),
      "date_end" = max(t, na.rm = TRUE),
      "date_peak" = t[which.max(intensity)],
      across(all_of(target_col), funcs, .names = "{.col}_{.fn}")
    ) |>
    drop_na()
  
}


get_intersections <- function(data, x, y1, y2) {
  data |>
    mutate(
      "diff" = get({{y1}}) - get({{y2}}),
      "sign_change" = diff * lag(diff, default = 0) < 0,
      "prev_diff" = lag(diff),
      "prev_x" = lag(get({{x}})),
      "prev_y1" = lag(get({{y1}}))
    ) |>
    filter(sign_change) |>
    rowwise() |>
    mutate(
      "x_intersect" = approx(c(prev_diff, diff), c(prev_x, get({{x}})), xout = 0)$y,
      "y" = approx(c(prev_x, get({{x}})), c(prev_y1, get({{y1}})), xout = x_intersect)$y,
    ) |>
    select(
      !!x := x_intersect, 
      !!y1 := y,
      !!y2 := y,
    )
}

