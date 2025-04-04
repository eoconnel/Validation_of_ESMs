library(tidyverse)

plot_event <- function(clm, event_num, spread = NULL, 
                       start_date = NULL, end_date = NULL) {
  
  # Getting start and end dates
  if (!is.null(spread)) {
    start_date <- min(clm$t[clm$event_no == event_num], na.rm = TRUE) - days(spread)
    end_date <- max(clm$t[clm$event_no == event_num], na.rm = TRUE) + days(spread)
  }
  
  # Finding intersection points of temp and thresh
  intersections <- clm |>
    mutate("t" = as.numeric(as_datetime(t))) |>
    get_intersections("t", "temp", "thresh") |>
    mutate("t" = as_datetime(t))
  
  # Constructing dataframe for geom_ribbon
  area_clm <- clm |>
    mutate(t = as_datetime(t)) |>
    bind_rows(intersections) |>
    arrange(t) |>
    mutate(
      "event_no" = pmax(lag(event_no), lead(event_no), na.rm = TRUE)
    ) |>
    filter(event_no == event_num) |>
    filter(temp >= thresh)
  
  # Filtering climatology
  event_clm <- clm |>
    mutate(t = as_datetime(t)) |>
    filter(between(t, start_date, end_date))
  
  # Pivot longer for plotting
  clm_long <- event_clm |>
    pivot_longer(
      c(temp, seas, thresh),
      names_to = "type",
      values_to = "value"
    )
  
  # Constructing MHW plot
  clm_long |>
    ggplot(aes(x = t)) +
    geom_ribbon(data = area_clm, aes(ymin = thresh, ymax = temp), 
                lwd = 0.7, fill = "red", alpha = 0.9) +
    geom_line(aes(y = value, col = type), lwd = 0.7) +
    scale_color_manual(
      name = "",
      values = c("black", "blue", "darkgreen"),
      breaks = c("temp", "seas", "thresh"),
      labels = c("Temperature", "Climatology", "Threshold")
    ) +
    scale_x_datetime(date_labels = "%b %Y") +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black",
                                       linewidth = 0.5, linetype = "solid")
    )
}

plot_all_events <- function(clm, start_date = NULL, end_date = NULL) {
  
  # Getting start and end dates for plotting
  if (is.null(start_date)) {
    start_date <- min(clm$t)
  }
  if (is.null(end_date)) {
    end_date <- max(clm$t)
  }
  
  # Finding intersection points of temp and thresh
  intersections <- clm |>
    mutate("t" = as.numeric(as_datetime(t))) |>
    get_intersections("t", "temp", "thresh") |>
    mutate("t" = as_datetime(t))
  
  # Constructing dataframe for geom_ribbon
  area_clm <- clm |>
    mutate(t = as_datetime(t)) |>
    bind_rows(intersections) |>
    arrange(t) |>
    mutate(
      "event" = lag(event) | lead(event)
    ) |>
    filter(between(t, start_date, end_date)) |>
    filter(event)
  
  # Filtering climatology
  event_clm <- clm |>
    mutate(t = as_datetime(t)) |>
    filter(between(t, start_date, end_date))
  
  # Pivot longer for plotting
  clm_long <- event_clm |>
    pivot_longer(
      c(temp, seas, thresh),
      names_to = "type",
      values_to = "value"
    )
  
  # Constructing MHW plot
  clm_long |>
    ggplot(aes(x = t)) +
    geom_ribbon(data = area_clm, aes(ymin = thresh, ymax = temp), 
                lwd = 0.7, fill = "red", alpha = 0.9) +
    geom_line(aes(y = value, col = type), lwd = 0.7) +
    scale_color_manual(
      name = "",
      values = c("black", "blue", "darkgreen"),
      breaks = c("temp", "seas", "thresh"),
      labels = c("Temperature", "Climatology", "Threshold")
    ) +
    scale_x_datetime(date_labels = "%b %Y") +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black",
                                       linewidth = 0.5, linetype = "solid")
    )
}

plot_exceedances <- function(clm, start_date = NULL, end_date = NULL) {
  
  # Getting start and end dates for plotting
  if (is.null(start_date)) {
    start_date <- min(clm$t)
  }
  if (is.null(end_date)) {
    end_date <- max(clm$t)
  }
  
  # Finding intersection points of temp and thresh
  intersections <- clm |>
    mutate("t" = as.numeric(as_datetime(t))) |>
    get_intersections("t", "temp", "thresh") |>
    mutate("t" = as_datetime(t))
  
  # Constructing dataframe for geom_ribbon
  event_clm <- clm |>
    mutate(t = as_datetime(t)) |>
    bind_rows(intersections) |>
    arrange(t) |>
    mutate(
      "event" = event | lag(event) | lead(event),
      "exceedance" = if_else(temp >= thresh, temp - thresh, NA)
    ) |>
    filter(between(t, start_date, end_date))
  
  # Constructing MHW plot
  event_clm |>
    ggplot(aes(x = t)) +
    geom_ribbon(data = filter(event_clm, event), aes(ymin = 0, ymax = exceedance), 
                lwd = 0.7, fill = "red", alpha = 0.9) +
    geom_line(aes(y = exceedance), na.rm = TRUE, lwd = 0.7) +
    geom_hline(yintercept = 0, col = "darkgreen") +
    scale_x_datetime(date_labels = "%b %Y") +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black",
                                       linewidth = 0.5, linetype = "solid")
    )
}

plot_event_metric <- function(events, ind, metric, start_ind = NULL, end_ind = NULL) {
  
  # Getting start and end dates for plotting
  if (is.null(start_ind)) {
    start_ind <- min(events[[ind]])
  }
  if (is.null(end_ind)) {
    end_ind <- max(events[[ind]])
  }
  
  events_plot <- events |>
    filter(between(get({{ind}}), start_ind, end_ind)) |>
    add_row(!!ind := start_ind, !!metric := NA) |>
    add_row(!!ind := end_ind, !!metric := NA) |>
    ggplot(aes(x = .data[[ind]], y = .data[[metric]])) +
    geom_segment(aes(xend = .data[[ind]], yend = 0), na.rm = TRUE) +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw()
  
  if (is.Date(events[[ind]][1])) {
    events_plot <- events_plot + scale_x_date(date_labels = "%b %Y")
  } else if (is.POSIXct(events[[ind]][1])) {
    events_plot <- events_plot + scale_x_datetime(date_labels = "%b %Y")
  }
  
  events_plot
}

compare_event <- function(clms, nums, labels, start_date, end_date) {
  
  n <- length(clms)
  clm <- NULL
  area_clm <- NULL
  for (i in 1:n) {
    
    intersections <- clms[[i]] |>
      mutate("t" = as.numeric(as_datetime(t))) |>
      get_intersections("t", "temp", "thresh") |>
      mutate("t" = as_datetime(t))
    
    area_clm <- clms[[i]] |>
      mutate("t" = as_datetime(t)) |>
      bind_rows(intersections) |>
      arrange(t) |>
      mutate(
        "event_no" = pmax(lag(event_no), lead(event_no), na.rm = TRUE),
        "label" = labels[[i]]
      ) |>
      filter(event_no == nums[[i]]) |>
      bind_rows(area_clm) |>
      arrange(t)
    
    clm <- clms[[i]] |>
      mutate("t" = as_datetime(t)) |>
      filter(between(t, start_date, end_date)) |>
      mutate(
        "label" = labels[[i]]
      ) |>
      bind_rows(clm) |>
      arrange(t)
  }
  
  # Pivot climatologies longer for plotting
  clm_long <- clm |>
    pivot_longer(
      c(temp, seas, thresh),
      names_to = "type",
      values_to = "value"
    ) |>
    mutate("label" = factor(label, levels = labels))
  
  area_clm <- area_clm |>
    mutate("label" = factor(label, levels = labels))
  
  # Constructing MHW plot
  clm_long |> 
    ggplot(aes(x = t)) +
    geom_ribbon(data = area_clm, aes(ymin = thresh, ymax = temp), 
                lwd = 0.7, fill = "red", alpha = 0.9) +
    geom_line(aes(y = value, col = type), lwd = 0.7) +
    facet_wrap(~label, ncol = 1) +
    scale_color_manual(
      name = "",
      values = c("black", "blue", "darkgreen"),
      breaks = c("temp", "seas", "thresh"),
      labels = c("Temperature", "Climatology", "Threshold")
    ) +
    scale_x_datetime(date_labels = "%b %Y") +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black",
                                       linewidth = 0.5, linetype = "solid")
    )
}

compare_all_events <- function(clms, labels, start_date, end_date) {
  
  n <- length(clms)
  clm <- NULL
  area_clm <- NULL
  
  for (i in 1:n) {
    
    # Finding intersection points of temp and thresh
    intersections <- clms[[i]] |>
      mutate("t" = as.numeric(as_datetime(t))) |>
      get_intersections("t", "temp", "thresh") |>
      mutate("t" = as_datetime(t))
    
    area_clm <- clms[[i]] |>
      mutate("t" = as_datetime(t)) |>
      bind_rows(intersections) |>
      arrange(t) |>
      mutate(
        "event" = event | lag(event) | lead(event),
        "label" = labels[[i]]
      ) |>
      filter(between(t, start_date, end_date)) |>
      filter(event) |>
      bind_rows(area_clm) |>
      arrange(t)
    
    clm <- clms[[i]] |>
      mutate("t" = as_datetime(t)) |>
      filter(between(t, start_date, end_date)) |>
      mutate(
        "label" = labels[[i]]
      ) |>
      bind_rows(clm) |>
      arrange(t)
  }
  
  # Pivot climatologies longer for plotting
  clm_long <- clm |>
    pivot_longer(
      c(temp, seas, thresh),
      names_to = "type",
      values_to = "value"
    ) |>
    mutate("label" = factor(label, levels = labels))
  
  area_clm <- area_clm |>
    mutate("label" = factor(label, levels = labels))
  
  # Constructing MHW plot
  clm_long |>
    ggplot(aes(x = t)) +
    geom_ribbon(data = area_clm, aes(ymin = thresh, ymax = temp), 
                lwd = 0.7, fill = "red", alpha = 0.9) +
    geom_line(aes(y = value, col = type), lwd = 0.7) +
    facet_wrap(~label, ncol = 1) +
    scale_color_manual(
      name = "",
      values = c("black", "blue", "darkgreen"),
      breaks = c("temp", "seas", "thresh"),
      labels = c("Temperature", "Climatology", "Threshold")
    ) +
    scale_x_datetime(date_labels = "%b %Y") +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black",
                                       linewidth = 0.5, linetype = "solid")
    )
}

compare_exceedances <- function(clms, labels, start_date, end_date) {
  
  n <- length(clms)
  event_clm <- NULL
  
  for (i in 1:n) {
    
    intersections <- clms[[i]] |>
      mutate("t" = as.numeric(as_datetime(t))) |>
      get_intersections("t", "temp", "thresh") |>
      mutate("t" = as_datetime(t))
    
    event_clm <- clms[[i]] |>
      mutate("t" = as_datetime(t)) |>
      bind_rows(intersections) |>
      arrange(t) |>
      mutate(
        "event" = event | lag(event) | lead(event),
        "exceedance" = if_else(temp >= thresh, temp - thresh, NA),
        "label" = labels[[i]]
      ) |>
      filter(between(t, start_date, end_date)) |>
      bind_rows(event_clm) |>
      arrange(t)
  }
  
  event_clm <- event_clm |>
    mutate("label" = factor(label, levels = labels))
  
  # Constructing MHW plot
  event_clm |>
    ggplot(aes(x = t)) +
    geom_ribbon(data = filter(event_clm, event), aes(ymin = 0, ymax = exceedance), 
                lwd = 0.7, fill = "red", alpha = 0.9) +
    geom_line(aes(y = exceedance), na.rm = TRUE, lwd = 0.7) +
    geom_hline(yintercept = 0, col = "darkgreen") +
    scale_x_datetime(date_labels = "%b %Y") +
    facet_wrap(~label, ncol = 1) +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black",
                                       linewidth = 0.5, linetype = "solid")
    )
}

compare_event_metric <- function(events1, events2, ind, metric, 
                                 label1, label2, start_ind = NULL, end_ind = NULL) {
  
  # Getting start and end dates for plotting
  if (is.null(start_ind)) {
    start_ind <- min(events1[[ind]], events2[[ind]])
  }
  if (is.null(end_ind)) {
    end_ind <- max(events1[[ind]], events2[[ind]])
  }
  
  # Combining events
  events <- bind_rows(
    mutate(events1, "label" = label1),
    mutate(events2, "label" = label2)
  )
  
  events_plot <- events |>
    filter(between(get({{ind}}), start_ind, end_ind)) |>
    add_row(!!ind := start_ind, !!metric := NA) |>
    add_row(!!ind := end_ind, !!metric := NA) |>
    ggplot(aes(x = .data[[ind]], y = .data[[metric]])) +
    geom_segment(aes(xend = .data[[ind]], yend = 0)) +
    facet_wrap(~label, ncol = 1) +
    labs(
      x = "",
      y = "Temperature [\u00B0C]"
    ) +
    theme_bw()
  
  if (is.Date(events[[ind]][1])) {
    events_plot <- events_plot + scale_x_date(date_labels = "%b %Y")
  } else if (is.POSIXct(events[[ind]][1])) {
    events_plot <- events_plot + scale_x_datetime(date_labels = "%b %Y")
  }
  
  events_plot
}

taylor_diagram <- function(x, y, sd_max, labels, print_results = FALSE) {
  
  sd_vec <- 1
  cor_vec <- 1
  
  for (i in 1:ncol(y)) {
    sd_vec[i+1] <- sd(y[,i]) / sd(x)
    cor_vec[i+1] <- cor(x, y[,i])
  }
  
  if (print_results) {
    print(sd_vec)
    print(cor_vec)
  }
  
  taylor_data <- tibble(
    "sd" = sd_vec,
    "cor" = cor_vec,
    "label" = factor(labels, levels = labels)
  ) 
  
  taylor_data |>
    ggplot(aes(x = cor, y = sd, col = label, shape = label)) +
    geom_point(size = 3) +
    coord_radial(start = 0, end = pi/2, expand = FALSE) +
    scale_x_continuous(
      name = "Normalized standard deviation",
      limits = c(0,1),
      breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
      labels = c(" R = 0", " R = 0.1", " R = 0.2", " R = 0.3", " R = 0.4",
                 " R = 0.5", " R = 0.6", " R = 0.7", " R = 0.8", " R = 0.9",
                 " R = 1")
    ) +
    scale_y_continuous(
      name = "",
      breaks = c(0, 0.5, 0.75, 1, 1.25),
      limits = c(0, sd_max),
      sec.axis = dup_axis()
    ) +
    labs(
      col = "",
      shape = ""
    )
}


taylor_diagram2 <- function(x, y, sd_max, sites, types, print_results = FALSE) {
  
  sd_vec <- 1
  cor_vec <- 1
  n_site <- length(sites)
  n_type <- length(types)
  
  sd_mat <- matrix(NA, nrow = n_site, ncol = n_type)
  cor_mat <- matrix(NA, nrow = n_site, ncol = n_type)
  
  sd_mat[,1] <- 1
  cor_mat[,1] <- 1
  
  for (i in 1:n_site) {
    for (j in 1:(n_type-1)) {
      x_temp <- x[,i]
      y_temp <- y[, (i-1)*(n_type-1) + j]
      
      sd_mat[i, (j+1)] <- sd(y_temp) / sd(x_temp)
      cor_mat[i, (j+1)] <- cor(x_temp, y_temp)
    }
  }
  
  if (print_results) {
    print(sd_mat)
    print(cor_mat)
  }
  
  sd_vec <- as.numeric(t(sd_mat))
  cor_vec <- as.numeric(t(cor_mat))
  
  taylor_data <- tibble(
    "sd" = sd_vec,
    "cor" = cor_vec,
    "type" = factor(rep(types, n_site), levels = types),
    "site" = factor(rep(sites, each = n_type), levels = sites)
  ) 
  
  taylor_data |>
    ggplot(aes(x = cor, y = sd, col = type, shape = type)) +
    geom_point(size = 3) +
    coord_radial(start = 0, end = pi/2, expand = FALSE) +
    scale_x_continuous(
      name = "Normalized standard deviation",
      limits = c(0,1),
      breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
      labels = c(" R = 0", " R = 0.1", " R = 0.2", " R = 0.3", " R = 0.4",
                 " R = 0.5", " R = 0.6", " R = 0.7", " R = 0.8", " R = 0.9",
                 " R = 1")
    ) +
    scale_y_continuous(
      name = "",
      breaks = c(0, 0.5, 0.75, 1, 1.25),
      limits = c(0, sd_max),
      sec.axis = dup_axis()
    ) +
    facet_wrap(~site, ncol = 2) +
    labs(
      col = "",
      shape = ""
    )
}
