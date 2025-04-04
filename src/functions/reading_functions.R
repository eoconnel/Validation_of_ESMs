# Name: reading_functions.R
# Author: Ethan O'Connell
# Date: April 19th, 2024
# Description: Function definitions for reading and processing 
#              FVCOM, GLORYS and temperature observation data.

# Loading libraries
library(ncdf4)
library(R.matlab)
library(readr)
library(lubridate)
library(dplyr)


##################################################################
# Defining function for reading in FVCOM data
##################################################################
readFVCOM <- function(filename, process = "depth_average") {
  
  # Opening netCDF File
  fvcom <- nc_open(filename)
  
  # Extracting Variables
  temp <- ncvar_get(fvcom, "temp")   # temperature
  lon <- ncvar_get(fvcom, "lon")     # longitude
  lat <- ncvar_get(fvcom, "lat")     # latitude
  h <- ncvar_get(fvcom, "h")         # depth
  depth <- fvcom$dim$siglay$vals
  
  # Getting Time Data
  time <- ncvar_get(fvcom, "time")   # time
  time_units <- ncatt_get(fvcom, "time", "units") # time units
  time_str <- strsplit(time_units$value, " ")
  origin <- unlist(time_str)[3]
  datetime <- as.POSIXct(as.vector(time*86400), origin=origin, tz="UTC")
  
  # Closing file
  nc_close(fvcom)
  
  # Processing data (default: depth average)
  if (process == "SST") {
    temp <- temp[1, ]
  } else if (process == "bottom") {
    temp <- temp[length(depth), ]
  } else {
    temp <- colMeans(temp)
  }
  
  # Constructing result
  list(
    data = tibble(
      "time" = datetime,
      "temp" = temp,
    ),
    latitude = lat,
    longitude = lon,
    max_depth = h
  )
}


##################################################################
# Defining function for reading in FVCOM 2.0 data
##################################################################
readFVCOM2 <- function(filename, ph_ind = 1, th_ind = 17) {
  
  # Opening matlab data file
  fvcom <- readMat(filename)
  
  # Defining site strings
  sites <- c("Port l'Hebert", "Taylor's Head Shallow")
  
  # Extracting times
  times <- fvcom$time[,]
  time <- as.POSIXct((times - 1.0)*86400, origin="0000-01-01", tz = "UTC")
  
  # Extracting Port l'Hebert and Taylor Head data
  ph_temp <- fvcom$sst[ph_ind,]
  th_temp <- fvcom$sst[th_ind,]
  
  # Constructing result
  list(
    data = tibble(
      "time" = rep(time, 2),
      "temp" = c(ph_temp, th_temp),
      "site" = rep(c("Port l'Hebert", "Taylor's Head Shallow"), 
                   each = length(times))
    ),
    sites = sites,
    latitude = fvcom$Lat[c(ph_ind, th_ind)],
    longitude = fvcom$Lon[c(ph_ind, th_ind)]
  )
}


##################################################################
# Defining function for reading in GLORYS data
##################################################################
readGLORYS <- function(filename, process = "depth_average") {
  
  # Opening netCDF file
  glorys <- nc_open(filename)
  
  # Extracting variables
  temp <- ncvar_get(glorys, "thetao")         # temperature
  temp <- na.omit(temp)                       # removing missing data
  depth <- glorys$dim$depth$vals[nrow(temp)]  # depth
  lat <- glorys$dim$latitude$vals             # latitude
  lon <- glorys$dim$longitude$vals            # longitude

  # Getting time data
  time <- glorys$dim$time$vals
  time_units <- glorys$dim$time$units
  time_str <- strsplit(time_units, " ")
  origin <- unlist(time_str)[3]
  datetime <- as.POSIXct(as.vector(time*3600), origin=origin, tz="UTC")
  
  # Closing file
  nc_close(glorys)
  
  # Processing data (default: depth average)
  if (process == "SST") {
    temp <- temp[1, ]
  } else if (process == "bottom") {
    temp <- temp[length(depth), ]
  } else {
    temp <- colMeans(temp)
  }
  
  # Constructing result
  list(
    data = tibble(
      "time" = datetime,
      "temp" = temp,
    ),
    latitude = lat,
    longitude = lon,
    max_depth = depth[length(depth)]
  )
}


##################################################################
# Defining function to read an CMIP6 data file
##################################################################
readCMIP6 <- function(filename, process = "sst") {
  
  # Reading matlab file
  cmip6 <- readMat(filename)$cmip6
  cmip6_dim <- dim(cmip6)
  
  all_time <- NULL
  all_temp <- NULL
  
  # Looping through the file dimensions
  for (i in 1:cmip6_dim[3]) {
    
    # Extracting time
    time <- as_vector(as.POSIXct((cmip6[[1 + cmip6_dim[1]*(i-1)]] - 719529)*86400, 
                                 origin = "1970-01-01", tz = "UTC"))
    
    # Appending times
    all_time <- c(time, all_time)
    all_time <- sort(all_time)
    
    # Extracting temperatures
    temp <- cmip6[[2 + cmip6_dim[1]*(i-1)]]
    
    depth_ind <- 1
    # Processing temperatures
    if (process == "bottom") {
      depth_ind <- nrow(temp)
    } else if (process == "depth_average") {
      depth_ind <- 1:nrow(temp)
    }
    temp <- colMeans(temp[depth_ind, , drop = FALSE], na.rm = TRUE)
    
    # Appending temperatures
    all_temp <- c(all_temp, temp)
  }
  
  tibble(time = all_time, temp = all_temp)
}


##################################################################
# Wrapper function for reading observation data files
##################################################################
readOBS <- function(directory = NULL, filename = NULL, sites = "all", 
                    hourly = TRUE) {
  
  # Reading single file is filename is provided
  if (!is.null(filename)) {
    return(readOBS_file(filename, sites, hourly = hourly))
  }
  
  # Checking if directory is provided
  if (is.null(directory)) {
    return(NULL)
  }
  
  # Extracting files from directory
  files <- list.files(path = directory)
  if (length(files) == 0) {
    return(NULL)
  }
  filenames <- paste(directory, files, sep = "")
  
  # Looping through files
  obs_data <- NULL
  for (file in filenames) {
    file_data <- readOBS_file(file, sites, hourly = hourly)
    obs_data <- bind_rows(obs_data, file_data) |>
      arrange(time)
  }
  
  # Returning data
  obs_data
}


##################################################################
# Defining function to read an observation data file
##################################################################
readOBS_file <- function(filename, sites, hourly = TRUE) {
  
  # Reading data
  data <- read_csv(filename, show_col_types = FALSE)
  
  # Extracting site and temperature
  data_site <- data |>
    select("time" = time.UTC, "temp" = temp.C, "site" = Site)
  
  if (all(sites != "all")) {
    data_site <- data_site |>
      filter(site %in% sites)
  }
  
  # Converting times if needed
  if (!is.POSIXct(data_site$time)) {
    data_site <- data_site |>
      mutate(
        time = mdy_hm(time)
      )
  }
  
  # Checking for averaging over hourly time grid
  if (hourly) {
    
    # Extracting hour and dates
    data_site <- data_site |>
      mutate(
        date = date(time),
        hour = hour(time)
      )
    
    # Grouping by hour and averaging
    data_site <- data_site |>
      group_by(site, date, hour) |>
      summarise(
        hourly_temp = mean(temp),
        .groups = "drop"
      )
    
    # Reconstructing time variable
    data_site <- data_site |>
      mutate(
        datetime = make_datetime(year = year(date), month = month(date), 
                                 day = day(date), hour)
      ) |>
      select("time" = datetime, temp = "hourly_temp", "site" = site)
  }
  
  # Returning time series
  data_site
}


##################################################################
# Defining function to read all data files together
##################################################################
readALL <- function(sites, process, observation_path,
                    fvcom_files, glorys_files, fvcom2_file) {
  
  # Reading observation data
  obs_data <- readOBS(
    directory = observation_path, 
    sites = sites,
    hourly = TRUE
  )
  
  # Initializing FVCOM and GLORYS data
  fvcom_data <- NULL
  glorys_data <- NULL
  
  # Reading FVCOM and GLORYS data
  for (i in 1:length(sites)) {
    
    # Reading data
    fvcom_temp <- readFVCOM(fvcom_files[i], process = process)$data
    glorys_temp <- readGLORYS(glorys_files[i], process = process)$data
    
    # Adding FVCOM data to data frame
    fvcom_temp <- mutate(fvcom_temp, "site" = sites[i])
    fvcom_data <- bind_rows(fvcom_data, fvcom_temp)
    
    # Adding GLORYS data to data frame
    glorys_temp <- mutate(glorys_temp, "site" = sites[i])
    glorys_data <- bind_rows(glorys_data, glorys_temp)
  }
  
  # Reading FVCOM2 data
  fvcom2_data <- readFVCOM2(fvcom2_file)$data
  
  # Determining minimum starting time
  start_time <- min(
    min(obs_data$time),
    min(fvcom_data$time),
    min(fvcom2_data$time),
    min(glorys_data$time)
  )
  
  # Determining maximum end time
  end_time <- max(
    max(obs_data$time),
    max(fvcom_data$time),
    max(fvcom2_data$time),
    max(glorys_data$time)
  )
  
  # Contructing common time sequence
  common_time <- seq(start_time, end_time, by = "hours")
  
  # Initilizing full sequence data
  obs_full <- NULL
  fvcom_full <- NULL
  fvcom2_full <- NULL
  glorys_full <- NULL
  all_data <- NULL
  
  # Matching data to common time sequence
  for (i in 1:length(sites)) {
    
    # Initializing common temperature vector
    obs_full_vec <- rep(NA, length(common_time))
    fvcom_full_vec <- rep(NA, length(common_time))
    fvcom2_full_vec <- rep(NA, length(common_time))
    glorys_full_vec <- rep(NA, length(common_time))
    
    # Extracting data for specific site
    obs_temp <- obs_data |>
      filter(site == sites[i])
    fvcom_temp <- fvcom_data |>
      filter(site == sites[i])
    fvcom2_temp <- fvcom2_data |>
      filter(site == sites[i])
    glorys_temp <- glorys_data |>
      filter(site == sites[i])
    
    # Matching data to common time sequence
    obs_full_vec[common_time %in% obs_temp$time] <- obs_temp$temp
    fvcom_full_vec[common_time %in% fvcom_temp$time] <- fvcom_temp$temp
    fvcom2_full_vec[common_time %in% fvcom2_temp$time] <- fvcom2_temp$temp
    glorys_full_vec[common_time %in% glorys_temp$time] <- glorys_temp$temp
    
    # Combining data to single data frame
    all_data_temp <- tibble(
      "time" = common_time,
      "site" = sites[i],
      "OBS" = obs_full_vec,
      "FVCOM" = fvcom_full_vec,
      "FVCOM2" = fvcom2_full_vec,
      "GLORYS" = glorys_full_vec,
    )
    
    # Appending data to single data frame
    all_data <- bind_rows(all_data, all_data_temp)
  }
  
  # Returning data frame
  all_data
}


##################################################################
# Defining function to interpolate over missing observations
##################################################################
interpolate_obs <- function(filename, gap_size = 12) {
  
  # Reading raw hourly data
  hourly_data <- read_csv(filename, show_col_types = FALSE)
  
  # Extracting sites
  sites <- unique(hourly_data$site)
  
  # Initializing new data
  interpolated_obs <- NULL
  
  # Looping over each site
  for (i in 1:length(sites)) {
    
    # Extracting observations from specific site
    site_obs <- hourly_data |>
      filter(site == sites[i])
    
    # Getting start and end times of observations
    st_time <- min(site_obs$time[!is.na(site_obs$OBS)])
    en_time <- max(site_obs$time[!is.na(site_obs$OBS)])
    
    # Filtering between start and end times
    filtered_site <- site_obs |>
      filter(between(time, st_time, en_time))
    
    # Getting indices of missing values
    missing_ind <- which(is.na(filtered_site$OBS))
    gaps <- 1
    
    # Identifying continuous gaps
    for (j in 2:length(missing_ind)) {
      if (missing_ind[j] == (missing_ind[j-1] + 1)) {
        gaps[j] <- gaps[j-1]
      } else {
        gaps[j] <- gaps[j-1] + 1
      }
    }
    
    # Eliminating gaps of more than 12 hours
    gap_numbers <- tibble(gaps) |>
      group_by(gaps) |>
      summarize(count = n()) |>
      filter(count <= gap_size) |>
      pull(gaps)
    
    # Linearizing over gaps of less than 12 hours
    for (j in 1:length(gap_numbers)) {
      ind <- missing_ind[gaps == gap_numbers[j]]
      slope <- (filtered_site$OBS[ind[length(ind)]+1] - filtered_site$OBS[ind[1]-1]) / (length(ind) + 1)
      filtered_site$OBS[ind] <- filtered_site$OBS[ind[1]-1] + slope*(1:length(ind))
    }
    
    # Replacing data with interpolated version
    site_obs$OBS[site_obs$time %in% filtered_site$time] <- filtered_site$OBS
    
    # Removing linearized data at Port l'Hebert from 2020-05-15 to 2020-05-25
    if (sites[i] == "Port l'Hebert") {
      site_obs$OBS[between(site_obs$time, 
                           ymd_hms("2020-05-15 00:00:00"), 
                           ymd_hms("2020-05-25 23:00:00"))] <- NA
    }
    
    interpolated_obs <- bind_rows(interpolated_obs, site_obs)
  }
  
  # Returning interpolated data
  interpolated_obs
}


##################################################################
# Defining function to average to daily data
##################################################################
average_daily <- function(filename) {
  
  # Reading hourly data
  hourly <- read_csv(filename, show_col_types = FALSE)
  
  # Extracting sites
  sites <- unique(hourly_data$site)
  
  # Initializing daily data
  daily_data <- NULL
  
  # Looping over each site
  for (i in 1:length(sites)) {
    
    # Extracting site data
    hourly_site <- hourly |>
      filter(site == sites[i])
    
    # Identifying observations with incomplete days
    obs_missing <- hourly_site |>
      select(time, site, OBS) |>
      drop_na() |>
      mutate("day" = date(time)) |>
      group_by(day) |>
      summarise(hours = n(), .groups = "drop") |>
      filter(hours < 24)
    
    # Identifying FVCOM with incomplete days
    fvcom_missing <- hourly_site |>
      select(time, site, FVCOM) |>
      drop_na() |>
      mutate("day" = date(time)) |>
      group_by(day) |>
      summarise(hours = n(), .groups = "drop") |>
      filter(hours < 12)
    
    # Identifying FVCOM2 with incomplete days
    fvcom2_missing <- hourly_site |>
      select(time, site, FVCOM2) |>
      drop_na() |>
      mutate("day" = date(time)) |>
      group_by(day) |>
      summarise(hours = n(), .groups = "drop") |>
      filter(hours < 12)
    
    # Averaging over daily observations
    obs_daily <- hourly_site |>
      select(time, site, OBS) |>
      drop_na() |>
      filter(!(date(time) %in% obs_missing$day)) |>
      mutate("day" = date(time)) |>
      group_by(day) |>
      summarize(
        temp = mean(OBS)
      ) |>
      rename("time" = day)
    
    # Averaging over daily FVCOM data
    fvcom_daily <- hourly_site |>
      select(time, site, FVCOM) |>
      drop_na() |>
      filter(!(date(time) %in% fvcom_missing$day)) |>
      mutate("day" = date(time)) |>
      group_by(day) |>
      summarise(
        temp = mean(FVCOM)
      ) |>
      rename("time" = day)
    
    # Averaging over daily FVCOM2 data
    fvcom2_daily <- hourly_site |>
      select(time, site, FVCOM2) |>
      drop_na() |>
      filter(!(date(time) %in% fvcom2_missing$day)) |>
      mutate("day" = date(time)) |>
      group_by(day) |>
      summarise(
        temp = mean(FVCOM2)
      ) |>
      rename("time" = day)
    
    # Extracting daily GLORYS data
    glorys_daily <- hourly_site |>
      select(time, site, GLORYS) |>
      drop_na() |>
      mutate("day" = date(time)) |>
      select("time" = day, site, "temp" = GLORYS)
    
    # Getting overall start time
    st_time <- min(
      min(obs_daily$time),
      min(fvcom_daily$time),
      min(fvcom2_daily$time),
      min(glorys_daily$time)
    )
    
    # Getting overall end time
    en_time <- max(
      max(obs_daily$time),
      max(fvcom_daily$time),
      max(fvcom2_daily$time),
      max(glorys_daily$time)
    )
    
    # Constructing common time sequence
    common_time <- seq(st_time, en_time, by = "days")
    
    # Initializing common temperature vector
    obs_full_vec <- rep(NA, length(common_time))
    fvcom_full_vec <- rep(NA, length(common_time))
    fvcom2_full_vec <- rep(NA, length(common_time))
    glorys_full_vec <- rep(NA, length(common_time))
    
    # Matching data to common time sequence
    obs_full_vec[common_time %in% obs_daily$time] <- obs_daily$temp
    fvcom_full_vec[common_time %in% fvcom_daily$time] <- fvcom_daily$temp
    fvcom2_full_vec[common_time %in% fvcom2_daily$time] <- fvcom2_daily$temp
    glorys_full_vec[common_time %in% glorys_daily$time] <- glorys_daily$temp
    
    # Combining data to single data frame
    daily_data_temp <- tibble(
      "time" = common_time,
      "site" = sites[i],
      "OBS" = obs_full_vec,
      "FVCOM" = fvcom_full_vec,
      "FVCOM2" = fvcom2_full_vec,
      "GLORYS" = glorys_full_vec,
    )
    
    # Appending daily data
    daily_data <- bind_rows(daily_data, daily_data_temp)
  }
  
  # Returning daily data
  daily_data
}

