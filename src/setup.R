
# Loading functions
source("functions/mhw_functions.R")
source("functions/evt_functions.R")
source("functions/plotting_functions.R")
source("functions/utility_functions.R")
source("functions/reading_functions.R")

# Defining file names
files <- list(
  observation_path = "../data/observations_all_sites_2017_2020/",
  fvcom = c("../data/Port_LHebert/FVCOM_Temperature_PortLHebert_2017_to_2022.nc",
            "../data/TaylorHead/FVCOM_Temperature_TaylorHead_2017_to_2022.nc"),
  glorys = c("../data/Port_LHebert/GLORYS_Temperature_PortLHebert_2000to2020.nc",
             "../data/TaylorHead/GLORYS_Temperature_TaylorHead_2000to2020.nc"),
  fvcom2 = "../data/fvcom2/eshore_v12_md_step02_v00_new.mat",
  hourly = "../data/all_data_hourly.csv",
  hourly_clean = "../data/all_data_hourly_clean.csv",
  daily_clean = "../data/all_data_daily_clean.csv" 
)

# Defining data availability periods
dates <- list(
  ph_obs = c(ymd("2017-05-27"), ymd("2021-10-23")),
  th_obs = c(ymd("2017-06-30"), ymd("2021-11-16")),
  fvcom = c(ymd("2017-01-01"), ymd("2022-12-31")),
  fvcom2 = c(ymd("2018-01-01"), ymd("2020-02-28")),
  glorys = c(ymd("2000-01-01"), ymd("2020-12-31")),
  freq = c(ymd("2018-05-28"), ymd("2020-03-20")),
  extreme = c(ymd("2018-01-01"), ymd("2021-10-23"))
)

# Reading raw data files if applicable
if (read_data) {
    
  # Reading raw hourly data
  hourly_data <- readALL(sites, process, files$observation_path,
                         files$fvcom, files$glorys, files$fvcom2)
    
  # Writing hourly data to single csv
  write_csv(hourly_data, file = files$hourly)
    
  # Interpolating hourly data
  hourly_clean <- interpolate_obs(files$hourly)
    
  # Writing cleaned hourly data to single csv
  write_csv(hourly_clean, file = files$hourly_clean)
    
  # Averaging over daily data
  daily_clean <- average_daily(files$hourly_clean)
    
  # Writing daily data to single csv
  write_csv(daily_clean, file = files$daily_clean)
  
  # Removing data variables
  rm("daily_clean", "hourly_clean", "hourly_data")
}

# Reading daily and hourly data
daily_dataframe <- read_csv(files$daily_clean, show_col_types = FALSE)
hourly_dataframe <- read_csv(files$hourly_clean, show_col_types = FALSE)

# Replacing site labels
for (i in 1:length(sites)) {
  hourly_dataframe$site[hourly_dataframe$site == sites[i]] <- site_labels[i]
  daily_dataframe$site[daily_dataframe$site == sites[i]] <- site_labels[i]
}

daily_data <- list(
  glorys = daily_dataframe |>
    select("t" = time, "temp" = GLORYS, site) |>
    filter(between(t, dates$glorys[1], dates$glorys[2])),
  fvcom = daily_dataframe |>
    select("t" = time, "temp" = FVCOM, site) |>
    filter(between(t, dates$fvcom[1], dates$fvcom[2])),
  fvcom2 = daily_dataframe |>
    select("t" = time, "temp" = FVCOM2, site) |>
    filter(between(t, dates$fvcom2[1], dates$fvcom2[2])),
  obs = daily_dataframe |>
    select("t" = time, "temp" = OBS, site) |>
    filter(between(t, dates$ph_obs[1], dates$th_obs[2]))
)

hourly_data <- list(
  glorys = hourly_dataframe |>
    select("t" = time, "temp" = GLORYS, site) |>
    filter(between(t, dates$glorys[1], dates$glorys[2])),
  fvcom = hourly_dataframe |>
    select("t" = time, "temp" = FVCOM, site) |>
    filter(between(t, dates$fvcom[1], dates$fvcom[2])),
  fvcom2 = hourly_dataframe |>
    select("t" = time, "temp" = FVCOM2, site) |>
    filter(between(t, dates$fvcom2[1], dates$fvcom2[2])),
  obs = hourly_dataframe |>
    select("t" = time, "temp" = OBS, site) |>
    filter(between(t, dates$ph_obs[1], dates$th_obs[2]))
)
