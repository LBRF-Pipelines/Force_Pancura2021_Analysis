#######################################
### MEP preprocessing import script ###
#######################################

# Author: Austin Hurst


### Import required packages ###

library(readr)
library(purrr)
library(tibble)
library(dplyr)
library(tidyr)
library(arrow)

source("./_Scripts/_settings.R")
source("./_Scripts/utils.R")

options(readr.show_progress = FALSE)



### Import trial metadata ###

# Get file lists

mep_trial_files <- list.files(
  "./_Data/frames/", pattern = "*.csv",
  full.names = TRUE
)

signal_files <- list.files(
  "./_Data/signals/", pattern = "*.arrow",
  full.names = TRUE
)


# Import participant master list

master_path <- "./_Data/Master List.csv"

participant_dat <- read_csv(master_path, skip = 5, col_types = cols()) %>%
  filter(`Participant ID` != "PILOT") %>%
  mutate(id = as.numeric(gsub("P([0-9]+).*$", "\\1", `Participant ID`))) %>%
  rename(
    kviq = `KVIQ (kinesthetic)`
  ) %>%
  select(c(id, Age, Gender, Handedness, Order, kviq, RMT))

names(participant_dat) <- tolower(names(participant_dat))


# Import Signal frame metadata

trialdat <- map_df(mep_trial_files, function(f) {
  pid <- as.integer(gsub("[Pp]_([0-9]+).*$", "\\1", basename(f)))
  df <- read_csv(f, col_types = cols(), progress = FALSE)
  add_column(df, id = pid, .before = 1)
})


# Import & extract the force target for each trial

targetdat <- map_df(signal_files, function(f) {

  # Import signal data and grab the force target channel
  df <- read_feather(f, col_select = c("frame", "time", "target"))

  # Extract the trial's target force from the signal
  df <- df %>%
    filter(time > 0.5 & time < 1.0) %>%
    group_by(frame) %>%
    summarize(force_target = round(mean(target) * 10) * 10)

  # Append the participant ID and return the data
  pid <- as.integer(gsub("[Pp]_([0-9]+).*$", "\\1", basename(f)))
  add_column(df, id = pid, .before = 1)

})


# Wrangle metadata to be easier to work with

trialdat <- trialdat %>%
  select(-c(pulse_interval, pwr_b)) %>%
  left_join(targetdat, by = c("id", "frame")) %>%
  left_join(select(participant_dat, c(id, order)), by = "id") %>%
  mutate(
    order = as.factor(order),
    state = as.factor(state),
    force_target = as.factor(force_target)
  ) %>%
  relocate(order, .after = id)


# Get frame numbers of rest trials for each ID

rest_trials <- trialdat %>%
  filter(state == "HOLD") %>%
  select(c(id, frame))



### Import EMG & force transducer data ###

signaldat <- map_df(signal_files, function(f) {

  # Import the EMG and force transducer signal data
  df <- read_feather(f, col_select = c("frame", "time", "fds", "force"))

  # Append the participant ID
  pid <- as.integer(gsub("[Pp]_([0-9]+).*$", "\\1", basename(f)))
  df <- add_column(df, id = pid, .before = 1)

  # Drop rest trials
  rest_for_id <- subset(rest_trials, id == pid)
  df <- df %>%
    anti_join(rest_for_id, by = "frame")

  df
})

signaldat <- signaldat %>%
  rename(emg = fds)



### Import labview data ###

labviewfiles <- list.files(
  "./_Data/labview/", pattern = "*.csv",
  full.names = TRUE
)

# Set column types for subsequent import

labview_names <- c("maxGrip", "condition", "block", "trial", "target")

col_overrides <- cols(
  maxGrip = col_double(),
  condition = col_double(),
  block = col_double(),
  trial = col_double(),
  target = col_double(),
)

# Import all transducer data into a single data frame

labviewdat <- map_df(labviewfiles, function(f) {
  id_num <- as.numeric(gsub("^p(\\d+).*", "\\1", basename(f)))
  data <-  read_csv(f, col_names = labview_names, col_types = col_overrides)
  data <- add_column(data, id = id_num, .before = 1)
  data
})

labviewdat <- labviewdat %>%
  relocate(trial, .before = 2)
