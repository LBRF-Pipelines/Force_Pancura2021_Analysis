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
library(tidyverse)

source("./_Scripts/_settings.R")
source("./_Scripts/_importcfs.R")
source("./_Scripts/utils.R")



### Import trial metadata files ###

# Get file list

cfs_files <- list.files(
  "./_Data/cfs/", pattern = "*.cfs",
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


# Import trial metadata

trialvars <- c('StateL', 'Power A')

verbose <- TRUE

optlog("\n - Importing trial metadata")

trialdat <- map_df(cfs_files, function(f) {

  # Import CFS file into R
  cfs_dat <- import_cfs(f, c("Target"), trialvars)
  optlog(".")

  # Append participant ID
  pid <- as.integer(gsub("[Pp]_([0-9]+).*$", "\\1", basename(f)))
  df <- cfs_dat$trial %>%
    add_column(id = pid, .before = 1)

  # Summarize force targets from signal data
  force_targets <- cfs_dat$signal %>%
    group_by(frame) %>%
    filter(time > 0.5 & time < 1.0) %>%
    summarize(force_target = round(mean(Target) * 10) * 10)

  # Join the summarized force data to the rest of the trial data
  df %>%
    left_join(force_targets, by = "frame")

})


# Wrangle metadata to be easier to work with

trialdat <- trialdat %>%
  rename(
    state = StateL,
    pwr_a = `Power A`
  ) %>%
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



### Import signal data ###

optlog("\n - Importing signal data")
signaldat <- map_df(cfs_files, function(f) {

  # Import CFS file into R
  cfs_dat <- import_cfs(f, c("FDS", "Force"), c())
  optlog(".")

  # Append participant ID
  pid <- as.integer(gsub("[Pp]_([0-9]+).*$", "\\1", basename(f)))
  df <- cfs_dat$signal %>%
    add_column(id = pid, .before = 1)

  # Drop rest trials
  rest_for_id <- subset(rest_trials, id == pid)
  df <- df %>%
    anti_join(rest_for_id, by = "frame")

  df
})

signaldat <- signaldat %>%
  rename(emg = FDS, force = Force)



### Import labview data ###

labviewfiles <- list.files(
  "./_Data/csv/", pattern = "*.csv",
  full.names = TRUE
)

# Set column types for subsequent import

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
  data = read_csv(f,col_names = c("maxGrip", "condition", "block", "trial", "target"), col_types = col_overrides)
  data <- add_column(data, id = id_num, .before = 1)
  data
})

labviewdat<- labviewdat %>%
  relocate(trial, .before = 2)


### Write out imported/wrangled data to an RDS ###

# Note: This doesn't save much time, but it can avoid issues with R being
# unstable after loading the reticulate package

optlog("\n - Saving data to .Rdata...\n")
save(
  participant_dat, trialdat, signaldat,
  file = "./imported_data.RData"
)

