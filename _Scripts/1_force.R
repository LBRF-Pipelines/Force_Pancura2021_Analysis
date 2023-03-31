###########################################
### Force Data Import & Analysis script ###
###########################################

# written by Austin Hurst, Jack Solomon & Devan Pancura
# correspondence: jack.solomon@dal.ca devanpancura@dal.ca

### Import required packages ###

library(dplyr)
library(tidyr)

source("./_Scripts/_settings.R")
source("./_Scripts/utils.R")


### Run after import script ###


### Additional cleaning of the trial metadata ###

# Check for and remove any extra trials per session

trialdat <- trialdat %>%
  filter(frame <= 248)


# Remove bad participants

bad_ids <- c(4, 34, 48)

dropped_participants <- participant_dat %>%
  filter(id %in% bad_ids | is.na(rmt))

participant_dat <- anti_join(participant_dat, dropped_participants)



### Preprocess signal data ###

# Discard unnecessary signal data following MEP period

signaldat <- signaldat %>%
  subset(time * 1000 <= max_time)



### Remove unnecessary data to shorten dataset ###

# Get & downsample force transducer data for window of interest
# (2.2 seconds of data after GO!)

forcedat <- signaldat %>%
  select(-c(emg)) %>%
  subset(time > 3 & time < 5) %>%
  group_by(id, frame) %>%
  filter(row_number() %% 10 == 0) %>%
  ungroup()


# Add trial number to Signal Data

trialmap <- forcedat %>%
  group_by(id, frame) %>%
  summarize(trial = 1) %>%
  mutate(trial = cumsum(trial)) %>%
  ungroup()

forcedat <- forcedat %>%
  left_join(trialmap, by = c("id", "frame")) %>%
  select(c(id, frame, trial, time, force))


# Add target and max grip to signal data

forcedat <- forcedat %>%
  group_by(id, trial) %>%
  right_join(
    select(labviewdat, c(id, trial, target, maxGrip)),
    by = c("id", "trial")
  ) %>%
  select(id, frame, trial, time, force, target, maxGrip) %>%
  ungroup()


# Add order to signal data

forcedat <- forcedat %>%
  right_join(
    select(trialdat, c(id, frame, order)),
    by = c("id", "frame")
  ) %>%
  select(frame, id, trial, time, force, target, maxGrip, order)


# Change force to current grip based on max grip

forcedat <- forcedat %>%
  group_by(id, trial) %>%
  mutate(force = round(((force - 1) / (1 - maxGrip)) * -100, 2))


# Separate by order (ME first vs MI first) & drop MI trials from force analysis

order1force <- forcedat %>%
  subset(order == 1 & time > 4.8 & time < 5 & trial < 81)

order2force <- forcedat %>%
  subset(order == 2 & time > 4.8 & time < 5 & trial > 80)

order1force <- order1force %>%
  group_by(id, trial) %>%
  summarize(sumforce = mean(force))

order2force <- order2force %>%
  group_by(id, trial) %>%
  summarize(sumforce = mean(force))

order1force <- order1force %>%
  right_join(
    select(labviewdat, c(id, trial, target)),
    by = c("id", "trial")
  ) %>%
  select(id, trial, target, sumforce)

order2force <- order2force %>%
  right_join(
    select(labviewdat, c(id, trial, target)),
    by = c("id", "trial")
  ) %>%
  select(id, trial, target, sumforce)


#Determine if target force was achieved

order1force <- order1force %>%
  mutate(target_achieved = abs(sumforce - target) < 3)

order2force <- order2force %>%
  mutate(target_achieved = abs(sumforce - target) < 3)


#Remove bad trials

force_data_1 <- order1force %>%
  subset(target_achieved != FALSE)

force_data_2 <- order2force %>%
  subset(target_achieved != FALSE)

drop_trials_1 <- order1force %>%
  subset(target_achieved != TRUE)

drop_trials_2 <- order2force %>%
  subset(target_achieved != TRUE)

drop_trials <- full_join(drop_trials_1, drop_trials_2)


#Rejoin force data

force_data <- full_join(force_data_1, force_data_2) %>%
  filter(!(id %in% bad_ids))


#Discard extra data frames (saves space, big files)

rm(drop_trials_1, drop_trials_2, force_data_1, force_data_2, order1force, order2force, trialmap, col_overrides)



#Visualizing participant accuracy across levels of effector load

all_data <- full_join(force_data, drop_trials) %>%
  filter(!(id %in% bad_ids))

accuracy <- all_data %>%
  group_by(id, trial, target) %>%
  summarize(error_rate = (abs(sumforce - target) / target) * 100)

accuracy <- accuracy %>%
  group_by(id, trial, target) %>%
  summarize(accuracy = 100 - error_rate)

accuracy <- accuracy %>%
  group_by(id, target) %>%
  mutate(
    accuracy,
    sd = sd(accuracy)
  )
  
avg_accuracy <- accuracy %>%
  group_by(target) %>%
  summarize(
    avg_acc = mean(accuracy),
    sd = sd(accuracy)
  )
