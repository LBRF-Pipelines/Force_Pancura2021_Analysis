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


# Remove participants where we were unable to establish RMT

no_rmt_ids <- subset(participant_dat, is.na(rmt))$id


# Remove participants with equipment issues during collection
# (force transducer not working properly, 0% grip accuracy)

bad_ids <- c(4, 48)


# Drop all bad IDs from participant data

dropped_participants <- participant_dat %>%
  filter(id %in% bad_ids | id %in% no_rmt_ids)

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


# Get the mean force in the 200 ms pre-pulse for each trial

force_means <- forcedat %>%
  subset(!(id %in% bad_ids)) %>%
  subset(time > 4.8 & time < 5) %>%
  group_by(id, trial) %>%
  summarize(
    order = order[1],
    imagery = ((order[1] == 1) == (trial > 80))[1],
    target = target[1],
    sumforce = mean(force)
  )


# Check whether each physical trial was within 3% of its target force

force_means_pp <- force_means %>%
  subset(!imagery) %>%
  mutate(
    target_achieved = abs(sumforce - target) < 3
  )

bad_by_force <- subset(force_means_pp, !target_achieved) %>%
  mutate(bad_by_force = TRUE)



#Visualizing participant accuracy across levels of effector load

accuracy <- force_means_pp %>%
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
