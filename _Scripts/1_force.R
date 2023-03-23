###########################################
### Force Data Import & Analysis script ###
###########################################

# written by Austin Hurst, Jack Solomon & Devan Pancura
# correspondence: jack.solomon@dal.ca devanpancura@dal.ca

### Import required packages ###

library(dplyr)
library(tidyr)
library(ggplot2)

source("./_Scripts/_settings.R")
source("./_Scripts/utils.R")


### Run after import script ###


### Additional cleaning of the trial metadata ###

# Check for and remove any extra trials per session

trialdat <- trialdat %>%
  filter(frame <= 248)

### Remove bad participants ###

participant_dat <- participant_dat %>%
  na.omit()

bad_id1 <- participant_dat %>%
  subset(id == 4)

bad_id2 <- participant_dat %>%
  subset(id == 34)

bad_id3 <- participant_dat %>%
  subset(id == 48)

dropped_participants <- full_join(bad_id1, bad_id2)
dropped_participants <- full_join(dropped_participants, bad_id3)

participant_dat <- anti_join(participant_dat, dropped_participants)

### Preprocess signal data ###

# Discard unnecessary signal data following MEP period

signaldat <- signaldat %>%
  subset(time * 1000 <= max_time)


# Discard rest trials from signal data 

signaldat <- signaldat %>%
  subset(frame!= "21")


### Remove unnecessary data to shorten dataset ###

#Separating countdown data (first 3 seconds) from the rest of trial

countdown <- signaldat %>%
  subset(time < 3)

#Separating GO! data (2 seconds of data following countdown) from the rest of trial

signaldat <- signaldat %>%
  subset(time > 3 & time < 5.2)


#Downsampling Force data

forcedat <- signaldat %>%
  group_by(id, frame) %>%
  filter(row_number() %% 10 == 0) %>%
  ungroup()

#Add trial number to Signal Data

trialmap <- signaldat %>%
  group_by(id, frame) %>%
  summarize(trial = 1) %>%
  mutate(trial = cumsum(trial)) %>%
  ungroup()

forcedat <- forcedat %>%
  left_join(trialmap, by = c("id", "frame")) %>%
  select(c(id, frame, trial, time, force))

#Add target and max grip to signal data

forcedat <- forcedat %>%
  group_by(id, trial) %>%
  right_join(
    select(labviewdat, c(id, trial, target, maxGrip)),
    by = c("id", "trial")
  )%>%
  select(id, frame, trial, time, force, target, maxGrip) %>%
  ungroup()

#Add order to signal data

forcedat <- forcedat %>%
  right_join(
    select(trialdat, c(id, frame, order)),
    by = c("id", "frame")
  )%>%
  select(frame, id, trial, time, force, target, maxGrip, order)

#Change force to current grip based on max grip

forcedat <- forcedat %>%
  mutate(force = round(((force-1)/(1-maxGrip))*-100,2))

#Separate by order (ME first vs MI first) and drop MI trials from force analysis

order1force <- forcedat %>%
  subset(order == 1 & time > 4.8 & time < 5 & trial < 81)

order2force <- forcedat %>%
  subset(order == 2 & time > 4.8 & time < 5 & trial > 80)

order1force <- order1force %>%
  group_by(id, trial) %>%
  summarize (sumforce = mean(force))

order2force <- order2force %>%
  group_by(id, trial) %>%
  summarize (sumforce = mean(force))

order1force <- order1force %>%
  right_join(
    select(labviewdat, c(id, trial, target)),
    by = c("id", "trial")
  )%>%
  select(id, trial, target, sumforce)

order2force <- order2force %>%
  right_join(
    select(labviewdat, c(id, trial, target)),
    by = c("id", "trial")
  )%>%
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

force_data <- full_join(force_data_1, force_data_2)

bad1 <- force_data %>%
  subset(id == 4)
bad2 <- force_data %>%
  subset(id == 34)
bad = full_join(bad1, bad2)
bad3 <- force_data %>%
  subset(id == 48)
bad = full_join(bad, bad3)

force_data <- anti_join(force_data, bad)

#Discard extra data frames (saves space, big files)

rm(drop_trials_1, drop_trials_2, force_data_1, force_data_2, order1force, order2force, trialmap, col_overrides)


#Some descriptives about force trial data

ggplot(force_data, aes(x = target)) +
  geom_bar(color = "navyblue", fill = "navyblue", width = 3) +
  labs(x = "Force (% MVC)", y = "Number of Successful Trials") +
  geom_text(aes(label = ..count..), stat = "count", size = 8, vjust = 1.5, colour = "white") +
  plot_theme


#Visualizing participant accuracy across levels of effector load

all_data <- full_join(force_data, drop_trials)

bad1 <- all_data %>%
  subset(id == 4)
bad2 <- all_data %>%
  subset(id == 34)
bad = full_join(bad1, bad2)
bad3 <- all_data %>%
  subset(id == 48)
bad = full_join(bad, bad3)

all_data <- anti_join(all_data, bad)

accuracy <- all_data %>%
  group_by(id, trial, target) %>%
  summarize(error_rate = (abs(sumforce-target)/target)*100)

accuracy <- accuracy %>%
  group_by(id, trial, target) %>%
  summarize(accuracy = 100-error_rate)

accuracy <- accuracy %>%
  group_by(id, target) %>%
  summarize(accuracy,
            sd = sd(accuracy))

ggplot(accuracy, aes(x = target, y = accuracy)) +
  geom_point(color = "navyblue", size = 2) +
  labs(x = "Force (% MVC)", y = "Accuracy (%)") +
  geom_errorbar(aes(
    ymin = accuracy - sd, 
    ymax = accuracy + sd), size = 0.3, width = 1, col = "navyblue") +
  plot_theme
  
avg_accuracy <- accuracy %>%
  group_by(target) %>%
  summarize(avg_acc = mean(accuracy),
            sd = sd(accuracy))

plot_theme <- theme_classic() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    legend.key = element_blank(),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 22),
    axis.text = element_text(color = "black", size = 18),
    axis.title.y = element_text(margin = margin(r = 0.5, unit = "cm")),
    axis.title.x = element_text(margin = margin(t = 0.5, unit = "cm")),
    legend.title = element_text(size = 22),
    axis.ticks = element_line(color = "black", size = 1),
    legend.text = element_text(size = 22),
    legend.spacing.x = unit(0.5, "cm"),
    legend.box.margin = margin(l = 1, unit = "cm"),
  )


ggplot(avg_accuracy, aes(x = target, y = avg_acc)) +
  geom_point(size = 4.5, col = "navyblue") +
  geom_errorbar(aes(
    ymin = avg_acc - sd, 
    ymax = avg_acc + sd), size = 0.5, col = "navyblue") +
  xlab("Force (% MVC)") +
  ylab("Accuracy (%)") +
  plot_theme
