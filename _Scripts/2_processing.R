################################
### MEP summarization script ###
################################

# Author: Austin Hurst & Devan Pancura


### Run after force script ###

# Grab 1000 ms of EMG at rest (0.5 to 1.5s) for muscle activity comparison

emg_rest <- signaldat %>%
  subset(time >= 0.5 & time < 1.5) %>%
  select(-force)


# Isolate time windows for MEPs (1000 ms pre-pulse to 200 ms post-pulse)

signaldat <- signaldat %>%
  subset(time > 4 & time < 5.2)


# Add trial number and order to EMG data

trialnums <- signaldat %>%
  group_by(id, frame) %>%
  summarize(trial = 1) %>%
  mutate(trial = cumsum(trial))

signaldat <- signaldat %>%
  left_join(trialnums, by = c("id", "frame")) %>%
  select(c(id, trial, time, emg, force))

signaldat <- signaldat %>%
  left_join(
    select(participant_dat, c(id, order)),
    by = c("id")
  ) %>%
  select(c(id, trial, order, time, emg, force))

emg_rest <- emg_rest %>%
  left_join(trialnums, by = c("id", "frame")) %>%
  select(c(id, trial, time, emg))


# Drop bad IDs & ME trials in which target force was not achieved

signaldat <- signaldat %>%
  subset(!(id %in% bad_ids)) %>%
  anti_join(bad_by_force, by = c("id", "trial"))


# Demean and filter line noise from the EMG signal for each trial

signaldat <- signaldat %>%
  group_by(id, trial) %>%
  mutate(
    emg = emg - median(emg),
    emg_filt = dft_filt(emg, line_hz, srate, window = 1000)
  )

emg_rest <- emg_rest %>%
  group_by(id, trial) %>%
  mutate(
    emg = emg - median(emg),
    emg_filt = dft_filt(emg, line_hz, srate)
  )


# Summarize EMG activity during rest for reference

resting_emg <- emg_rest %>%
  group_by(id, trial) %>%
  summarize(
    rms_pre = sqrt(mean(emg_filt ** 2)),
    madmed_pre = mad(emg_filt)
  )


# Summarize EMG activity in pre-pulse period for later checking

prepulse_emg <- signaldat %>%
  filter(time * 1000 < pulse_on - 1) %>%
  group_by(id, trial) %>%
  summarize(
    rms = sqrt(mean(emg_filt ** 2)),
    madmed = mad(emg_filt),
    max_deviation = max(abs(emg_filt - median(emg_filt)) / mad(emg_filt))
  ) %>%
  left_join(resting_emg, by = c("id", "trial"))


# Summarize TMS pulse artifact amplitudes in EMG signal

pulse_tmin <- (pulse_on - 5) / 1000
pulse_tmax <- (max_pulse_on + 5) / 1000

pulse_onset <- 5000 # ms
pulse_window <- 4 # ms
pulse <- 1

emg_pulses <- signaldat %>%
  filter(time >= pulse_tmin & time <= pulse_tmax) %>%
  filter(time * 1000 >= (pulse_onset)) %>%
  filter(time * 1000 <= (pulse_onset + pulse_window)) %>%
  summarize(
    pulse_amp = max(abs(emg_filt))
  )


# Extract post-pulse MEP windows following each supra-threshold pulse

# NOTE: MEP windowing could be refined further using average MEPs

mep_windows <- signaldat %>%
  filter(time * 1000 >= pulse_on) %>%
  filter(time * 1000 <= (max_pulse_on + 100)) %>%
  filter(time * 1000 <= (pulse_onset + mep_window[2])) %>%
  filter(time * 1000 >= (pulse_onset + mep_window[1])) %>%
  select(c(id, trial, time, emg, emg_filt))


# Detect and correct for recordings with inverted electrodes

mean_meps <- mep_windows %>%
  mutate(n = 1:n()) %>%
  group_by(id, n) %>%
  mutate(
    mean_emg = mean(emg_filt)
  )

peak_orders <- mean_meps %>%
  group_by(id) %>%
  summarize(
    nmin = min(n[mean_emg == min(mean_emg)]),
    nmax = min(n[mean_emg == max(mean_emg)]),
    max_first = nmax < nmin
  ) %>%
  select(c(id, max_first))

signaldat <- signaldat %>%
  left_join(peak_orders, by = "id") %>%
  mutate(
    emg = ifelse(max_first, emg, -emg),
    emg_filt = ifelse(max_first, emg_filt, -emg_filt)
  ) %>%
  select(-max_first)

mep_windows <- mep_windows %>%
  left_join(peak_orders, by = "id") %>%
  mutate(
    emg = ifelse(max_first, emg, -emg),
    emg_filt = ifelse(max_first, emg_filt, -emg_filt)
  ) %>%
  select(-max_first)



### Summarize MEPs using EMG data ###

# Detect peaks/amplitudes of MEPs

mep_peaks <- mep_windows %>%
  group_by(id, trial) %>%
  group_modify(
    ~ detect_peaks(.x$emg_filt, .x$time)
  ) %>%
  mutate(
    ampl = vmax - vmin,
    wonky = tmin < tmax
  )


# Detect onsets of MEPs

mep_onsets <- mep_windows %>%
  left_join(
    select(mep_peaks, c(id, trial, tmax, wonky)),
    by = c("id", "trial")
  ) %>%
  filter(time <= tmax & !wonky) %>%
  left_join(
    select(prepulse_emg, c(id, trial, madmed)),
    by = c("id", "trial")
  ) %>%
  group_by(id, trial) %>%
  summarize(
    onset = detect_onset(emg_filt, time, madmed)
  )

mep_peaks <- mep_peaks %>%
  left_join(mep_onsets, by = c("id", "trial"))


# Join summarized MEP data to trial metadata

sr_meps <- signaldat %>%
  right_join(
    subset(mep_peaks),
    by = c("id", "trial")
  ) %>%
  left_join(
    prepulse_emg, by = c("id", "trial")
  )



### Removing bad MEPs from data ###

# Find MI trials where pre-pulse EMG activity is > 3x the activity during rest
# or 4x the overall median pre-pulse EMG during imagery

prepulse_emg <- prepulse_emg %>%
  subset(!(id %in% bad_ids)) %>%
  left_join(select(participant_dat, c(id, order)), by = "id") %>%
  mutate(
    imagery = ((order == 1) == (trial > 80)),
    madm_ratio = madmed / madmed_pre
  )

median_madm_mi <- median(subset(prepulse_emg, imagery)$madmed)

bad_by_emg <- prepulse_emg %>%
  ungroup() %>%
  mutate(
    bad_by_ratio = imagery & (madm_ratio > 3),
    bad_by_madm = imagery & (madmed > 4 * median_madm_mi)
  ) %>%
  subset(bad_by_ratio | bad_by_madm)


# Find MEPs with small amplitudes (< 25 mV) and weird MEP peak times

bad_by_ampl <- mep_peaks %>%
  subset(!(id %in% bad_ids) & !wonky) %>%
  subset(ampl < 25)

bad_by_tmax <- mep_peaks %>%
  subset(!(id %in% bad_ids) & !wonky & ampl > 25) %>%
  group_by(id) %>%
  mutate(
    med_tmax = median(tmax),
    tmax_diff = tmax - med_tmax
  ) %>%
  subset(abs(tmax_diff) > 0.0071)


# Mark 'missed trigger' trials

missedtrigger <- mep_peaks %>%
  subset(id %in% c(1, 15, 40, 42) & trial == 81)


# Remove all bad trials from MEP data

trimmed_meps <- mep_peaks %>%
  subset(!(id %in% bad_ids)) %>%
  subset(!wonky) %>%
  anti_join(bad_by_emg, by = c("id", "trial")) %>%
  anti_join(bad_by_ampl, by = c("id", "trial")) %>%
  anti_join(bad_by_tmax, by = c("id", "trial")) %>%
  anti_join(missedtrigger, by = c("id", "trial"))


#Separate by order (ME first vs MI first) and isolate MI trials for filtering

order1emg <- sr_meps %>%
  subset(order == 1 & trial > 80)

order2emg <- sr_meps %>%
  subset(order == 2 & trial < 81)

order1emg_ME <- sr_meps %>%
  subset(order == 1 & trial < 81)

order2emg_ME <- sr_meps %>%
  subset(order == 2 & trial > 80)

#Drop wonky trials

order1emg <- order1emg %>% #229 trials dropped
  filter(wonky != TRUE)

order2emg <- order2emg %>% #181 trials dropped
  filter(wonky != TRUE)

order1emg_ME <- order1emg_ME %>% #24 trials dropped
  filter(wonky != TRUE)

order2emg_ME <- order2emg_ME %>% #39 trials dropped
  filter(wonky != TRUE)

#Drop trials with EMG > 2 SDs from mean

order1emg <- order1emg %>% #38 trials dropped
  subset(time < 5) %>%
  summarize(id, trial, order, time, emg, force, emg_filt, vmax, vmin, tmax, tmin, ampl, wonky, onset, rms, madmed, max_deviation, z_madmed,
            mean = mean(emg)) %>%
  filter(z_madmed < 2 && z_madmed > -2)

order2emg <- order2emg %>% #1 trial dropped
  subset(time < 5) %>%
  summarize(id, trial, order, time, emg, force, emg_filt, vmax, vmin, tmax, tmin, ampl, wonky, onset, rms, madmed, max_deviation, z_madmed,
            mean = mean(emg)) %>%
  filter(z_madmed < 2 && z_madmed > -2)

#Rejoin data frames once bad EMG MI trials are dropped

MI_meps <- full_join(order1emg, order2emg)

ME_meps <- full_join(order1emg_ME, order2emg_ME)

sr_meps <- full_join(MI_meps, ME_meps)

rm(order1emg, order1emg_ME, order2emg, order2emg_ME, MI_meps, ME_meps) #don't need dfs takes up space


#Cut extra columns of data for visualization

trimmed_meps <- sr_meps %>%
  subset(select = -c(emg, emg_filt, vmax, vmin, tmax, tmin, wonky, rms, madmed, max_deviation))

#Add target to trimmed MEP data

trimmed_meps <- trimmed_meps %>%
  right_join(
    select(labviewdat, c(id, trial, target)),
    by = c("id", "trial")
  )%>%
  select(id, trial, order, time, force, ampl, onset, target)

#Drop missed trigger trials

bad_trial1 <-trimmed_meps %>%
  subset(id == 1 & trial == 81)

bad_trial2 <- trimmed_meps %>%
  subset(id == 15 & trial == 81)

bad_trial3 <- trimmed_meps %>%
  subset(id == 40 & trial == 81)

bad_trial4 <- trimmed_meps %/%
  subset(id == 42 & trial == 81)

missedtrigger <- full_join(bad_trial1, bad_trial2)
missedtrigger <- full_join(missedtrigger, bad_trial3)
missedtrigger <- full_join(missedtrigger, bad_trial4)

trimmed_meps <- anti_join(trimmed_meps, missedtrigger)

#Drop Bad Participants
#These participants have multiple levels of the DV effector load in which there were no successful trials
#and therefore no MEPs to use for analysis. >50% of trials were dropped, therefore the participants are
#dropped from further analysis

bad_id1 <- trimmed_meps %>%
  subset(id == 4)

bad_id2 <- trimmed_meps %>%
  subset(id == 34)

bad_id3 <- trimmed_meps %>%
  subset(id == 48)

bad = full_join(bad_id1, bad_id2)

dropped_participants <- full_join(bad_id3, bad)

trimmed_meps <- anti_join(trimmed_meps, dropped_participants)

wonkyones <- sr_meps %>%
  subset(wonky != FALSE)

trial_ex <- sr_meps %>%
  subset(id == 60 & trial == 149)

ggplot(trial_ex, aes(x = time, y = emg_filt)) +
  geom_line() +
  ylab("Filtered EMG") +
  xlab("Time") +
  plot_theme

#Separate trimmed data by block and condition

trimmed_MI_1 <- trimmed_meps %>%
  subset(order == 1 & trial > 80) #MI first

trimmed_MI_2 <- trimmed_meps %>% #MI second
  subset(order == 2 & trial < 81)

trimmed_ME_2 <- trimmed_meps %>%
  subset(order == 1 & trial < 81) #ME second

trimmed_ME_1 <- trimmed_meps %>%
  subset(order == 2 & trial > 80) #ME first

trimmed_MI <- full_join(trimmed_MI_1, trimmed_MI_2) #rejoin conditions for H1/2

trimmed_ME <- full_join(trimmed_ME_1, trimmed_ME_2) #rejoin conditions for H1/2


#Create one average MEP value per target in each df

trimmed_MI <- trimmed_MI %>%
  group_by(id, target) %>%
  summarize (mean_amp = mean(ampl))

trimmed_ME <- trimmed_ME %>%
  group_by(id, target) %>%
  summarize (mean_amp = mean(ampl))

trimmed_MI_1 <- trimmed_MI_1 %>%
  group_by(id, target) %>%
  summarize (mean_amp = mean(ampl))

trimmed_MI_2 <- trimmed_MI_2 %>%
  group_by(id, target) %>%
  summarize (mean_amp = mean(ampl))

trimmed_ME_1 <- trimmed_ME_1 %>%
  group_by(id, target) %>%
  summarize (mean_amp = mean(ampl))

trimmed_ME_2 <- trimmed_ME_2 %>%
  group_by(id, target) %>%
  summarize (mean_amp = mean(ampl))


#Playing around with plotting participant data

trimmed_ME2_test <- trimmed_ME_2 %>%
  group_by(target) %>%
  summarize(id, 
            mean = mean_amp,
            sd = sd(mean_amp),
            label = "ME",
            order = 2)

trimmed_ME1_test <- trimmed_ME_1 %>%
  group_by(target) %>%
  summarize(id,
            mean = mean_amp,
            sd = sd(mean_amp),
            label = "ME",
            order = 1)

trimmed_MI1_test <- trimmed_MI_1 %>%
  group_by(target) %>%
  summarize(id,
            mean = mean_amp,
            sd = sd(mean_amp),
            label = "MI",
            order = 1)

trimmed_MI2_test <- trimmed_MI_2 %>%
  group_by(target) %>%
  summarize(id, mean = mean_amp,
            sd = sd(mean_amp),
            label = "MI",
            order = 2)

trimmed_data_test <-trimmed_ME2_test %>%
  group_by(target) %>%
  summarize(id, mean, sd, label, order) %>%
  full_join(trimmed_ME1_test) %>%
  full_join(trimmed_MI1_test) %>%
  full_join(trimmed_MI2_test) %>%
  mutate(order = as.factor(order),
         label = as.factor(label),
         id = as.factor(id))

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

pd = position_dodge(5)
plot1_colors <- c("lightsalmon", "indianred1", "tomato", "firebrick1", "firebrick2",
                  "firebrick3", "firebrick","darkred", "darkorange4",
                  "darkorange3", "darkorange2", "darkorange", "orange", "tan1",
                  "tan2", "tan3", "tan4", "burlywood4", "peachpuff4", "seashell4", "slategray",
                  "steelblue", "steelblue4", "steelblue3", "royalblue1", "royalblue3",
                  "royalblue4", "navyblue", "midnightblue", "black", "grey36",
                  "lightsalmon", "indianred1", "tomato", "firebrick1", "firebrick2",
                  "firebrick3", "firebrick","darkred", "darkorange4", "darkorange3",
                  "darkorange2", "darkorange", "orange", "tan1")

#Displays individual participants raw data (for observing trends)

participant <- trimmed_data_test %>%
  subset(id == 60)

ggplot(participant, aes(x = target, y = mean, col = order, group = order)) +
  geom_point(size = 1) +
  geom_pointrange(aes(
    ymin = mean - sd, 
    ymax = mean + sd), position = pd, size = 0.5) +
  scale_color_manual(values = plot_colors) +
  ylab("Mean MEP Amplitude") +
  xlab("Force") +
  facet_wrap(~label, nrow = 1) +
  plot_theme + 
  labs(group = "Order", color = "Order")

#Displays the raw data for all participants in one graph

trimmed_data_test <- trimmed_data_test %>%
  mutate(label = recode(label, ME = "Overt", MI = "Imagery"))

ggplot(trimmed_data_test, aes(x = id, y = mean, col = id, group = id)) +
  geom_point(size = 1) +
  geom_pointrange(aes(
    ymin = mean - sd, 
    ymax = mean + sd), position = pd, size = 0.5) +
  scale_color_manual(values = plot1_colors) + #, labels = c("1", "2", "3", "4", "5", "6", "7")) +
  ylab("Mean MEP Amplitude (uV)") +
  xlab("Participant ID") +
  facet_wrap(~label*target, nrow = 2) +
  plot_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(group = " Participant ID", color = "Participant ID")

#Create 1 total MEP value per target across all participants

trimmed_MI <- trimmed_MI %>%
  group_by(target) %>%
  summarize (mean_amp = mean(mean_amp))

trimmed_ME <- trimmed_ME %>%
  group_by(target) %>%
  summarize(mean_amp = mean(mean_amp))

trimmed_MI_1test <- trimmed_MI_1 %>%
  group_by(target) %>%
  summarize(mean = mean(mean_amp),
            sd = sd(mean_amp),
            label = "MI",
            order = 1)

trimmed_MI_2test <- trimmed_MI_2 %>%
  group_by(target) %>%
  summarize(mean = mean(mean_amp),
            sd = sd(mean_amp),
            label = "MI",
            order = 2)

trimmed_ME_1test <- trimmed_ME_1 %>%
  group_by(target) %>%
  summarize(mean = mean(mean_amp),
            sd = sd(mean_amp),
            label = "ME",
            order = 1)

trimmed_ME_2test <- trimmed_ME_2 %>%
  group_by(target) %>%
  summarize(mean = mean(mean_amp),
            sd = sd(mean_amp),
            label = "ME",
            order = 2) %>%
  full_join(trimmed_ME_1test) %>%
  full_join(trimmed_MI_1test) %>%
  full_join(trimmed_MI_2test) %>%
  mutate(order = as.factor(order),
         label = as.factor(label))


#Code samples from Austin

ggplot(trimmed_MI, aes(x = target, y = mean_amp)) + 
  labs(x = "Target", y = "MEP Amplitude") + 
  geom_col(color = "red", width = 3)
    
ggplot(trimmed_ME, aes(x = target, y = mean_amp)) + 
  labs(x = "Target", y = "MEP Amplitude") + 
  geom_col(color = "green", width = 3)


ggplot (trimmed_MI_1, aes(x = target, y = mean_amp)) +
  labs(x = "Target", y = "MEP Amplitude") +
  geom_col(color = "purple", width = 3)

ggplot(trimmed_MI_2, aes(x = target, y = mean_amp)) + 
  labs(x = "Target", y = "MEP Amplitude") + 
  geom_col(color = "orange", width = 3)

  
pd <- position_dodge(3)
plot_colors <- c("navyblue", "darkorange3")

trimmed_ME_2test <- trimmed_ME_2test %>%
  mutate(label = recode(label, ME = "Overt", MI = "Imagery"))

ggplot(trimmed_ME_2test, aes(x = target, y = mean, col = order, group = order)) + #look for outlier with id as col and group
  geom_point(position = pd, size = 4) +
  geom_pointrange(aes(
    ymin = mean - sd, 
    ymax = mean + sd), position = pd, size = 1) +
  scale_color_manual(values = plot_colors, labels = c("1", "2")) +
  xlab("Force (% MVC)") +
  ylab("Mean MEP Amplitude (uV)") +
  facet_wrap(~label, nrow = 1) +
  plot_theme +
  labs(group = "Order", color = "Order")

