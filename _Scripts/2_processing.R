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
    max_deviation = max(abs(emg_filt - median(emg_filt)) / mad(emg_filt)),
    minmax = max(emg_filt) - min(emg_filt)
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
  group_by(id, trial) %>%
  mutate(time = time - min(time)) %>%
  group_by(id, time) %>%
  summarize(
    mean_emg = mean(emg_filt)
  )

peak_orders <- mean_meps %>%
  group_by(id) %>%
  summarize(
    tmin_avg = min(time[mean_emg == min(mean_emg)]),
    tmax_avg = min(time[mean_emg == max(mean_emg)]),
    max_first = tmax_avg < tmin_avg
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
