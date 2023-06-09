##################################
### Force MEP modelling script ###
##################################

# Author: Austin Hurst & Devan Pancura



### Removing bad MEPs from data ###

# NOTE: For this study, electrodes are on tricep area, which has a lot more
# motor units (groups of motor neurons + muscle fibres) next to one another,
# greatly increasing the likelihood of eliciting MEPs from other nearby muscles
# with a pulse, which can seemingly often happen in the inverted direction

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


# Find MEPs with very small amplitudes (< 25 mV) or amplitudes smaller than
# the pre-pulse noise for that trial

bad_by_ampl <- mep_peaks %>%
  subset(!(id %in% bad_ids)) %>%
  subset(ampl < 25) %>%
  mutate(bad_by_ampl = TRUE)

bad_by_snr <- mep_peaks %>%
  left_join(prepulse_emg, by = c("id", "trial")) %>%
  subset(!(id %in% bad_ids) & ampl > 25) %>%
  subset((ampl / minmax) <= 1) %>%
  mutate(bad_by_snr = TRUE)


# Find MEPs with implausibly early peak onsets (< 15 ms post-pulse)

bad_by_onset <- mep_peaks %>%
  left_join(prepulse_emg, by = c("id", "trial")) %>%
  subset(!(id %in% bad_ids) & ampl > 25 & (ampl / minmax) > 1) %>%
  mutate(
    peak_time = ifelse(tmin < tmax, tmin, tmax)
  ) %>%
  subset(peak_time < 5.015) %>%
  mutate(bad_by_onset = TRUE)


# Mark 'missed trigger' trials

missedtrigger <- mep_peaks %>%
  subset(id %in% c(1, 15, 40, 42) & trial == 81) %>%
  mutate(missed_trigger = TRUE)


# Remove all bad trials from MEP data

trimmed_meps <- mep_peaks %>%
  subset(!(id %in% bad_ids)) %>%
  anti_join(bad_by_emg, by = c("id", "trial")) %>%
  anti_join(bad_by_ampl, by = c("id", "trial")) %>%
  anti_join(bad_by_snr, by = c("id", "trial")) %>%
  anti_join(bad_by_onset, by = c("id", "trial")) %>%
  anti_join(missedtrigger, by = c("id", "trial"))

all_bads <- mep_peaks %>%
  anti_join(trimmed_meps, by = c("id", "trial"))



### Gather / prepare data for models ###

# Join trial metadata to mep peak data

mep_dat <- trimmed_meps %>%
  left_join(select(participant_dat, c(id, order)), by = "id") %>%
  mutate(
    imagery = ((order == 1) == (trial > 80))
  ) %>%
  left_join(
    select(labviewdat, c(id, trial, target)), by = c("id", "trial")
  )


# Check for participants with no MEPs in a given cell or that are missing
# more than 50% of all trials in a block

cell_counts <- mep_dat %>%
  mutate(
    imagery = as.factor(ifelse(imagery, "Imagery", "Physical")),
    target = as.factor(target)
  ) %>%
  group_by(id, imagery, target, .drop = FALSE) %>%
  summarize(
    count = n()
  ) %>%
  mutate(
    target = as.factor(as.numeric(as.character(target)))
  )

too_many_missing <- cell_counts %>%
  group_by(id, imagery) %>%
  summarize(
    pct_usable = mean(count) * 10
  ) %>%
  subset(pct_usable < 50)

missing_cells <- cell_counts %>%
  subset(count == 0)

mep_dat <- mep_dat %>%
  anti_join(missing_cells, by = "id") %>%
  anti_join(too_many_missing, by = "id")


# Check for participants with excessively low imagery abilities (per the KVIQ)

low_kviq <- subset(participant_dat, kviq < 10)

mep_dat <- mep_dat %>%
  anti_join(low_kviq, by = "id")


# Z-score the MEP amplitudes separately for each participant

mep_dat_z <- mep_dat %>%
  group_by(id, imagery) %>%
  mutate(
    ampl_z = (ampl - mean(ampl)) / sd(ampl),
    vmax_z = (vmax - mean(vmax)) / sd(vmax),
  ) %>%
  mutate(
    imagery = as.factor(imagery),
    target = as.factor(target),
    order = as.factor(ifelse(order == 2, "MI->ME", "ME->MI"))
  )
