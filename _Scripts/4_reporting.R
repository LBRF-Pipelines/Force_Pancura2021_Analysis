### Import required libraries ###

library(dplyr)
library(tidyr)
library(emmeans)
library(bayestestR)


### Get demographics info ###

demographics <- participant_dat %>%
  semi_join(mep_dat, by = "id") %>%
  mutate(
    gender = as.factor(gender),
    handedness = as.factor(handedness),
    order = as.factor(order)
  )

summary(demographics)


# Get KVIQ means per group

demographics %>%
  group_by(order) %>%
  summarize(kviq = mean(kviq))



### Summarize dropped trial info ###

bad_trials <- labviewdat %>%
  select(c(id, trial, target)) %>%
  left_join(select(participant_dat, c(id, order)), by = "id") %>%
  mutate(
    imagery = ((order == 1) == (trial > 80))
  ) %>%
  left_join(
    select(bad_by_force, c(id, trial, bad_by_force)),
    by = c("id", "trial")
  ) %>%
  left_join(
    select(bad_by_emg, c(id, trial, bad_by_ratio, bad_by_madm)),
    by = c("id", "trial")
  ) %>%
  left_join(
    select(bad_by_ampl, c(id, trial, bad_by_ampl)),
    by = c("id", "trial")
  ) %>%
  left_join(
    select(bad_by_snr, c(id, trial, bad_by_snr)),
    by = c("id", "trial")
  ) %>%
  left_join(
    select(bad_by_onset, c(id, trial, bad_by_onset)),
    by = c("id", "trial")
  ) %>%
  # Replace NAs from joins with FALSE
  mutate(across(where(is.logical), ~ replace_na(.x, FALSE))) %>%
  # Exclude dropped participants
  semi_join(mep_dat, by = "id")

bad_trials_summary <- bad_trials %>%
  mutate(
    bad_by_madm = bad_by_madm & !bad_by_ratio,
    bad_by_ampl = bad_by_ampl & !(bad_by_madm | bad_by_ratio),
    bad_by_snr = bad_by_snr & !(bad_by_ampl | bad_by_madm | bad_by_ratio),
    bad_by_onset = (
      bad_by_onset & !(bad_by_snr | bad_by_ampl | bad_by_madm | bad_by_ratio)
    ),
    bad_by_any = (
      bad_by_force | bad_by_onset | bad_by_snr | bad_by_ampl |
      bad_by_madm | bad_by_ratio
    )
  ) %>%
  group_by(imagery, target) %>%
  summarize(
    bad_by_force = mean(bad_by_force) * 100,
    bad_by_ratio = mean(bad_by_ratio) * 100,
    bad_by_madm = mean(bad_by_madm) * 100,
    bad_by_ampl = mean(bad_by_ampl) * 100,
    bad_by_snr = mean(bad_by_snr) * 100,
    bad_by_onset = mean(bad_by_onset) * 100,
    usable = mean(!bad_by_any) * 100
  )



### Summarize models ###

# Summarize model for overt execution trials

emm_phys <- emmeans(mod_phys, ~ target * order, level = 0.9)

emm_phys_rel <- contrast(
  emmeans(mod_phys, ~ target), interaction = "trt.vs.ctrl"
)
describe_posterior(emm_phys_rel, ci = 0.9, rope_ci = 1)

emm_phys_seq <- contrast(
  emmeans(mod_phys, ~ target), interaction = "consec"
)
describe_posterior(emm_phys_seq, ci = 0.9, rope_ci = 1)

emm_phys_diffs <- contrast(emm_phys, interaction = "trt.vs.ctrl")
describe_posterior(emm_phys_diffs, ci = 0.9, rope_ci = 1)


# Summarize model for imagery trials

emm_mi <- emmeans(mod_mi, ~ target * order, level = 0.9)

emm_mi_rel <- contrast(
  emmeans(mod_mi, ~ target), interaction = "trt.vs.ctrl"
)
describe_posterior(emm_mi_rel, ci = 0.9, rope_ci = 1)

emm_mi_seq <- contrast(
  emmeans(mod_mi, ~ target), interaction = "consec"
)
describe_posterior(emm_mi_seq, ci = 0.9, rope_ci = 1)

emm_mi_diffs <- contrast(emm_mi, interaction = "trt.vs.ctrl")
describe_posterior(emm_mi_diffs, ci = 0.9, rope_ci = 1)
