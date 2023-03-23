################################
### MEP Standardizing Script ###
################################

# Author: Devan Pancura


### Run after processing script ###

#Standardizing to MI10 Condition

target10_MI1 <- trimmed_MI_1 %>%
  subset(target == 10)

target10_MI1 <- target10_MI1 %>%  
  summarize(id,
            MI10 = mean_amp)

trimmed_MI_1 <- trimmed_MI_1 %>%
  group_by(id) %>%
  right_join(target10_MI1, by = c("id")) %>%
  select(c(id, target, mean_amp, MI10))

target10_MI2 <- trimmed_MI_2 %>%
  subset(target == 10)

target10_MI2 <- target10_MI2 %>%  
  summarize(id, 
            MI10 = mean_amp)

trimmed_MI_2 <- trimmed_MI_2 %>%
  group_by(id) %>%
  right_join(target10_MI2, by = c("id")) %>%
  select(c(id, target, mean_amp, MI10))

target10_ME1 <- trimmed_ME_1 %>%
  subset(target == 10)

target10_ME1 <- target10_ME1 %>%  
  summarize(id,
            ME10 = mean_amp)

trimmed_ME_1 <- trimmed_ME_1 %>%
  group_by(id) %>%
  right_join(target10_ME1, by = c("id")) %>%
  select(c(id, target, mean_amp, ME10))

target10_ME2 <- trimmed_ME_2 %>%
  subset(target == 10)

target10_ME2 <- target10_ME2 %>%  
  summarize(id,
            ME10 = mean_amp)

trimmed_ME_2 <- trimmed_ME_2 %>%
  group_by(id) %>%
  right_join(target10_ME2, by = c("id")) %>%
  select(c(id, target, mean_amp, ME10)) %>%
  na.omit()

testME1 <- trimmed_ME_1 %>%
  summarize(id, target, mean_amp, ME10,
            standardized_value = (mean_amp - ME10)/ME10,
            condition = "ME",
            order = 1)

testME2 <- trimmed_ME_2 %>%
  summarize(id, target, mean_amp, ME10,
            standardized_value = (mean_amp - ME10)/ME10,
            condition = "ME",
            order = 2)

testMI1 <- trimmed_MI_1 %>%
  summarize(id, target, mean_amp, MI10,
            standardized_value = (mean_amp - MI10)/MI10,
            condition = "MI",
            order = 1)

testMI2 <- trimmed_MI_2 %>%
  summarize(id, target, mean_amp, MI10,
            standardized_value = (mean_amp - MI10)/MI10,
            condition = "MI",
            order = 2)

ind_standard_data <- full_join(testME1, testME2)
ind_standard_data <- full_join(ind_standard_data, testMI1)
ind_standard_data <- full_join(ind_standard_data, testMI2)

testME1 <- testME1 %>%
  group_by(target) %>%
  summarize(mean = mean(standardized_value),
            sd = sd(standardized_value))

testME2 <- testME2 %>%
  group_by(target) %>%
  summarize(mean = mean(standardized_value),
            sd = sd(standardized_value))

testMI1 <- testMI1 %>%
  group_by(target) %>%
  summarize(mean = mean(standardized_value),
            sd = sd(standardized_value))

testMI2 <- testMI2 %>%
  group_by(target) %>%
  summarize(mean = mean(standardized_value),
            sd = sd(standardized_value))

standard_data <- testMI2 %>%
  full_join(testMI1) %>%
  full_join(testME2) %>%
  full_join(testME1) %>%
  mutate(order = as.factor(order),
         condition = as.factor(condition))

ggplot(standard_data, aes(x = target, y = mean, col = order, group = order)) +
  geom_point(position = pd, size = 4) +
  geom_pointrange(aes(
    ymin = mean - sd, 
    ymax = mean + sd), position = pd, size = 1) +
  scale_color_manual(values = plot_colors, labels = c("1", "2")) +
  xlab("Effector Load") +
  ylab("Standardized Means") +
  facet_wrap(~condition, nrow = 1) +
  plot_theme +
  labs(group = "Order", color = "Order")


### Cannot run an ANOVA - data does not meet the assumptions ###


#Z-Score Transformation

tMI1 <- trimmed_meps %>%
  subset(order == 1 & trial > 80) #MI first

tMI2 <- trimmed_meps %>% #MI second
  subset(order == 2 & trial < 81)

tME2 <- trimmed_meps %>%
  subset(order == 1 & trial < 81) #ME second

tME1 <- trimmed_meps %>%
  subset(order == 2 & trial > 80) #ME first

tMI1 <- tMI1 %>%
  group_by(id, trial) %>%
  summarize(order, ampl, target, condition = "MI")

tMI2 <- tMI2 %>%
  group_by(id, trial) %>%
  summarize(order, ampl, target, condition = "MI")

tME1 <- tME1 %>%
  group_by(id, trial) %>%
  summarize(order, ampl, target, condition = "ME")

tME2 <- tME2 %>%
  group_by(id, trial) %>%
  summarize(order, ampl, target, condition = "ME")

tMI = full_join(tMI1, tMI2) %>%
  ungroup() %>%
  distinct() %>%
  group_by(id) %>%
  mutate(z_ampl = (ampl - mean(ampl))/sd(ampl))

tME = full_join(tME1, tME2) %>%
  ungroup() %>%
  distinct() %>%
  group_by(id) %>%
  mutate(z_ampl = (ampl - mean(ampl))/sd(ampl))

dat <- full_join(tMI, tME) %>%
  ungroup() %>%
  mutate(order = as.factor(order),
         condition = as.factor(condition),
         id = as.factor(id))

z_transform <- dat %>%
  group_by(target, order, condition) %>%
  summarize(mean = mean(z_ampl),
            sd = sd(z_ampl))

plot_colors <- c("navyblue", "darkorange3")

z_transform <- z_transform %>%
  mutate(condition = recode(condition, ME = "Overt", MI = "Imagery"))

plot_colors <- c("black", "black")
plot_shapes <- c(1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 
                 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19, 1, 19)

ggplot(z_transform, aes(x = target, y = mean, col = order, group = order)) +
  geom_pointrange(aes(
    ymin = mean - sd, 
    ymax = mean + sd), size = 1, position = position_dodge(4)) +
  scale_color_manual(values = plot_colors, labels = c("1","2")) +
  xlab("Force (% MVC)") +
  ylab("Z-Transformed Mean MEP Amplitude") +
  facet_wrap(~condition, nrow = 1) +
  plot_theme +
  labs(group = "Order", color = "Order")


