################################
### MEP Standardizing Script ###
################################

# Author: Devan Pancura

install.packages("ez")
install.packages("pastecs")
install.packages("lme4")
library(ez)
library(pastecs)
library(lme4)

#Mixed Model ANOVA on Z-Score Transformed Data

with(dat, table(id, condition, target))

dat <- dat %>%
  mutate(condition = recode(condition, ME = "Overt", MI = "Imagery"))

dat <- dat %>%
  mutate(target = as.factor(target)) %>%
  rename(ID = id, Trial = trial, Order = order, Force = target, Modality = condition, Mean = z_ampl) %>%
  summarize(ID, Order, Force, Modality, Mean)


model <- ezANOVA(data = dat, dv = .(Mean), wid = .(ID), between = .(Order, Modality), 
                  within = .(Force), type = 2, detailed = TRUE)
model


ezPlot(data = dat, dv = Mean,  wid = ID, within = .(Force, Modality),
       x = Force, split = Modality) +
      plot_theme +
      scale_color_manual(values = plot_colors) +
      xlab("Force (% MVC)") +
      ylab("Z-Transformed Mean MEP Amplitude")


#Bonferroni Post Hoccing

imagery <- dat %>%
  subset(Modality == "Imagery")

overt <- dat %>%
  subset(Modality == "Overt")

overt80 <- overt %>%
  subset(Force == 80)

overt70 <- overt %>%
  subset(Force == 70)

t.test(overt80$Mean, overt70$Mean, p.adjust.method = "bonferroni")

##Some descriptive stats about trials total

ggplot(dat, aes(x = target)) +
  geom_bar(color = "navyblue", fill = "navyblue", width = 3) +
  labs(x = "Force (% MVC)", y = "Number of Successful Trials") +
  facet_wrap(~condition, nrow = 2) +
  geom_text(aes(label = ..count..), stat = "count", size = 8, vjust = 1.5, colour = "white") +
  plot_theme

##Testing differences between groups

t.test(kviq ~ order, data = participant_dat)

##Post-hoc power analysis

install.packages("pwr")
library(pwr)

pwr.anova.test(k = 2, n = 23, power = 0.8)
