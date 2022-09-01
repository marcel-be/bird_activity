library(dplyr)
library(ggplot2)
library(DHARMa)
library(glmmTMB)
install.packages("glmmTMB")

path<- "J:/rts/rts_activity/"


df<- fread(paste0(path,"bird_data_storage/activity_characteristics/activity_characteristics.csv"))

boxplot(df$act_at_sunrise_mean ~ df$species_en)
boxplot(df$act_at_sunrise_mean ~ df$ring_ID)


########################################################################
#### plotting larissa

# sunrise
df %>% 
  group_by(ring_ID,species_en) %>% 
  summarise_each(funs(mean)) %>%  
  ggplot(aes(y = fct_reorder(ring_ID, act_at_sunrise_mean), x = act_at_sunrise_mean, group=species_en, color=species_en)) +
  geom_point() +
  geom_pointrange(aes(xmin = act_at_sunrise_lowerCI, xmax = act_at_sunrise_upperCI)) +
  theme_bw() +
  xlim(0, 1) +
  ylab("ID") +
  xlab("Mean probability of activity at sunrise")


# peak activity
data_new_indiv %>% 
  group_by(time_to_rise_std) %>% 
  summarise_each(funs(mean)) %>% 
  ggplot(aes(x = time_to_rise_std, y = mu))+
  geom_ribbon(aes(ymin = ci_lower ,
                  ymax = ci_upper), 
              fill = "grey", color = "grey") +
  geom_line(size = .8) + 
  geom_point(data = max_activity, 
             aes(x = time_to_rise_std, y = mu),
             col = "red", size = 3) +
  geom_vline(xintercept = max_activity$time_to_rise_std, linetype = "dotted", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_bw(14) +
  xlab("Time since sunrise") + 
  ylab("Activity probability \n") + 
  ylim(0, 1) +
  ggtitle("150007_0_10_1")

data_new_indiv_NEW %>% 
  # group_by(time_to_rise_std) %>% 
  # summarise_each(funs(mean)) %>% 
  ggplot(aes(x = time_to_rise_std, y = mu)) +
  geom_ribbon(aes(ymin = ci_lower ,
                  ymax = ci_upper), 
              fill = "grey", color = "grey") +
  geom_line(aes(group = 1, color = act_over_05), size = .8) + 
  geom_point(data = max_activity, 
             aes(x = time_to_rise_std, y = mu),
             col = "red", size = 3) +
  geom_vline(xintercept = max_activity$time_to_rise_std, linetype = "dotted", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_bw(14) +
  xlab("Time since sunrise") + 
  ylab("Activity probability \n") + 
  ylim(0, 1) +
  ggtitle("150007_0_10_1")


