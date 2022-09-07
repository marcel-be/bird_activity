library(tidyverse)
library(lubridate)
library(data.table)
library(DHARMa)
library(glmmTMB)


path<- "J:/rts/rts_activity/"

df<- fread(paste0(path,"bird_data_storage/activity_characteristics/activity_characteristics.csv"), stringsAsFactors = T)
df$date_f<- as.factor(as.character(df$date_f))
str(df)

df_agg<- df %>% 
  group_by(ring_ID,species_en) %>% 
  summarise_each(funs=mean)


#####################################################################################
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
  xlab("Mean probability of activity at sunrise") # problem: how to deal with uncertainty?? Does it make sense to use the CI??

df %>% 
  group_by(ring_ID,species_en) %>% 
  summarise(mean = mean(act_at_sunrise_mean),
            sd = sd(act_at_sunrise_mean)) %>%  
  ggplot(data=., aes(y = fct_reorder(ring_ID, mean), x = mean, group=species_en, color=species_en)) +
  geom_point() +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd)) +
  theme_light() +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"),
         legend.position = c(0.85, 0.3),
         legend.title =  element_text(size=14),
         legend.text = element_text(size=12))+
  scale_colour_discrete("Species")+
  xlim(0, 1) +
  ylab("Bird Individual") +
  xlab("Mean probability of activity at sunrise") 

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_at_sunrise_individuals" , ".png"),
       width = 8, height = 10)



df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(act_at_sunrise_mean),
            sd = sd(act_at_sunrise_mean)) %>%  
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  geom_point(size=2) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1) +
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")+
  xlim(0, 1) +
  ylab("Species") +
  xlab("Mean probability of activity at sunrise") 

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_at_sunrise_species" , ".png"),
       width = 8, height = 10)



mod1<- glmmTMB(act_at_sunrise_mean ~  species_en + (1|ring_ID), data=df)
summary(mod1)
car::Anova(mod1)


mod2<- glmmTMB(act_at_sunrise_mean ~  ring_ID, data=df)
summary(mod2)

AIC(mod1, mod2)


df$E<- residuals(mod1, type="pearson")
df$fit<- fitted(mod1)

ggplot(data=df, aes(y=E, x=fit))+
  geom_point(alpha=0.2)+
  xlab("Fitted values") +
  ylab("Pearson residuals")+
  theme(text = element_text(size=15))+
  geom_hline(yintercept=0)+
  #geom_smooth(se=F)+
  facet_wrap(~species_en)

plot(df$E~df$species_en)


############################################################################################
# activity onset

df %>% 
  filter(ring_ID != "90850007") %>% 
  group_by(ring_ID,species_en) %>% 
  summarise(mean = mean(steepest_ascend),
            sd = sd(steepest_ascend)) %>%  
  ggplot(data=., aes(y = fct_reorder(ring_ID, mean), x = mean, group=species_en, color=species_en)) +
  geom_point() +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd)) +
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.85, 0.3),
        legend.title =  element_text(size=14),
        legend.text = element_text(size=12))+
  scale_colour_discrete("Species")+
  #xlim(0, 1) +
  ylab("Bird Individual") +
  xlab("Mean time of activity onset \n (centered around time of sunrise)") 

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_onset_individuals" , ".png"),
       width = 8, height = 10)

df %>% 
  filter(ring_ID != "90850007") %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(steepest_ascend),
            sd = sd(steepest_ascend)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  geom_point(size=2) +
  geom_point(data=df_agg, aes(y= species_en,x = steepest_ascend), col="black")+
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1) +
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")+
  #xlim(0, 1) +
  ylab("Species") +
  xlab("Mean time of activity onset \n (centered around time of sunrise)")

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_onset_species" , ".png"),
       width = 8, height = 10)


mod3<- glmmTMB(steepest_ascend ~  species_en + (1|ring_ID), data=df)
summary(mod3)
car::Anova(mod3)

df$E<- residuals(mod3, type="pearson")
df$fit<- fitted(mod3)

ggplot(data=df, aes(y=E, x=fit))+
  geom_point(alpha=0.2)+
  geom_hline(yintercept=0)+
  facet_wrap(~species_en)

plot(df$E~df$species_en)



############################################################################################
# activity end

df %>% 
  filter(ring_ID != "90850007") %>% 
  group_by(ring_ID,species_en) %>% 
  summarise(mean = mean(steepest_descend),
            sd = sd(steepest_descend)) %>%  
  ggplot(data=., aes(y = fct_reorder(ring_ID, mean), x = mean, group=species_en, color=species_en)) +
  geom_point() +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd)) +
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.85, 0.3),
        legend.title =  element_text(size=14),
        legend.text = element_text(size=12))+
  scale_colour_discrete("Species")+
  #xlim(0, 1) +
  ylab("Bird Individual") +
  xlab("Mean time of activity end \n (centered around time of sunrise)") 

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_end_individuals" , ".png"),
       width = 8, height = 10)

df %>% 
  filter(ring_ID != "90850007") %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(steepest_descend),
            sd = sd(steepest_descend)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  geom_point(size=2) +
  geom_point(data=df_agg, aes(y= species_en,x = steepest_descend), col="black")+
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1) +
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")+
  #xlim(0, 1) +
  ylab("Species") +
  xlab("Mean time of activity end \n (centered around time of sunrise)")

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_end_species" , ".png"),
       width = 8, height = 10)


mod3<- glmmTMB(steepest_ascend ~  species_en + (1|ring_ID), data=df)
summary(mod3)
car::Anova(mod3)

df$E<- residuals(mod3, type="pearson")
df$fit<- fitted(mod3)

ggplot(data=df, aes(y=E, x=fit))+
  geom_point(alpha=0.2)+
  geom_hline(yintercept=0)+
  facet_wrap(~species_en)

plot(df$E~df$species_en)




########################################################################
#### plotting larissa

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


