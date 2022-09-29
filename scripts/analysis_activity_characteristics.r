library(tidyverse)
library(lubridate)
library(data.table)
library(DHARMa)
library(glmmTMB)
library(viridis)


path<- "J:/rts/rts_activity/"

# Data preperation ####

df<- fread(paste0(path,"bird_data_storage/activity_characteristics/activity_characteristics.csv"), stringsAsFactors = T)
df$date_f<- as.factor(as.character(df$date_f))
str(df)

df<- df %>%   
  filter(ring_ID != "90850007")


df$ring_ID_num <- NA
df_vca<- data.frame()
num_vec<- seq(1:15)
for(i in 1:nlevels(df$species_en)){
  
  df_sub<- df %>% 
    filter(species_en==levels(species_en)[i]) %>% 
    droplevels()
  
  for(j in 1:nlevels(df_sub$ring_ID)){
    df_num<- df %>% 
      filter(species_en==levels(df$species_en)[i]) %>% 
      droplevels() %>% 
      filter(ring_ID==levels(df_sub$ring_ID)[j]) %>% 
      mutate(ring_ID_num = num_vec[j])
    df_vca<- rbind(df_vca, df_num)
  }
}


df_agg_ind <- df_vca %>% 
  select(ring_ID, species_en, act_at_sunrise_mean, steepest_ascend, steepest_descend , ring_ID_num) %>% 
  group_by(ring_ID,species_en) %>% 
  summarise(mean_start=mean(steepest_ascend),
            mean_end=mean(steepest_descend),
            sd_start=sd(steepest_ascend),
            sd_end=sd(steepest_descend),
            ring_ID_num= mean(ring_ID_num),
            mean_act_at_sunrise = mean(act_at_sunrise_mean),
            sd_act_at_sunrise = sd(act_at_sunrise_mean))

df_agg_spec <- df %>% 
  select(ring_ID, species_en, act_at_sunrise_mean, steepest_ascend, steepest_descend ) %>% 
  group_by(species_en) %>% 
  summarise(mean_start=mean(steepest_ascend),
            mean_end=mean(steepest_descend),
            sd_start=sd(steepest_ascend),
            sd_end=sd(steepest_descend)) %>% 
  mutate(ring_ID = "spec_mean")

df_plot <- rbind(df_agg_ind, df_agg_spec)
df_plot$species_ID <-  as.numeric(df_plot$species_en)




##############################################################################
# plotting activity onset ####

df_plot %>%
  ggplot(data=., aes(y = fct_reorder(ring_ID, species_ID), x = mean_start, group=species_en, color=species_en)) +
  geom_point() +
  geom_pointrange(aes(xmin = mean_start-sd_start, xmax = mean_start+sd_start)) +
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
  #facet_wrap(~species_en)


df %>% 
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
  group_by(species_en) %>% 
  summarise(mean = mean(steepest_ascend),
            sd = sd(steepest_ascend)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, ring_ID_num), x = mean_start, xmin = mean_start-sd_start, xmax = mean_start+sd_start), size=0.3, orientation="y", position=position_jitter(width=0, height=0.5)) +
  scale_color_brewer(palette="Dark2")+
 #scale_color_viridis(discrete = TRUE, option = "C")+
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



##############################################################################
# plotting activity end ####

df %>% 
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




##############################################################################
# plotting sunrise ####

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(act_at_sunrise_mean),
            sd = sd(act_at_sunrise_mean)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, ring_ID_num), x = mean_act_at_sunrise, xmin = mean_act_at_sunrise-sd_act_at_sunrise, xmax = mean_act_at_sunrise+sd_act_at_sunrise), size=0.3, orientation="y", position=position_jitter(width=0, height=0.5)) +
  scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")+
  #xlim(-1, 2) +
  ylab("Species") +
  xlab("Mean activity at sunrise \n (centered around time of sunrise)")

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_sunrise_species" , ".png"),
       width = 8, height = 10)

##############################################################################
# plotting larissa ####

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



##############################################################################
# Variance Component Analysis ####

## ANOVA
nest <- aov(df$steepest_ascend ~ df$species_en +  df$species_en:df$ring_ID)
summary(nest)
car::Anova(nest)
car::Anova(lm(steepest_ascend ~ species_en +  species_en:ring_ID, df), type='3')
#https://www.flutterbys.com.au/stats/tut/tut9.2a.html

## VCA Package 

library(VCA)
df<- as.data.frame(df)

vca <-  fitVCA(steepest_ascend ~ species_en/ring_ID, df)
vca

inf <- VCAinference(vca, VarVC=TRUE)
inf
############ CHECK


## Factor "species_en" and "ring_ID" as numbers (as in package-example)
# makes no difference!


# for species

#df_vca$ring_ID_num<- as.factor(as.character(df_vca$ring_ID_num))
df_vca <- as.data.frame(df_vca)

vca2 <-  fitVCA(steepest_ascend ~ species_en_num/ring_ID_num, df_vca, method = c("anova"))
vca2