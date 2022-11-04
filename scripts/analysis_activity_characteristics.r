library(tidyverse)
library(lubridate)
library(data.table)
library(DHARMa)
library(glmmTMB)
library(viridis)
library(scico)


path<- "J:/rts/rts_activity/"
path<- "G:/rts_activity/"
path<- "F:/Uni_Arbeit/rts_activity/"

# Data preperation ####

df<- fread(paste0(path,"bird_data_storage/activity_characteristics/activity_characteristics_individual.csv"), stringsAsFactors = T)
df$date_f<- as.factor(as.character(df$date_f))
str(df)


df_agg_ind <- df %>% 
  group_by(ring_ID,species_en) %>% 
  summarise(mean_start=mean(steepest_ascend,na.rm = TRUE),
            mean_end=mean(steepest_descend,na.rm = TRUE),
            sd_start=sd(steepest_ascend,na.rm = TRUE),
            sd_end=sd(steepest_descend,na.rm = TRUE),
            mean_act_at_sunrise = mean(act_at_sunrise_rel,na.rm = TRUE),
            sd_act_at_sunrise = sd(act_at_sunrise_rel,na.rm = TRUE),
            mean_auc = mean(auc,na.rm = TRUE),
            sd_auc = sd(auc,na.rm = TRUE),
            n=n())

df_agg_spec <- df %>% 
  select(ring_ID, species_en, act_at_sunrise_rel, act_at_sunset_rel, steepest_ascend, steepest_descend ) %>% 
  group_by(species_en) %>% 
  summarise(mean_start=mean(steepest_ascend,na.rm = TRUE),
            mean_end=mean(steepest_descend,na.rm = TRUE),
            sd_start=sd(steepest_ascend,na.rm = TRUE),
            sd_end=sd(steepest_descend,na.rm = TRUE),
            mean_sunrise = mean(act_at_sunrise_rel,na.rm = TRUE),
            sd_sunrise = sd(act_at_sunrise_rel,na.rm = TRUE),
            mean_sunset = mean(act_at_sunset_rel,na.rm = TRUE),
            sd_sunset = sd(act_at_sunset_rel,na.rm = TRUE),
            ) %>% 
  mutate(ring_ID = "spec_mean")

df_plot <- rbind(df_agg_ind, df_agg_spec)
df_plot$species_ID <-  as.numeric(df_plot$species_en)

##############################################################################
# plotting some diagnistics
ggplot(df_agg_ind, aes(y=sd_start, x=n, group=species_en, color=species_en))+
  geom_point()+
  facet_wrap(~species_en)
summary(lm(sd_start~n*species_en, data=df_agg_ind))

ggplot(df_agg_ind, aes(y=sd_end, x=n, group=species_en, color=species_en))+
  geom_point()+
  facet_wrap(~species_en)
summary(lm(sd_end~n*species_en, data=df_agg_ind))




cor(df_agg_ind$sd_start, df_agg_ind$n)


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
  summarise(mean = mean(steepest_ascend,na.rm = TRUE),
            sd = sd(steepest_ascend,na.rm = TRUE)) %>%  
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
  summarise(mean = mean(steepest_ascend,na.rm = TRUE),
            sd = sd(steepest_ascend,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_start), x = mean_start, xmin = mean_start-sd_start, xmax = mean_start+sd_start), size=0.3, orientation="y", position=position_jitter(width=0, height=0.5)) +
 #scale_color_brewer(palette="Dark2")+
 #scale_color_viridis(discrete = TRUE, option = "C")+
 #scale_colour_scico_d(palette = 'vikO')+
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

E<- residuals(mod3, type="pearson")
fit<- fitted(mod3)

ggplot(data=data.frame(cbind(E,fit)), aes(y=E, x=fit))+
  geom_point(alpha=0.2)+
  geom_hline(yintercept=0)+
  facet_wrap(~species_en)

plot(df$E~df$species_en)



##############################################################################
# plotting activity end ####

df %>% 
  group_by(ring_ID,species_en) %>% 
  summarise(mean = mean(steepest_descend),na.rm = TRUE,
            sd = sd(steepest_descend,na.rm = TRUE)) %>%  
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
  summarise(mean = mean(steepest_descend,na.rm = TRUE),
            sd = sd(steepest_descend,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_end), x = mean_end, xmin = mean_end-sd_end, xmax = mean_end+sd_end), size=0.3, orientation="y", position=position_jitter(width=0, height=0.5)) +
  #scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  #scale_colour_scico_d(palette = 'vikO')+
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")+
  #xlim(0, 1) +
  ylab("Species") +
  xlab("Mean time of activity end \n (centered around time of sunrise)")

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_end_species" , ".png"),
       width = 8, height = 10)


mod3<- glmmTMB(steepest_descend ~  species_en + (1|ring_ID), data=df)
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
# plotting relative activity  at sunrise ####

# Activity at sunrise is divided by total activity (AUC) to allow for comparability across individuals

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(act_at_sunrise_rel,na.rm = TRUE),
            sd = sd(act_at_sunrise_rel,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_act_at_sunrise), x = mean_act_at_sunrise, xmin = mean_act_at_sunrise-sd_act_at_sunrise, xmax = mean_act_at_sunrise+sd_act_at_sunrise), size=0.3, orientation="y", position=position_jitter(width=0, height=0.5)) +
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
# plotting relative activity  at sunset ####

# Activity at sunrise is divided by total activity (AUC) to allow for comparability across individuals

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(act_at_sunrise_rel,na.rm = TRUE),
            sd = sd(act_at_sunrise_rel,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_act_at_sunrise), x = mean_act_at_sunrise, xmin = mean_act_at_sunrise-sd_act_at_sunrise, xmax = mean_act_at_sunrise+sd_act_at_sunrise), size=0.3, orientation="y", position=position_jitter(width=0, height=0.5)) +
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
# plotting AUC ####

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(auc,na.rm = TRUE),
            sd = sd(auc,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, ring_ID_num), x = mean_auc, xmin = mean_auc-sd_auc, xmax = mean_auc+sd_auc), size=0.3, orientation="y", position=position_jitter(width=0, height=0.5)) +
  scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  theme_light() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none")+
  #xlim(-1, 2) +
  ylab("Species") +
  xlab("Mean activity [area under curve]")

ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "activity_AUC_species" , ".png"),
       width = 8, height = 10)







##############################################################################
# Variance Component Analysis ####

## ANOVA
act_start <- aov(df$steepest_ascend ~ df$species_en / df$ring_ID)
summary(act_start)
car::Anova(act_start)

summary(glmmTMB(steepest_ascend ~ species_en + (1|ring_ID), data=df))

## VCA Package 

library(VCA)
df_vca<- as.data.frame(df_vca)

vca <-  fitVCA(steepest_ascend ~ species_en/ring_ID, df_vca , method="anova")
vca
inf <- VCAinference(vca, VarVC=TRUE)
inf

plotRandVar(vca, term="species_en:ring_ID", mode="student", pick=T) 
abline(h=c(-3, 3), lty=2, col="red", lwd=2)
mtext(side=4, at=c(-3, 3), col="red", line=.25, las=1, text=c(-3, 3))
############ CHECK



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

