########################################################################################################################## 6. Plotting of all activity characteristics (not the activity curves)   ####
#### Graphical inference of variation in data between species, individuals and within individuals (days)


library(tidyverse)
library(data.table)
library(viridis)


path<- "/Users/pandadamda/rts_activity/"
path<- "E:/Uni_Arbeit/rts_activity/"

## load data
df<- fread(paste0(path,"data/bird_data_storage/activity_characteristics/activity_characteristics_individual.csv"), stringsAsFactors = T)
df$date_f<- as.factor(as.character(df$date_f))


## load metadata
meta<- read.csv(paste0(path,"data/bird_data_storage/tags_overview.csv"), stringsAsFactors = T)
meta<- meta %>% 
  filter(ring_ID %in% levels(df$ring_ID)) %>% 
  select(ring_ID, ID, year, sex, age, brood_patch, weight_g, wing_mm, recapture) %>% 
  distinct()

df<- df %>% 
  left_join(meta, by=c("ring_ID", "ID"))

## change bird names
df$species_en<- gsub("_", " ", df$species_en)



##############################################################################
## data aggregation for plotting (inidvidual level)
df_agg_ind <- df %>% 
  group_by(ring_ID, species_en, sex) %>% 
  summarise(mean_start=mean(steepest_ascend,na.rm = TRUE),
            mean_end=mean(steepest_descend,na.rm = TRUE),
            sd_start=sd(steepest_ascend,na.rm = TRUE),
            sd_end=sd(steepest_descend,na.rm = TRUE),
            mean_act_at_sunrise = mean(act_at_sunrise_rel,na.rm = TRUE),
            sd_act_at_sunrise = sd(act_at_sunrise_rel,na.rm = TRUE),
            mean_act_at_sunset = mean(act_at_sunset_rel,na.rm = TRUE),
            sd_act_at_sunset = sd(act_at_sunset_rel,na.rm = TRUE),
            mean_auc = mean(auc,na.rm = TRUE),
            sd_auc = sd(auc,na.rm = TRUE),
            mean_auc_sunrise = mean(auc_sunrise,na.rm = TRUE),
            sd_auc_sunrise = sd(auc_sunrise,na.rm = TRUE),
            mean_length_act = mean(length_act, na.rm=T),
            sd_length_act = sd(length_act, na.rm=T),
            n=n()) %>% 
  mutate(sex=as.character(sex))
df_agg_ind[is.na(df_agg_ind$sex),]$sex <- "NA"

## data aggregation for plotting (species level)
df_agg_spec <- df %>% 
  select(ring_ID, species_en, act_at_sunrise_rel, act_at_sunset_rel, steepest_ascend, steepest_descend, length_act ) %>% 
  group_by(species_en) %>% 
  summarise(mean_start=mean(steepest_ascend,na.rm = TRUE),
            mean_end=mean(steepest_descend,na.rm = TRUE),
            sd_start=sd(steepest_ascend,na.rm = TRUE),
            sd_end=sd(steepest_descend,na.rm = TRUE),
            mean_sunrise = mean(act_at_sunrise_rel,na.rm = TRUE),
            sd_sunrise = sd(act_at_sunrise_rel,na.rm = TRUE),
            mean_sunset = mean(act_at_sunset_rel,na.rm = TRUE),
            sd_sunset = sd(act_at_sunset_rel,na.rm = TRUE),
            mean_length_act = mean(length_act, na.rm=T),
            sd_length_act = sd(length_act, na.rm=T)
            ) %>% 
  mutate(ring_ID = "spec_mean")

df_plot <- rbind(df_agg_ind, df_agg_spec)
df_plot$species_ID <-  as.numeric(df_plot$species_en)




##############################################################################
# plotting activity start ####

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(steepest_ascend,na.rm = TRUE),
            sd = sd(steepest_ascend,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.5) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_start), x = mean_start, xmin = mean_start-sd_start, xmax = mean_start+sd_start), size=0.6, orientation="y",alpha=0.4, position=position_jitter(width=0, height=0.5)) +
 #scale_color_brewer(palette="Dark2")+
 #scale_color_viridis(discrete = TRUE, option = "C")+
 #scale_colour_scico_d(palette = 'vikO')+
  theme_light() +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) +
  #xlim(0, 1) +
  ylab("Species") +
  xlab("Mean time of activity start \n (centered around sunrise)")

ggsave(filename = paste0(path, "manuscript/figures/" , "activity_start_species" , ".png"),
       width = 8, height = 10)


lm<- glmmTMB(steepest_ascend ~ sex + (1|species_en) + (1|species_en:ring_ID), data=df)
lm<- glmmTMB(steepest_ascend ~ sex *species_en + (1|ring_ID), data=df)
summary(lm)
car::Anova(lm)



##############################################################################
# plotting activity end ####

df %>% 
  group_by(species_en) %>% 
  summarise(order = mean(steepest_ascend,na.rm = TRUE),
            mean = mean(steepest_descend,na.rm = TRUE),
            sd = sd(steepest_descend,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en,order), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd),size=1.5) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_end), x = mean_end, xmin = mean_end-sd_end, xmax = mean_end+sd_end), size=0.6, orientation="y",  alpha=0.4, position=position_jitter(width=0, height=0.5)) +
  #scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  #scale_colour_scico_d(palette = 'vikO')+
  theme_light() +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) +
  #xlim(0, 1) +
  ylab("Species") +
  xlab("Mean time of activity end \n (centered around sunrise)")

ggsave(filename = paste0(path, "manuscript/figures/" , "activity_end_species" , ".png"),
       width = 8, height = 10)
lm<- glmmTMB(steepest_descend ~ sex + (1|species_en) + (1|species_en:ring_ID), data=df)
summary(lm)



##############################################################################
# plotting relative activity at sunrise ####

# Activity at sunrise is divided by total activity (AUC) to allow for comparability across individuals

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(act_at_sunrise_rel,na.rm = TRUE),
            sd = sd(act_at_sunrise_rel,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.5) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_act_at_sunrise), x = mean_act_at_sunrise, xmin = mean_act_at_sunrise-sd_act_at_sunrise,alpha=0.4, xmax = mean_act_at_sunrise+sd_act_at_sunrise),alpha=0.4, size=0.6, orientation="y", position=position_jitter(width=0, height=0.5)) +
  #scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  theme_light() +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) +
  #xlim(-1, 2) +
  ylab("Species") +
  xlab("Mean relative activity at sunrise \n (centered around sunrise)")

ggsave(filename = paste0(path, "manuscript/figures/" , "activity_sunrise_species" , ".png"),
       width = 8, height = 10)

lm<- glmmTMB(act_at_sunrise_rel ~ sex + (1|species_en) + (1|species_en:ring_ID), data=df)
summary(lm)



##############################################################################
# plotting relative activity  at sunset ####
# Activity at sunrise is divided by total activity (AUC) to allow for comparability across individuals

df %>% 
  group_by(species_en) %>% 
  summarise(order = mean(act_at_sunrise_rel,na.rm = TRUE),
            mean = mean(act_at_sunset_rel,na.rm = TRUE),
            sd = sd(act_at_sunset_rel,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, order), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.5) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_act_at_sunset), x = mean_act_at_sunset, xmin = mean_act_at_sunset-sd_act_at_sunset, xmax = mean_act_at_sunset+sd_act_at_sunset),alpha=0.4, size=0.6, orientation="y", position=position_jitter(width=0, height=0.5)) +
 # scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  theme_light() +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) +
  #xlim(-1, 2) +
  ylab("Species") +
  xlab("Mean relative activity at sunset \n (centered around sunrise)")

ggsave(filename = paste0(path, "manuscript/figures/" , "activity_sunset_species" , ".png"),
       width = 8, height = 10)
lm<- glmmTMB(act_at_sunset_rel ~ sex + (1|species_en) + (1|species_en:ring_ID), data=df)
summary(lm)




##############################################################################
# plotting AUC ####

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(auc,na.rm = TRUE),
            sd = sd(auc,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.5) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_auc), x = mean_auc, xmin = mean_auc-sd_auc, xmax = mean_auc+sd_auc), size=0.6, orientation="y",alpha=0.4, position=position_jitter(width=0, height=0.5)) +
 # scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  theme_light() +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) +
  #xlim(-1, 2) +
  ylab("Species") +
  xlab("Mean activity [area under curve]")

ggsave(filename = paste0(path, "manuscript/figures/" , "activity_AUC_species" , ".png"),
       width = 8, height = 10)

lm<- glmmTMB(auc ~ sex + (1|species_en) + (1|species_en:ring_ID), data=df)
summary(lm)


##############################################################################
# plotting AUC around sunrise ####

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(auc_sunrise,na.rm = TRUE),
            sd = sd(auc_sunrise,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.5) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_auc), x = mean_auc_sunrise, xmin = mean_auc_sunrise-sd_auc_sunrise, xmax = mean_auc_sunrise+sd_auc_sunrise), size=0.6, orientation="y",alpha=0.4, position=position_jitter(width=0, height=0.5)) +
  # scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  theme_light() +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) +
  #xlim(-1, 2) +
  ylab("Species") +
  xlab("Mean activity around sunrise [area under curve +- 1 hour]")

ggsave(filename = paste0(path, "manuscript/figures/" , "activity_AUCsunrise_species" , ".png"),
       width = 8, height = 10)

lm<- glmmTMB(auc_sunrise ~ sex  + (1|species_en) + (1|species_en:ring_ID), data=df)
summary(lm)


##############################################################################
# plotting Activity length ####

df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(length_act,na.rm = TRUE),
            sd = sd(length_act,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean), x = mean, group=species_en, color=species_en)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2) +
  #geom_point(data=df_agg_ind, aes(y = species_en,x = mean_start, group=species_en, color=species_en), position=position_jitter(width=0, height=0.5))+
  geom_pointrange(data= df_agg_ind, aes(y = fct_reorder(species_en, mean_length_act), x = mean_length_act, xmin = mean_length_act-sd_length_act, xmax = mean_length_act+sd_length_act), size=0.3, orientation="y",alpha=0.5, position=position_jitter(width=0, height=0.5)) +
  scale_color_brewer(palette="Dark2")+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  theme_light() +
  theme(legend.position = "none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) +
  #xlim(-1, 2) +
  ylab("Species") +
  xlab("Activity length")

ggsave(filename = paste0(path, "manuscript/figures/" , "activity_length_species" , ".png"),
       width = 8, height = 10)

lm<- glmmTMB(length_act ~ sex + (1|species_en) + (1|species_en:ring_ID), data=df)
summary(lm)





################################################################################################
######################## influence of "brood_patch" 

df %>% 
  filter(species_en!="Common_Blackbird") %>% 
  filter(species_en!="Wood_Warbler") %>% 
  group_by(species_en,brood_patch) %>% 
  summarise(mean = mean(act_at_sunset_rel,na.rm = TRUE),
            sd = sd(act_at_sunset_rel,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y=species_en, x = mean, group=brood_patch, color=brood_patch)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("rel. Activity at sunrise")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) 
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/brood_patch/" , "activity_sunrise_brood_patch" , ".png"),
       width = 8, height = 10)
xtabs(~df_meta$species_en+df_meta$brood_patch)


df %>% 
  filter(species_en!="Common_Blackbird") %>% 
  filter(species_en!="Wood_Warbler") %>% 
  group_by(species_en,brood_patch) %>% 
  summarise(mean = mean(auc,na.rm = TRUE),
            sd = sd(auc,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean),  x = mean, group=brood_patch, color=brood_patch)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("Mean Activity (Area under Curve)")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/brood_patch/" , "activity_brood_patch" , ".png"),
       width = 8, height = 10)
xtabs(~df_meta$species_en+df_meta$brood_patch)


df %>% 
  filter(species_en!="Common_Blackbird") %>% 
  filter(species_en!="Wood_Warbler") %>% 
  group_by(species_en,brood_patch) %>% 
  summarise(mean = mean(length_act,na.rm = TRUE),
            sd = sd(length_act,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y = fct_reorder(species_en, mean),  x = mean, group=brood_patch, color=brood_patch)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("Activity length")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/brood_patch/" , "activity_length_brood_patch" , ".png"),
       width = 8, height = 10)


################################################################################################
######################## influence of "sex" 

df %>% 
  group_by(species_en,sex) %>% 
  summarise(mean = mean(act_at_sunset_rel,na.rm = TRUE),
            sd = sd(act_at_sunset_rel,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y=species_en, x = mean, group=sex, color=sex)) +
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("rel. Activity at sunrise")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold")) 
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/sex/" , "activity_sunrise_sex" , ".png"),
       width = 8, height = 10)

df %>% 
  group_by(species_en,sex) %>% 
  summarise(mean = mean(auc,na.rm = TRUE),
            sd = sd(auc,na.rm = TRUE)) %>%   
  ggplot(data=., aes(y=species_en, x = mean, group=sex, color=sex))+
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("Mean Activity (Area under Curve)")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/sex/" , "activity_auc_sex" , ".png"),
       width = 8, height = 10)


df %>% 
  group_by(species_en,sex) %>% 
  summarise(mean = mean(length_act,na.rm = TRUE),
            sd = sd(length_act,na.rm = TRUE)) %>%      
  ggplot(data=., aes(y=species_en, x = mean, group=sex, color=sex))+
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("Activity length")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/sex/" , "activity_length_sex" , ".png"),
       width = 8, height = 10)

df %>% 
  group_by(species_en,sex) %>% 
  summarise(mean = mean(steepest_ascend,na.rm = TRUE),
            sd = sd(steepest_ascend,na.rm = TRUE)) %>%      
  ggplot(data=., aes(y=species_en, x = mean, group=sex, color=sex))+
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("Activity start")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/sex/" , "activity_start_sex" , ".png"),
       width = 8, height = 10)

df %>% 
  group_by(species_en,sex) %>% 
  summarise(mean = mean(steepest_descend,na.rm = TRUE),
            sd = sd(steepest_descend,na.rm = TRUE)) %>%      
  ggplot(data=., aes(y=species_en, x = mean, group=sex, color=sex))+
  #geom_point(size=3) +
  geom_pointrange(aes(xmin = mean-sd, xmax = mean+sd), size=1.2)+
  ggtitle("Activity end")+
  ylab("Species")+
  theme_light() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24,face="bold"))
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/sex/" , "activity_end_sex" , ".png"),
       width = 8, height = 10)



##############################################################################
###### Correlation plots
traits <- df[, c("auc","auc_sunrise", "act_at_sunset_rel","act_at_sunrise_rel","steepest_descend","steepest_ascend", "length_act")]
cor_mat <- cor(traits, use="pairwise.complete.obs")

corrplot::corrplot(cor_mat, method = "color", type = "upper",
         addCoef.col = "black", # show correlation values
         tl.cex = 0.8)


# plotting "mean activity" ~ "Activity length" ####

ggplot(df, aes(x=log(weight_g), y=auc))+
  geom_point(alpha=0.5)+
  theme(legend.position="none")+
  #facet_wrap(~species_en)+
  geom_smooth(method="lm", se=F)
ggsave(filename = paste0(path, "plots/model_output/analysis_activity/" , "weight_vs_AUC_species" , ".png"),
       width = 16, height = 10)

ggplot(df, aes(x=length_act, y=mean_act, color=ring_ID, group=ring_ID))+
  geom_point()+
  theme(legend.position="none")+
  geom_smooth(method="lm", se=F)+
  facet_wrap(~species_en, scales = "free")


### PCA for all traits
traits_pca<- traits %>% filter(!is.na(length_act))
df_pca<- df %>% filter(!is.na(length_act))
p <- prcomp(scale(traits_pca), center = TRUE, scale. = TRUE)
summary(p)

df_pca$PC1 <- p$x[,1]
df_pca$PC2 <- p$x[,2]

m_PC1 <- glmmTMB(PC1 ~ sex + weight_g + age + (1|species_en) + (1|ring_ID), data=df_pca)
m_PC2 <- glmmTMB(PC2 ~ sex + weight_g + age +  (1|species_en) + (1|ring_ID), data=df_pca)
summary(m_PC1)
summary(m_PC2)

