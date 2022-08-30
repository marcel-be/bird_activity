library(tidyverse);library(lattice); library(lubridate); library(hms)
library(lme4); library(viridis); library(effects); library(performance)
library(insight); library(ggthemes); library(stringr); library(data.table)
library(damr); library(behavr); library(ggetho); library(zeitgebr)
library(sleepr)

# Installing suite of rethomics packages
#devtools::install_github("rethomics/behavr")
#devtools::install_github("rethomics/ggetho")
#devtools::install_github("rethomics/damr")
#devtools::install_github("rethomics/scopr")
#devtools::install_github("rethomics/sleepr")
#devtools::install_github("rethomics/zeitgebr")

rm(list=ls())
path<- "J:/rts/rts_activity/"

df <- fread(paste0(path, "bird_data_storage/tags_overview.csv"), key = "ring_ID")

################################################################################
################################################################################
####   1.  Overview Plots of all tags ####

## Overview: Number of species/individuals; Duration of tags; Coverage per week

plot_list<- list()

(p1.1<- df %>% 
    count(species_en, year) %>%
    as.data.frame(.) %>% 
    mutate(year=as.factor(year)) %>% 
    mutate(species_en = fct_reorder(species_en, n)) %>% 
    ggplot(. ,aes(y=species_en, x=n, fill=year,label = round(n,digits=1))) +
    geom_bar(stat="identity")+
    #geom_text(aes(y=species_en, label=year), vjust=0, 
    #         color="white", size=3.5)+
    scale_fill_brewer(palette="Paired")+
    ggtitle("Count of tagged birds \n(including repeated captures of Individuals)")+
    geom_text(size = 5, position = position_stack(vjust = 1.025))+
    xlab("Count")+
    ylab("Species")+
    theme_minimal()+
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"),
          legend.title=element_blank(),
          legend.position = c(0.9, 0.1),
          legend.text = element_text(size=15))
)

(p1.2<- df %>%
    filter(!duplicated(.[["ring_ID"]])) %>% 
    count(species_en, year) %>%
    as.data.frame(.) %>% 
    mutate(year=as.factor(year)) %>% 
    mutate(species_en = fct_reorder(species_en, n)) %>% 
    ggplot(. ,aes(y=species_en, x=n, fill=year,label = round(n,digits=1))) +
    geom_bar(stat="identity")+
    #geom_text(aes(y=species_en, label=year), vjust=0, 
    #         color="white", size=3.5)+
    scale_fill_brewer(palette="Paired")+
    ggtitle("Count of tagged Individuals \n(some individuals were tagged two times)")+
    geom_text(size = 5, position = position_stack(vjust = 1.025))+
    xlab("Count")+
    ylab("Species")+
    theme_minimal()+
    theme(text = element_text(size=15),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"),
          legend.title=element_blank(),
          legend.position = c(0.9, 0.1),
          legend.text = element_text(size=15))
)

(p2<- df %>%
  count(family, year) %>%
  as.data.frame(.) %>% 
  mutate(year=as.factor(year)) %>% 
  mutate(family = fct_reorder(family, n)) %>% 
  ggplot(. ,aes(y=family, x=n, fill=year,label = round(n,digits=1))) +
  geom_bar(stat="identity")+
  #geom_text(aes(y=species_en, label=year), vjust=0, 
  #         color="white", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  geom_text(size = 5, position = position_stack(vjust = 1.025))+
  xlab("Count")+
  ylab("Family")+
  theme_minimal()+
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"),
        legend.title=element_blank(),
        legend.position = c(0.9, 0.1),
        legend.text = element_text(size=15))
)
(p3<- df %>% 
  group_by(species_en) %>% 
  summarise(mean = mean(time_total),
            sd = sd(time_total)) %>% 
  mutate(species_en = fct_reorder(species_en, mean)) %>% 
  ggplot(., aes(x=mean, y=species_en, label = round(mean,digits=1)))+
  geom_bar(stat="identity", col="white", fill="steelblue")+
  geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd), width=.2,
                position=position_dodge(.9), col="grey")+
  geom_text(size = 3, position = position_stack(vjust = 1.025))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Duration of transmitted signal [days]") + 
  ylab("Species")
)
(p4<- df %>% 
    filter(year==2019) %>%
       mutate(ID = fct_reorder(ID, time_total)) %>% 
       ggplot(. , aes(x=time_total, y=ID, label = round(time_total,digits=1)))+
       geom_bar(stat="identity", col="white", fill="steelblue")+
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA))+
       geom_text(size = 3)+
       theme_bw()+
       theme(text = element_text(size=15),
             panel.border = element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
             axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"))+
       xlab("Duration of transmitted signal [days]") + 
       ylab("Tags")+
      ggtitle("Tag duration 2019")
    )
(p5<- df %>% 
    filter(year==2020) %>%
    mutate(ID = fct_reorder(ID, time_total)) %>% 
    ggplot(. , aes(x=time_total, y=ID, label = round(time_total,digits=1)))+
    geom_bar(stat="identity", col="white", fill="steelblue")+
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA))+
    geom_text(size = 3)+
    theme_bw()+
    theme(text = element_text(size=15),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"))+
    xlab("Duration of transmitted signal [days]") + 
    ylab("Tags")+
    ggtitle("Tag duration 2020")
)
(p6<- df %>% 
    filter(year==2021) %>%
    mutate(ID = fct_reorder(ID, time_total)) %>% 
    ggplot(. , aes(x=time_total, y=ID, label = round(time_total,digits=1)))+
    geom_bar(stat="identity", col="white", fill="steelblue")+
    scale_x_continuous(expand = c(0, 0), limits = c(0, NA))+
    geom_text(size = 3)+
    theme_bw()+
    theme(text = element_text(size=15),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"))+
    xlab("Duration of transmitted signal [days]") + 
    ylab("Tags")+
    ggtitle("Tag duration 2021")
)

plot_list[[1]]<- p1.1
plot_list[[2]]<- p1.2
plot_list[[3]]<- p2
plot_list[[4]]<- p3
plot_list[[5]]<- p4
plot_list[[6]]<- p5
plot_list[[7]]<- p6
plot_list[[8]]<- pweek # see further down in script!!!


ggsave(filename = paste0(path, "plots/data_overview/" , "summary" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)


################################################################################
################################################################################
################################################################################
####  2.  Descriptive Analysis of raw data   ####

## Activity Plots for daily rhythm (activity ~ time of day) on different resolutions (per species, per tag, weekly/daily base)

df_1min_meta<- fread(paste0(path, "bird_data_storage/tags_1min_withmeta.csv"), key = "ID")

df_1min_meta$hour <- hour(df_1min_meta$timestamp) # Hour of the day
df_1min_meta$minute <- minute(df_1min_meta$timestamp)
df_1min_meta$month_f <- factor(month(df_1min_meta$date))
df_1min_meta$year_f <- factor(df_1min_meta$year)
df_1min_meta$date_f <- as.factor(as.character(df_1min_meta$date_CET))
df_1min_meta$week <- week(df_1min_meta$timestamp)
df_1min_meta$time_of_day <- as.numeric(as_hms(df_1min_meta$timestamp))/3600 # numeric time of day

########################################################
### 2.1. Rough Overview across all years

## Plot Species means by hour of the day (all years)
df_1min_meta %>% 
  group_by(species_en,hour) %>% 
  summarise(mean_activity=mean(activity, na.rm=T)) %>%
  ggplot(., aes(x=hour, y=mean_activity, group=species_en, fill=species_en, color=species_en)) +
  geom_point(size=3, alpha=.9) + 
  facet_wrap(~species_en) + 
  ylim(0,1) +
  ylab("% Activity") + 
  xlab("hour") + 
  scale_color_viridis(discrete = T) +
  theme_bw(22) + 
  theme(legend.position = "none")


## 1h-plots for all tags for each year (one plot per tag). 
(df_1min_meta %>% 
    filter(year==2021) %>%    # change year if wanted
    #coverage > .8) %>% 
    group_by(species_en,hour,ID) %>% 
    summarise(mean_activity=mean(activity)) %>%
    ggplot(., aes(x=hour, y=mean_activity, group=species_en, fill=species_en, color=species_en)) +
    geom_point(size=1, alpha=.9) + facet_wrap(~ID) + ylim(0,1) +
    ylab("% Activity") + xlab("hour") + scale_color_viridis(discrete = T) +
    theme_bw(14) + 
    theme(legend.position = "none",
          panel.spacing.y=unit(0, "lines"),
          strip.text.y = element_text(size = 7, angle = 0, hjust = 0))) %>% 
  ggsave(filename = paste0(path, "plots/tag_check/" , "2021_all_tags_1h" , ".png"), . , width = 15, height = 9)


#### create data with 10min-intervals for better overview (tag-check-plots)
  df_10min <- df_1min_meta %>% 
  mutate(species_en = as.factor(species_en),
         interval   = as_hms(floor_date(timestamp_CET, unit="10minutes"))) %>% 
  group_by(ID,ring_ID,species_en,year_f,month,week,date_f,ydate,hour,interval) %>% 
  summarise(n_intervals=length(activity),
            n_active=length(activity[activity==1]),
            n_passive=length(activity[activity==0]),
            time_of_day=mean(time_of_day),
            time_to_rise_std=mean(time_to_rise_std),
            ring_ID=as.factor(ring_ID),
            ID=as.factor(ID),
            coverage_daily=as.factor(mean(coverage_daily))) 


# plott for weekly coverage of data:

pweek<- df_10min %>% 
    ggplot(., aes(y = n_active/n_intervals, x = ydate,
                  group = species_en, color = year_f)) +
    geom_point(alpha = .25, size = 2) + 
    # geom_contour_filled() + 
    # geom_density_2d(colour = "black", alpha = .5) +
   # scale_color_viridis(discrete = T) +
    facet_wrap(~ species_en) + theme_bw(14) +
    ggtitle("coverage of data by week")+
    xlab("Day of the Year")+
    ylab("Activity")
ggsave(filename = paste0(path, "plots/data_overview/" , "coverage_week" , ".png"),
       width = 15, height = 9)


########################################################
#### 2.2. Plots for verifying all tags:

## One Curve for each species per week
plot_list<- list()
for(i in 1:nlevels(df_10min$species_en)){
p <- df_10min %>% 
  filter(species_en == levels(species_en)[i]) %>% 
  ggplot(., aes(y = n_active/n_intervals, x = time_of_day, 
                group = ID, color  =  ID)) +
  geom_point(alpha = .15) + 
  geom_smooth(method = 'gam', 
              formula = y ~ s(x, bs = "tp", k = 9)) +
  #scale_color_viridis(discrete = T) +
  scale_fill_brewer(palette="Set1")+
  facet_wrap(~ week) + 
  theme_bw()+
  ggtitle(levels(df_10min$species_en)[i])

plot_list[[i]] <-p
}

ggsave(filename = paste0(path, "plots/tag_check/" , "species_by_week" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)


## One Plot for each tag AND week:
plot_list<- list()
for(i in 1:nlevels(df_10min$ID)){
  p <- df_10min %>% 
    filter(ID == levels(ID)[i]) %>% 
    ggplot(., aes(y = n_active/n_intervals, x = time_of_day)) +
    geom_point(alpha = .20) + 
    geom_smooth(method = 'gam', 
                formula = y ~ s(x, bs = "tp", k = 15)) +
    #scale_color_viridis(discrete = T) + 
    facet_wrap(~ week) +
    theme_bw()+
    ggtitle(levels(df_10min$ID)[i])
  
  plot_list[[i]] <-p
}

ggsave(filename = paste0(path, "plots/tag_check/" , "tags_by_week" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)

## for each tag AND day

plot_list<- list()
for(i in 1:nlevels(df_10min$ID)){
  df <- df_10min %>%  
    filter(ID == levels(ID)[i]) 
  
  name<- paste0(as.character(df$species_en[1]),
                " | brood patch=", as.character(df$brood_patch[1]),
                " | year=", as.character(df$year_f[1]),
                " | sex=", as.character(df$sex[1]))
  
  p<-  ggplot(df, aes(y = n_active/n_intervals, x = time_to_rise_std)) +
    geom_point(alpha = .20) + 
    geom_smooth(method = 'gam', 
                formula = y ~ s(x, bs = "tp", k = 15)) +
    #scale_color_viridis(discrete = T) + 
    facet_wrap(~ date_f) +
    #facet_wrap(~ ydate + coverage_daily, labeller = label_both) +
    theme_bw()+
    ggtitle(paste0(levels(df_10min$ID)[i]," ", name))
  
  plot_list[[i]] <-p
}

ggsave(filename = paste0(path, "plots/tag_check/" , "tags_by_day" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)


# all tags that look odd (after first visual inspection of plots above):
df_10min2<- df_10min %>% 
  filter(ID=="150037_20" |ID=="150114_1_10_2" |ID=="150156_4_30_1" |ID=="210513_150087_40" |ID=="210518_150077_40" |ID=="210610_150099_20" |ID=="210617_150007_10" | ID=="210701_150039_40")# list of fishy tags (from visual a/p-inspection)
df_10min$ID<- as.factor(df_10min$ID)
plot_list<- list()
df_10min2$ID<- droplevels(df_10min2$ID)



################################################################################
################################################################################
####  3. Rethomics Plots ####

# Load and prepare Data 
df_1min <- fread(paste0(path, "bird_data_storage/tags_1min_nometa.csv"), key = "ID")
#df_1min$t <- df_1min$t*60 # As seconds
metadata <- fread(paste0(path, "bird_data_storage/tags_overview.csv"), key = "ID")
metadata$species<- metadata$species_de


#df_1min <- data.table(df_1min, key = "ID", stringsAsFactors = T)
#metadata <- data.table(metadata, key = "ID", stringsAsFactors = T)

dt_raw <- behavr(df_1min,metadata)

summary(dt_raw)
summary(dt_raw, detailed = TRUE)

# to remove animals with long inactive periods (="dead") (if wanted; not used in further analysis)
dt_curated <- curate_dead_animals(dt_raw, moving_var = activity)
nrow(dt_curated)/nrow(dt_raw)
summary(dt_curated)
#dt <- dt_curated # first check for issued datasets and remove them 

dt<- dt_raw
dt[meta = TRUE]
summary(dt)
summary(dt, detailed = TRUE)


################################################################################
####  Plot data with ggetho (rethomics) ####

## in case, select by year:
plot_list<- list()
(p2019<- dt %>% 
    filter(year==2019) %>% 
    ggetho(., aes(x=t, y=ID, z=activity)) + stat_tile_etho() + theme_bw(14) + ggtitle("2019"))
(p2020<- dt %>% 
    filter(year==2020) %>% 
    ggetho(., aes(x=t, y=ID, z=activity)) + stat_tile_etho() + theme_bw(14) + ggtitle("2020"))
(p2021<- dt %>% 
    filter(year==2021) %>% 
    ggetho(., aes(x=t, y=ID, z=activity)) + stat_tile_etho() + theme_bw(14) + ggtitle("2021"))

plot_list[[1]]<- p2019
plot_list[[2]]<- p2020
plot_list[[3]]<- p2021

ggsave(filename = paste0(path, "plots/tag_check/" , "summary_ggetho" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)

#########################################################################################
# further rethomic-plots (playaround...):

ggetho(dt, aes(x=t, y=ID, z=activity)) + stat_bar_tile_etho()
ggetho(dt, aes(x=t, y=ID, z=activity)) + stat_tile_etho() + theme_bw(14) #best 
ggetho(dt, aes(x=t, y=species, z=activity)) + stat_tile_etho() + theme_bw(14) 
ggetho(dt, aes(x=t, y=ID, z=activity)) + geom_tile() + theme_bw(14) # black&white

ggetho(dt, aes(x=t, y=activity)) + 
  stat_pop_etho() +
  facet_wrap(~species) +
  theme_bw(14)

ggetho(dt, aes(x=t, y=activity), time_wrap = hours(24)) + 
  stat_pop_etho()

ggetho(dt, aes(x=t, y=activity), time_wrap = hours(24)) + 
  stat_pop_etho() +
  facet_wrap(~ species) +
  theme_bw(14)

ggetho(dt, aes(x=t, y=activity), time_wrap = hours(24)) + 
  # time_offset = hours(12)) + 
  stat_pop_etho() +
  facet_wrap(~ ID) +
  theme_bw(14)

ggetho(dt, aes(x=t, z=activity), multiplot = 2) + stat_bar_tile_etho()

ggetho(dt, aes(x=t, z=activity), multiplot = 2) + 
  stat_bar_tile_etho() +
  facet_wrap( ~ ID) + theme_bw(14)


################################################################################
#### periodogram (rethomics) ####

per_dt <- periodogram(activity, dt, FUN = chi_sq_periodogram, 
                      resample_rate = 1/mins(10))
per_dt

ggperio(per_dt, aes(period, power, colour=species)) + 
  stat_pop_etho() + 
  scale_color_wsj() + 
  theme_bw(14)

ggperio(per_dt, aes(period, power, colour=species)) + 
  geom_line() + 
  facet_wrap( ~ ID) + 
  #scale_color_wsj() + 
  theme_bw(14)


## Circadian rythm analysis:

per_xsq_dt <- periodogram(activity, dt, FUN = chi_sq_periodogram)
per_xsq_dt

per_xsq_dt <- find_peaks(per_xsq_dt)
per_xsq_dt

ggperio(per_xsq_dt) + 
  geom_line(aes(group = ID, colour=species)) + theme_bw(14)

ggperio(per_xsq_dt, aes(y = power - signif_threshold, colour=species)) + 
  stat_pop_etho() + theme_bw(14)

summary_dt <- rejoin(per_xsq_dt[peak==1])
summary_dt

ggplot(summary_dt, aes(species, period, fill= species)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(aes(size=power -  signif_threshold), alpha=.5) +
  scale_y_hours(name = "Period") 

################################################################################
#### Bout analysis (rethomics) ####

bout_dt <- bout_analysis(activity, dt)
ggetho(bout_dt, aes(y=duration / 60, colour=species), time_wrap = hours(24)) + 
  stat_pop_etho() + 
  facet_grid(year ~ .) +
  scale_y_continuous(name= "Bout length (min)")

