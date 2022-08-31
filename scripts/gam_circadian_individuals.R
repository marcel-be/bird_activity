library(tidyverse)
library(lubridate)
library(data.table)
library(hms)
library(mgcv)
library(mgcViz)
library(gratia)
library(viridis)
library(tRackIT)
library(DHARMa)
library(itsadug)

path<- "J:/rts/rts_activity/"
df_1min<- fread(paste0(path, "bird_data_storage/tags_1min_withmeta.csv")) # A/P data with metadata (see data_preperation script)


##################################################################################################################
#### 1. Data preperation:

## pool woodpecker species
df_1min$species_en<- gsub("Black_Woodpecker", "woodpecker",
                          gsub("Great_Spotted_Woodpecker", "woodpecker",
                               gsub("Middle_Spotted_Woodpecker", "woodpecker", df_1min$species_en)))

## data format
df_1min$hour           <- hour(df_1min$timestamp) # Hour of the day
df_1min$minute         <- minute(df_1min$timestamp)
df_1min$month_f        <- factor(month(df_1min$date))
df_1min$year_f         <- factor(df_1min$year)
df_1min$ydate_f        <- as.factor(df_1min$ydate)
df_1min$species_en     <- as.factor(df_1min$species_en)
df_1min$ring_ID        <- as.factor(df_1min$ring_ID)
df_1min$date_CET       <- date(df_1min$date_CET)
df_1min$date_f         <- as.factor(df_1min$date_CET)
df_1min$ID             <- as.factor(df_1min$ID)
df_1min$week           <- week(df_1min$timestamp)
df_1min$brood_patch    <- as.factor(df_1min$brood_patch)
df_1min$time_of_day    <- as.numeric(as_hms(df_1min$timestamp_CET))/3600 # numeric time of day (not used)
df_1min$timestamp_CET  <- fasttime::fastPOSIXct(df_1min$timestamp_CET, tz="CET")
#df_1min$start_datetime<- fasttime::fastPOSIXct(df_1min$start_datetime, tz="CET") # only of needed in further analysis
#df_1min$stop_datetime <- fasttime::fastPOSIXct(df_1min$stop_datetime, tz="CET") # only of needed in further analysis

# nrow= 2128539

## set range of "time_to_rise" as -5 to 18. Not all Individuals cover the whole range from -8.7 to 22.1 (all data). This max and min values will be used for cc-smoother. -5 to 18 is the range that all bird individuals fall into (see script "range_time_to_rise.R")
nrow(df_1min[df_1min$time_to_rise_std<=18 & df_1min$time_to_rise_std>=-5,])/nrow(df_1min) # 95.4%
df_1min<- df_1min %>%
  filter(time_to_rise_std >= -5) %>% 
  filter(time_to_rise_std <= 18)

## Subset of coverage > 75%:
nrow(df_1min[df_1min$coverage_daily>=0.75,])/nrow(df_1min) # 94.1 %
df_1min<- df_1min %>% 
  filter(coverage_daily >= 0.75) %>% 
  droplevels()

## Subset of Tags with > 3 days of data (date of capture was removed already)
nrow(df_1min[df_1min$time_total >= 3 & df_1min$ID != "210408_150113_40" & df_1min$ID != "210622_150007_40",])/nrow(df_1min) # 0.99 %
df_1min<- df_1min %>% 
  filter(time_total >= 3) %>% 
  filter(ID != "210408_150113_40") %>% 
  filter(ID != "210622_150007_40")

## Subset of species with at least 4 individuals:
df_1min  %>% 
  count(species_en, ring_ID) %>%
  as.data.frame(.) %>% 
  count(species_en) %>%
  mutate(species_en = fct_reorder(species_en, n)) %>% 
  ggplot(. ,aes(y=species_en, x=n, label = round(n,digits=1))) +
  geom_bar(stat="identity", color="steelblue", fill="steelblue")+
  ggtitle("Count of tagged individuals \n(some individuals were tagged two times)")+
  geom_text(size = 5, position = position_stack(vjust = 1.025))+
  xlab("Count")+
  ylab("Species")+
  theme_minimal()+
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"),
        legend.title=element_blank(),
        legend.position = c(0.9, 0.1),
        legend.text = element_text(size=15))


df_1min<- df_1min %>% 
  filter(species_en=="Eurasian_Blackcap"|
           species_en=="Great_Tit"|
           species_en=="European_Robin" |
           species_en=="Common_Blackbird" |
           species_en=="Eurasian_Jay" |
           species_en=="Eurasian_Blue_Tit"| 
           species_en=="Common_Chaffinch"|
           species_en=="Wood_Warbler"
  ) %>% 
  droplevels() # subset of species with most individuals


## set maximum and minimum sunset time (for cc-smoother)
min_set  <- min(df_1min$time_to_rise_std)
max_set  <- max(df_1min$time_to_rise_std)
min_ydate<- min(df_1min$ydate)
max_ydate<- max(df_1min$ydate)

## marking the start of each time series (start of each Bird ID) as "TRUE". Will be incorporated as autocorrelation structure. Each ID is treated as autocorrelated time series.
df_1min<- start_event(df_1min, column="timestamp_CET", event="ID") # Package "itsadug"

str(df_1min)

## create 10-min-intervals for faster plotting:
df_10min <- df_1min %>% 
  mutate(species_en = as.factor(species_en),
         interval   = as_hms(floor_date(timestamp_CET, unit="5minutes"))) %>% 
  group_by(ID, ring_ID, species_en,year_f,brood_patch,sex, month,week,ydate, date_f,hour,interval) %>% 
  summarise(n_intervals=length(activity),
            n_active=length(activity[activity==1]),
            n_passive=length(activity[activity==0]),
            time_of_day=mean(time_of_day),
            time_to_rise_std = mean(time_to_rise_std)) %>% 
  mutate(active_prop =  n_active/n_intervals)

df_10min <- df_10min[seq(1, nrow(df_10min),5), ] # reduce data for faster plotting

df_1min_short <- df_1min[seq(1, nrow(df_1min),1), ] # reduce data for faster modelling


##########################################################################################################
## 2. Activity values for all Individuals
# possible to loop over all IDs and write all value into one table for further analysis

plot_list<- list()
n_ID<- nlevels(df_1min$ID)
a<- seq(1,(2*n_ID)-1, length.out=n_ID)
b<- seq(2, 2*n_ID   , length.out=n_ID)

for(i in 1:nlevels(df_1min$ID)){
  
  print(levels(df_1min$ID)[i])
  
  df<-  df_1min %>% 
    filter(ID==levels(ID)[i]) %>% 
    droplevels()
  
  min_set  <- min(df$time_to_rise_std)
  max_set  <- max(df$time_to_rise_std)
  
  gam <-df %>% 
    bam(activity ~ 
          s(time_to_rise_std, m=2,  bs="cc", k=50) +
          s(date_f, bs="re"),  # k equals the number of levels of grouping variable
        method ="fREML", 
        family ="binomial",
        discrete = T, 
        knots=list(time_to_rise_std=c(min_set, max_set)),
        data = .)
  
  rho<- start_value_rho(gam, plot=F)
  
  gam <- df %>% 
    bam(activity ~ 
          s(time_to_rise_std, m=2,  bs="cc", k=50) +
          s(date_f, bs="re"),  # k equals the number of levels of grouping variable
        method ="fREML", 
        family ="binomial",
        discrete = T, 
        knots=list(time_to_rise_std=c(min_set, max_set)),
        rho= rho,
        data = .)
  
  time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
  data_new <- df%>% 
    expand(ID, time_to_rise_std_seq) %>% 
    rename(time_to_rise_std = time_to_rise_std_seq) %>% 
    mutate(date_f = NA)
  
  ## get predicted values
  pred <- predict.gam(gam, 
                      newdata = data_new, 
                      exclude = "s(date_f)",
                      newdata.guaranteed=TRUE,
                      se = TRUE, 
                      type="link")
  data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
  data_new$se_min <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
  data_new$se_max <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))
  
  ## get derivatives
  # find steepest slope of curve as a measure of "start" and "end" of activity
  # redo! https://rdrr.io/cran/gratia/man/derivatives.html
  
  
  # deri<- derivatives(gam, type = "central",term = "s(time_to_rise_std)", newdata=data_new)
  # slope_max<- deri %>% 
  #  group_by(data) %>% 
  #  summarise_each(funs(mean)) %>% 
  #  filter(derivative == max(derivative)|
  #          derivative == min(derivative)) %>% 
  #  rename(time_to_rise_std = data) %>% 
  #  select(-var, -smooth) 
  
  name<- paste0(as.character(df$ID[1]),
                " | ", as.character(df$species_en[1]),
                " | brood patch=", as.character(df$brood_patch[1]),
                " | year=", as.character(df$year_f[1]),
                " | sex=", as.character(df$sex[1]),
                " | ringID=", as.character(df$ring_ID[1]))
  
  p<- data_new %>% 
    group_by(time_to_rise_std) %>% 
    summarise_each(funs(mean)) %>% 
    ggplot(data = ., 
           aes(x = time_to_rise_std, y = mu))+
    geom_ribbon(aes(ymin = se_min ,
                    ymax = se_max), 
                fill = "grey", color = "grey") +
    geom_line(size = .8) + 
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    #geom_vline(xintercept = slope_max$time_to_rise_std[1],  color = "blue", size=0.5)+
    #geom_vline(xintercept = slope_max$time_to_rise_std[2],  color = "red", size=0.5)+
    theme_bw(14) +
    xlab("Time since sunrise") + 
    ylab("Activity probability \n") + 
    ylim(0, 1)+
    ggtitle(name)
  
  plot_list[[a[i]]] <-p
  
  df <- df_10min %>%  
    filter(ID == levels(ID)[i]) 
  
  name<- paste0(as.character(df$species_en[1]),
                " | brood patch=", as.character(df$brood_patch[1]),
                " | year=", as.character(df$year_f[1]),
                " | sex=", as.character(df$sex[1]))
  
  p<-  ggplot(df, aes(y = n_active/n_intervals, x = time_to_rise_std)) +
    geom_point(alpha = .20) + 
    geom_smooth(method = 'gam', 
                formula = y ~ s(x, bs = "tp", k = 25)) +
    #scale_color_viridis(discrete = T) + 
    facet_wrap(~ date_f) +
    #facet_wrap(~ ydate + coverage_daily, labeller = label_both) +
    theme_bw()+
    ggtitle(paste0(levels(df_10min$ID)[i]," ", name))
  
  plot_list[[b[i]]] <-p
}

ggsave(filename = paste0(path, "plots/model_output/" , "curve_by_tag" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)
