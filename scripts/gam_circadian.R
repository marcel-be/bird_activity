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

# useful links
# https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/

path<- "J:/rts/rts_activity/"
df_1min<- fread(paste0(path, "bird_data_storage/tags_1min_withmeta.csv")) # A/P data with metadata (see data_preperation script)


##################################################################################################################
#### Data preperation:

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


## Subset of coverage > 75%:
nrow(df_1min[df_1min$coverage_daily>=0.75,])/nrow(df_1min) # 6.06%
df_1min<- df_1min %>% 
  filter(coverage_daily >= 0.75) %>% 
  droplevels()

## Subset of Tags with > 3 days of data (date of capture was removed already)
nrow(df_1min[df_1min$time_total >= 3,])/nrow(df_1min) # 0.8 %
df_1min<- df_1min %>% 
  filter(time_total >= 3)


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
         species_en=="European_Jay" |
         species_en=="Eurasian_Blue_Tit"| 
         species_en=="Common_Chaffinch"|
         species_en=="woodpecker"|
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
  group_by(ID, ring_ID, species_en,year_f,month,week,ydate,hour,interval) %>% 
  summarise(n_intervals=length(activity),
            n_active=length(activity[activity==1]),
            n_passive=length(activity[activity==0]),
            time_of_day=mean(time_of_day),
            time_to_rise_std = mean(time_to_rise_std)) %>% 
  mutate(active_prop =  n_active/n_intervals)

df_10min <- df_10min[seq(1, nrow(df_10min),1), ] # reduce data for faster plotting

df_1min_short <- df_1min[seq(1, nrow(df_1min),1), ] # reduce data for faster modelling


##################################################################################################################
#### Some descriptives /collinearity

summary(df_1min) 

df<- df_1min  %>% 
  count(species_en, ring_ID) %>%
  as.data.frame(.) %>% 
  count(ring_ID) %>%
  as.data.frame(.)
sum(df$n) # number of individuals

table(df_1min$activity) # okay
table(df_1min$species_en) # unbalanced? Difference of n of almost factor 10
hist(table(df_1min$date_f)) # okayish

boxplot(time_to_rise_std ~ date_f, 
        data = df_1min,
        xlab = "date", ylab = "Time to Sunrise") # two dates (first days of the season) show higher mean time_to_sunrise. But looks good overall
boxplot(time_to_rise_std ~ species_en, 
        data = df_1min,
        xlab = "date", ylab = "Time to Sunrise") # perfect

hist(df_1min$coverage_daily)
boxplot(coverage_daily ~ species_en, 
        data = df_1min,
        xlab = "date", ylab = "Time to Sunrise") # including "coverage" does not improve the models



##################################################################################################################
#### Models for circadian rhythm of bird activity (Daily rhythm)

## activity ~ time of day (centered around sunrise)
## Group level trends with different smoothness for each group 
## Bird species as grouping factor
## Corrected for Inidividual ID (RingID) and date --> activity pattern for typical Individual on a typical day/date

ggplot(df_10min, aes(y = n_active/n_intervals, x = time_to_rise_std)) + # what we aim to model
  geom_point(alpha = .2) + 
  geom_smooth(se = T) +
  facet_wrap(~species_en)

## model GI: global smoother plus group-level smoothers with differing wiggliness:
gam_GI<- bam(activity ~ species_en +
                s(time_to_rise_std, m=2,  bs="cc", k=50) +
                s(time_to_rise_std, by=species_en, m=1,  bs="cc", k=50) + # no intercept-per-species when using the "by"-argument --> specify separately (fixed or random)
                s(ring_ID, bs="re") + # random intercept for each individual
               #s(species_en, bs="re")+ # random intercept for each species
                s(ring_ID, time_to_rise_std, bs="re") + # random slope for each bird individual
                s(date_f, bs="re"),  # k equals the number of levels of grouping variable
              method ="fREML", 
              family ="binomial",
              discrete = T, 
              knots=list(time_to_rise_std=c(min_set, max_set)),
              rho= 0.27,
              AR.start = df_1min_short$start.event,
              data = df_1min_short)
# specify m=1 instead of m=2 for the group-level smoothers. This reduces collinearity between the global smoother and the group-specific terms which occasionally leads to high uncertainty around the global smoother.
# explicitly include a random effect ( s(species, bs="re")) for the intercept (or fixed factor in our case), as group-specific intercepts are not incorporated into factor "by" variable smoothers
summary(gam_GI)
anova(gam_GI)


## model I: Group-specific smoothers with different wiggliness (no global smoother)
## underlying assumption is that group-level smooth terms do not share or deviate from a common form
## Species-curves are not drawn towards a global trend (more "individualistic" than model GI)
gam_I<- bam(activity ~ species_en +
              s(time_to_rise_std, by=species_en, m=2,  bs="cc", k=50) + # no intercept-per-ID when using the "by"-argument --> specify separately (fixed or random)
              s(ring_ID, bs="re") + # random intercept for each individual
             #s(species_en, bs="re")+ # random intercept for each species
              s(ring_ID, time_to_rise_std, bs="re") + # random slope for each bird individual
              s(date_f, bs="re"),  # k equals the number of levels of grouping variable
            method ="fREML", 
            family ="binomial",
            discrete = T, 
            knots=list(time_to_rise_std=c(min_set, max_set)),
            rho= 0.54,
            AR.start = df_1min_short$start.event,
            data = df_1min_short)
#summary(gam_I)

AIC(gam_GI, gam_I) # go for gam_I

#################################################################################################################
#### Diagnostics

## check residuals
#appraise(gam_I)
gam.check(gam_I)

E<- residuals(gam_I, type="pearson")
fit<- fitted(gam_I)

p<- ggplot(data=df_1min_short, aes(y=E, x=fit))+
  geom_point(alpha=0.2)+
  xlab("Fitted values") +
  ylab("Pearson residuals")+
  theme(text = element_text(size=15))+
  geom_hline(yintercept=0)+
  geom_smooth()+
  ylim(-15, 70)+
  facet_wrap(~species_en)
ggsave(filename = paste0(path, "plots/model_output/diagnostics/" , "residuals_fitted" , ".png"),
       plot=p, width = 15, height = 9)

p<- ggplot(data=df_1min_short, aes(y=E, x=time_to_rise_std))+
  geom_point(alpha=0.2)+
  xlab("time_to_rise_std") +
  ylab("Pearson residuals")+
  theme(text = element_text(size=15))+
  geom_hline(yintercept=0)+
  ylim(-15, 70)+
  facet_wrap(~species_en)
ggsave(filename = paste0(path, "plots/model_output/diagnostics/" , "residuals_covariate" , ".png"),
       plot=p, width = 15, height = 9)
  
plot(E~df_1min_short$species_en, ylim=c(-10,60))
plot(E~df_1min_short$date_f)


# DHARMa
#simulationOutput <- simulateResiduals(fittedModel = gam_I, plot = F)
#plot(simulationOutput)
#testDispersion(simulationOutput) # underdispersion problem


## temporal autocorrelation
# https://mran.microsoft.com/snapshot/2016-10-12/web/packages/itsadug/vignettes/acf.html
acf_GI<- acf_resid(gam_I)
start_value_rho(gam_GI, plot=TRUE) #0.27
start_value_rho(gam_I, plot=TRUE) #0.54



#############################################################################################
#### Prediction: get fitted values

#data_new <- expand.grid(time_to_rise_std = seq(min_set,max_set, length.out=10),
#                        species_en = levels(df_1min$species_en),
#                        date_f = "2020-06-01",
#                        ID = levels(df_1min$ID)
#                        )

# approach one: predict for a random date (used as random factor in model)
date_fixed = "2020-06-01"
time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
data_new <- df_1min %>%  
  expand(nesting(species_en, ring_ID),date_fixed, time_to_rise_std_seq) %>% 
  rename(date_f       = date_fixed,
         time_to_rise_std = time_to_rise_std_seq)

# approach 2: predict for all dates and aggregate mean activity for all dates before plotting. Prediction takes way longer.
time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
data_new <- df_1min %>%  
  expand(nesting(species_en, ring_ID),date_f, time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq)

data_new$time_to_rise_std[data_new$time_to_rise_std == min(abs(data_new$time_to_rise_std))] <- 0 # get activity predictions for time of sunrise

pred <- predict(gam_I, 
                newdata = data_new, 
                #exclude_terms = c(s(ID),s(ID, time_to_rise_std), s(date_f)),
                se = TRUE, 
                type="link")

data_new<- data_new %>% 
  filter(date_f!= "2019-08-15") %>% 
  filter(date_f!="2020-04-05")

data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
data_new$se_min <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
data_new$se_max <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))


## plot results
p<- data_new %>% 
  group_by(species_en, ring_ID, time_to_rise_std) %>% 
  summarise_each(funs(mean)) %>% 
ggplot(data = ., 
       aes(x = time_to_rise_std, y = mu, 
           color = ring_ID, group = ring_ID))+
  #geom_point(data=df_10min, alpha = .1, 
  #           aes(x = time_to_rise_std, y = n_active/n_intervals)) +
  geom_ribbon(aes(ymin = se_min ,
                  ymax = se_max), 
              fill = "grey", color = "grey") +
  geom_line(size = .8) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #scale_color_wsj() +
  theme_bw(14) +
  xlab("Time since sunrise") + 
  ylab("Activity probability \n") + 
  ylim(0, 1) +
  #theme(legend.position = c(0.85,0.88))+
  theme(legend.position = "none") +
  facet_wrap(~species_en)
ggsave(filename = paste0(path, "plots/model_output/" , "circadian_ID" , ".png"),
       plot=p, width = 15, height = 9)


p<- data_new %>% 
  group_by(species_en, time_to_rise_std) %>% 
  summarise_each(funs(mean)) %>% 
ggplot(data = ., 
       aes(x = time_to_rise_std, y = mu, 
           color = species_en, group = species_en))+
  #geom_point(data=df_10min, alpha = .1, 
  #           aes(x = time_to_rise_std, y = n_active/n_intervals)) +
  geom_ribbon(aes(ymin = se_min ,
                  ymax = se_max), 
              fill = "grey", color = "grey") +
  geom_line(size = .8) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #scale_color_wsj() +
  theme_bw(14) +
  xlab("Time since sunrise") + 
  ylab("Activity probability \n") + 
  ylim(0, 1) +
  theme(legend.position = c(0.85,0.88))
ggsave(filename = paste0(path, "plots/model_output/" , "circadian_species" , ".png"),
       plot=p, width = 15, height = 9)



#####################################################################################################
#### Extract timing of activity values (from final curve) 

# either get all values from model that includes all data (for species level), or calculate one curve for each ID respectively (get activity values for each individual)
# Curve-Values for IDs cannot be drawn from global model, as ID is only used as random factor (curves drawn towards a global pattern)

## 1. activity values from global model (species):

# activity at time_to_rise_std = ~0
data_new %>% 
  filter(time_to_rise_std == 0) %>% 
  group_by(species_en) %>% 
  summarise(activity = mean(mu))

# activity > 0.5
data_new %>% 
  filter(mu > 0.5) %>% 
  group_by(species_en) %>% 
  summarise(a.onset = as_hms(min(time_to_rise_std)*60),
            a.end = as_hms(max(time_to_rise_std)*60))

# Peak activity: what are the two highest values for p(activity)?
data_new %>% 
  group_by(species_en) %>% 
  filter(mu == max(mu)) %>% 
  summarise(peak.a = max(mu),
            peak.a.low = se_min,
            peak.a.up = se_max)

# Peak activity timing: time of day with maximum p(activity)
data_new %>% 
  group_by(species_en) %>% 
  filter(mu == max(mu)) %>%  
  summarise(t.peak = as_hms(time_to_rise_std*60))

# finding maximum and minimum slope (potential measures for "start" and "end" of activity)
df<- derivatives(gam_I, , type = "central", term = "s(time_to_rise_std)", partial_match = TRUE)
df$species<- str_split_fixed(df$smooth, "species_en",2)[,2]

ggplot(data = df, aes(x = data, y = derivative,
                      group=species, color=species))+
  geom_line(alpha = .8)

df %>% 
  group_by(species) %>% 
  summarise(max(derivative))



## 2. Activity values for certain IDs
# possible to loop over all IDs and write all value into one table for further analysis

plot_list<- list()
for(i in 1:nlevels(df_1min$ID)){
  
  df<-  df_1min %>% 
    filter(ID==levels(ID)[i])
  
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
  
  rho<- start_value_rho(gam, plot=TRUE)
  
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
  data_new <- df_1min %>% 
    filter(ID==ID_sub) %>% 
    droplevels() %>% 
    expand(ID,date_f, time_to_rise_std_seq) %>% 
    rename(time_to_rise_std = time_to_rise_std_seq)
  
  ## get predicted values
  pred <- predict(gam, 
                  newdata = data_new, 
                  #exclude_terms = c(s(ID),s(ID, time_to_rise_std), s(date_f)),
                  se = TRUE, 
                  type="link")
  data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
  data_new$se_min <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
  data_new$se_max <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))
  
  ## get derivatives
  # find steepest slope of curve as a measure of "start" and "end" of activity
  # redo! https://rdrr.io/cran/gratia/man/derivatives.html
  df<- derivatives(gam, type = "central", newdata=data_new)
  slope_max<- df %>% 
    filter(derivative == max(derivative)|
             derivative == min(derivative)) %>% 
    rename(time_to_rise_std = data) %>% 
    select(-var, -smooth) 
  
  #ggplot(data = df, aes(x = data, y = derivative))+
  # geom_line(alpha = .8)+
  # geom_point(data=slope_max, aes(x=time_to_rise_std, y=derivative))
  
  
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
    geom_vline(xintercept = slope_max$time_to_rise_std[1],  color = "blue", size=0.5)+
    geom_vline(xintercept = slope_max$time_to_rise_std[2],  color = "red", size=0.5)+
    theme_bw(14) +
    xlab("Time since sunrise") + 
    ylab("Activity probability \n") + 
    ylim(0, 1)
  
  
  
  plot_list[[i]] <-p
}

ggsave(filename = paste0(path, "plots/model_output/" , "curve_by_tag" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)