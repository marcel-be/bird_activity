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

mean(df_1min$daylength)
##################################################################################################################
#### 1. Data preperation: ####

## pool woodpecker species
df_1min$species_en<- gsub("Black_Woodpecker", "woodpecker",
                          gsub("Great_Spotted_Woodpecker", "woodpecker",
                               gsub("Middle_Spotted_Woodpecker", "woodpecker", df_1min$species_en)))

## change data format
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


## marking the start of each time series (start of each Bird ID) as "TRUE". Will be incorporated as autocorrelation structure. Each ID is treated as autocorrelated time series.
df_1min<- start_event(df_1min, column="timestamp_CET", event="ID") # Package "itsadug"


## save final dataframe
fwrite(df_1min, paste0(path, "bird_data_storage/tags_1_min_for_analysis.csv"))
#df_1min<- fread(paste0(path, "bird_data_storage/tags_1_min_for_analysis.csv"), stringsAsFactors = T)

df_1min_short <- df_1min[seq(1, nrow(df_1min),1), ] # reduce data for faster modelling

## set maximum and minimum sunset time (for cc-smoother in model)
min_set  <- min(df_1min$time_to_rise_std)
max_set  <- max(df_1min$time_to_rise_std)
min_ydate<- min(df_1min$ydate)
max_ydate<- max(df_1min$ydate)


## create 10-min-intervals for faster plotting:
df_10min <- df_1min %>% 
  mutate(species_en = as.factor(species_en),
         interval   = as_hms(floor_date(timestamp_CET, unit="10minutes"))) %>% 
  group_by(ID, ring_ID, species_en,year_f,brood_patch,sex, month,week,ydate, date_f,hour,interval) %>% 
  summarise(n_intervals=length(activity),
            n_active=length(activity[activity==1]),
            n_passive=length(activity[activity==0]),
            time_of_day=mean(time_of_day),
            time_to_rise_std = mean(time_to_rise_std)) %>% 
  mutate(active_prop =  n_active/n_intervals)

df_10min <- df_10min[seq(1, nrow(df_10min),1), ] # reduce data for faster plotting
fwrite(df_10min, paste0(path, "bird_data_storage/tags_10_min_for_plotting.csv"))


##################################################################################################################
#### 2. Some descriptives / check for collinearity ####

summary(df_1min) 

df<- df_1min  %>% 
  count(species_en, ring_ID) %>%
  as.data.frame(.) %>% 
  count(ring_ID) %>%
  as.data.frame(.)
sum(df$n) # number of individuals

df<- df_1min  %>% 
  count(species_en, ID) %>%
  as.data.frame(.) %>% 
  count(ID) %>%
  as.data.frame(.)
sum(df$n) # number of Tags

table(df_1min$activity) # okay
table(df_1min$species_en) # unbalanced? Difference of n of almost factor 10
table(df_1min$ring_ID)
hist(table(df_1min$date_f)) # okayish

boxplot(time_to_rise_std ~ ring_ID, 
        data = df_1min,
        xlab = "Individual", ylab = "Time to Sunrise")

boxplot(time_to_rise_std ~ date_f, 
        data = df_1min,
        xlab = "date", ylab = "Time to Sunrise") #  looks good overall --> standardizing of time-to-sunrise seems to work fine!

boxplot(time_to_rise_std ~ species_en, 
        data = df_1min,
        xlab = "date", ylab = "Time to Sunrise") # perfect

hist(df_1min$coverage_daily)
boxplot(coverage_daily ~ species_en, 
        data = df_1min,
        xlab = "date", ylab = "Time to Sunrise") # including "coverage" into the models does not improve the models
# coverage of most tags > 90%

# data coverage across year for selected species
p<- df_1min %>% 
  mutate(activity=1) %>% 
  ggplot(., aes(y = activity, x = ydate,
                group = ring_ID, color = ring_ID)) +
  geom_point(alpha = .25, size = 2) + 
  facet_wrap(~ species_en, ncol=1) + 
  theme_bw() +
  theme(legend.position = "none")+
  ggtitle("coverage of data over year")+
  ylab("Activity")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(filename = paste0(path, "plots/model_output/diagnostics/" , "coverage_year" , ".png"),
       plot=p, width = 15, height = 9)




##################################################################################################################
#### 3. Models for circadian rhythm of bird activity (Daily rhythm) ####

## activity ~ time of day (centered around sunrise)
## Group level trends with different smoothness for each group 
## Bird species as grouping factor
## Corrected for Inidividual ID (RingID) and date --> activity pattern for typical Individual on a typical day/date

ggplot(df_10min, aes(y = n_active/n_intervals, x = time_to_rise_std)) + # what we aim to model
  geom_point(alpha = .2) + 
  geom_smooth(se = T) +
  facet_wrap(~species_en)

## model GI: global smoother plus group-level smoothers with differing wiggliness (not used in further steps):
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



## model I: Group-specific smoothers with different wiggliness (no global smoother - used for further steps)
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
            rho= 0.53,
            AR.start = df_1min_short$start.event,
            data = df_1min_short)
#summary(gam_I)

## same model for bird individuals instead of species (draw activity-curve-values from this model)
gam_I<- bam(activity ~ ring_ID +
              s(time_to_rise_std, by=ring_ID, m=2,  bs="cc", k=75) + # no intercept-per-ID when using the "by"-argument --> specify separately (fixed or random)
              s(species_en, bs="re") + # random intercept for each species
              s(species_en, time_to_rise_std, bs="re") + # random slope for each species
              s(date_f, bs="re"),  # k equals the number of levels of grouping variable
            method ="fREML", 
            family ="binomial",
            discrete = T, 
            knots=list(time_to_rise_std=c(min_set, max_set)),
            rho= 0.52,
            AR.start = df_1min_short$start.event,
            data = df_1min_short)


AIC(gam_GI, gam_I) # go for gam_I

## save model

saveRDS(gam_I, file=paste0(path,"bird_data_storage/models/gam_I_individuals.rda"))
gam_I<- readRDS(paste0(path,"bird_data_storage/models/gam_I_individuals.rda"))

#################################################################################################################
#### 5. Diagnostics ####

## check residuals
#appraise(gam_I)
par(mfrow=c(2,2))
gam.check(gam_I)
par(mfrow=c(1,1))

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
  facet_wrap(~ring_ID)
ggsave(filename = paste0(path, "plots/model_output/diagnostics/" , "residuals_fitted_individual" , ".png"),
       plot=p, width = 15, height = 9)

p<- ggplot(data=df_1min_short, aes(y=E, x=fit))+
  geom_point(alpha=0.2)+
  xlab("Fitted values") +
  ylab("Pearson residuals")+
  theme(text = element_text(size=15))+
  geom_hline(yintercept=0)+
  geom_smooth()+
  ylim(-15, 70)+
  facet_wrap(~species_en)
ggsave(filename = paste0(path, "plots/model_output/diagnostics/" , "residuals_fitted_species" , ".png"),
       plot=p, width = 15, height = 9)

p<- ggplot(data=df_1min_short, aes(y=E, x=time_to_rise_std))+
  geom_point(alpha=0.2)+
  xlab("time_to_rise_std") +
  ylab("Pearson residuals")+
  theme(text = element_text(size=15))+
  geom_hline(yintercept=0)+
  ylim(-15, 70)+
  facet_wrap(~ring_ID)
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
#acf_GI<- acf_resid(gam_I)
start_value_rho(gam_GI, plot=TRUE) #0.54
start_value_rho(gam_I, plot=TRUE) #0.52



#############################################################################################
#### 6. Prediction: get fitted values ####

#data_new <- expand.grid(time_to_rise_std = seq(min_set,max_set, length.out=10),
#                        species_en = levels(df_1min$species_en),
#                        date_f = "2020-06-01",
#                        ID = levels(df_1min$ID)
#                        )

## approach one: predict for a random date (used as random factor in model) - don´t use because of variation of outcomes for different dates!
date_fixed = "2020-06-01"
time_to_rise_std_seq<- seq(min_set, max_set, length.out=100)
data_new <- df_1min %>%  
  expand(nesting(species_en, ring_ID),date_fixed, time_to_rise_std_seq) %>% 
  rename(date_f       = date_fixed,
         time_to_rise_std = time_to_rise_std_seq)

## approach 2: predict for all possible dates (random factor "date_f") and then aggregate mean activity of all dates ("typical date") before plotting --> Mean activity-prediction (across all dates) for each datapoint of "time_to_rise" will be calculated. Prediction takes long! --> Similar to approach 3
# Problem: for each bird all possible dates (2019-2021) are calculated
time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
data_new <- df_1min %>%  
  expand(nesting(species_en, ring_ID),date_f, time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq)

## approach 3: Set "date_f" in new x-dataframe as "NA". Exclude "s(date_f)" in predict-function and make sure that newdata.guaranteed=TRUE. # Results are very similar to approach 2 (mean value for activity for a "typical" date). Much faster than approach 2, but no date-specific results. For comparison of different approaches see script "gam_prediction_comparison.R"
time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
data_new <- df_1min %>%  
  expand(nesting(species_en, ring_ID),time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq) %>% 
  mutate(date_f=NA)

## approach 4: calculate predictions for all dates where birds were transmitting signals. Question: can this dates be compared with each other, as the data comes from different month/years?
time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
data_new <- df_1min %>%  
  expand(nesting(species_en, ring_ID, date_f), time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq)



data_new$time_to_rise_std[data_new$time_to_rise_std == -min(abs(data_new$time_to_rise_std))] <- 0 # get activity predictions for time of sunrise

pred <- predict.gam(gam_I, 
                    newdata = data_new, 
                    #exclude = "s(date_f)",
                    #newdata.guaranteed=TRUE,
                    se = TRUE, 
                    type="link")
# this prediction excludes random effect date --> very similar to predictions of mean activity of all dates ("typical date")


#data_new<- data_new %>% 
#  filter(date_f!= "2019-08-15") %>% 
#  filter(date_f!="2020-04-05")

data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
data_new$ci_lower <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
data_new$ci_upper <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))


## save model predictions
fwrite(data_new, paste0(path,"bird_data_storage/models/model_prediction.csv"))


#######################################################################################################
#### 7. plot results ####

## at individual level (one curve per individual - use individual model)
p<- data_new %>% 
  group_by(species_en, ring_ID, time_to_rise_std) %>% 
  summarise_each(funs(mean)) %>% 
ggplot(data = ., 
       aes(x = time_to_rise_std, y = mu, 
           color = ring_ID, group = ring_ID))+
  #geom_point(data=df_10min, alpha = .1, 
  #           aes(x = time_to_rise_std, y = n_active/n_intervals)) +
  geom_ribbon(aes(ymin = ci_lower ,
                  ymax = ci_upper), 
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

## at species level (one curve per species - use species model)
p<- data_new %>% 
  group_by(species_en, time_to_rise_std) %>% 
  summarise_each(funs(mean)) %>% 
ggplot(data = ., 
       aes(x = time_to_rise_std, y = mu, 
           color = species_en, group = species_en))+
  #geom_point(data=df_10min, alpha = .1, 
  #           aes(x = time_to_rise_std, y = n_active/n_intervals)) +
  geom_ribbon(aes(ymin = ci_lower ,
                  ymax = ci_upper), 
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

## activity across year (only working if "date_f" is defined in "data_new":
df<- data_new %>% 
  group_by(species_en, ring_ID, date_f) %>%   
  filter(mu == max(mu)) %>% 
  summarise(peak = max(mu))  %>% 
  mutate(ydate = yday(date_f))
ggplot(data = df,  aes(x = ydate, y = peak, group=ring_ID, color=ring_ID))+
  geom_line(size = .8) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  facet_wrap(~species_en)+
  theme(legend.position = "none")
rm(df)


#####################################################################################################
#### 8. Extract timing of activity values (from final curve) ####

# either get all values from model that includes all data (for species level), or calculate one curve for each ID respectively (get activity values for each individual)
# Curve-Values for IDs cannot be drawn from global model, as ID is only used as random factor (curves drawn towards a global pattern)

df_10min<- fread(paste0(path, "bird_data_storage/tags_10_min_for_plotting.csv"), stringsAsFactors = T)
data_new<- fread(paste0(path,"bird_data_storage/models/model_prediction.csv"), stringsAsFactors=T)
data_new$date_f<- as.factor(as.character(data_new$date_f)) 
df_1min<- fread(paste0(path, "bird_data_storage/tags_1_min_for_analysis.csv"), stringsAsFactors = T)
df_1min$date_f <- as.factor(as.character(df_1min$date_f))

## activity values from global model (Individuals):

## 1. activity at time_to_rise_std = 0
df<- data_new %>% 
  filter(time_to_rise_std == 0) %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  summarise(act_at_sunrise_mean = mu,
            act_at_sunrise_lowerCI = ci_lower, 
            act_at_sunrise_upperCI = ci_upper)

## 2. activity > 0.5 (doesn´t make sence since some individuals don´t reach 0.5 activity margin)
df2<- data_new %>% 
  filter(mu > 0.5) %>% 
  group_by(ring_ID,species_en, date_f) %>% 
  summarise(#act_0.5_onset = as_hms(min(time_to_rise_std)*60), 
            #act_0.5_end = as_hms(max(time_to_rise_std)*60)) %>%   # output in minutes
            act_0.5_onset = min(time_to_rise_std), 
            act_0.5_end = max(time_to_rise_std)) %>%   # output in hours (decimal numbers)
  as.data.frame() %>% 
  select(-species_en) %>% 
  left_join(df,., by=c("ring_ID","date_f"))


## 3. Peak activity: what is the highest value for p(activity)?
df3 <- data_new %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(mu == max(mu)) %>% 
  summarise(peak_act = max(mu),
            peak_act_lowerCI = ci_lower,
            peak_act_upperCI = ci_upper) %>% 
  as.data.frame() %>% 
  select(-species_en) %>% 
  left_join(df2,., by=c("ring_ID","date_f"))

## 4. mean activity 
df4<- data_new %>% 
  group_by(ring_ID, species_en,  date_f) %>% 
  summarise(mean_act = mean(mu)) %>% 
  as.data.frame() %>% 
  select(-species_en) %>% 
  left_join(df3,., by=c("ring_ID","date_f"))

## 5. Peak activity timing: time of day with maximum p(activity) --> makes little sense as the curves are very wobbly
df5<- data_new %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(mu == max(mu)) %>%  
  summarise(t.peak = time_to_rise_std) %>% 
  as.data.frame() %>% 
  select(-species_en) %>% 
  left_join(df4,., by=c("ring_ID","date_f"))

## 6. area under the curve?
# ...

## 7. finding maximum and minimum slope (potential measures for "start" and "end" of activity)

#  calculate differences (increase) in mu for each datapoint for each bird and date (in comparison to the respective, former datapoint) --> highest value for this difference represents the strongest increase in mu ergo the highest slope of the activity slope
# "diff"-function calculates the difference to the former datapoint
# calculating "diff" for the whole dataset results in big diff-values between day-shifts (steepest as/descent at first or last timestamp of the day) --> needs to be calculated for each day seperately
# new dataframe "data_new_diff" will have less datapoints than "data_new", because for the first datapoint of each bird individual and date, no differnece in mu can be calculated. This applies only for the very first datapoint of each individual and date (out of thousands) and can be neglected
# no idea if calculating the CI makes any sense here...

data_new_diff<- data.frame()

for(i in 1:nlevels(data_new$ring_ID)){ # find a nicer, vectorized approach!
  
diff_1 <- data_new %>% 
  filter(ring_ID==levels(ring_ID)[i]) %>% 
  droplevels()

for(j in 1:nlevels(diff_1$date_f)){

diff_2 <- diff_1 %>% 
  filter(date_f==levels(date_f)[j]) %>% 
           droplevels()
  
diff_3 <- slice_tail(diff_2, n=(nrow(diff_2)-1)) %>% 
  mutate(diff=diff(diff_2$mu),
         diff_lowerCI = diff(diff_2$ci_lower),
         diff_upperCI = diff(diff_2$ci_upper))

data_new_diff<- rbind(diff_3,data_new_diff)
}
}

# steepest_ascend:
df6<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(time_to_rise_std < 5) %>% 
  filter(diff == max(diff,na.rm=TRUE)) %>% 
  summarise(steepest_ascend = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df5,., by=c("ring_ID","date_f"))


# steepest_ascend lower CI: 
df7<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(diff_lowerCI == max(diff_lowerCI,na.rm=TRUE)) %>% 
  summarise(steepest_ascend_lowerCI = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df6,., by=c("ring_ID","date_f")) #does this make sense???

# steepest_ascend upper  CI: 
df8<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(diff_upperCI == max(diff_upperCI,na.rm=TRUE)) %>% 
  summarise(steepest_ascend_upperCI = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df7,., by=c("ring_ID","date_f")) #does this make sense???

# steepest_descend
df9<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(time_to_rise_std > 10) %>% 
  filter(diff == min(diff,na.rm=TRUE)) %>% 
  summarise(steepest_descend = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df8,., by=c("ring_ID","date_f"))

# steepest_descend lower CI: 
df10<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(diff_lowerCI == min(diff_lowerCI,na.rm=TRUE)) %>% 
  summarise(steepest_descend_lowerCI = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df9,., by=c("ring_ID","date_f")) #does this make sense???

# steepest_descend upper  CI: 
df11<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(diff_upperCI == min(diff_upperCI,na.rm=TRUE)) %>% 
  summarise(steepest_descend_upperCI = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df10,., by=c("ring_ID","date_f")) #does this make sense???


# 8. Time of sunset of each date
# can´t be used as the data was corrected for the mean daylength
dfXX<- df_1min %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  summarise(time_of_sunset = mean(daylength))%>% 
  as.data.frame %>% 
  select(-species_en) %>% 
  left_join(df11,., by=c("ring_ID","date_f")) 


# 9. Activity at sunset
# the mean daylength that was used to standardize "time to rise" is 15.51328
sunset<- data_new[which(abs(data_new$time_to_rise_std - 15.51328) == min(abs(data_new$time_to_rise_std - 15.51328)))]$time_to_rise_std[1]
df12<- data_new %>% 
  filter(time_to_rise_std == sunset) %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  summarise(act_at_sunset_mean = mu,
            act_at_sunset_lowerCI = ci_lower, 
            act_at_sunset_upperCI = ci_upper)%>% 
  as.data.frame %>% 
  select(-species_en) %>% 
  left_join(df11,., by=c("ring_ID","date_f"))

df_act_charac <- df12

## safe file
fwrite(df_act_charac, paste0(path,"bird_data_storage/activity_characteristics/activity_characteristics.csv"))



###### plot in one pdf

plot_list<- list()

for(i in 1:nlevels(data_new$ring_ID)){ #nlevels(data_new$ring_ID)
  
  print(levels(data_new$ring_ID)[i])
  
  df<-  data_new %>% 
    filter(ring_ID==levels(ring_ID)[i]) %>% 
    droplevels()
  
  df_10min_sub<- df_10min %>% 
    filter(ring_ID==levels(data_new$ring_ID)[i]) %>% 
    droplevels()
  
  df_act_charac_sub<- df_act_charac %>% 
    filter(ring_ID==levels(data_new$ring_ID)[i]) %>% 
    droplevels() %>% 
    summarise_each(funs(mean))
  
  name<- paste0(as.character(df$ring_ID[1]),
                " | ", as.character(df$species_en[1]),
                " | brood patch=", as.character(df_10min_sub$brood_patch[1]),
                " | date=", as.character(df_10min_sub$date_f[1]),
                #" | year=", as.character(df_10min_sub$year_f[1]),
                " | sex=", as.character(df_10min_sub$sex[1]))
  
p<- df %>% 
    group_by(time_to_rise_std) %>% 
    summarise_each(funs(mean)) %>% 
    ggplot(data = ., 
           aes(x = time_to_rise_std, y = mu))+
    geom_ribbon(aes(ymin = ci_lower ,
                    ymax = ci_upper), 
                fill = "grey", color = "grey") +
    geom_line(size = .8) + 
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_hline(yintercept = df_act_charac_sub$peak_act, linetype = "dashed", color="steelblue", size=1) +
      geom_text(data= df_act_charac_sub, aes(x=-3.5, label="peak activity", y=peak_act), colour="steelblue", angle=0, vjust = 1.1, text=element_text(size=11))+  
    geom_vline(xintercept = 0, linetype = "dashed", color="orange", size=1)+
      geom_text(aes(x=0, label="sunrise", y=0), colour="orange", angle=0, hjust = -0.1, text=element_text(size=11))+
    geom_vline(xintercept = 15.51328, linetype = "dashed", color="navy",  size=1) +
      geom_text(aes(x=15.51328, label="sunset", y= 0), colour="navy", angle=0, hjust = -0.1, text=element_text(size=11))+
    geom_vline(xintercept = df_act_charac_sub$steepest_ascend, linetype = "dashed", color="darkgreen",  size=1) +
      geom_text(data= df_act_charac_sub, aes(x=steepest_ascend, label="activity \n onset", y=0.25), colour="darkgreen", angle=0, hjust = -0.1, text=element_text(size=11))+
    geom_vline(xintercept = df_act_charac_sub$steepest_descend, linetype = "dashed", color="red", size=1) +
      geom_text(data= df_act_charac_sub, aes(x=steepest_descend, label="activity \n offset", y=0.25), colour="red", angle=0, hjust = -0.1, text=element_text(size=11))+
    theme_bw(14) +
    xlab("Time since sunrise") + 
    ylab("Activity probability \n") + 
    ylim(0, 1)+
    ggtitle(name)

  plot_list[[i]] <-p
}
  
ggsave(filename = paste0(path, "plots/model_output/" , "activity_characteristics" , ".pdf"),
         plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
         width = 15, height = 9)
  
  



