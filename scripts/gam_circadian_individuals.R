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

##################################################################################################################
#### 1. Data import:

df_1min<- fread(paste0(path, "bird_data_storage/tags_1_min_for_analysis.csv"))

df_1min$ydate_f        <- as.factor(df_1min$ydate)
df_1min$species_en     <- as.factor(df_1min$species_en)
df_1min$ring_ID        <- as.factor(df_1min$ring_ID)
df_1min$date_CET       <- date(df_1min$date_CET)
df_1min$date_f         <- as.factor(df_1min$date_CET)
df_1min$ID             <- as.factor(df_1min$ID)
df_1min$brood_patch    <- as.factor(df_1min$brood_patch)
df_1min$timestamp_CET  <- fasttime::fastPOSIXct(df_1min$timestamp_CET, tz="CET")

## set maximum and minimum sunset time (for cc-smoother in model)
min_set  <- min(df_1min$time_to_rise_std)
max_set  <- max(df_1min$time_to_rise_std)
min_ydate<- min(df_1min$ydate)
max_ydate<- max(df_1min$ydate)


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


##########################################################################################################
## 2. model for one individual

df<-  df_1min %>% 
  filter(ring_ID==levels(df_1min$ring_ID)[1]) %>% 
  droplevels()

min_set  <- min(df$time_to_rise_std)
max_set  <- max(df$time_to_rise_std)

gam <- bam(activity ~ 
        s(time_to_rise_std, by=date_f, m=2,  bs="cc", k=50)+
        s(date_f, bs="re"),  
      method ="fREML", 
      family ="binomial",
      discrete = T, 
      knots=list(time_to_rise_std=c(min_set, max_set)),
      data = df)

rho<- start_value_rho(gam, plot=F)

gam <- df %>% 
  bam(activity ~ 
        s(time_to_rise_std, by=date_f, m=2,  bs="cc", k=30)+
        s(date_f, bs="re"),   
      method ="fREML", 
      family ="binomial",
      discrete = T, 
      knots=list(time_to_rise_std=c(min_set, max_set)),
      rho= rho,
      data = .)

par(mfrow=c(2,2))
gam.check(gam)
par(mfrow=c(1,1))

time_to_rise_std_seq<- seq(min_set,max_set, length.out=1000)
data_new <- df%>% 
  expand(nesting(ring_ID, date_f), time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq)

## get predicted values
pred <- predict.gam(gam, 
                    newdata = data_new, 
                    #exclude = "s(date_f)",
                    #newdata.guaranteed=TRUE,
                    se = TRUE, 
                    type="link")
data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
data_new$se_min <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
data_new$se_max <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))


data_new %>% 
  ggplot(data = ., 
         aes(x = time_to_rise_std, y = mu, group=date_f, color=date_f))+
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
  ylim(0, 1)





##########################################################################################################
## 3. Daily Activity curves for all Tags/Individuals
# the loop creates one ativity curve for each Tag and attaches the daily activity-plots. The output is used to check for raw-data-issues

plot_list<- list()
n_ID<- nlevels(df_1min$ID)
a<- seq(1,(2*n_ID)-1, length.out=n_ID)
b<- seq(2, 2*n_ID   , length.out=n_ID)
data_new_final<- data.frame()

for(i in 63:nlevels(df_1min$ring_ID)){
  
  print(levels(df_1min$ring_ID)[i])
  
  df<-  df_1min %>% 
    filter(ring_ID==levels(ring_ID)[i]) %>% 
    droplevels()
  
  min_set  <- min(df$time_to_rise_std)
  max_set  <- max(df$time_to_rise_std)
  
  gam <- bam(activity ~ 
               s(time_to_rise_std, by=date_f, m=2,  bs="cc", k=35)+
               s(date_f, bs="re"),  
             method ="fREML", 
             family ="binomial",
             discrete = T, 
             knots=list(time_to_rise_std=c(min_set, max_set)),
             data = df)
  
  rho<- start_value_rho(gam, plot=F)
  
  gam <- df %>% 
    bam(activity ~ 
          s(time_to_rise_std, by=date_f, m=2,  bs="cc", k=30)+
          s(date_f, bs="re"),   
        method ="fREML", 
        family ="binomial",
        discrete = T, 
        knots=list(time_to_rise_std=c(min_set, max_set)),
        rho= rho,
        data = .)
  
  time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
  data_new <- df%>% 
    expand(nesting(ring_ID, date_f), time_to_rise_std_seq) %>% 
    rename(time_to_rise_std = time_to_rise_std_seq)
  
  ## get predicted values
  pred <- predict.gam(gam, 
                      newdata = data_new, 
                      #exclude = "s(date_f)",
                      #newdata.guaranteed=TRUE,
                      se = TRUE, 
                      type="link")
  data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
  data_new$se_min <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
  data_new$se_max <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))
  
  data_new_final<- rbind(data_new_final, data_new)

  name<- paste0(as.character(df$ring_ID[1]),
                " | ", as.character(df$species_en[1]),
                #" | brood patch=", as.character(df$brood_patch[1]),
                #" | year=", as.character(df$year_f[1]),
                " | sex=", as.character(df$sex[1]))
  
  p<- data_new %>% 
    ggplot(data = ., 
           aes(x = time_to_rise_std, y = mu, group=date_f, color=date_f))+
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
    filter(ring_ID == levels(ring_ID)[i]) 
  
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

ggsave(filename = paste0(path, "plots/model_output/" , "curve_by_individual" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)

## save model predictions
fwrite(data_new_final, paste0(path,"bird_data_storage/models/model_prediction_individuals.csv"))


#################################################################################################################################
#### 3. Extract timing of activity values (from final curve) ####

df_10min<- fread(paste0(path, "bird_data_storage/tags_10_min_for_plotting.csv"), stringsAsFactors = T)
data_new<- fread(paste0(path,"bird_data_storage/models/model_prediction_individuals.csv"), stringsAsFactors=T)
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
    
    diff_3<- slice_tail(diff_2, n=(nrow(diff_2)-1)) %>% 
      mutate(diff=diff(diff_2$mu))
    
    data_new_diff<- rbind(diff_3,data_new_diff)
  }
}

# steepest_ascend:
df6<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f) %>% 
  filter(time_since_sunrise_std >= -1.5 & time_since_sunrise_std <= 16.5)%>% 
  filter(time_since_sunrise_std <=  1.5 & time_since_sunrise_std >= 13.5) %>% 
  filter(diff == max(diff,na.rm=TRUE)) %>% 
  summarise(steepest_ascend = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df5,., by=c("ring_ID","date_f"))

# steepest_descend
df9<- data_new_diff %>% 
  group_by(ring_ID, species_en, date_f)  %>% 
  filter(time_since_sunrise_std >= -1.5 & time_since_sunrise_std <= 16.5)%>% 
  filter(time_since_sunrise_std <=  1.5 & time_since_sunrise_std >= 13.5) %>% 
  filter(diff == min(diff,na.rm=TRUE)) %>% 
  summarise(steepest_descend = time_to_rise_std) %>% 
  as.data.frame() %>%  
  select(-species_en) %>% 
  left_join(df8,., by=c("ring_ID","date_f"))


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


