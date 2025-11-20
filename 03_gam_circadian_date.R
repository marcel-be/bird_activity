library(tidyverse)
library(data.table)
library(mgcv)
library(gratia)


path<- "/Users/pandadamda/rts_activity/"
path<- "E:/Uni_Arbeit/rts_activity/"

##################################################################################################################
#### 1. Data import and preperation:

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
         interval   = as_hms(floor_date(timestamp_CET, unit="10minutes"))) %>% 
  group_by(ID, ring_ID, species_en,year_f,brood_patch,sex, month,week,ydate, date_f,hour,interval) %>% 
  summarise(n_intervals=length(activity),
            n_active=length(activity[activity==1]),
            n_passive=length(activity[activity==0]),
            time_of_day=mean(time_of_day),
            time_to_rise_std = mean(time_to_rise_std)) %>% 
  mutate(active_prop =  n_active/n_intervals)



##################################################################################################################
#### 2. Models for circadian rhythm of bird activity (Daily rhythm) ####
# the loop creates one activity-curve for each Tag (from GAM) and attaches the daily activity-plots (from raw data). The output is used to check for raw-data-issues (seeing some weird curve pattern for certain days that should be removed)

plot_list<- list()
n_ID<- nlevels(df_1min$ID)
a<- seq(1,(2*n_ID)-1, length.out=n_ID) # used as a counter for the plot-list
b<- seq(2, 2*n_ID   , length.out=n_ID) # used as a counter for the plot-list
data_new_final<- data.frame() # new dataframe for the predictions from the model 

for(i in 1:nlevels(df_1min$ring_ID)){
  
  print(i)
  print(levels(df_1min$ring_ID)[i])
  
  df<-  df_1min %>% 
    filter(ring_ID==levels(ring_ID)[i]) %>% 
    droplevels()
  
  min_set  <- -5
  max_set  <- 18
  #min_set  <- min(df$time_to_rise_std)
  #max_set  <- max(df$time_to_rise_std)
  
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
  
  time_to_rise_std_seq<- seq(min_set,max_set, length.out=1000)
  time_to_rise_std_seq[1] <- 0  #  set minimal value as Zero. Important to get activity predictions for the exact time of sunrise
  data_new <- df %>% 
    expand(nesting(species_en,ring_ID, ID, date_f), time_to_rise_std_seq) %>% 
    rename(time_to_rise_std = time_to_rise_std_seq)

  ## get predicted values
  pred <- predict.gam(gam, 
                      newdata = data_new, 
                      #exclude = "s(date_f)",
                      #newdata.guaranteed=TRUE,
                      se = TRUE, 
                      type="link")
  data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
  data_new$ci_lower <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
  data_new$ci_upper <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))
  
  data_new_final<- rbind(data_new_final, data_new)

  name<- paste0(as.character(df$ring_ID[1]), # for plotting
                " | ", as.character(df$species_en[1]),
                " | sex=", as.character(df$sex[1]))
  
  p<- data_new %>% 
    ggplot(data = ., 
           aes(x = time_to_rise_std, y = mu, group=date_f, color=date_f))+
    geom_ribbon(aes(ymin = ci_lower ,
                    ymax = ci_upper), 
                fill = "grey", color = "grey") +
    geom_line(size = .8) + 
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    theme_bw(14) +
    xlab("Time since sunrise") + 
    ylab("Activity probability \n") + 
    ylim(0, 1)+
    ggtitle(name)
  
  plot_list[[a[i]]] <-p
  
  df <- df_10min %>%  
    filter(ring_ID == levels(ring_ID)[i]) %>% 
    droplevels()
  
  name<- paste0(as.character(df$species_en[1]),
                " | brood patch=", as.character(df$brood_patch[1]),
                " | year=", as.character(df$year_f[1]),
                " | sex=", as.character(df$sex[1]))
  
  p<-  ggplot(df, aes(y = n_active/n_intervals, x = time_to_rise_std)) +
    geom_point(alpha = .20) + 
    geom_smooth(method = 'gam', 
                formula = y ~ s(x, bs = "tp", k = 25)) +
    facet_wrap(~ date_f) +
    theme_bw()+
    ggtitle(paste0(levels(df$ID)[1]," ", name))
  
  plot_list[[b[i]]] <-p
}

ggsave(filename = paste0(path, "plots/model_output/diagnostics/" , "curve_by_individual" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)
# these plots can be used for further checking of the data. Fishy curves can be removed or taken care of in the next scripts where the activity characteristics (like peak activity etc) will be calculated

data_new_final<- data_new_final %>% 
  filter(ring_ID!="90850086") %>% 
  droplevels() # remove

## save model predictions
fwrite(data_new_final, paste0(path,"bird_data_storage/models/model_prediction_dates.csv"))

