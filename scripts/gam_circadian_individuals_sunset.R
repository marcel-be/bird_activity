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

#path<- "J:/rts/rts_activity/"
path<- "G:/rts_activity/"

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
min_set  <- min(df_1min$time_to_set_std)
max_set  <- max(df_1min$time_to_set_std)
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
            time_to_set_std = mean(time_to_set_std)) %>% 
  mutate(active_prop =  n_active/n_intervals)

df_10min <- df_10min[seq(1, nrow(df_10min),5), ] # reduce data for faster plotting


##########################################################################################################
## 2. model for one individual

df<-  df_1min %>% 
  filter(ring_ID==levels(df_1min$ring_ID)[1]) %>% 
  droplevels()

min_set  <- min(df$time_to_set_std)
max_set  <- max(df$time_to_set_std)

gam <- bam(activity ~ 
        s(time_to_set_std, by=date_f, m=2,  bs="cc", k=50)+
        s(date_f, bs="re"),  
      method ="fREML", 
      family ="binomial",
      discrete = T, 
      knots=list(time_to_set_std=c(min_set, max_set)),
      data = df)

rho<- start_value_rho(gam, plot=F)

gam <- df %>% 
  bam(activity ~ 
        s(time_to_set_std, by=date_f, m=2,  bs="cc", k=30)+
        s(date_f, bs="re"),   
      method ="fREML", 
      family ="binomial",
      discrete = T, 
      knots=list(time_to_set_std=c(min_set, max_set)),
      rho= rho,
      data = .)

par(mfrow=c(2,2))
gam.check(gam)
par(mfrow=c(1,1))

time_to_set_std_seq<- seq(min_set,max_set, length.out=1000)
data_new <- df%>% 
  expand(nesting(ring_ID, date_f), time_to_set_std_seq) %>% 
  rename(time_to_set_std = time_to_set_std_seq)

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
         aes(x = time_to_set_std, y = mu, group=date_f, color=date_f))+
  geom_ribbon(aes(ymin = se_min ,
                  ymax = se_max), 
              fill = "grey", color = "grey") +
  geom_line(size = .8) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #geom_vline(xintercept = slope_max$time_to_set_std[1],  color = "blue", size=0.5)+
  #geom_vline(xintercept = slope_max$time_to_set_std[2],  color = "red", size=0.5)+
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

for(i in 1:nlevels(df_1min$ring_ID)){
  
  print(levels(df_1min$ring_ID)[i])
  
  df<-  df_1min %>% 
    filter(ring_ID==levels(ring_ID)[i]) %>% 
    droplevels()
  
  #min_set  <- -5
  #max_set  <- 18
  min_set  <- min(df$time_to_set_std)
  max_set  <- max(df$time_to_set_std)
  
  gam <- bam(activity ~ 
               s(time_to_set_std, by=date_f, m=2,  bs="cc", k=35)+
               s(date_f, bs="re"),  
             method ="fREML", 
             family ="binomial",
             discrete = T, 
             knots=list(time_to_set_std=c(min_set, max_set)),
             data = df)
  
  rho<- start_value_rho(gam, plot=F)
  
  gam <- df %>% 
    bam(activity ~ 
          s(time_to_set_std, by=date_f, m=2,  bs="cc", k=30)+
          s(date_f, bs="re"),   
        method ="fREML", 
        family ="binomial",
        discrete = T, 
        knots=list(time_to_set_std=c(min_set, max_set)),
        rho= rho,
        data = .)
  
  time_to_set_std_seq<- seq(min_set,max_set, length.out=1000)
  time_to_set_std_seq[1] <- 0  # get activity predictions for time of sunrise
  data_new <- df %>% 
    expand(nesting(species_en,ring_ID, date_f), time_to_set_std_seq) %>% 
    rename(time_to_set_std = time_to_set_std_seq)

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

  name<- paste0(as.character(df$ring_ID[1]),
                " | ", as.character(df$species_en[1]),
                #" | brood patch=", as.character(df$brood_patch[1]),
                #" | year=", as.character(df$year_f[1]),
                " | sex=", as.character(df$sex[1]))
  
  p<- data_new %>% 
    ggplot(data = ., 
           aes(x = time_to_set_std, y = mu, group=date_f, color=date_f))+
    geom_ribbon(aes(ymin = ci_lower ,
                    ymax = ci_upper), 
                fill = "grey", color = "grey") +
    geom_line(size = .8) + 
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    #geom_vline(xintercept = slope_max$time_to_set_std[1],  color = "blue", size=0.5)+
    #geom_vline(xintercept = slope_max$time_to_set_std[2],  color = "red", size=0.5)+
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
  
  p<-  ggplot(df, aes(y = n_active/n_intervals, x = time_to_set_std)) +
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

ggsave(filename = paste0(path, "plots/model_output/diagnostics/" , "curve_by_individual_sunset" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)

## save model predictions
fwrite(data_new_final, paste0(path,"bird_data_storage/models/model_prediction_individuals_sunset.csv"))



#################################################################################################################################
#### 4. Plot Model Results (one curve per individual)

# Exclude: 
# 90619209 & 2019−07−08
# 6415106 & 2019−07−28
# 90619323 & 2021−04−21
# V188907 & 2020−06−07

df_10min<- fread(paste0(path, "bird_data_storage/tags_10_min_for_plotting.csv"), stringsAsFactors = T)
data_new<- fread(paste0(path,"bird_data_storage/models/model_prediction_individuals_sunset.csv"), stringsAsFactors=T)
data_new$date_f<- as.factor(as.character(data_new$date_f)) 
df_1min<- fread(paste0(path, "bird_data_storage/tags_1_min_for_analysis.csv"), stringsAsFactors = T)
df_1min$date_f <- as.factor(as.character(df_1min$date_f))

data_new_agg<- data_new %>% 
  group_by(species_en, time_to_set_std) %>% 
  summarise(sd = sd(mu),
            mu = mean(mu),
            n  = n()) %>% 
  mutate(se= sd / sqrt(n),
         ci_lower = mu - 1.96 * se,
         ci_upper = mu + 1.96 * se,
         ring_ID=NA)

data_new_agg<- data_new %>% 
  group_by(species_en, time_to_set_std) %>% 
  summarise_each(funs(mean))
# which "data_new_agg" to use?
  
  
# at individual level (one curve per individual)
p<- data_new %>% 
  group_by(species_en, ring_ID, time_to_set_std) %>% 
  summarise_each(funs(mean)) %>% 
  ggplot(data = ., 
         aes(x = time_to_set_std, y = mu, 
             color = ring_ID, group = ring_ID))+
  #geom_ribbon(data=data_new_agg, 
 #             aes(ymin = ci_lower ,
 #                 ymax = ci_upper), 
#            fill = "grey", color = "grey", alpha=0.6)+
  #geom_ribbon(aes(ymin = ci_lower ,
  #                ymax = ci_upper), 
  #            fill = "grey", color = "grey", alpha=0.6) +
  #geom_point(data=df_10min, alpha = .1, 
  #           aes(x = time_to_set_std, y = n_active/n_intervals)) +
  geom_line(size = .7) +
 # geom_line(data =data_new_agg, 
#            aes(x = time_to_set_std, y = mu), 
 #           color="black", size=1.1, alpha=0.5)+
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #scale_color_wsj() +
  theme_bw(14) +
  xlab("Time since sunrise") + 
  ylab("Activity probability \n") + 
  ylim(0, 1) +
  #theme(legend.position = c(0.85,0.88))+
  theme(legend.position = "none") +
  facet_wrap(~species_en)
ggsave(filename = paste0(path, "plots/model_output/date_specific_models/" , "circadian_ID_sunset" , ".png"),
       plot=p, width = 15, height = 9)


