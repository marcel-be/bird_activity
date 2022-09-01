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
        s(time_to_rise_std, m=2,  bs="cc", k=50) +
        s(date_f, bs="re"),  # k equals the number of levels of grouping variable
      method ="fREML", 
      family ="binomial",
      discrete = T, 
      knots=list(time_to_rise_std=c(min_set, max_set)),
      data = df)

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

time_to_rise_std_seq<- seq(min_set,max_set, length.out=1000)
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

steepest_ascend  <- data_new[which(diff(data_new$mu) == max(diff(data_new$mu))),]
steepest_descend <- data_new[which(diff(data_new$mu) == min(diff(data_new$mu))),]

data_new %>% 
  group_by(time_to_rise_std) %>% 
  summarise_each(funs(mean)) %>% 
  ggplot(data = ., 
         aes(x = time_to_rise_std, y = mu))+
  geom_ribbon(aes(ymin = se_min ,
                  ymax = se_max), 
              fill = "grey", color = "grey") +
  geom_line(size = .8) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = steepest_ascend$time_to_rise_std, linetype = "solid", col="green") +
  geom_vline(xintercept = steepest_descend$time_to_rise_std, linetype = "solid", col="red") +
  #geom_vline(xintercept = slope_max$time_to_rise_std[1],  color = "blue", size=0.5)+
  #geom_vline(xintercept = slope_max$time_to_rise_std[2],  color = "red", size=0.5)+
  theme_bw(14) +
  xlab("Time since sunrise") + 
  ylab("Activity probability \n") + 
  ylim(0, 1)





##########################################################################################################
## 3. Activity curves for all Tags/Individuals
# the loop creates one ativity curve for each Tag and attaches the daily activity-plots. The output is used to check for raw-data-issues

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


##########################################################################################################
## 3. Determine activity characterstics of all Individuals

