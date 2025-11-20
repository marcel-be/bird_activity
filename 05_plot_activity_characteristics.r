#################################################################################################################################
#### Plotting of Activity Characteristics and Curves for each Individual and day, to check if it makes sense to keep/change values

library(tidyverse)
library(data.table)

path<- "J:/rts/rts_activity/"
path<- "G:/rts_activity/"
path<- "/Users/pandadamda/rts_activity/"

## load data
data_new<- fread(paste0(path,"data/bird_data_storage/models/model_prediction_dates.csv"), stringsAsFactors=T)
data_new$date_f<- as.factor(as.character(data_new$date_f)) 
df_act_charac <- fread(paste0(path,"bird_data_storage/activity_characteristics/activity_characteristics_individual.csv"), stringsAsFactors = T)

## change bird names
data_new$species_en<- gsub("_", " ", data_new$species_en)


###################################################
## Plotting one example (one day of one individual)

df_act_charac_sub<- df_act_charac %>% 
  filter(ring_ID=="6415104" & date_f=="2020-06-11") %>% 
  droplevels() 

data_new %>%
  filter(ring_ID=="6415104" & date_f=="2020-06-11") %>% 
  ggplot(data = ., 
         aes(x = time_to_rise_std, y = mu))+
  geom_ribbon(aes(ymin = ci_lower ,
                  ymax = ci_upper), 
              fill = "grey", color = "grey") +
  geom_line(size = .8) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  #geom_hline(yintercept = df_act_charac_sub$peak_act, linetype = "dashed", color="steelblue", size=1) +
  #geom_text(data= df_act_charac_sub, aes(x=-3.5, label="peak activity", y=peak_act), colour="steelblue", angle=0, vjust = 1.1, size=7)+  
  geom_vline(xintercept = 0, linetype = "dashed", color="orange", size=1)+
  geom_text(aes(x=0, label="sunrise", y=0.2), colour="orange", angle=0, hjust = -0.1, size=7)+
  geom_vline(xintercept = 15.51328, linetype = "dashed", color="navy",  size=1) +
  geom_text(aes(x=15.51328, label="sunset", y= 0.2), colour="navy", angle=0, hjust = -0.1, size=7)+
  geom_vline(xintercept = df_act_charac_sub$steepest_ascend, linetype = "dashed", color="darkgreen",  size=1) +
  geom_text(data= df_act_charac_sub, aes(x=steepest_ascend, label="activity\n onset", y=0.6), colour="darkgreen", angle=0, hjust = 1.1, size=7)+
  geom_vline(xintercept = df_act_charac_sub$steepest_descend, linetype = "dashed", color="red", size=1) +
  geom_text(data= df_act_charac_sub, aes(x=steepest_descend, label="activity\n offset", y=0.25), colour="red", angle=0, hjust = 1.1, size=7)+
  theme_bw(14) +
  xlab("Time since sunrise") + 
  ylab("Activity probability") + 
  ylim(0, 1)+
  ggtitle("European Jay on 11.06.2021 ")+
  theme_bw()+
  theme(text = element_text(size=30),
        axis.text = element_text(face="bold"))
  

ggsave(filename = paste0(path, "output/model_output/diagnostics/" , "jay_example" , ".png"),
       width = 10, height = 7)





#######################################################
## One plot per Individual (several daily curves in one plot) in one pdf
## faster, can be used for presentation/publication 

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

ggsave(filename = paste0(path, "plots/model_output/diagnostics/activity_charcteristics/" , "activity_characteristics_individual" , ".pdf"),
       plot = gridExtra::marrangeGrob(plot_list, nrow=1, ncol=1), 
       width = 15, height = 9)

