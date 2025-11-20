########################################################################################################################## 2. Descriptive statistics ####

library(tidyverse)
library(data.table)


path<- "/Users/pandadamda/rts_activity/"
path<- "E:/Uni_Arbeit/rts_activity/"


df<- fread(paste0(path, "data/bird_data_storage/tags_1_min_for_analysis.csv"), key = "ID")


### Some descriptives
mean(df$length_act, na.rm=T)
max(df$length_act, na.rm=T)
min(df$length_act, na.rm=T)

mean(df$steepest_ascend, na.rm=T)
max(df$steepest_descend, na.rm=T)
min(df$steepest_ascend, na.rm=T)

#...