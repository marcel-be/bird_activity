## approach one: predict for a random date (used as random factor in model) - don´t use because of variation of outcomes for different dates!
date_fixed = "2020-06-01"
time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
data_new0 <- df_1min %>%  
  expand(nesting(species_en, ring_ID),date_fixed, time_to_rise_std_seq) %>% 
  rename(date_f       = date_fixed,
         time_to_rise_std = time_to_rise_std_seq)

## approach 2: predict for all possible dates (random factor "date_f") and then aggregate mean activity of all dates ("typical date") before plotting --> Mean activity-prediction for each date will be calculated. Prediction takes long! --> Similar to approach 3
time_to_rise_std_seq<- seq(min_set,max_set, length.out=10)
data_new1 <- df_1min %>%  
  expand(nesting(species_en, ring_ID),date_f, time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq)

data_new2 <- data_new1

## approach 3: Set "date_f" in new x-dataframe as "NA". Exclude "s(date_f)" in predict-function and make sure that newdata.guaranteed=TRUE. # Results are very similar to approach 2 
time_to_rise_std_seq<- seq(min_set,max_set, length.out=10)
data_new3 <- df_1min %>%  
  expand(nesting(species_en, ring_ID),time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq) %>% 
  mutate(date_f=NA)



pred <- predict.gam(gam_I, 
                    newdata = data_new1, 
                    #exclude = "s(date_f)",
                    #newdata.guaranteed=TRUE,
                    se = TRUE, 
                    type="link")

pred2 <- predict.gam(gam_I, 
                     newdata = data_new1, 
                     exclude = "s(date_f)",
                     newdata.guaranteed=TRUE,
                     se = TRUE, 
                     type="link")

pred3 <- predict.gam(gam_I, 
                     newdata = data_new3, 
                     exclude = "s(date_f)",
                     newdata.guaranteed=TRUE,
                     se = TRUE, 
                     type="link")


data_new1$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
data_new1$se_min <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
data_new1$se_max <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))

data_new2$mu     <- exp(pred2$fit)/(1+exp(pred2$fit)) # inverse link function (logit scale)
data_new2$se_min <- exp(pred2$fit + 1.96 * pred2$se.fit) / (1+exp(pred2$fit + 1.96 * pred2$se.fit)) # 95% CV
data_new2$se_max <- exp(pred2$fit - 1.96 * pred2$se.fit) / (1+exp(pred2$fit - 1.96 * pred2$se.fit))

data_new3$mu     <- exp(pred3$fit)/(1+exp(pred3$fit)) # inverse link function (logit scale)
data_new3$se_min <- exp(pred3$fit + 1.96 * pred3$se.fit) / (1+exp(pred3$fit + 1.96 * pred3$se.fit)) # 95% CV
data_new3$se_max <- exp(pred3$fit - 1.96 * pred3$se.fit) / (1+exp(pred3$fit - 1.96 * pred3$se.fit))



df<- data_new1 %>% 
  group_by(species_en, ring_ID, time_to_rise_std) %>% 
  summarise_each(funs(mean)) 
df2<- data_new2 %>% 
  group_by(species_en, ring_ID, time_to_rise_std) %>% 
  summarise_each(funs(mean)) 
df3<- data_new3


plot(df$mu ~ df2$mu)
plot(df$mu ~ df3$mu)
plot(df2$mu ~ df3$mu)

plot(df$se_min ~ df2$se_min)
plot(df$se_min ~ df3$se_min)
plot(df2$se_min ~ df3$se_min)

plot(df$se_max ~ df2$se_max)
plot(df$se_max ~ df3$se_max)
plot(df2$se_max ~ df3$se_max)

cor(df$mu , df2$mu)
cor(df$mu , df3$mu)
cor(df2$mu , df3$mu)

cor(df$se_max, df2$se_max)
cor(df$se_max, df3$se_max)
cor(df2$se_max, df3$se_max)

cor(df$se_min, df2$se_min)
cor(df$se_min, df3$se_min)
cor(df2$se_min, df3$se_min)