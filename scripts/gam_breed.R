table(df_1min$brood_patch)

df_1min_short<- df_1min %>% 
    filter(brood_patch=="yes" | brood_patch=="no") %>% 
  droplevels()

gam_breed<- bam(activity ~ brood_patch +
              s(time_to_rise_std, by=brood_patch, m=2,  bs="cc", k=50) + # no intercept-per-ID when using the "by"-argument --> specify separately (fixed or random)
              s(species_en, bs="re") +
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
summary(gam_breed)


time_to_rise_std_seq<- seq(min_set,max_set, length.out=100)
data_new <- df_1min_short %>%  
  expand(nesting(species_en,ring_ID,brood_patch),date_f, time_to_rise_std_seq) %>% 
  rename(time_to_rise_std = time_to_rise_std_seq)

pred <- predict(gam_breed, 
                newdata = data_new, 
                #exclude_terms = c(s(ID),s(ID, time_to_rise_std), s(date_f)),
                se = TRUE, 
                type="link")

data_new$mu     <- exp(pred$fit)/(1+exp(pred$fit)) # inverse link function (logit scale)
data_new$se_min <- exp(pred$fit + 1.96 * pred$se.fit) / (1+exp(pred$fit + 1.96 * pred$se.fit)) # 95% CV
data_new$se_max <- exp(pred$fit - 1.96 * pred$se.fit) / (1+exp(pred$fit - 1.96 * pred$se.fit))


data_new %>% 
  group_by(ring_ID, brood_patch, time_to_rise_std) %>% 
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
  facet_wrap(~brood_patch)+
  theme(legend.position = "none")
ggsave(filename = paste0(path, "plots/model_output/" , "circadian_breeding_state_1" , ".png"),
       width = 15, height = 9)

data_new %>% 
  group_by(species_en, brood_patch, time_to_rise_std) %>% 
  summarise_each(funs(mean)) %>% 
  ggplot(data = ., 
         aes(x = time_to_rise_std, y = mu, 
             color = brood_patch, group = brood_patch))+
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
  facet_wrap(~species_en)
ggsave(filename = paste0(path, "plots/model_output/" , "circadian_breeding_state_2" , ".png"),
       width = 15, height = 9)
