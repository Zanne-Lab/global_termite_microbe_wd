library(tidyverse)

mm <- read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T", "C") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara")


# Microbial decomposition wood blocks
mm %>% 
  filter(termite_exposure == 0) %>% 
  group_by(site) %>% 
  summarise(mean_k_value = mean(k_value), 
            mean_temp = mean(temp), 
            termite_exposure = 0, 
            n = n()) -> df_m

# Termide discoverd wood blocks
mm %>% 
  filter(termite_exposure == 1) %>% 
  group_by(site) %>% 
  summarise(mean_k_value = mean(k_value), 
            mean_temp = mean(temp), 
            termite_exposure = 1, 
            n = n()) -> df_t

df <- bind_rows(df_m, df_t)

ggplot(df,aes(x=mean_temp,y=mean_k_value,col=as.factor(termite_exposure)))+geom_point()+geom_smooth(method="lm") + scale_y_log10()+theme_bw()
#Microbial Q10
mm<-lm(log(mean_k_value)~mean_temp, data=df_m, weights=n)
#slope and CI
exp(coef(mm)[2]*10)
exp(confint(mm)[2, ]*10)

#Termite Q10 value
mt<-lm(log(mean_k_value)~mean_temp, data=df_t, weights=n)
exp(coef(mt)[2]*10)
exp(confint(mt)[2, ]*10)

#Average Q10 all wood blocks
ma<-lm(log(mean_k_value)~mean_temp, data=df, weights=n)
exp(coef(ma)[2]*10)
exp(confint(ma)[2, ]*10)

#comparison of termites and microbial q10
exp(coef(mt)[2]*10)/exp(coef(mm)[2]*10)
