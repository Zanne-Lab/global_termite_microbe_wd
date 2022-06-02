library(tidyverse)
library(easystats)
library(ggeffects)
library(patchwork)
#read data and subset
mm <- read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T", "C") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara") %>%
  mutate(date_diff_years = date_diff/365)


#get site level decay
mm %>% 
  group_by(site, termite_exposure) %>% 
  summarise(temp = temp[1], 
            prec = prec[1], 
            forestcover100 = forestcover100[1], 
            alt = alt[1],
            date_diff_years=date_diff_years[1],
            lat = lat[1],
            N_pc = N_pc[1],
            C_pc = C_pc[1],
            n = n(), 
            k_value= mean(k_value),
            #kval = log(k_value),
            discovered = mean(termite_exposure)) %>%
  ungroup() %>% 
  mutate(termite_exposure = factor(termite_exposure), 
         abs_lat = abs(lat), 
         prec_m = prec/1000, 
         k_value_tr = log(k_value))->mmm


#Termite decay across absolute latitude
modlatt<-lm(k_value_tr ~abs_lat,
            data=filter(mmm,termite_exposure ==1)
)

#Microbial decay across absolute latitude
modlatm<-lm(k_value_tr ~abs_lat,
            data=filter(mmm,termite_exposure ==0)
)

#Termite decay across elevation
modaltt<-lm(k_value_tr ~alt,
            data=filter(mmm,termite_exposure ==1)
)

#Microbial decay across elevation
modaltm<-lm(k_value_tr ~alt,
            data=filter(mmm,termite_exposure ==0)
)

#Termite decay across temperature
modtempt<-lm(k_value_tr ~temp,
             data=filter(mmm,termite_exposure ==1)
)

#Microbial decay across temperature
modtempm<-lm(k_value_tr ~temp,
             data=filter(mmm,termite_exposure ==0)
)

#Termite decay across precipitation
modprect<-lm(k_value_tr ~prec,
             data=filter(mmm,termite_exposure ==1)
)
#Microbial decay across precipitation
modprecm<-lm(k_value_tr ~prec,
             data=filter(mmm,termite_exposure ==0)
)


table <-bind_rows(model_parameters(modlatt),
                  model_parameters(modlatm),
                  model_parameters(modaltt),
                  model_parameters(modaltm),
                  model_parameters(modtempt),
                  model_parameters(modtempm),
                  model_parameters(modprect),
                  model_parameters(modprecm))

print_md(table)


bind_rows(r2(modlatt, exponentiate = TRUE)[[2]],
          r2(modlatm, exponentiate = TRUE)[[2]],
          r2(modaltt, exponentiate = TRUE)[[2]],
          r2(modaltm, exponentiate = TRUE)[[2]],
          r2(modtempt,exponentiate = TRUE)[[2]],
          r2(modtempm, exponentiate = TRUE)[[2]],
          r2(modprect, exponentiate = TRUE)[[2]],
          r2(modprecm,exponentiate = TRUE)[[2]])

#Supplementary table S4 
#multivariate model for microbe (termite undiscovered) wood decay (k) 
#versus climatic sensitivities

m1w <- lm(k_value_tr ~ temp*prec_m ,
          data=filter(mmm,termite_exposure ==0), weights=n)
model_parameters(m1w)
r2(m1w)
summary(m1w)

#multivariate model for termite discovered wood decay (k) 
#versus climatic sensitivities
m2w <- lm(k_value_tr ~ temp*prec_m ,
          data=filter(mmm,termite_exposure ==1), weights=n)
model_parameters(m2w)
r2(m2w)
summary(m2w)

#comparison with termite treatment and temperature
mod_temp_trt <- lm(k_value_tr ~temp*termite_exposure,
                   data=mmm)

model_parameters(mod_temp_trt)


### bivariate models for decay (k) versus key spatial versus climatic variables
#and wood chemistry


#Termite decay across absolute latitude
modlatt_NC<-lm(k_value_tr ~abs_lat+N_pc+C_pc,
            data=filter(mmm,termite_exposure ==1)
)

#Microbial decay across absolute latitude
modlatm_NC<-lm(k_value_tr ~abs_lat+N_pc+C_pc,
            data=filter(mmm,termite_exposure ==0)
)

#Termite decay across elevation
modaltt_NC<-lm(k_value_tr ~alt+N_pc+C_pc,
            data=filter(mmm,termite_exposure ==1)
)

#Microbial decay across elevation
modaltm_NC<-lm(k_value_tr ~alt+N_pc+C_pc,
            data=filter(mmm,termite_exposure ==0)
)

#Termite decay across temperature
modtempt_NC<-lm(k_value_tr ~temp+N_pc+C_pc,
             data=filter(mmm,termite_exposure ==1)
)

#Microbial decay across temperature
modtempm_NC<-lm(k_value_tr ~temp+N_pc+C_pc,
             data=filter(mmm,termite_exposure ==0)
)

#Termite decay across precipitation
modprect_NC<-lm(k_value_tr ~prec+N_pc+C_pc,
             data=filter(mmm,termite_exposure ==1)
)
#Microbial decay across precipitation
modprecm_NC<-lm(k_value_tr ~prec+N_pc+C_pc,
             data=filter(mmm,termite_exposure ==0)
)

#parameter values
table <-bind_rows(model_parameters(modlatt_NC),
                  model_parameters(modlatm_NC),
                  model_parameters(modaltt_NC),
                  model_parameters(modaltm_NC),
                  model_parameters(modtempt_NC),
                  model_parameters(modtempm_NC),
                  model_parameters(modprect_NC),
                  model_parameters(modprecm_NC))

print_md(table)
#r2 values

bind_rows(r2(modlatt_NC, exponentiate = TRUE)[[2]],
          r2(modlatm_NC, exponentiate = TRUE)[[2]],
          r2(modaltt_NC, exponentiate = TRUE)[[2]],
          r2(modaltm_NC, exponentiate = TRUE)[[2]],
          r2(modtempt_NC,exponentiate = TRUE)[[2]],
          r2(modtempm_NC, exponentiate = TRUE)[[2]],
          r2(modprect_NC, exponentiate = TRUE)[[2]],
          r2(modprecm_NC,exponentiate = TRUE)[[2]])

#multivariate model for microbe (termite undiscovered) wood decay (k) 
#versus climatic sensitivities and wood chemistry

m1w_NC <- lm(k_value_tr ~ temp*prec_m+N_pc+C_pc ,
          data=filter(mmm,termite_exposure ==0), weights=n)
model_parameters(m1w_NC)
r2(m1w_NC)
summary(m1w_NC)

#multivariate model for microbe (termite undiscovered) wood decay (k) versus climatic sensitivities and wood chemistry
m2w_NC <- lm(k_value_tr ~ temp*prec_m+N_pc+C_pc ,
             data=filter(mmm,termite_exposure ==1), weights=n)
model_parameters(m2w_NC)
r2(m2w_NC)
summary(m2w_NC)
nrow(filter(mmm,termite_exposure ==1))


## check for mesh hole effect on undiscovered baits (no effect observed)
mm %>% 
  filter(termite_exposure == 0) %>% 
  group_by(site, trt) %>% 
  summarise(temp = temp[1], 
            prec = prec[1], 
            forestcover100 = forestcover100[1], 
            alt = alt[1],
            date_diff_years=date_diff_years[1],
            lat = lat[1],
            N_pc = N_pc[1],
            C_pc = C_pc[1],
            n = n(), 
            k_value= mean(k_value),
            #kval = log(k_value),
            discovered = mean(termite_exposure)) %>%
  ungroup() %>% 
  mutate(trt = factor(trt), 
         abs_lat = abs(lat), 
         prec_m = prec/1000, 
         k_value_tr = log(k_value+abs(min(k_value))+0.001))->site_means_u
hist(site_means_u$k_value)
hist(site_means_u$k_value_tr)

m1w_u <- lm(k_value_tr ~ trt*temp*prec +
              trt*alt + trt*abs_lat, data=site_means_u, weights=n)
summary(m1w_u)[['r.squared']]
car::qqPlot(m1w_u)
car::residualPlot(m1w_u)
car::Anova(m1w_u, type='II')
