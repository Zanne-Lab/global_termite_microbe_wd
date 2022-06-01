library(tidyverse)
library(ggeffects)
library(lme4)
library(easystats)
library(boot)

mm <- read_csv("data/decomp_with_covar.csv")  %>%
  mutate(abs_lat = abs(lat),
         date_diff_years = date_diff/365) %>%
  filter(trt %in% c("T") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara")
####Discovery models Table S1
#discovery over temperature model
mod_temp<-glm(termite_exposure ~ temp,
            data=filter(mm, trt == 'T'), offset=date_diff_years, family=binomial(link = "logit"),
            control=list(epsilon = 1e-10, maxit = 100))
temp <- model_parameters(mod_temp,exponentiate = TRUE)



#discovery over precipitation model
modpre<-glm(termite_exposure~prec,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
pre <- model_parameters(modpre,exponentiate = TRUE)

#discovery~absolute absolute latitude model
modlat<-glm(termite_exposure~abs_lat,
             data=mm,offset=date_diff_years,family = binomial(link = "logit"),
             control = list(epsilon = 1e-10,maxit = 100))
lat <- model_parameters(modlat,exponentiate = TRUE)

#discovery~elevation model
modalt<-glm(termite_exposure~alt,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
alt <- model_parameters(modalt,exponentiate = TRUE)

print_md(bind_rows(temp,pre, lat, alt))
r2_mcfadden(mod_temp)
r2_mcfadden(modpre)
r2_mcfadden(modlat)
r2_mcfadden(modalt)
kableExtra::kable(bind_rows(temp,pre, lat, alt),"pipe")
####

#multivariate model for probably of termite discovery versus temperature and precipitation
modinteraction<-glm(termite_exposure ~ prec*temp,
                    data=filter(mm, trt == 'T'), offset=date_diff_years, family=binomial(link = "logit"),
                    control=list(epsilon = 1e-100, maxit = 10000))
nrow(modinteraction$data)
summary(modinteraction)
print_md(model_parameters(modinteraction,exponentiate = TRUE))

ggpredict(model = modinteraction, "temp [21.3]")

v <- seq(from = 0,to = 30,by = 0.5)
v2 <- c(250, 2000, 2700)
v3 <- c(7, 25)
l <- ggpredict(modinteraction,terms = c("temp [v3]", "prec [v2]"))
print_md(l)
r2_mcfadden(modinteraction)


#### Models with wood chemistry as covariables 

#Bivariate models
####Discovery models Table S6
#discovery over temperature model
mod_temp_NC<-glm(termite_exposure ~ temp+N_pc+C_pc,
              data=filter(mm, trt == 'T'), offset=date_diff_years, family=binomial(link = "logit"),
              control=list(epsilon = 1e-10, maxit = 100))
temp <- model_parameters(mod_temp_NC,exponentiate = TRUE)


#discovery over precipitation model
modpre_NC<-glm(termite_exposure~prec+N_pc+C_pc,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
pre <- model_parameters(modpre_NC,exponentiate = TRUE)

#discovery~absolute absolute latitude model
modlat_NP<-glm(termite_exposure~abs_lat+N_pc+C_pc,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
lat <- model_parameters(modlat_NP,exponentiate = TRUE)

#discovery~elevation model
modalt_NC<-glm(termite_exposure~alt+N_pc+C_pc,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
alt <- model_parameters(modalt_NC,exponentiate = TRUE)

print_md(bind_rows(temp,pre, lat, alt))
r2_mcfadden(mod_temp_NC)
r2_mcfadden(modpre_NC)
r2_mcfadden(modlat_NP)
r2_mcfadden(modalt_NC)
kableExtra::kable(bind_rows(temp,pre, lat, alt),"pipe")


#multivariate model for probably of termite discovery versus climatic sensitivities and wood chemistry
modinteraction_stoic<-glm(termite_exposure ~ prec*temp+N_pc+C_pc,
                          data=filter(mm, trt == 'T'), offset=date_diff_years, family=binomial(link = "logit"),
                          control=list(epsilon = 1e-100, maxit = 10000))

print_md(model_parameters(modinteraction_stoic,exponentiate = TRUE))
check_model(modinteraction_stoic)

nrow(modinteraction_stoic$data)
r2_mcfadden(modinteraction_stoic)
