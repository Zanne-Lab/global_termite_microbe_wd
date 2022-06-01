library(tidyverse)
library(easystats)
library(ggeffects)
library(patchwork)
library(ggnewscale)
#load data
a_dat <- read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T", "C") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara") -> mm

mm$date_diff_years<-mm$date_diff/365
#sumarise discovery data at the site level
mm %>% 
  group_by(site, termite_exposure) %>% 
  summarise(temp = temp[1], 
            prec = prec[1], 
            forestcover100 = forestcover100[1], 
            alt = alt[1],
            date_diff_years=date_diff_years[1],
            lat = lat[1], 
            n = n(), 
            k_value= mean(k_value),
            #kval = log(k_value),
            discovered = mean(termite_exposure)) %>%
  ungroup() %>% 
  mutate(termite_exposure = factor(termite_exposure), 
         abs_lat = abs(lat), 
         prec_m = prec/1000, 
         k_value_tr = log(k_value))->mmm
#generate absolute latitude model for termite exposed wood blocks
modlatt<-lm(k_value_tr ~abs_lat,
            data=filter(mmm,termite_exposure ==1)
)
#generate model prediction 
v <- c(0.6,2,5, 10, 12, 15, 17,20,22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50, 52)
modlatt_pred <- ggpredict(modlatt, c("abs_lat [v]"))
modlatm1<-lm(k_value_tr ~1,
             data=filter(mmm,termite_exposure ==0)
)

#generate absolute latitude model for undiscovered wood blocks
modlatm<-lm(k_value_tr ~abs_lat,
            data=filter(mmm,termite_exposure ==0)
)


#generate prediction for undiscovered wood blocks
v <- c(0.6,2,5, 10, 12, 15, 17,20,22, 25, 27, 30, 32, 35, 37, 40, 42, 45, 47, 50, 52)
modlatm_pred <- ggpredict(modlatm, c("abs_lat [v]"))

#Figure S2 b absolute laitude and decay 
b <-ggplot(mmm,aes(x=abs_lat,y=k_value, colour=termite_exposure))+
  geom_point()+
  geom_line(aes(x=x, y=exp(predicted)), colour='#317A22',data=modlatm_pred, 
            inherit.aes=FALSE) + 
  geom_ribbon(aes(x=x, ymin=exp(conf.low), ymax=exp(conf.high)), data=modlatm_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#317A22')+
  scale_color_manual(name = "Termite discovery",
                     breaks = c(0,1),
                     values = c("#317A22", "#DCBB50"),
                     labels = c("Undiscovered", "Discovered"),
                     guide = "legend")+
  labs(x = expression(`Absolute Latitude `(degree)), y = 'k (per year)')+
  theme_classic()+ theme(legend.position = "none", axis.title.y = element_blank(),
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))


#Generate elevation model for discovered blocks
modaltt<-lm(k_value_tr ~alt,
            data=filter(mmm,termite_exposure ==1)
)

modaltt2<-lm(k_value_tr ~poly(alt,2),
             data=filter(mmm,termite_exposure ==1)
)

#generate model predtion
v <- c(0,400,800,1200,1600,2000,2400,3200)
modaltt_pred <- ggpredict(modaltt, c("alt [v]"))

#Generate elevation model for undiscovered blocks
modaltm<-lm(k_value_tr ~alt,
            data=filter(mmm,termite_exposure ==0)
)
#generate model predictions
modaltm_pred <- ggpredict(modaltm, c("alt [v]"))

#Figure S2 A elevation and decay 
a <- ggplot(mmm,aes(x=alt,y=k_value, colour = termite_exposure))+
  geom_point()+
  geom_line(aes(x=x, y=exp(predicted)), colour='#317A22',data=modaltm_pred, 
            inherit.aes=FALSE) + 
  geom_ribbon(aes(x=x, ymin=exp(conf.low), ymax=exp(conf.high)), data=modaltm_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#317A22')+
  geom_line(aes(x=x, y=exp(predicted)), colour='#DCBB50', data=modaltt_pred, 
            inherit.aes=FALSE) + 
  geom_ribbon(aes(x=x, ymin=exp(conf.low), ymax=exp(conf.high)), data=modaltt_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#DCBB50')+
  scale_color_manual(name = "Termite discovery",
                     breaks = c(0,1),
                     values = c("#317A22", "#DCBB50"),
                     labels = c("Undiscovered", "Discovered"),
                     guide = "legend")+
  labs(x = "Elevation (m)", y = 'k (per year)')+
  theme_classic()+ theme(legend.position = "none",
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))




#Generate temperature model for discovered blocks
modtempt<-lm(k_value_tr ~temp,
             data=filter(mmm,termite_exposure ==1))
#generate model predictions
v <- c(0,4,6,10,12,14,16,18,22,24,26.5)
modtempt_pred <- ggpredict(modtempt, c("temp [v]"))

#Generate temperature model for undiscovered blocks
modtempm<-lm(k_value_tr ~temp,
             data=filter(mmm,termite_exposure ==0)
)
#generate model predictions
modtempm_pred <- ggpredict(modtempm, c("temp [v]"))


#Figure S2 C elevation and decay 
c <- ggplot(mmm,aes(x=temp,y=k_value, colour = termite_exposure))+
  geom_point()+
  geom_line(aes(x=x, y=exp(predicted)), colour='#317A22', data=modtempm_pred, 
            inherit.aes=FALSE) + 
  geom_ribbon(aes(x=x, ymin=exp(conf.low), ymax=exp(conf.high)), data=modtempm_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#317A22')+
  geom_line(aes(x=x, y=exp(predicted)), colour='#DCBB50', data=modtempt_pred, 
            inherit.aes=FALSE) + 
  geom_ribbon(aes(x=x, ymin=exp(conf.low), ymax=exp(conf.high)), data=modtempt_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#DCBB50')+
  scale_color_manual(name = "Termite discovery",
                     breaks = c(0,1),
                     values = c("#317A22", "#DCBB50"),
                     labels = c("Undiscovered", "Discovered"),
                     guide = "legend")+
  labs(x = expression(`MAT `(degree * C)), y = 'k (per year)')+
  theme_classic()+ theme(legend.position = "none",
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))


##Generate precipitation model for discovered blocks
modprect<-lm(k_value_tr ~prec,
             data=filter(mmm,termite_exposure ==1)
)
#model prediciton
v <- c(0,4,6,10,12,14,16,18,22,24,26.5)
modprect_pred <- ggpredict(modprect, c("prec"))

##Generate precipitation model for undiscovered blocks
modprecm<-lm(k_value_tr ~prec,
             data=filter(mmm,termite_exposure ==0)
)
#model prediction
modprecm_pred <- ggpredict(modprecm, c("prec"))

#Figure S2 D precipitation and decay 
d <- ggplot(mmm,aes(x=prec,y=k_value, colour = termite_exposure))+
  geom_point()+
  geom_line(aes(x=x, y=exp(predicted)), colour='#317A22', data=modprecm_pred, 
            inherit.aes=FALSE) + 
  geom_ribbon(aes(x=x, ymin=exp(conf.low), ymax=exp(conf.high)), data=modprecm_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#317A22')+
  #geom_line(aes(x=x, y=exp(predicted), colour='#DCBB50'), data=modprect_pred, 
  #         inherit.aes=FALSE) + 
  #eom_ribbon(aes(x=x, ymin=exp(conf.low), ymax=exp(conf.high)), data=modprect_pred, 
  #           inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#DCBB50')+
  scale_color_manual(name = "Termite discovery",
                     breaks = c(0,1),
                     values = c("#317A22", "#DCBB50"),
                     labels = c("Undiscovered", "Discovered"),
                     guide = "legend")+
  labs(x = "MAP (mm)", y = 'k (per year)')+
  theme_classic()+ theme(axis.title.y=element_blank(),
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))

#join panels
(a+b)/(c+d)+patchwork::plot_annotation(tag_levels = "A")+theme(title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                                                               legend.key.size = unit(8,units = "points"))
#save plots
ggsave("figures/figure_S2.png")
ggsave("figures/figure_S2.pdf", dpi = 300)
