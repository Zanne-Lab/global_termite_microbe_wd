library(tidyverse)
library(ggeffects)
library(patchwork)
library(easystats)
#read in data
a <- read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T", "C") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara")

#filter only termite records
mm<-a %>%
  mutate(abs_lat = abs(lat)) %>%
  filter(trt %in% c("T")) %>%
  filter(!is.na(termite_exposure))
names(mm)
mm$date_diff_years<-mm$date_diff/365

#generate site mean DB
mm %>% 
  filter(trt == 'T') %>% 
  group_by(site) %>% 
  summarise(temp = temp[1], 
            prec = prec[1], 
            forestcover100 = forestcover100[1], 
            alt = alt[1],
            date_diff_years=date_diff_years[1],
            lat = lat[1], 
            n = n(), 
            discovered = mean(termite_exposure),
            abs_lat = abs(lat))->mmm
#Generate termite discovery across precpitation model
modpre<-glm(termite_exposure~prec,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
#generate model predictions
modpre_pred <- ggpredict(modpre, c("prec"))


#Figure S1 D discovery across precipitation
b<-ggplot(mmm,aes(x=prec,y=discovered*100))+
  geom_point()+
  geom_line(aes(x=x, y=predicted*100), data=modpre_pred, 
            inherit.aes=FALSE, colour='#000000') + 
  geom_ribbon(aes(x=x, ymin=conf.low*100, ymax=conf.high*100), data=modpre_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#000000')+
  labs(x = "MAP (mm)", y = "% discovered")+
  theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))
#Generate termite discovery across temperature model
modtemp<-glm(termite_exposure~temp,
             data=mm,offset=date_diff_years,family = binomial(link = "logit"),
             control = list(epsilon = 1e-10,maxit = 100))
#generate model predictions
modtemp_pred <- ggpredict(modtemp, c("temp [all]")) %>% filter(x > 1.5)


#Figure S1 C discovery model across temperature
a<-ggplot(mmm,aes(x=temp,y=discovered*100)) +
  geom_point()+
  geom_line(aes(x=x, y=predicted*100), data=modtemp_pred, 
            inherit.aes=FALSE, colour='#000000') + 
  geom_ribbon(aes(x=x, ymin=conf.low*100, ymax=conf.high*100), data=modtemp_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#000000')+
  labs(x = expression(`MAT `(degree * C)), y = "% discovered")+
  theme_classic()+ theme(legend.position = "none",
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))

#Generate termite discovery across elevation model
modalt<-glm(termite_exposure~alt,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
#generate model predictions
modalt_pred <- ggpredict(modalt, c("alt"))

#Figure S1 A discovery across elevation
e<-ggplot(mmm,aes(x=alt,y=discovered*100))+ 
  geom_point()+
  geom_line(aes(x=x, y=predicted*100), data=modalt_pred, 
            inherit.aes=FALSE, colour='#000000') + 
  geom_ribbon(aes(x=x, ymin=conf.low*100, ymax=conf.high*100), data=modalt_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#000000')+
  labs(x = "Elevation (m)", y = "% discovered")+
  theme_classic()+ theme(legend.position = "none",
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))

#Generate termite discovery across absolute latitude model
modlat<-glm(termite_exposure~abs_lat,
            data=mm,offset=date_diff_years,family = binomial(link = "logit"),
            control = list(epsilon = 1e-10,maxit = 100))
#generate model predition
modlat_pred <- ggpredict(modlat, c("abs_lat"))

#Figure S1 B discovery across absolute latitude
f<-ggplot(mmm,aes(x=abs_lat,y=discovered*100))+
  geom_point()+
  geom_line(aes(x=x, y=predicted*100), data=modlat_pred, 
            inherit.aes=FALSE, colour='#000000') + 
  geom_ribbon(aes(x=x, ymin=conf.low*100, ymax=conf.high*100), data=modlat_pred, 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#000000')+
  labs(x =expression(`Absolute Latitude `(degree)), y = "% discovered")+
  theme_classic()+ theme(legend.position = "none",axis.title.y=element_blank(),
                         title = element_text(size=10), text = element_text(size = 9, family = "Helvetica"),
                         legend.key.size = unit(8,units = "points"))


(e+f) / (a + b ) + plot_annotation(tag_levels = "A")
ggsave("figures/figure_S1.png",width=9.2,height=4)
ggsave("figures/figure_S1.pdf",dpi = 300)


#Figure S1 caption, termite discovery quantiles
mmm %>% summarise(quantile(discovered, c(0.05, 0.5, 0.95)))
