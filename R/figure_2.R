library(tidyverse)
library(ggeffects)
library(patchwork)
library(jpeg)
#load data
mm <- read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T", "C") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara")
#calculate date since deployment using years as units
mm$date_diff_years<-mm$date_diff/365

#termite discovery with temp and precip model
mod2<-glm(termite_exposure ~ temp * prec,
          data=filter(mm, trt == 'T'), offset=date_diff_years, family=binomial(link = "logit"),
          control=list(epsilon = 1e-10, maxit = 100))
#calculate model prediction over temperature holding precipitation constant on three values

#precipitation values
x1 <- 250
y1 <- 2000
z1 <- 2700
v <- c(x1, y1, z1)
#temperature values
v2 <- seq(1.8, 26.5, by = 0.5)
#predict model over climate space
mod2_pred <- ggpredict(mod2, c("temp [v2]","prec [v]"))

#add images for termite and fungi
macrotermite_d <- readJPEG("images/soldier.jpg",native = TRUE)
fungi <- readJPEG("images/fungi.jpg",native = TRUE)
termite_fungi <- readJPEG("images/termiteFungi.jpg",native = TRUE)

#figure 2 A termite discovery over climate (temperature and rainfall)
mm %>% 
  filter(trt == 'T') %>% 
  group_by(site) %>% 
  summarise(temp = temp[1], 
            prec = prec[1], 
            alt = alt[1], 
            lat = lat[1], 
            n = n(), 
            discovered = mean(termite_exposure)) %>% 
  ggplot(aes(temp, discovered*100, color=prec, size = n)) + 
  geom_point() + 
  scale_x_log10()+
  #scale_colour_viridis_c(option = "E")+ #317A22
  scale_colour_gradient2(low='#DCBB50',midpoint = 1700,mid = '#727272', high='#317A22', name='MAP (mm)') +
  geom_line(aes(x=x, y=predicted*100), data=filter(mod2_pred, group==x1), 
            inherit.aes=FALSE, colour='#DCBB50') + 
  geom_ribbon(aes(x=x, ymin=conf.low*100, ymax=conf.high*100), data=filter(mod2_pred, group==x1), 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#DCBB50') + 
  geom_line(aes(x=x, y=predicted*100), data=filter(mod2_pred, group==y1), 
            inherit.aes=FALSE, colour='#727272') + 
  geom_ribbon(aes(x=x, ymin=conf.low*100, ymax=conf.high*100), data=filter(mod2_pred, group==y1), 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#727272') + 
  geom_line(aes(x=x, y=predicted*100), data=filter(mod2_pred, group==z1), 
            inherit.aes=FALSE, colour='#317A22') + 
  geom_ribbon(aes(x=x, ymin=conf.low*100, ymax=conf.high*100), data=filter(mod2_pred, group==z1), 
              inherit.aes=FALSE, alpha=0, linetype='dashed', colour='#317A22') + 
  xlab(expression(`MAT `(degree * C))) + ylab('% discovered') + 
  labs(title = "A")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 9, family = "Helvetica"), 
        title = element_text(size=10),
        axis.title.x = element_blank(),axis.text.x=element_blank())->a_plot







#Model for decay value at site level across MAT and MAP for termite discovered and undiscovered wood

#mean site values 
site_means <- mm %>% 
  group_by(site, termite_exposure) %>% 
  summarise(k_value = mean(k_value), 
            temp = temp[1], 
            prec = prec[1], 
            alt = alt[1], 
            lat = lat[1], 
            n = n()) %>% 
  ungroup() %>% 
  mutate(termite_exposure = factor(termite_exposure), 
         abs_lat = abs(lat), 
         prec_m = prec/1000, 
         k_value_tr = log(k_value))
#decay model by termite exposure and across precip and temperature
m2w <- lm(k_value_tr ~ termite_exposure*temp, data=site_means, weights=n)

x1 <- 250
y1 <- 2000
z1 <- 2700
v <- c(x1, y1, z1)
#Predict model over temperature and precipitation
mod2_pred <- ggpredict(m2w, c("temp"))

#generate DB with prediction for plotting
pred_m <- bind_rows(data.frame(termite_exposure=factor(0, levels=c('0', '1')), 
                               alt=mean(site_means$alt), 
                               abs_lat=mean(site_means$abs_lat), 
                               temp=seq(min(site_means$temp), max(site_means$temp), 1)))
pred_m$k_value_tr <- predict(m2w, newdata=pred_m, se.fit=FALSE)
pred_m$k_value_tr_se <- predict(m2w, newdata=pred_m, se.fit=TRUE)$se.fit
pred_m <- pred_m %>% 
  mutate(k_value = exp(k_value_tr), 
         k_value_upper = exp(k_value_tr + k_value_tr_se), 
         k_value_lower = exp(k_value_tr - k_value_tr_se))

#figure 2b microbial decay across temperature and precipitation using prediction from model
b_plot <- ggplot(filter(site_means, termite_exposure=='0'), aes(x=temp, y=k_value, size=n)) + 
  geom_point(alpha=0.75) + 
  labs(x=expression(`MAT `(degree * C)), y='k (per year)') + 
  scale_size(name='# wood blocks') + 
  scale_x_log10()+
  stat_smooth(aes(x=temp, y=k_value), data=filter(pred_m, termite_exposure=='0'), 
              method='glm', formula = y~x,method.args = list(family = gaussian(link = 'log')), 
              se=FALSE, inherit.aes=FALSE, colour='black') + 
  stat_smooth(aes(x=temp, y=k_value_upper), data=filter(pred_m, termite_exposure=='0'), 
              method='glm', formula = y~x,method.args = list(family = gaussian(link = 'log')), 
              se=FALSE, inherit.aes=FALSE, colour='black', linetype='dashed') + 
  stat_smooth(aes(x=temp, y=k_value_lower), data=filter(pred_m, termite_exposure=='0'), 
              method='glm', formula = y~x,method.args = list(family = gaussian(link = 'log')), 
              se=FALSE, inherit.aes=FALSE, colour='black', linetype='dashed') + 
  guides(size = "none")+
  labs(title = "B")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        text = element_text(size = 9, family = "Helvetica"), 
        title = element_text(size=10),
        axis.text.x=element_blank())



### plot k predictions for termite discoverd wood blocks
pred_t <- data.frame(termite_exposure=factor(1, levels=c('0', '1')), 
                     alt=mean(site_means$alt), 
                     abs_lat=mean(site_means$abs_lat), 
                     prec=mean(site_means$prec), 
                     temp=seq(min(site_means$temp), max(site_means$temp), 1))
pred_t$k_value_tr <- predict(m2w, newdata=pred_t, se.fit=FALSE)
pred_t$k_value_tr_se <- predict(m2w, newdata=pred_t, se.fit=TRUE)$se.fit
pred_t <- pred_t %>% 
  mutate(k_value = exp(k_value_tr), 
         k_value_upper = exp(k_value_tr + k_value_tr_se), 
         k_value_lower = exp(k_value_tr - k_value_tr_se))
#Figure 2c decay of termite discovered blocks
c_plot <- ggplot(filter(site_means, termite_exposure=='1'), aes(x=temp, y=k_value, size=n)) + 
  geom_point(alpha=0.75) + 
  labs(x=expression(`MAT `(degree * C)), y='k (per year)') + 
  scale_x_log10()+
  scale_size(name='# wood blocks') + 
  stat_smooth(aes(x=temp, y=k_value), data=filter(pred_t, termite_exposure=='1'), 
              method='glm', formula = y~x,method.args = list(family = gaussian(link = 'log')), 
              se=FALSE, inherit.aes=FALSE, colour='black') + 
  stat_smooth(aes(x=temp, y=k_value_upper), data=filter(pred_t, termite_exposure=='1'), 
              method='glm', formula = y~x,method.args = list(family = gaussian(link = 'log')), 
              se=FALSE, inherit.aes=FALSE, colour='black', linetype='dashed') + 
  stat_smooth(aes(x=temp, y=k_value_lower), data=filter(pred_t, termite_exposure=='1'), 
              method='glm', formula = y~x,method.args = list(family = gaussian(link = 'log')), 
              se=FALSE, inherit.aes=FALSE, colour='black', linetype='dashed') + 
  labs(title = "C")+
  theme_classic()+
  theme(text = element_text(size = 9, family = "Helvetica"), title = element_text(size=10))

# Put panel a, b and c together
fig_2 <-a_plot+inset_element(macrotermite_d,left = 0.0,
                             bottom = 0.65,
                             right = 0.5,
                             top = 0.95)+
  b_plot+inset_element(fungi,left = 0.0,
                       bottom = 0.65,
                       right = 0.5,
                       top = 0.95)+plot_layout(ncol = 1)+
  c_plot+inset_element(termite_fungi,left = 0.0,
                       bottom = 0.65,
                       right = 0.5,
                       top = 0.95)+plot_layout(ncol = 1)

fig_2 &  xlim(1, 30)
ggsave(fig_2,filename = "figures/fig_2.pdf",height = 13,width = 8, dpi = 300)
ggsave(fig_2,filename = "figures/figure_2.png",height = 13,width = 8)
