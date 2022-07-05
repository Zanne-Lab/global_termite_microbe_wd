library(tidyverse)
library(ggeffects)
library(khroma)

# discovery effects on decay within biomes (supplemental analysis)
read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara") %>% 
  group_by(site, biome) %>%
  summarise(across(c(k_value,temp, prec,alt, lat,forestcover100, aridity),mean), 
            sample_size=n(), 
            discovery=mean(termite_exposure)*100) %>% 
  mutate(k_value_tr = log(k_value), 
         biome = factor(biome,levels = c("Boreal forest",
                                         "Temperate rain forest",
                                         "Temperate seasonal forest",
                                         "Woodland/shrubland",
                                         "Temperate grassland/desert",
                                         "Tropical rain forest",
                                         "Tropical seasonal forest/savanna",
                                         "Subtropical desert"))) %>% 
  ungroup() %>% 
  group_by(biome) %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  ungroup() %>% 
  mutate(biome = fct_drop(biome)) -> mm
table(mm$biome)
mm %>% group_by(biome) %>% summarise(n = n(), discovery = max(discovery))

sunset <- colour("sunset")
Biomes <- sunset(9)
names(Biomes) <- c("Tundra" ,"Boreal forest" ,"Temperate seasonal forest",
                   "Temperate rain forest" ,"Tropical rain forest", 
                   "Tropical seasonal forest/savanna","Subtropical desert" ,
                   "Temperate grassland/desert" ,"Woodland/shrubland")
Biomes <- Biomes[names(Biomes) %in% unique(mm$biome)]
Biomes <- c(Biomes, All='black')

m1 <- lm(k_value_tr ~ discovery * biome, data=mm, weights=sample_size)
car::Anova(m1)
summary(m1)
pred_d <- ggeffect(m1, terms='discovery') %>% 
  mutate(group = 'All')
pred_db <- ggpredict(m1, terms=c('discovery [0:100 by=1]', 'biome')) %>% 
  left_join(mm %>% 
              rename(group = biome) %>% 
              group_by(group) %>% 
              summarise(min_k_tr = min(k_value_tr), 
                        max_k_tr = max(k_value_tr), 
                        min_x = floor(min(discovery)), 
                        max_x = ceiling(max(discovery)))) %>% 
  group_by(group) %>% 
  filter(x < max_x)

k_tr <- c(`0.01`=log(0.01), 
          `0.1`=log(0.1), 
          `1`=log(1), 
          `5`=log(5))

ggplot(pred_db, aes(x=x, y=predicted, colour=group, fill=group)) + 
  geom_line(data=pred_d, mapping=aes(x=x, y=predicted), 
            size=1.5) +
  geom_ribbon(data=pred_d, mapping=aes(ymin=conf.low, ymax=conf.high, colour=NULL), alpha=0.5) +
  geom_line(size=0.75, linetype='longdash') + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,
                  colour=NULL), alpha=0.5) +
  scale_colour_manual(breaks = names(Biomes),
                      labels = names(Biomes),
                      values = Biomes) +
  scale_fill_manual(breaks = names(Biomes),
                      labels = names(Biomes),
                      values = Biomes) +
  geom_point(data=mm, mapping=aes(x=discovery, y=k_value_tr, colour=biome), 
           size=0.75, inherit.aes=FALSE) +
  scale_y_continuous(breaks = k_tr, labels=names(k_tr)) + 
  labs(x='% discovered', y='k (per year)', colour='Whittaker biomes', fill='Whittaker biomes') +
  theme_classic()+
  theme(axis.text = element_text(size = 9),
        title = element_text(size=10),
        legend.key.size = unit(8,units = "points")) -> fig_s4
fig_s4
ggsave(fig_s4,filename = "figures/figure_s4.png",height=4,width=8)

