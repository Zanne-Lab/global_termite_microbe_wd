library(tidyverse)
library(maps)
library(ggforce)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgdal)
library(plotbiomes)
library(patchwork)
library(ggeffects)
#read data
a <- read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T", "C") &
           !is.na(termite_exposure) &
           species != "Simarouba_amara")


#discovery model
mm<-a %>%
  filter(trt %in% c("T")) %>%
  filter(!is.na(termite_exposure))
names(mm)
mm$date_diff_years<-mm$date_diff/365
mod1<-glm(termite_exposure~temp*prec,
          data=mm,offset=date_diff_years,family = binomial(link = "logit"),
          control = list(epsilon = 1e-10,maxit = 100))
ggpredict(mod1, c("temp[all]","prec[all]")) ->l

l <- as.data.frame(l)
l$group <- as.numeric(as.character(l$group))

#calculate site level values

mm %>%
  filter(!is.na(lat) & !is.na(realms)) %>%
  group_by(country, site,realms) %>%
  summarise(lat=unique(lat),
            long=unique(long),
            k_mean =mean(k_value, na.rm=TRUE),
            prec = mean(prec, na.rm = TRUE),
            temp=mean(temp, na.rm = TRUE)) %>%
  ungroup()%>%
  group_by(realms) %>%
  add_count(realms) %>%
  mutate(realms_n = paste0(realms, " (n = ",n,")"))->site_mean
#join discovery model prediction to site level variables
site_mean <- left_join(site_mean,l, by =c("prec"="group", "temp"="x"))

#map locations Figure 1 A
world <- ne_countries(scale = "medium", returnclass = "sf") 

#site locations
sites <- ggplot(world)+
  geom_sf(fill= "white")+
  geom_point(data=site_mean,
             aes(x=long, y=lat,colour=predicted*100), size = 3.5)+
  scale_colour_viridis_c(option = "B", name = "% discovered", 
                         breaks = c(0,10,30,50,70,80,90,100))+
  xlab(expression(paste("Longitude (",degree,")"))) +
  ylab(expression(paste("Latitude (",degree,")")))+
  labs(title = "A")+
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE)+
  theme_classic()+
  theme(text = element_text(size = 9, family = "Helvetica"), title = element_text(size=10))


#Map discovery across biomes Figure 1 B
Whittaker_biomes$precp_mm <- Whittaker_biomes$precp_cm*10
biome <-ggplot()+geom_polygon(data = Whittaker_biomes,aes(x = temp_c, y = precp_mm, fill = biome), 
                              colour = "gray98", size = 1) +
  scale_fill_manual(name = "Whittaker biomes", 
                    breaks = names(Ricklefs_colors), labels = names(Ricklefs_colors),
                    values = Ricklefs_colors) +
  xlab(expression(paste("MAT (",degree," C)"))) +
  scale_y_continuous("MAP (mm)")+
  geom_point(data = site_mean,
             aes(x = temp,
                 y = prec,
                 colour = predicted*100), # map `status` variable to colour
             shape = 16,
             size  = 2.5)+
  scale_color_viridis_c(option = "B", name = "% discovered", 
                        breaks = c(0,10,30,50,70,80,90,100), guide = "none")+
  labs(title = "B")+
  theme_classic()+
  theme(text = element_text(size = 9, family = "Helvetica"), title = element_text(size=10),
        legend.position = "none")

#biomes boxplots by termite discovery figure 1 C
a_dat <- read_csv("data/decomp_with_covar.csv") %>%
  filter(trt %in% c("T", "C") &
           species != "Simarouba_amara")


a %>%
  group_by(site, termite_exposure, biome) %>%
  summarise(across(c(k_value,temp, prec,alt, lat,forestcover100, aridity),mean),sample_size=n())->mm

mm %>%
  group_by(termite_exposure, biome) %>%
  summarise(sample_size=n(), k_value= max(k_value))->ss

mm %>%
  group_by( biome) %>%
  summarise(sample_size=n(), k_value= min(k_value), termite_exposure= NA)->sss

mm$biome <-factor(mm$biome,levels = c("Boreal forest",
                                      "Temperate rain forest",
                                      "Temperate seasonal forest",
                                      "Woodland/shrubland",
                                      "Temperate grassland/desert",
                                      "Tropical rain forest",
                                      "Tropical seasonal forest/savanna",
                                      "Subtropical desert"))

boxplot_biomes <- ggplot(mm, aes(x = biome, y =k_value, fill = biome, linetype =as.factor(termite_exposure) ))+
  geom_boxplot()+
  scale_linetype_discrete(labels = c("Undiscovered","Discovered"),name   = "Termite discovery")+
  scale_y_log10()+
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors)+
  geom_text(data = filter(ss,termite_exposure==0 &
                            !(biome %in% c("Boreal forest", "Temperate grassland/desert"))),
            aes(y = k_value,x = biome,label = sample_size),
            vjust= -1, hjust=2,size=5*.36)+
  geom_text(data = filter(ss,biome %in% c("Boreal forest", "Temperate grassland/desert")),
            aes(y = k_value,x = biome,label = sample_size),
            vjust= -1,size=5*.36)+
  geom_text(data = filter(ss,termite_exposure==1),
            aes(y = k_value-0.23,x = biome,label = sample_size),
            vjust= -.25,hjust=-2,size=5*.36)+
  geom_text_repel(data = sss,
                  aes(y = k_value,x = biome,label = biome),direction = "y",
                  nudge_y = -.5,
                  max.overlaps = 50,
                  #vjust = .5,
                  arrow = arrow(angle = 10,length = unit(5, "points")),size=5*.36)+
  labs(y='k', x = "Biome",title = "C")+
  theme_classic()+
  theme(axis.text.x = element_blank(),axis.text = element_text(size = 9),
        title = element_text(size=10),
        legend.key.size = unit(8,units = "points"))
#join panels
#discovery by biome and decay by biomes-termite discovery
p<-biome+coord_fixed(1/100)+boxplot_biomes+  plot_layout(guides = "collect",nrow = 1)
#join all panels
l <- sites/p+plot_layout(heights = c(2,1))
ggsave("figures/fig_1.pdf",plot = l,width = 10,height = 7,dpi = 300)
ggsave("figures/fig_1.png",plot = l,width = 10,height = 7,dpi = 300)
