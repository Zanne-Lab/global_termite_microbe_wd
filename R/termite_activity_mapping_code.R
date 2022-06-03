library(raster)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)

## Mapping where shifts in the climate make termite activity likely to increase

#download 1970-2000 "climate norms" for tmean
tmean <- mean(getData('worldclim', var = "tmean", res = 2.5)) / 10

#manual download of CMIP6 layers since automatic seems to be bugged in raster
rast_name_list <-
  as.list(list.files("future_clim/", full.names = TRUE))
rast_list <- stack(lapply(rast_name_list, raster))

#get tmeans
rast_list <- stack(lapply(rast_name_list, raster, band = 1))
#get precips
precip_rast_list <- stack(lapply(rast_name_list, raster, band = 12))

#record keeping
term_increase_area_vec <- NA
mean_warming_vec <- NA
#str(fut)

#find experimental range of precips
biomes <-
  readr::read_csv("processed_data/forestcover_atsitelevel.csv")
biomes$tmean <- raster::extract(tmean, cbind(biomes$long, biomes$lat))
MAP <- sum(getData('worldclim', var = "prec", res = 2.5))
#plot(MAP)
biomes$MAP <- raster::extract(MAP, cbind(biomes$long, biomes$lat))
min_precip <- min(biomes$MAP)
max_precip <- max(biomes$MAP)


#for graphical purposes
EqualEarthProj <-
  "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#reprojecting climate norms
MAP_ee <- raster::projectRaster(MAP, crs = EqualEarthProj)
tmean_ee <- raster::projectRaster(tmean, crs = EqualEarthProj)
#read_in_model_from_table_s3
dm <- readRDS("discovery.rds")

#
#looping through alternate futures
#

for (i in 1:length(rast_name_list)) {
  print(i)
  # this is a check to exclude misnamed files coming from worldclim
  if (max(getValues(rast_list[[i]]), na.rm = T) > 50) {
    term_increase_area_vec[i] <- NA
    mean_warming_vec[i] <- NA
    next
  }
  #data frame for calculating areas
  df_current <-
    data.frame(
      prec = getValues(MAP),
      temp = getValues(tmean),
      area = getValues(terra::area(tmean)),
      date_diff_years = 1
    )
  # excluding Antarctica from the futures to match the climate norms
  df_future <-
    data.frame(
      prec = getValues(crop(precip_rast_list[[i]], tmean)),
      temp = getValues(crop(rast_list[[i]], tmean)),
      area = getValues(terra::area(tmean)),
      date_diff_years = 1
    )
  
  #estimate from model
  df_current$current <- round(predict(dm, newdata = df_current, type = "response"), 2)
  df_future$future <- round(predict(dm, newdata = df_future, type = "response"), 2)
  
  #calculate term activity areas for clim norms
  df_current %>% mutate(high_term=current>0.5) %>%
    group_by(high_term) %>%
    summarize(sum(area))
  
  #average change weighted by terra areas
  mean_temp_current <- weighted.mean(x=df_current$temp,df_current$area,na.rm=TRUE)
  mean_temp_future <- weighted.mean(x=df_future$temp,df_future$area,na.rm=TRUE)
  
  mean_warming_vec[i] <-
    mean_temp_future - mean_temp_current
  
  current_discovery <- setValues(tmean, values = df_current$current)
  future_discovery <- setValues(tmean, values = df_future$future)

  
 # get into the right projection for plotting
  cd_ee <-
   raster::projectRaster(current_discovery, crs = EqualEarthProj)
  fd_ee <-
   raster::projectRaster(future_discovery, crs = EqualEarthProj)
  map_future_ee <-
   raster::projectRaster(crop(precip_rast_list[[i]], MAP), crs = EqualEarthProj)
  mat_future_ee <-
   raster::projectRaster(crop(rast_list[[i]], MAP), crs = EqualEarthProj)
  
  
#setting up for geom_tile  
  y <- data.frame(
    long = xFromCell(MAP_ee, 1:length(MAP_ee)),
    MAP = getValues(crop(map_future_ee, MAP_ee)),
    MAT = getValues(crop(mat_future_ee, MAP_ee))
  )
  y$lat <- yFromCell(MAP_ee, 1:length(MAP_ee))
  y$fd <- getValues(fd_ee)
  y$cd <- getValues(cd_ee)

  yy <- filter(y, !is.na(MAP)) #remove water

  #note that there are other edge cases of transitions, but this only focuses on <50% to >50%
  yy$Discovery <- dplyr::case_when(
    yy$MAP < min_precip * 0.9 |
      yy$MAP > max_precip * 1.1 ~ "Unable to predict",
    yy$cd > 0.5 ~ "Current >50%",
    yy$fd > 0.5 ~ "Mid century expansion >50%",
    yy$cd > 0.05 ~ "Continuing >5 & <50%",
    yy$cd <= 0.05 ~ "Continuing <5%",
    is.na(yy$MAP) ~ "water"
  )
  # 
  # table(yy$Discovery)
  # 
   z2 <- dplyr::filter(yy, !is.na(Discovery) & Discovery != "water")
  
  #term increase
  term_incr <-
    sum(getValues(
      terra::area(future_discovery > 0.5 &
                    current_discovery < 0.5) * future_discovery > 0.5 &
        current_discovery < 0.5
    ),
    na.rm = T)
  term_increase_area_vec[i] <- term_incr

  #let geom_tile do the work
   p1 <- ggplot() +
     geom_tile(data = z2 , aes(x = long, y = lat,
                               fill = Discovery)) +
     scale_fill_manual(values = c("grey", "#97B669", "#A09700", "#DCBB50", "#fffdd0")) +
     theme_bw() + theme_void() +
     theme(legend.position = "bottom")
  
   ggsave(paste0("another_view_v3", i, ".pdf"),
          p1,
          width = 11,
          height = 4)
   ggsave(paste0("another_view_v3.", i, ".jpg"),
          p1,
          width = 11,
          height = 4)
}


zz <-
  data.frame(
    term_inc = term_increase_area_vec,
    mean_warming = mean_warming_vec,
    scenar = unlist(rast_name_list)
  )
write_csv(zz, "summary_output.csv")

library(tidyverse)
read_csv("summary_output.csv") %>%
  mutate(`Shared \nsocioeconomic \npathway` = case_when(
    grepl("ssp126", scenar) ~ "ssp126",
    grepl("ssp585", scenar) ~ "ssp585"
  )) -> zzz

zzz %>%
  ggplot(aes(x = mean_warming, y = term_inc / 1000000, col = `Shared \nsocioeconomic \npathway`)) +
  geom_point(size = 3) + theme_classic() +
  xlab("Global mean terrestrial warming in 2041-2060 relative to 1970-2000 baseline") +
  ylab("Increase in terrestrial area with termite discovery >50% (km2*106)") +
  scale_color_manual(values = c("#97B669", "#DCBB50"))
ggsave("all_scenarios.pdf", width = 7, height = 5)


# get only legend for figure assembly
legend <- cowplot::get_legend(p1)
grid.newpage()
grid.draw(legend)
