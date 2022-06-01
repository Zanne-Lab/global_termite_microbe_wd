#download the data from repo place it in data folder

if(!file.exists("data/global_termite_microbe_wd.xlsx")) {
  download.file("https://figshare.com/ndownloader/files/35496569",
                "data/global_termite_microbe_wd.xlsx")
}

library(readxl)
library(readr)
library(dplyr)
# library(tidyr)
library(magrittr)
read_excel("data/global_termite_microbe_wd.xlsx", sheet=1) %>% 
  left_join(read_csv("data/covars_by_site.csv")) %>% 
  rename(lat = latitude, 
         long = longitude, 
         trt = treatment, 
         species = wood_used, 
         termite_exposure = termite_discovery) %>% 
  mutate(prec_m = prec/1000) %>% 
  write_csv("data/decomp_with_covar.csv")
  
  
