library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(stringr)

genera<-c("Azemiops","Agkistrodon","Atropoides","Bothriechis","Bothrocophias","Bothrops","Cerrophidion","Crotalus","Lachesis","Metlapilcoatlus","Mixcoatlus","Ophryacus","Porthidium","Sistrurus","Calloselasma","Deinagkistrodon","Garthius","Gloydius","Hypnale","Ovophis","Protobothrops","Trimeresurus","Tropidolaemus","Atheris","Bitis","Causus","Cerastes","Daboia","Echis","Eristicophis","Macrovipera","Montatheris","Montivipera","Proatheris","Pseudocerastes","Vipera")
nwcrotes<-c("Agkistrodon","Atropoides","Bothriechis","Bothrocophias","Bothrops","Cerrophidion","Crotalus","Lachesis","Metlapilcoatlus","Mixcoatlus","Ophryacus","Porthidium","Sistrurus")

`%notin%` <- Negate(`%in%`)

roll<-st_read("~/Dropbox/Projects/2020_MacroCharDisp/Occurence_RangeMaps/shp_databases/Roll2017_distributions/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
roll <- roll %>% filter(Group=="snake") %>% 
  separate(Binomial, c("genus","species")) %>% 
  mutate(Species=paste(genus,species,sep="_"), Subspecies="") %>% 
  filter(genus %in% genera) %>% 
  filter(genus %notin% nwcrotes) %>%
  select(Species,Subspecies)

load("~/Dropbox/GitHub/VenomMaps/data/shiny-data3.RData")

combined_distribution<-combined_distribution %>% select(!admin2)
combined_distribution<-rbind(combined_distribution, roll)

rm(genera, nwcrotes, `%notin%`, roll)

save.image("~/Dropbox/GitHub/VenomMaps/data/shiny-data.RData")

countries<-ne_countries(returnclass = "sf")

tmp<-st_intersection(st_buffer(combined_distribution,0), countries)

tmp2<-st_drop_geometry(tmp) %>% select(Species,admin) %>% group_by(Species) %>% distinct() %>% summarise(admin=paste(admin,collapse='; '))

write.csv(tmp2,"tmp_country_overlap.csv", quote=F, row.names = F)
