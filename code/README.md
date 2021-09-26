# Venom Maps

### Rhett M. Rautsaw, Gustavo Jiménez, Christoph Grünwald, Erich P. Hofmann, Laura Alencar, Marcio Martins, Paola Carrasco, Tiffany Doan, & Christopher L. Parkinson

2021 September 26

***

This repository contains the code used to

1.  Update occurrence record taxonomy for New World Crotalinae from several sources.
2.  Construct updated digital distributions all 159 species of New World Crotalinae.
3.  Use distributions to further update occurrence record taxonomy.

# Table of Contents

  - [Prep](#prep)
      - [Synonyms List](#synonyms-list)
  - [Collecting Occurrence Records](#collecting-occurrence-records)
      - [Filtering Databases](#filtering-databases)
      - [Combining Databases](#combining-databases)
      - [Cleaning `subspecies` Column](#cleaning-subspecies-column)
      - [Updating Taxonomy](#updating-taxonomy)
  - [Constructing Distribution Maps](#constructing-distribution-maps)
      - [Combining Databases](#combining-databases-1)
      - [qGIS](#qgis)
      - [Smooth & Clip Distributions](#smooth-clip-distributions)
      - [Extract Country Information for Each
        Species](#extract-country-information-for-each-species)
  - [Spatial Join Correction](#spatial-join-correction)
      - [Reference Phylogeny
        (bind.tip2.R)](#reference-phylogeny-bind.tip2.r)
      - [Occurence Point Cleaner
        (occ\_cleaner.R)](#occurence-point-cleaner-occ_cleaner.r)
      - [Manual Check](#manual-check)
  - [Final Occurrence Counts](#final-occurrence-counts)
  - [Summarizing ENM Results](#summarizing-enm-results)
      - [Averaging Selected Final
        Models](#averaging-selected-final-models)
      - [Gathering Final Model
        Statistics](#gathering-final-model-statistics)
      - [Summary Statistics](#summary-statistics)
      - [Plotting AUC and Variable
        Contributions](#plotting-auc-and-variable-contributions)
      - [Plotting Example Niche Models](#plotting-example-niche-models)
      - [Plotting Example Bioclim vs. Combo
        Model](#plotting-example-bioclim-vs.-combo-model)
      - [Summarizing Mean Combo Models (Supplemental Table
        1)](#summarizing-mean-combo-models-supplemental-table-1)
  - [References](#references)


# Prep

``` r
library(readr)
library(readxl)
library(tidyverse)
library(data.table)
library(rgdal)
library(raster)
library(sp)
library(sf)
library(smoothr)
library(ape)
library(phytools)
library(kableExtra)
library(ggspatial)
library(ggnewscale)
library(ggpubr)
library(patchwork)
```

## Synonyms List

I obtained a list of all synonyms for species in “Crotalinae” from the
[Reptile Database (May 2020)](http://www.reptile-database.org/) (Thanks
Peter Uetz\!) and combined this with a full list of viper species to
create a long list of synonyms.

``` r
synonyms<-read.csv("../OtherData/reptile_database/Crotalinae_synonyms_sorted2.csv",colClasses = rep("character",4))
ReptileDB <- read_xlsx("../OtherData/reptile_database/reptile_checklist_2020_04.xlsx")
Viperidae <- ReptileDB[ReptileDB$Familyetc %like% "Viperidae", ]$Species
synonyms2<-c(synonyms$Species,synonyms$Synonym, Viperidae)
synonyms2<-unique(synonyms2)

rm(ReptileDB, Viperidae)
```

# Collecting Occurrence Records

Data for Viperidae was downloaded from [GBIF](https://www.gbif.org/),
[HerpMapper](https://www.herpmapper.org/),
[iDigBio](https://www.idigbio.org/), [Bison](https://bison.usgs.gov/),
[Brazil Snake
Atlas](https://bioone.org/journals/South-American-Journal-of-Herpetology/volume-14/issue-sp1/SAJH-D-19-00120.1/Atlas-of-Brazilian-Snakes--Verified-Point-Locality-Maps-to/10.2994/SAJH-D-19-00120.1.short),
[BioWeb Ecuador](https://bioweb.bio/), and custom
databases/georeferencing.

``` r
gbif<-read_tsv("occ_databases/gbif/2021-08-19/0000311-210819072339941.csv")
herpmapper<-read.csv("occ_databases/herpmapper/2021-08-19_herpmapper.csv")
idigbio<-read.csv("occ_databases/idigbio/2021-08-19_d8076d32-8a81-4d08-a458-c243c5be7ab5/occurrence.csv")
bison<-read_csv("occ_databases/bison/2021-08-19/bison-Viperidae.csv")
brazil_atlas<-read_xlsx("occ_databases/Nogueira2019_BrazilAtlas/Table_S2_Revised_FINAL_OCC_RECORDS.xlsx")
bioweb<-read_tsv("occ_databases/bioweb/2021-07-07_BioWeb_Viperidae.txt")
custom<-read_xlsx("occ_databases/custom/combined_db.xlsx")
```

## Filtering Databases

Databases were cleaned to include the following columns

`ID, species, subspecies, locality, latitude, longitude, accuracy,
source, voucher, issues/notes`

``` r
gbif2<- gbif %>% #[gbif$species %in% synonyms2,] %>% 
  mutate(ID=paste0("gbif_",gbifID), locality=paste(countryCode, stateProvince, locality, sep="; "), voucher=paste(institutionCode, collectionCode, catalogNumber, sep=":")) %>% 
  filter(taxonRank!="GENUS" & taxonRank!="FAMILY" & basisOfRecord!="FOSSIL_SPECIMEN") %>%
  dplyr::select(ID,
         species=species,
         subspecies=infraspecificEpithet,
         locality,
         latitude=decimalLatitude,
         longitude=decimalLongitude,
         accuracy_m=coordinateUncertaintyInMeters,
         source = basisOfRecord,
         voucher = voucher,
         issue)

herpmapper2<- herpmapper %>% #[herpmapper$Taxon %in% synonyms2,] %>% 
  mutate(locality=paste(Country, Level.1, Level.2, sep="; "), source="observation") %>% 
  dplyr::select(ID,
         species=Taxon,
         locality,
         latitude=Latitude,
         longitude=Longitude,
         accuracy_m=Accuracy,
         source,
         voucher=Vouchers)

idigbio2<- idigbio %>% #[str_to_sentence(idigbio$gbif.canonicalName) %in% synonyms2,] %>%
  separate(idigbio.geoPoint, into = c("latitude", "longitude"), sep = ", \"lon\": ") %>%
  mutate(latitude=as.numeric(gsub("\\{\"lat\": ","",latitude))) %>% mutate(longitude=as.numeric(gsub("\\}","",longitude))) %>%
  mutate(gbif.canonicalName=str_to_sentence(gbif.canonicalName, locale = "en")) %>%
  mutate(ID=paste0("idigbio_", coreid), locality=paste(idigbio.isoCountryCode,dwc.stateProvince,dwc.county,dwc.locality,sep="; ")) %>%
  mutate(voucher=paste(dwc.institutionCode, idigbio.collectionName, dwc.catalogNumber, sep=":")) %>%
  dplyr::select(ID,
         species=gbif.canonicalName,
         subspecies=dwc.infraspecificEpithet,
         locality,
         latitude,
         longitude,
         accuracy_m=dwc.coordinateUncertaintyInMeters,
         source=dwc.basisOfRecord,
         voucher,
         issue=idigbio.flags)

bison2<- bison %>% #[bison$scientificName %in% synonyms2,] %>% 
  mutate(ID=paste0("bison_",bisonID), locality=paste(countryCode, calculatedState, calculatedCounty, generalComments, sep="; ")) %>% 
  mutate(voucher=paste(ownerInstitutionCollectionCode, catalogNumber, sep=":")) %>%
  dplyr::select(ID,
         species=scientificName,
         locality,
         latitude=decimalLatitude,
         longitude=decimalLongitude,
         accuracy_m=coordinateUncertaintyInMeters,
         source = basisOfRecord,
         voucher)

brazil_atlas2 <- brazil_atlas[brazil_atlas$Species %in% synonyms2,] %>% 
  dplyr::select(ID,
         species=Species,
         latitude=Latitude,
         longitude=Longitude,
         locality="Locality detailed",
         county="Locality General",
         state="state/Ctry",
         voucher = Voucher,
         source = Source) %>%
  mutate(locality=paste(state, county, locality, sep="; ")) %>% dplyr::select(-state, -county)

bioweb2 <- bioweb %>% 
  filter(is.na(estadoDeterminacion)) %>%
  mutate(locality=paste(pais, provincia, canton, parroquia, localidad, localidadVerbatim, sep="; ")) %>%
  dplyr::select(ID,
         species,
         locality,
         latitude,
         longitude,
         source=fuenteCoordenadas,
         voucher)

custom2 <- custom %>% #[custom$species %in% synonyms2,] %>%
  mutate(subspecies=sub("^\\w+ \\w+","",subspecies)) %>%
  mutate(subspecies=sub("^ ","",subspecies))
```

## Combining Databases

The 7 different databases were then combined into one database.

``` r
combined_data<-Reduce(function(x, y) merge(x, y, all=TRUE), list(gbif2, brazil_atlas2, herpmapper2, idigbio2, bison2, bioweb2, custom2))
rm(gbif, brazil_atlas, herpmapper, idigbio, bison, bioweb, custom, gbif2, brazil_atlas2, herpmapper2, idigbio2, bison2, bioweb2, custom2)
```

## Cleaning `subspecies` Column

Given that each database had different standards for recording
subspecies, we first cleaned up this column to turn it into full
trinomial taxonomy.

``` r
combined_data2 <- combined_data %>% 
  filter(!is.na(species)) %>%
  # Remove subspecies from species column by separating, then combine genus & species back together
  separate(species, c("genus","species","subspecies2"), sep="\\s") %>% 
  mutate(species=paste0(genus," ",species)) %>%
  mutate_if(is.factor,as.character) %>%
  mutate_all(na_if,"") %>% 
  # Combine the two subspecies columns together
  unite(subspecies,c(subspecies,subspecies2), na.rm=T, remove=T) %>%
  mutate(subspecies=sapply(strsplit(subspecies, "_"),function(x) new_col = paste(unique(x), collapse = "_"))) %>%
  mutate_all(na_if,"")  %>%
  mutate(subspecies2=subspecies) %>%
  # Change to full trinomial instead of just subspecies name
  mutate(subspecies=ifelse(is.na(subspecies),species,paste0(species," ",subspecies)))
```

## Updating Taxonomy

Next, the `species` and `subspecies` columns were checked against list
of synonyms with current taxonomy. If the `species` and `subspecies`
columns matched the same taxon (or no subspecies was recorded), then the
species were left as recorded. If `species` and `subspecies` columns did
not match the same taxon (due to subspecies being elevated), the valid
species name that required the least amount of change to the name was
recorded.

``` r
least_diff_taxonomy<-function(name){
  if(name %in% synonyms$Synonym){
    as.character(synonyms[synonyms$Synonym==name,] %>% mutate(diff=adist(Species, Synonym)[,1]) %>% slice_min(diff, with_ties = F))[2]
  }
  else{
    NA
  }
}

combined_data3 <- combined_data2 %>% rowwise() %>%
  # Check species/subspecies ID by synonyms list
  mutate(species_match = least_diff_taxonomy(species), subspecies_match=least_diff_taxonomy(subspecies))

# If match, keep name. Else requires manual check. 
combined_data4 <- combined_data3 %>%
  mutate(species2=(ifelse(species_match==subspecies_match | is.na(subspecies_match), species_match, subspecies_match))) %>%
  mutate(match=(species==species2)) %>%
  mutate(subspecies2=(ifelse(match==TRUE,subspecies,paste0(species2," ",subspecies2)))) %>%
  mutate(subspecies2=sub(" NA","",subspecies2))

combined_data5 <- combined_data4 %>% 
  filter(!is.na(latitude)) %>%
  filter(!grepl(" NA", species)) %>%
  filter(!grepl("fossil",source))
  
write_csv(combined_data5,"occ_databases/z_combined_records/combined_records_v0.1.csv")
```

Taxonomy was then manually checked and further cleaned
(combined\_records\_v1).

# Constructing Distribution Maps

Preliminary distribution maps were obtained from
[IUCN](https://www.iucnredlist.org/resources/spatial-data-download),
[Heimes 2016 (Snakes of
Mexico)](https://www.chimaira.de/herpetofauna-mexicana-vol-1-snakes-of-mexico.html?___store=english&___from_store=default)
provided by Christoph Grünwald and Jason Jones of
[HerpMX](http://herp.mx/), [Roll et
al. 2017](https://www.nature.com/articles/s41559-017-0332-2), and
custom distribution maps made previously.

``` r
#IUCN<-readOGR("shp_databases/IUCN/REPTILES.shp")
#IUCN2 <- IUCN[IUCN$binomial %in% synonyms2,] 
#writeOGR(IUCN2, "shp_databases/IUCN/NW_Crotalinae.shp", driver="ESRI Shapefile", layer="NW_Crotalinae")
IUCN<-readOGR("shp_databases/IUCN/NW_Crotalinae.shp")

#Roll<-readOGR("shp_databases/Roll2017/GARD1.1_dissolved_ranges/modeled_reptiles.shp")
#Roll2 <- Roll[Roll$Binomial %in% synonyms2,] 
#writeOGR(Roll2, "shp_databases/Roll2017/GARD1.1_dissolved_ranges/NW_Crotalinae.shp", driver="ESRI Shapefile", layer="NW_Crotalinae")
Roll<-readOGR("shp_databases/Roll2017/GARD1.1_dissolved_ranges/NW_Crotalinae.shp")

Heimes<-readOGR("shp_databases/SnakesOfMexico/MEXICO_VIPERIDAE.shp")
Crotalus<-readOGR("shp_databases/custom/CrotalusAgkistrodon/CrotalusEdited_10kBuffer_Clip.shp")
Bothriechis<-readOGR("shp_databases/custom/Bothriechis/Bothriechis.shp")
```

## Combining Databases

These distributions were combined into one file and used as a reference
for the manual creation of the new distribution maps for each species of
New World Crotalinae.

``` r
IUCN <- IUCN %>% mutate(source_db="iucn") %>% dplyr::select(binomial, subspecies, source, source_db) %>% rename(species=binomial)
Roll <- Roll %>% mutate(source_db="roll") %>% dplyr::select(Binomial, source_db) %>% rename(species=Binomial)
Heimes <- Heimes %>% mutate(source_db="heimes") %>% rename(species=Species)
Crotalus <- Crotalus %>% mutate(source_db="custom") %>% rename(species=Species, subspecies=Subspecies)
Bothriechis <- Bothriechis %>% mutate(source_db="custom")
```

Underscores were added to the species names.

``` r
combined_distributions<-bind(IUCN,Roll,Heimes,Crotalus,Bothriechis)
combined_distributions$Species<-sub(" ","_",combined_distributions$species)
combined_distributions$Subspecies<-sub(" ","_",combined_distributions$subspecies)
combined_distributions$species<-NULL
combined_distributions$subspecies<-NULL
combined_distributions$id<-""
combined_distributions<-combined_distributions[,c(5,3,4,1,2)]
writeOGR(combined_distributions, "shp_databases/z_combined_distributions/combined_distributions2.shp", driver="ESRI Shapefile", layer="combined_distributions")
```

## qGIS

Nothing here but a lot of manual labor to create a distribution map for
each species/subspecies. To do this, we used the previous distribution
maps as reference and then carefully examined recent publications
alongside a detailed [relief map](https://maps-for-free.com/) and [The
Nature Conservancy Terrestrial
Ecoregions](https://geospatial.tnc.org/datasets/b1636d640ede4d6ca8f5e369f2dc368b/about).

## Smooth & Clip Distributions

Next, we took the distribution maps, smoothed them using `smoothr` and
clipped the distributions to the shoreline using a simplified
[GADM](https://gadm.org/) shapefile.

``` r
files<-list.files("working_map/Species_Ranges",pattern=".geojson",full.names = T)
names<-gsub("working_map/Species_Ranges/","",files) %>% gsub(".geojson","",.)

dir.create("Final/Distributions", recursive=T)

# world<-st_read(dsn="~/Dropbox/GIS_Shapefiles/Basemaps/GADM_Administrative/world.shp")
# world2<-st_simplify(world, dTolerance=0.001, preserveTopology = T)
# st_write(world2,"~/Dropbox/GIS_Shapefiles/Basemaps/GADM_Administrative/world_simplified.shp")
world2<-st_read(dsn="~/Dropbox/GIS_Shapefiles/Basemaps/GADM_Administrative/world_simplified.shp")

for(i in 1:length(files)){
  print(paste(i,names[i]))
  distribution<-st_read(dsn=files[i])
  distribution$id<-NULL
  distribution$NAME_0<-NULL

  distribution<-st_transform(st_buffer(st_transform(distribution,crs=3857),0),crs=4326)
  distribution_smooth<-smooth(distribution, method = "chaikin")
  distribution_clip<-try(st_intersection(x=distribution_smooth, y=world2))
  if(nrow(as.data.frame(distribution_clip))==0){
    distribution_clip<-distribution_smooth
  }
  
  st_write(distribution_clip, paste0("Final/Distributions/",names[i],".geojson"), layer=names[i], driver="GeoJSON")
}
```

## Extract Country Information for Each Species

Here I am performing a spatial join of the final distributions with a
layer of countries to obtain a list of countries in which each species
lives. This will be useful for later when trying to filter species in
the Shiny App by country.

``` r
files<-list.files("Final/Distributions",pattern=".geojson",full.names = T)
names<-gsub("Final/Distributions/","",files) %>% gsub(".geojson","",.)

countries<-st_read("~/Dropbox/GIS_Shapefiles/Basemaps/GADM_Administrative/countries.shp")

overlap<-matrix(ncol=3,nrow=1)
for(i in 1:length(files)){
  distribution<-st_read(files[i])
  tmp<-as.matrix(as.data.frame(st_join(distribution,countries)) %>% dplyr::select(Species, Subspecies,admin))
  overlap<-rbind(overlap,tmp)
  rm(tmp)
}

overlap <- as.data.frame(overlap) %>% dplyr::select(-Subspecies) %>% unique()
overlap <- overlap %>% arrange(Species, admin) %>% group_by(Species) %>% mutate(admin = paste0(admin, collapse="; ")) %>% unique()

write_csv(overlap, "Final/country_occupancy.csv")
```

# Spatial Join Correction

To further clean up the occurrence records, we next performed a spatial
join of occurrence records with the distributions. If records overlapped
with their expected distribution, they were left as recorded. However,
if points did not overlap with their expected distribution it is
possible that these points are erroneous or need to be updated. To check
this, we calculated the distance to all nearby distributions and also
the phylogenetic distance between the recorded species and species with
which the record overlapped.

## Reference Phylogeny (bind.tip2.R)

We used the tree published by [Zaher et
al. 2019](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0216148),
but this tree does not have complete sampling. Luckily, smaller (or more
focused) studies have shown us where in the phylogeny each species
belongs. Although we can’t get the branch lengths just right, we can at
least add them to a general position. First, lets see what species we
are missing in the tree based on our occurrence records.

``` r
tree<-read.nexus("../../2019_AnchorViperPhylo/PreviousPhylos/CollapsedTrees/2019_Zaher_pruned.tre")
data<-read_xlsx("occ_databases/z_combined_records/combined_records_v1.xlsx")

# Fix taxonomy in tree
tree$tip.label<-gsub("Bothrops_lojanus","Bothrocophias_lojanus", tree$tip.label) %>% gsub("Atropoides","Metlapilcoatlus", .) %>% gsub("Metlapilcoatlus_picadoi","Atropoides_picadoi",.) %>% gsub("Bothrops_bilineata","Bothrops_bilineatus",.) %>% gsub("Bothrops_matogrossensis","Bothrops_mattogrossensis",.) %>% gsub("Bothrops_pulchra","Bothrops_pulcher",.) %>% gsub("Bothrops_taeniata","Bothrops_taeniatus",.) %>% gsub("Macrovipera_lebetina","Macrovipera_lebetinus",.)

sort(setdiff(data$species2, tree$tip.label))
```

We are missing a lot of species. I’m going to only add the New World
taxa to the tree since that is the focus of my research. To add new tips
to the phylogeny, I developed `bind.tip2` which will place a new tip on
the phylogeny relative to other tips. Here, I use this function to place
new taxa sister to existing taxa or at the MRCA of two taxa in the tree.

``` r
source("bind.tip2.R")

# Agkistrodon_conanti (Burbrink et al. 2015)
tree<-bind.tip2(tree,"Agkistrodon_conanti", sp1="Agkistrodon_piscivorus", perc=0.5) 
# Agkistrodon_howardgloydi (Porras et al. 2013)
tree<-bind.tip2(tree,"Agkistrodon_howardgloydi", sp1="Agkistrodon_bilineatus", perc=0.5)
# Agkistrodon_laticinctus (Burbrink et al. 2015)
tree<-bind.tip2(tree,"Agkistrodon_laticinctus", sp1="Agkistrodon_contortrix", perc=0.5)
# Agkistrodon_russeolus (Porras et al. 2013)
tree<-bind.tip2(tree,"Agkistrodon_russeolus", sp1="Agkistrodon_howardgloydi", perc=0.5)
# Bothriechis_nubestris (Mason et al. 2019)
tree<-bind.tip2(tree,"Bothriechis_nubestris", sp1="Bothriechis_nigroviridis", perc=0.5)
# Bothrocophias_andianus (Timms et al. 2019)
tree<-bind.tip2(tree, tip.label="Bothrocophias_andianus", sp1="Bothrocophias_microphthalmus", sp2="Bothrocophias_hyoprora", perc=0.5)
# Bothrocophias_colombianus (unknown position)
##############
# Bothrocophias_myersi (Fenwick et al. 2009)
tree<-bind.tip2(tree,"Bothrocophias_myersi", sp1="Bothrocophias_campbelli", perc=0.5)
# Bothrops_ayerbei (Salazar-Valenzuela 2016)
tree<-bind.tip2(tree,"Bothrops_ayerbei", sp1="Bothrops_asper", perc=0.5)
# Bothrops_barnetti (Timms et al. 2019)
tree<-bind.tip2(tree,"Bothrops_barnetti", sp1="Bothrops_asper", sp2="Bothrops_lanceolatus", perc=0.5)
# Bothrops_jonathani (Carrasco et al. 2012)
tree<-bind.tip2(tree,"Bothrops_jonathani", sp1="Bothrops_alternatus", perc=0.5)
# Bothrops_lutzi (Carrasco et al. 2019)
tree<-bind.tip2(tree,"Bothrops_lutzi", sp1="Bothrops_erythromelas", perc=0.5)
# Bothrops_marmoratus  (Carrasco et al. 2019)
tree<-bind.tip2(tree,"Bothrops_marmoratus", sp1="Bothrops_neuwiedi", perc=0.5)
# Bothrops_medusa (unknown position)
##############
# Bothrops_monsignifer (Timms et al. 2019)
tree<-bind.tip2(tree,"Bothrops_monsignifer", sp1="Bothrops_jararaca", sp2="Bothrops_neuwiedi", perc=0.5)
# Bothrops_muriciensis (Rautsaw unpublished data)
tree<-bind.tip2(tree,"Bothrops_muriciensis", sp1="Bothrops_jararacussu", perc=0.5)
# Bothrops_oligolepis (Timms et al. 2019)
tree<-bind.tip2(tree,"Bothrops_oligolepis", sp1="Bothrops_chloromelas", perc=0.5)
# Bothrops_otavioi (Barbo et al. 2016)
tree<-bind.tip2(tree,"Bothrops_otavioi", sp1="Bothrops_insularis", perc=0.5)
# Bothrops_pirajai (Rautsaw unpublished data)
tree<-bind.tip2(tree,"Bothrops_pirajai", sp1="Bothrops_muriciensis", perc=0.5)
# Bothrops_sanctaecrucis (Timms et al. 2019)
tree<-bind.tip2(tree,"Bothrops_sanctaecrucis", sp1="Bothrops_monsignifer", perc=0.5)
# Bothrops_sazimai (Barbo et al. 2016)
tree<-bind.tip2(tree,"Bothrops_sazimai", sp1="Bothrops_insularis", perc=0.5)
# Bothrops_sonene (Carrasco et al. 2019)
tree<-bind.tip2(tree,"Bothrops_sonene", sp1="Bothrops_diporus", sp2="Bothrops_pauloensis", perc=0.5)
# Bothrops_venezuelensis (Timms et al. 2019)
tree<-bind.tip2(tree,"Bothrops_venezuelensis", sp1="Bothrops_barnetti", sp2="Bothrops_brazili", perc=0.5)
# Cerrophidion_sasai (Jadin et al. 2012)
tree<-bind.tip2(tree,"Cerrophidion_sasai", sp1="Cerrophidion_wilsoni", perc=0.5)
# Crotalus_angelensis (Meik et al. 2015)
tree<-bind.tip2(tree,"Crotalus_angelensis", sp1="Crotalus_mitchellii", perc=0.5)
# Crotalus_armstrongi (Blair et al. 2018)
tree<-bind.tip2(tree,"Crotalus_armstrongi", sp1="Crotalus_pusillus", perc=0.5)
# Crotalus_exiguus (Blair et al. 2018)
tree<-bind.tip2(tree,"Crotalus_exiguus", sp1="Crotalus_ravus", perc=0.5)
# Crotalus_brunneus (Blair et al. 2018)
tree<-bind.tip2(tree,"Crotalus_brunneus", sp1="Crotalus_exiguus", sp2="Crotalus_ravus", perc=0.5)
# Crotalus_campbelli (Blair et al. 2018)
tree<-bind.tip2(tree,"Crotalus_campbelli", sp1="Crotalus_armstrongi", perc=0.5)
# Crotalus_concolor (Schield et al. 2019)
tree<-bind.tip2(tree,"Crotalus_concolor", sp1="Crotalus_oreganus", perc=0.5)
# Crotalus_culminatus (Carbajal-Marquez et al. 2020)
tree<-bind.tip2(tree,"Crotalus_culminatus", sp1="Crotalus_simus", sp2="Crotalus_durissus", perc=0.5)
# Crotalus_ehecatl (Carbajal-Marquez et al. 2020)
tree<-bind.tip2(tree,"Crotalus_ehecatl", sp1="Crotalus_culminatus", perc=0.5)
# Crotalus_estebanensis (Ruiz-Sanchez et al. 2019)
tree<-bind.tip2(tree,"Crotalus_estebanensis", sp1="Crotalus_basiliscus", perc=0.5)
# Crotalus_helleri (Schield et al. 2019)
tree<-bind.tip2(tree,"Crotalus_helleri", sp1="Crotalus_concolor", perc=0.5)
# Crotalus_lorenzoensis (Ruiz-Sanchez et al. 2019)
tree<-bind.tip2(tree,"Crotalus_lorenzoensis", sp1="Crotalus_exsul", perc=0.5)
# Crotalus_lutosus (Schield et al. 2019)
tree<-bind.tip2(tree,"Crotalus_lutosus", sp1="Crotalus_helleri", sp2="Crotalus_concolor", perc=0.5)
# Crotalus_mictlantecuhtli (Carbajal-Marquez et al. 2020)
tree<-bind.tip2(tree,"Crotalus_mictlantecuhtli", sp1="Crotalus_ehecatl", sp2="Crotalus_culminatus", perc=0.5)
# Crotalus_morulus (Holding et al. 2021)
tree<-bind.tip2(tree,"Crotalus_morulus", sp1="Crotalus_lepidus", perc=0.5)
# Crotalus_ornatus (Holding et al. 2021)
tree<-bind.tip2(tree,"Crotalus_ornatus", sp1="Crotalus_molossus", perc=0.5)
# Crotalus_polisi (Meik et al. 2018)
tree<-bind.tip2(tree,"Crotalus_polisi", sp1="Crotalus_angelensis", perc=0.5)
# Crotalus_pyrrhus (Meik et al. 2015)
tree<-bind.tip2(tree,"Crotalus_pyrrhus", sp1="Crotalus_angelensis", sp2="Crotalus_polisi", perc=0.5)
# Crotalus_stephensi (Holding et al. 2021)
tree<-bind.tip2(tree,"Crotalus_stephensi", sp1="Crotalus_pyrrhus", sp2="Crotalus_mitchellii", perc=0.5)
# Crotalus_thalassoporus (Meik et al. 2018)
tree<-bind.tip2(tree,"Crotalus_thalassoporus", sp1="Crotalus_polisi", perc=0.5)
# Crotalus_tlaloci (Blair et al. 2018)
tree<-bind.tip2(tree,"Crotalus_tlaloci", sp1="Crotalus_pusillus", perc=0.5)
# Crotalus_tzabcan (Carbajal-Marquez et al. 2020)
tree<-bind.tip2(tree,"Crotalus_tzabcan", sp1="Crotalus_durissus", sp2="Crotalus_simus", perc=0.5)
# Lachesis_melanocephala (Zamudio & Greene 1997)
tree<-bind.tip2(tree,"Lachesis_melanocephala", sp1="Lachesis_stenophrys", perc=0.5)
# Mixcoatlus_browni (Jadin et al. 2011)
tree<-bind.tip2(tree,"Mixcoatlus_browni", sp1="Mixcoatlus_barbouri", perc=0.5)
# Ophryacus_smaragdinus (Grunwald et al. 2015)
tree<-bind.tip2(tree,"Ophryacus_smaragdinus", sp1="Ophryacus_undulatus", perc=0.5)
# Ophryacus_sphenophrys (Grunwald et al. 2015)
tree<-bind.tip2(tree,"Ophryacus_sphenophrys", sp1="Ophryacus_smaragdinus", sp2="Ophryacus_undulatus", perc=0.5)
# Porthidium_volcanicum
tree<-bind.tip2(tree,"Porthidium_volcanicum", sp1="Porthidium_porrasi", perc=0.5)
# Sistrurus_tergeminus (Holding et al. 2021)
tree<-bind.tip2(tree,"Sistrurus_tergeminus", sp1="Sistrurus_catenatus", perc=0.5)

# Ignoring Old World Taxa
# Atheris_acuminata, Atheris_broadleyi, Atheris_hirsuta, Atheris_katangensis, Atheris_mabuensis, Atheris_rungweensis, Atheris_subocularis, Bitis_heraldica, Bitis_inornata, Causus_bilineatus, Causus_maculatus, Daboia_siamensis, Gloydius_blomhoffi, Gloydius_caraganus, Gloydius_caucasicus, Gloydius_himalayanus, Gloydius_rickmersi, Macrovipera_lebetinus, Montatheris_hindii, Montivipera_kuhrangica, Ovophis_convictus, Ovophis_makazayazaya, Pseudocerastes_urarachnoides, Trimeresurus_andalasensis, Trimeresurus_brongersmai, Trimeresurus_cardamomensis, Trimeresurus_gunaleni, Trimeresurus_honsonensis, Trimeresurus_labialis, Trimeresurus_macrolepis, Trimeresurus_mcgregori, Trimeresurus_nebularis, Trimeresurus_phuketensis, Trimeresurus_rubeus, Trimeresurus_sabahi, Trimeresurus_salazar, Trimeresurus_sichuanensis, Trimeresurus_strigatus, Tropidolaemus_huttoni, Tropidolaemus_laticinctus, Tropidolaemus_philippensis, Tropidolaemus_subannulatus, Vipera_altaica, Vipera_anatolica, Vipera_darevskii, Vipera_graeca, Vipera_monticola, Vipera_olguni, Vipera_transcaucasiana, Vipera_walser

#tree$edge.length<-rep(1,length(tree$edge.length))

write.tree(tree,"Viperidae.newick")
```

## Occurence Point Cleaner (occ\_cleaner.R)

Next I developed `occ_cleaner` which will do all the work for me
including overlapping the points with the distributions, calculating
phylogenetic distance, and filtering to minimize phylogenetic and
geographic distance.

To run this function you need:

1.  Range maps in class `sf` in crs: 3857 or other projected CRS
      - `ranges<-st_read(dsn="file.geojson")`
      - `ranges<-st_transform(ranges,crs=3857)`
2.  Phylogeny in class `phylo`
      - `tree<-read.tree("tree.newick")`
3.  Occurrence records in class `sf` in crs:3857 or other projected CRS
      - `points<-read.csv("file.csv")`
      - `points<-st_as_sf(points,coords=c("longitude","latitude"),
        crs=4326)`
      - `points<-st_transform(points, crs=3857)`

I also took this opportunity to combine all my distribution maps into
one file (`Viper_Distributions.geojson`). *WARNING: occ\_cleaner can
take a very long time with many records. To speed this up I ran
`occ_cleaner` on a super computer with 80 cores and 1500 GB of RAM*

``` r
source("occ_cleaner.R")

files<-list.files("./Final/Distributions",pattern=".geojson",full.names = T)
ranges<-st_read(dsn=files[1])
ranges$NAME_0<-NULL
st_write(tmp, files[1],delete_dsn = T)
for(i in 2:length(files)){
  tmp<-st_read(dsn=files[i])
  tmp$NAME_0<-NULL
  st_write(tmp, files[i],delete_dsn = T)
  ranges<-rbind(ranges,tmp)
}
rm(i, tmp,files)
st_write(ranges, "Final/Viper_Distributions.geojson")

ranges<-st_read(dsn = "Final/Viper_Distributions.geojson")
ranges<-st_transform(ranges,crs=3857)

tree<-read.tree("Viperidae.newick")

occ<-read_xlsx("occ_databases/z_combined_records/combined_records_v1.xlsx")
occ_spdf<-occ[complete.cases(occ$latitude),] # Remove NA
occ_spdf_4326<-st_as_sf(occ_spdf, coords=c("longitude","latitude"), crs=4326)
points<-st_transform(occ_spdf_4326,crs=3857)

overlap<-occ_cleaner(points2, ranges, tree, maxdist=50000, k=20, parallel=8, maxphydist=5,
                     pt_id="ID", pt_species="species2", range_species="Species")

write.csv(overlap, "occ_databases/z_combined_records/overlap.csv")
```

## Manual Check

Removed all idigbio records because they were too messy.
gbif.canonicalNames did not match actual species ID and alternate
columns contained chaotic additional information (*e.g.*, (Author Year),
Author Year, Genus genus species species). Final version –\> v4

# Final Occurrence Counts

``` r
occ<-read_xlsx("Final/combined_records_v4.xlsx")
occ_clean <- occ %>% filter(!grepl("POTENTIAL_DUBIOUS_RECORD",flag_detailed))
occ_nodups <- occ %>% group_by(final_species) %>% distinct(latitude,longitude)

write.csv(occ_clean, "autokuenm/combined_records_v4_clean.csv", row.names = F)

#kable(table(occ$final_species), col.names = c("Species","Count"))
occ_count <- as.data.frame(table(occ$final_species)) %>% rename(Species=Var1, Total=Freq)
occ_clean_count <- as.data.frame(table(occ_clean$final_species))  %>% rename(Species=Var1, Clean=Freq)
occ_nodups_count <- as.data.frame(table(occ_nodups$final_species))  %>% rename(Species=Var1, NoDuplicates=Freq)

combined_table<-Reduce(function(x, y) merge(x, y, all=TRUE), list(occ_count, occ_clean_count, occ_nodups_count))
rm(occ_count, occ_clean_count, occ_nodups_count)

kable(combined_table)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

Species

</th>

<th style="text-align:right;">

Total

</th>

<th style="text-align:right;">

Clean

</th>

<th style="text-align:right;">

NoDuplicates

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Agkistrodon\_bilineatus

</td>

<td style="text-align:right;">

331

</td>

<td style="text-align:right;">

314

</td>

<td style="text-align:right;">

182

</td>

</tr>

<tr>

<td style="text-align:left;">

Agkistrodon\_conanti

</td>

<td style="text-align:right;">

3615

</td>

<td style="text-align:right;">

3612

</td>

<td style="text-align:right;">

2737

</td>

</tr>

<tr>

<td style="text-align:left;">

Agkistrodon\_contortrix

</td>

<td style="text-align:right;">

24545

</td>

<td style="text-align:right;">

23446

</td>

<td style="text-align:right;">

16686

</td>

</tr>

<tr>

<td style="text-align:left;">

Agkistrodon\_howardgloydi

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Agkistrodon\_laticinctus

</td>

<td style="text-align:right;">

1955

</td>

<td style="text-align:right;">

1938

</td>

<td style="text-align:right;">

1009

</td>

</tr>

<tr>

<td style="text-align:left;">

Agkistrodon\_piscivorus

</td>

<td style="text-align:right;">

19964

</td>

<td style="text-align:right;">

19498

</td>

<td style="text-align:right;">

13758

</td>

</tr>

<tr>

<td style="text-align:left;">

Agkistrodon\_russeolus

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

93

</td>

<td style="text-align:right;">

80

</td>

</tr>

<tr>

<td style="text-align:left;">

Agkistrodon\_taylori

</td>

<td style="text-align:right;">

68

</td>

<td style="text-align:right;">

66

</td>

<td style="text-align:right;">

50

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_acuminata

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_anisolepis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_barbouri

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

19

</td>

<td style="text-align:right;">

13

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_broadleyi

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_ceratophora

</td>

<td style="text-align:right;">

59

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:right;">

34

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_chlorechis

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

29

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_desaixi

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_hirsuta

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_hispida

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

11

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_katangensis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_mabuensis

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_nitschei

</td>

<td style="text-align:right;">

48

</td>

<td style="text-align:right;">

45

</td>

<td style="text-align:right;">

25

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_rungweensis

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

13

</td>

<td style="text-align:right;">

15

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_squamigera

</td>

<td style="text-align:right;">

204

</td>

<td style="text-align:right;">

160

</td>

<td style="text-align:right;">

132

</td>

</tr>

<tr>

<td style="text-align:left;">

Atheris\_subocularis

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Atropoides\_picadoi

</td>

<td style="text-align:right;">

39

</td>

<td style="text-align:right;">

36

</td>

<td style="text-align:right;">

26

</td>

</tr>

<tr>

<td style="text-align:left;">

Azemiops\_feae

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_arietans

</td>

<td style="text-align:right;">

1726

</td>

<td style="text-align:right;">

1585

</td>

<td style="text-align:right;">

1416

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_armata

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

14

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_atropos

</td>

<td style="text-align:right;">

179

</td>

<td style="text-align:right;">

65

</td>

<td style="text-align:right;">

157

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_caudalis

</td>

<td style="text-align:right;">

233

</td>

<td style="text-align:right;">

212

</td>

<td style="text-align:right;">

184

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_cornuta

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

18

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_gabonica

</td>

<td style="text-align:right;">

399

</td>

<td style="text-align:right;">

244

</td>

<td style="text-align:right;">

259

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_heraldica

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_inornata

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_nasicornis

</td>

<td style="text-align:right;">

394

</td>

<td style="text-align:right;">

263

</td>

<td style="text-align:right;">

226

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_parviocula

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_peringueyi

</td>

<td style="text-align:right;">

49

</td>

<td style="text-align:right;">

43

</td>

<td style="text-align:right;">

41

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_rhinoceros

</td>

<td style="text-align:right;">

66

</td>

<td style="text-align:right;">

54

</td>

<td style="text-align:right;">

49

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_rubida

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

20

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_schneideri

</td>

<td style="text-align:right;">

17

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

14

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_worthingtoni

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

</tr>

<tr>

<td style="text-align:left;">

Bitis\_xeropaga

</td>

<td style="text-align:right;">

12

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_aurifer

</td>

<td style="text-align:right;">

76

</td>

<td style="text-align:right;">

67

</td>

<td style="text-align:right;">

57

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_bicolor

</td>

<td style="text-align:right;">

66

</td>

<td style="text-align:right;">

49

</td>

<td style="text-align:right;">

54

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_guifarroi

</td>

<td style="text-align:right;">

36

</td>

<td style="text-align:right;">

36

</td>

<td style="text-align:right;">

24

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_lateralis

</td>

<td style="text-align:right;">

313

</td>

<td style="text-align:right;">

297

</td>

<td style="text-align:right;">

235

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_marchi

</td>

<td style="text-align:right;">

110

</td>

<td style="text-align:right;">

75

</td>

<td style="text-align:right;">

62

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_nigroviridis

</td>

<td style="text-align:right;">

152

</td>

<td style="text-align:right;">

125

</td>

<td style="text-align:right;">

86

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_nubestris

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

25

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_rowleyi

</td>

<td style="text-align:right;">

40

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

30

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_schlegelii

</td>

<td style="text-align:right;">

1656

</td>

<td style="text-align:right;">

1582

</td>

<td style="text-align:right;">

1148

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_supraciliaris

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

22

</td>

<td style="text-align:right;">

17

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothriechis\_thalassinus

</td>

<td style="text-align:right;">

53

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:right;">

40

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrocophias\_andianus

</td>

<td style="text-align:right;">

82

</td>

<td style="text-align:right;">

67

</td>

<td style="text-align:right;">

44

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrocophias\_campbelli

</td>

<td style="text-align:right;">

147

</td>

<td style="text-align:right;">

115

</td>

<td style="text-align:right;">

78

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrocophias\_colombianus

</td>

<td style="text-align:right;">

45

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrocophias\_hyoprora

</td>

<td style="text-align:right;">

372

</td>

<td style="text-align:right;">

369

</td>

<td style="text-align:right;">

206

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrocophias\_lojanus

</td>

<td style="text-align:right;">

44

</td>

<td style="text-align:right;">

44

</td>

<td style="text-align:right;">

18

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrocophias\_microphthalmus

</td>

<td style="text-align:right;">

231

</td>

<td style="text-align:right;">

220

</td>

<td style="text-align:right;">

162

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrocophias\_myersi

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:right;">

20

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_alcatraz

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_alternatus

</td>

<td style="text-align:right;">

4673

</td>

<td style="text-align:right;">

4642

</td>

<td style="text-align:right;">

2021

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_ammodytoides

</td>

<td style="text-align:right;">

502

</td>

<td style="text-align:right;">

465

</td>

<td style="text-align:right;">

257

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_asper

</td>

<td style="text-align:right;">

4529

</td>

<td style="text-align:right;">

4378

</td>

<td style="text-align:right;">

2308

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_atrox

</td>

<td style="text-align:right;">

4234

</td>

<td style="text-align:right;">

4158

</td>

<td style="text-align:right;">

1564

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_ayerbei

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_barnetti

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_bilineatus

</td>

<td style="text-align:right;">

1091

</td>

<td style="text-align:right;">

1074

</td>

<td style="text-align:right;">

433

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_brazili

</td>

<td style="text-align:right;">

320

</td>

<td style="text-align:right;">

305

</td>

<td style="text-align:right;">

179

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_caribbaeus

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

11

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_chloromelas

</td>

<td style="text-align:right;">

34

</td>

<td style="text-align:right;">

34

</td>

<td style="text-align:right;">

17

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_cotiara

</td>

<td style="text-align:right;">

448

</td>

<td style="text-align:right;">

438

</td>

<td style="text-align:right;">

159

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_diporus

</td>

<td style="text-align:right;">

1524

</td>

<td style="text-align:right;">

1472

</td>

<td style="text-align:right;">

889

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_erythromelas

</td>

<td style="text-align:right;">

889

</td>

<td style="text-align:right;">

884

</td>

<td style="text-align:right;">

119

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_fonsecai

</td>

<td style="text-align:right;">

372

</td>

<td style="text-align:right;">

372

</td>

<td style="text-align:right;">

100

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_insularis

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_itapetiningae

</td>

<td style="text-align:right;">

1118

</td>

<td style="text-align:right;">

1118

</td>

<td style="text-align:right;">

235

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_jararaca

</td>

<td style="text-align:right;">

5256

</td>

<td style="text-align:right;">

5229

</td>

<td style="text-align:right;">

1759

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_jararacussu

</td>

<td style="text-align:right;">

1688

</td>

<td style="text-align:right;">

1648

</td>

<td style="text-align:right;">

600

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_jonathani

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

20

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_lanceolatus

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

27

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_leucurus

</td>

<td style="text-align:right;">

3894

</td>

<td style="text-align:right;">

3874

</td>

<td style="text-align:right;">

312

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_lutzi

</td>

<td style="text-align:right;">

190

</td>

<td style="text-align:right;">

188

</td>

<td style="text-align:right;">

57

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_marajoensis

</td>

<td style="text-align:right;">

475

</td>

<td style="text-align:right;">

475

</td>

<td style="text-align:right;">

48

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_marmoratus

</td>

<td style="text-align:right;">

178

</td>

<td style="text-align:right;">

167

</td>

<td style="text-align:right;">

77

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_mattogrossensis

</td>

<td style="text-align:right;">

1226

</td>

<td style="text-align:right;">

1157

</td>

<td style="text-align:right;">

345

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_medusa

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_monsignifer

</td>

<td style="text-align:right;">

15

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:right;">

13

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_moojeni

</td>

<td style="text-align:right;">

2203

</td>

<td style="text-align:right;">

2187

</td>

<td style="text-align:right;">

605

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_muriciensis

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_neuwiedi

</td>

<td style="text-align:right;">

1493

</td>

<td style="text-align:right;">

1476

</td>

<td style="text-align:right;">

370

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_oligolepis

</td>

<td style="text-align:right;">

16

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_osbornei

</td>

<td style="text-align:right;">

28

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_otavioi

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_pauloensis

</td>

<td style="text-align:right;">

1709

</td>

<td style="text-align:right;">

1401

</td>

<td style="text-align:right;">

460

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_pictus

</td>

<td style="text-align:right;">

45

</td>

<td style="text-align:right;">

43

</td>

<td style="text-align:right;">

26

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_pirajai

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

18

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_pubescens

</td>

<td style="text-align:right;">

1736

</td>

<td style="text-align:right;">

1717

</td>

<td style="text-align:right;">

582

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_pulcher

</td>

<td style="text-align:right;">

69

</td>

<td style="text-align:right;">

65

</td>

<td style="text-align:right;">

36

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_punctatus

</td>

<td style="text-align:right;">

71

</td>

<td style="text-align:right;">

67

</td>

<td style="text-align:right;">

47

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_sanctaecrucis

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

42

</td>

<td style="text-align:right;">

25

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_sazimai

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_sonene

</td>

<td style="text-align:right;">

93

</td>

<td style="text-align:right;">

93

</td>

<td style="text-align:right;">

46

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_taeniatus

</td>

<td style="text-align:right;">

376

</td>

<td style="text-align:right;">

373

</td>

<td style="text-align:right;">

205

</td>

</tr>

<tr>

<td style="text-align:left;">

Bothrops\_venezuelensis

</td>

<td style="text-align:right;">

39

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

20

</td>

</tr>

<tr>

<td style="text-align:left;">

Calloselasma\_rhodostoma

</td>

<td style="text-align:right;">

174

</td>

<td style="text-align:right;">

156

</td>

<td style="text-align:right;">

89

</td>

</tr>

<tr>

<td style="text-align:left;">

Causus\_bilineatus

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

11

</td>

</tr>

<tr>

<td style="text-align:left;">

Causus\_defilippii

</td>

<td style="text-align:right;">

249

</td>

<td style="text-align:right;">

234

</td>

<td style="text-align:right;">

148

</td>

</tr>

<tr>

<td style="text-align:left;">

Causus\_lichtensteinii

</td>

<td style="text-align:right;">

190

</td>

<td style="text-align:right;">

136

</td>

<td style="text-align:right;">

139

</td>

</tr>

<tr>

<td style="text-align:left;">

Causus\_maculatus

</td>

<td style="text-align:right;">

381

</td>

<td style="text-align:right;">

374

</td>

<td style="text-align:right;">

229

</td>

</tr>

<tr>

<td style="text-align:left;">

Causus\_resimus

</td>

<td style="text-align:right;">

64

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

39

</td>

</tr>

<tr>

<td style="text-align:left;">

Causus\_rhombeatus

</td>

<td style="text-align:right;">

758

</td>

<td style="text-align:right;">

690

</td>

<td style="text-align:right;">

523

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerastes\_cerastes

</td>

<td style="text-align:right;">

174

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

156

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerastes\_gasperettii

</td>

<td style="text-align:right;">

67

</td>

<td style="text-align:right;">

64

</td>

<td style="text-align:right;">

55

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerastes\_vipera

</td>

<td style="text-align:right;">

114

</td>

<td style="text-align:right;">

109

</td>

<td style="text-align:right;">

75

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerrophidion\_godmani

</td>

<td style="text-align:right;">

431

</td>

<td style="text-align:right;">

397

</td>

<td style="text-align:right;">

281

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerrophidion\_petlalcalensis

</td>

<td style="text-align:right;">

32

</td>

<td style="text-align:right;">

15

</td>

<td style="text-align:right;">

17

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerrophidion\_sasai

</td>

<td style="text-align:right;">

77

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

44

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerrophidion\_tzotzilorum

</td>

<td style="text-align:right;">

128

</td>

<td style="text-align:right;">

126

</td>

<td style="text-align:right;">

90

</td>

</tr>

<tr>

<td style="text-align:left;">

Cerrophidion\_wilsoni

</td>

<td style="text-align:right;">

303

</td>

<td style="text-align:right;">

294

</td>

<td style="text-align:right;">

147

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_adamanteus

</td>

<td style="text-align:right;">

5463

</td>

<td style="text-align:right;">

5327

</td>

<td style="text-align:right;">

3677

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_angelensis

</td>

<td style="text-align:right;">

81

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:right;">

30

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_aquilus

</td>

<td style="text-align:right;">

468

</td>

<td style="text-align:right;">

389

</td>

<td style="text-align:right;">

275

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_armstrongi

</td>

<td style="text-align:right;">

56

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

36

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_atrox

</td>

<td style="text-align:right;">

31440

</td>

<td style="text-align:right;">

31123

</td>

<td style="text-align:right;">

20007

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_basiliscus

</td>

<td style="text-align:right;">

818

</td>

<td style="text-align:right;">

778

</td>

<td style="text-align:right;">

537

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_brunneus

</td>

<td style="text-align:right;">

49

</td>

<td style="text-align:right;">

48

</td>

<td style="text-align:right;">

34

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_campbelli

</td>

<td style="text-align:right;">

32

</td>

<td style="text-align:right;">

31

</td>

<td style="text-align:right;">

21

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_catalinensis

</td>

<td style="text-align:right;">

119

</td>

<td style="text-align:right;">

97

</td>

<td style="text-align:right;">

65

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_cerastes

</td>

<td style="text-align:right;">

14612

</td>

<td style="text-align:right;">

14280

</td>

<td style="text-align:right;">

8349

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_cerberus

</td>

<td style="text-align:right;">

855

</td>

<td style="text-align:right;">

788

</td>

<td style="text-align:right;">

594

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_concolor

</td>

<td style="text-align:right;">

597

</td>

<td style="text-align:right;">

589

</td>

<td style="text-align:right;">

388

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_culminatus

</td>

<td style="text-align:right;">

151

</td>

<td style="text-align:right;">

142

</td>

<td style="text-align:right;">

114

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_durissus

</td>

<td style="text-align:right;">

5713

</td>

<td style="text-align:right;">

5622

</td>

<td style="text-align:right;">

1855

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_ehecatl

</td>

<td style="text-align:right;">

85

</td>

<td style="text-align:right;">

81

</td>

<td style="text-align:right;">

65

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_enyo

</td>

<td style="text-align:right;">

889

</td>

<td style="text-align:right;">

846

</td>

<td style="text-align:right;">

601

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_ericsmithi

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_estebanensis

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_exiguus

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_helleri

</td>

<td style="text-align:right;">

13028

</td>

<td style="text-align:right;">

12776

</td>

<td style="text-align:right;">

8513

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_horridus

</td>

<td style="text-align:right;">

12994

</td>

<td style="text-align:right;">

12824

</td>

<td style="text-align:right;">

8632

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_intermedius

</td>

<td style="text-align:right;">

306

</td>

<td style="text-align:right;">

288

</td>

<td style="text-align:right;">

160

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_lannomi

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_lepidus

</td>

<td style="text-align:right;">

4044

</td>

<td style="text-align:right;">

3987

</td>

<td style="text-align:right;">

2428

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_lorenzoensis

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

16

</td>

<td style="text-align:right;">

8

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_lutosus

</td>

<td style="text-align:right;">

5281

</td>

<td style="text-align:right;">

5254

</td>

<td style="text-align:right;">

3706

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_mictlantecuhtli

</td>

<td style="text-align:right;">

36

</td>

<td style="text-align:right;">

36

</td>

<td style="text-align:right;">

28

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_mitchellii

</td>

<td style="text-align:right;">

451

</td>

<td style="text-align:right;">

407

</td>

<td style="text-align:right;">

255

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_molossus

</td>

<td style="text-align:right;">

5196

</td>

<td style="text-align:right;">

5122

</td>

<td style="text-align:right;">

3611

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_morulus

</td>

<td style="text-align:right;">

285

</td>

<td style="text-align:right;">

226

</td>

<td style="text-align:right;">

152

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_oreganus

</td>

<td style="text-align:right;">

13937

</td>

<td style="text-align:right;">

12977

</td>

<td style="text-align:right;">

8471

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_ornatus

</td>

<td style="text-align:right;">

2981

</td>

<td style="text-align:right;">

2962

</td>

<td style="text-align:right;">

2195

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_polisi

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_polystictus

</td>

<td style="text-align:right;">

260

</td>

<td style="text-align:right;">

253

</td>

<td style="text-align:right;">

127

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_pricei

</td>

<td style="text-align:right;">

1291

</td>

<td style="text-align:right;">

1019

</td>

<td style="text-align:right;">

539

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_pusillus

</td>

<td style="text-align:right;">

103

</td>

<td style="text-align:right;">

96

</td>

<td style="text-align:right;">

32

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_pyrrhus

</td>

<td style="text-align:right;">

3409

</td>

<td style="text-align:right;">

3307

</td>

<td style="text-align:right;">

2742

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_ravus

</td>

<td style="text-align:right;">

552

</td>

<td style="text-align:right;">

469

</td>

<td style="text-align:right;">

389

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_ruber

</td>

<td style="text-align:right;">

8996

</td>

<td style="text-align:right;">

8718

</td>

<td style="text-align:right;">

4890

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_scutulatus

</td>

<td style="text-align:right;">

12405

</td>

<td style="text-align:right;">

12341

</td>

<td style="text-align:right;">

7845

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_simus

</td>

<td style="text-align:right;">

297

</td>

<td style="text-align:right;">

266

</td>

<td style="text-align:right;">

157

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_stejnegeri

</td>

<td style="text-align:right;">

49

</td>

<td style="text-align:right;">

39

</td>

<td style="text-align:right;">

30

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_stephensi

</td>

<td style="text-align:right;">

875

</td>

<td style="text-align:right;">

845

</td>

<td style="text-align:right;">

686

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_tancitarensis

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_thalassoporus

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_tigris

</td>

<td style="text-align:right;">

1243

</td>

<td style="text-align:right;">

1218

</td>

<td style="text-align:right;">

835

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_tlaloci

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

46

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_totonacus

</td>

<td style="text-align:right;">

99

</td>

<td style="text-align:right;">

86

</td>

<td style="text-align:right;">

76

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_transversus

</td>

<td style="text-align:right;">

67

</td>

<td style="text-align:right;">

62

</td>

<td style="text-align:right;">

52

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_triseriatus

</td>

<td style="text-align:right;">

771

</td>

<td style="text-align:right;">

651

</td>

<td style="text-align:right;">

594

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_tzabcan

</td>

<td style="text-align:right;">

215

</td>

<td style="text-align:right;">

203

</td>

<td style="text-align:right;">

167

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_unicolor

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_vegrandis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_viridis

</td>

<td style="text-align:right;">

18198

</td>

<td style="text-align:right;">

17670

</td>

<td style="text-align:right;">

9396

</td>

</tr>

<tr>

<td style="text-align:left;">

Crotalus\_willardi

</td>

<td style="text-align:right;">

1072

</td>

<td style="text-align:right;">

912

</td>

<td style="text-align:right;">

493

</td>

</tr>

<tr>

<td style="text-align:left;">

Daboia\_mauritanica

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:right;">

32

</td>

<td style="text-align:right;">

35

</td>

</tr>

<tr>

<td style="text-align:left;">

Daboia\_palaestinae

</td>

<td style="text-align:right;">

240

</td>

<td style="text-align:right;">

224

</td>

<td style="text-align:right;">

212

</td>

</tr>

<tr>

<td style="text-align:left;">

Daboia\_russelii

</td>

<td style="text-align:right;">

151

</td>

<td style="text-align:right;">

138

</td>

<td style="text-align:right;">

134

</td>

</tr>

<tr>

<td style="text-align:left;">

Daboia\_siamensis

</td>

<td style="text-align:right;">

68

</td>

<td style="text-align:right;">

66

</td>

<td style="text-align:right;">

66

</td>

</tr>

<tr>

<td style="text-align:left;">

Deinagkistrodon\_acutus

</td>

<td style="text-align:right;">

67

</td>

<td style="text-align:right;">

38

</td>

<td style="text-align:right;">

63

</td>

</tr>

<tr>

<td style="text-align:left;">

Echis\_carinatus

</td>

<td style="text-align:right;">

265

</td>

<td style="text-align:right;">

184

</td>

<td style="text-align:right;">

218

</td>

</tr>

<tr>

<td style="text-align:left;">

Echis\_coloratus

</td>

<td style="text-align:right;">

156

</td>

<td style="text-align:right;">

136

</td>

<td style="text-align:right;">

125

</td>

</tr>

<tr>

<td style="text-align:left;">

Echis\_khosatzkii

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Echis\_leucogaster

</td>

<td style="text-align:right;">

47

</td>

<td style="text-align:right;">

34

</td>

<td style="text-align:right;">

30

</td>

</tr>

<tr>

<td style="text-align:left;">

Echis\_ocellatus

</td>

<td style="text-align:right;">

154

</td>

<td style="text-align:right;">

153

</td>

<td style="text-align:right;">

120

</td>

</tr>

<tr>

<td style="text-align:left;">

Echis\_omanensis

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

46

</td>

</tr>

<tr>

<td style="text-align:left;">

Echis\_pyramidum

</td>

<td style="text-align:right;">

92

</td>

<td style="text-align:right;">

79

</td>

<td style="text-align:right;">

46

</td>

</tr>

<tr>

<td style="text-align:left;">

Eristicophis\_macmahoni

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Garthius\_chaseni

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_blomhoffii

</td>

<td style="text-align:right;">

83

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:right;">

58

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_brevicaudus

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_caraganus

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_caucasicus

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_halys

</td>

<td style="text-align:right;">

68

</td>

<td style="text-align:right;">

63

</td>

<td style="text-align:right;">

62

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_himalayanus

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

19

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_intermedius

</td>

<td style="text-align:right;">

540

</td>

<td style="text-align:right;">

531

</td>

<td style="text-align:right;">

429

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_rickmersi

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_shedaoensis

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_strauchi

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_tsushimaensis

</td>

<td style="text-align:right;">

12

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

12

</td>

</tr>

<tr>

<td style="text-align:left;">

Gloydius\_ussuriensis

</td>

<td style="text-align:right;">

2614

</td>

<td style="text-align:right;">

2479

</td>

<td style="text-align:right;">

2569

</td>

</tr>

<tr>

<td style="text-align:left;">

Hypnale\_hypnale

</td>

<td style="text-align:right;">

206

</td>

<td style="text-align:right;">

164

</td>

<td style="text-align:right;">

121

</td>

</tr>

<tr>

<td style="text-align:left;">

Hypnale\_nepa

</td>

<td style="text-align:right;">

28

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:right;">

19

</td>

</tr>

<tr>

<td style="text-align:left;">

Hypnale\_zara

</td>

<td style="text-align:right;">

51

</td>

<td style="text-align:right;">

39

</td>

<td style="text-align:right;">

32

</td>

</tr>

<tr>

<td style="text-align:left;">

Lachesis\_acrochorda

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

62

</td>

</tr>

<tr>

<td style="text-align:left;">

Lachesis\_melanocephala

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:right;">

21

</td>

</tr>

<tr>

<td style="text-align:left;">

Lachesis\_muta

</td>

<td style="text-align:right;">

1321

</td>

<td style="text-align:right;">

1291

</td>

<td style="text-align:right;">

668

</td>

</tr>

<tr>

<td style="text-align:left;">

Lachesis\_stenophrys

</td>

<td style="text-align:right;">

117

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

57

</td>

</tr>

<tr>

<td style="text-align:left;">

Macrovipera\_lebetinus

</td>

<td style="text-align:right;">

189

</td>

<td style="text-align:right;">

146

</td>

<td style="text-align:right;">

136

</td>

</tr>

<tr>

<td style="text-align:left;">

Macrovipera\_schweizeri

</td>

<td style="text-align:right;">

31

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

24

</td>

</tr>

<tr>

<td style="text-align:left;">

Metlapilcoatlus\_indomitus

</td>

<td style="text-align:right;">

15

</td>

<td style="text-align:right;">

12

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Metlapilcoatlus\_mexicanus

</td>

<td style="text-align:right;">

351

</td>

<td style="text-align:right;">

330

</td>

<td style="text-align:right;">

219

</td>

</tr>

<tr>

<td style="text-align:left;">

Metlapilcoatlus\_nummifer

</td>

<td style="text-align:right;">

115

</td>

<td style="text-align:right;">

87

</td>

<td style="text-align:right;">

90

</td>

</tr>

<tr>

<td style="text-align:left;">

Metlapilcoatlus\_occiduus

</td>

<td style="text-align:right;">

82

</td>

<td style="text-align:right;">

66

</td>

<td style="text-align:right;">

60

</td>

</tr>

<tr>

<td style="text-align:left;">

Metlapilcoatlus\_olmec

</td>

<td style="text-align:right;">

93

</td>

<td style="text-align:right;">

78

</td>

<td style="text-align:right;">

67

</td>

</tr>

<tr>

<td style="text-align:left;">

Mixcoatlus\_barbouri

</td>

<td style="text-align:right;">

57

</td>

<td style="text-align:right;">

57

</td>

<td style="text-align:right;">

39

</td>

</tr>

<tr>

<td style="text-align:left;">

Mixcoatlus\_browni

</td>

<td style="text-align:right;">

42

</td>

<td style="text-align:right;">

38

</td>

<td style="text-align:right;">

22

</td>

</tr>

<tr>

<td style="text-align:left;">

Mixcoatlus\_melanurus

</td>

<td style="text-align:right;">

70

</td>

<td style="text-align:right;">

62

</td>

<td style="text-align:right;">

52

</td>

</tr>

<tr>

<td style="text-align:left;">

Montatheris\_hindii

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_albizona

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_bornmuelleri

</td>

<td style="text-align:right;">

31

</td>

<td style="text-align:right;">

28

</td>

<td style="text-align:right;">

17

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_bulgardaghica

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_kuhrangica

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_latifii

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_raddei

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

17

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_wagneri

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Montivipera\_xanthina

</td>

<td style="text-align:right;">

71

</td>

<td style="text-align:right;">

64

</td>

<td style="text-align:right;">

47

</td>

</tr>

<tr>

<td style="text-align:left;">

Ophryacus\_smaragdinus

</td>

<td style="text-align:right;">

29

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

26

</td>

</tr>

<tr>

<td style="text-align:left;">

Ophryacus\_sphenophrys

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Ophryacus\_undulatus

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

101

</td>

<td style="text-align:right;">

88

</td>

</tr>

<tr>

<td style="text-align:left;">

Ovophis\_convictus

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Ovophis\_makazayazaya

</td>

<td style="text-align:right;">

75

</td>

<td style="text-align:right;">

74

</td>

<td style="text-align:right;">

74

</td>

</tr>

<tr>

<td style="text-align:left;">

Ovophis\_monticola

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

19

</td>

<td style="text-align:right;">

24

</td>

</tr>

<tr>

<td style="text-align:left;">

Ovophis\_okinavensis

</td>

<td style="text-align:right;">

34

</td>

<td style="text-align:right;">

29

</td>

<td style="text-align:right;">

30

</td>

</tr>

<tr>

<td style="text-align:left;">

Ovophis\_tonkinensis

</td>

<td style="text-align:right;">

28

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

21

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_arcosae

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:right;">

57

</td>

<td style="text-align:right;">

39

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_dunni

</td>

<td style="text-align:right;">

382

</td>

<td style="text-align:right;">

355

</td>

<td style="text-align:right;">

261

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_hespere

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

15

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_lansbergii

</td>

<td style="text-align:right;">

457

</td>

<td style="text-align:right;">

420

</td>

<td style="text-align:right;">

192

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_nasutum

</td>

<td style="text-align:right;">

645

</td>

<td style="text-align:right;">

572

</td>

<td style="text-align:right;">

428

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_ophryomegas

</td>

<td style="text-align:right;">

235

</td>

<td style="text-align:right;">

195

</td>

<td style="text-align:right;">

152

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_porrasi

</td>

<td style="text-align:right;">

62

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

37

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_volcanicum

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Porthidium\_yucatanicum

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:right;">

147

</td>

<td style="text-align:right;">

98

</td>

</tr>

<tr>

<td style="text-align:left;">

Proatheris\_superciliaris

</td>

<td style="text-align:right;">

24

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

16

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_cornutus

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_elegans

</td>

<td style="text-align:right;">

29

</td>

<td style="text-align:right;">

12

</td>

<td style="text-align:right;">

22

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_flavoviridis

</td>

<td style="text-align:right;">

86

</td>

<td style="text-align:right;">

74

</td>

<td style="text-align:right;">

35

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_jerdonii

</td>

<td style="text-align:right;">

43

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

32

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_kaulbacki

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_mangshanensis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_mucrosquamatus

</td>

<td style="text-align:right;">

1673

</td>

<td style="text-align:right;">

1639

</td>

<td style="text-align:right;">

1576

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_sieversorum

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_tokarensis

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

8

</td>

</tr>

<tr>

<td style="text-align:left;">

Protobothrops\_xiangchengensis

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Pseudocerastes\_fieldi

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:right;">

25

</td>

<td style="text-align:right;">

23

</td>

</tr>

<tr>

<td style="text-align:left;">

Pseudocerastes\_persicus

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

10

</td>

</tr>

<tr>

<td style="text-align:left;">

Pseudocerastes\_urarachnoides

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Sistrurus\_catenatus

</td>

<td style="text-align:right;">

3337

</td>

<td style="text-align:right;">

2770

</td>

<td style="text-align:right;">

1877

</td>

</tr>

<tr>

<td style="text-align:left;">

Sistrurus\_miliarius

</td>

<td style="text-align:right;">

6681

</td>

<td style="text-align:right;">

6519

</td>

<td style="text-align:right;">

4568

</td>

</tr>

<tr>

<td style="text-align:left;">

Sistrurus\_tergeminus

</td>

<td style="text-align:right;">

3376

</td>

<td style="text-align:right;">

3262

</td>

<td style="text-align:right;">

2347

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_albolabris

</td>

<td style="text-align:right;">

594

</td>

<td style="text-align:right;">

506

</td>

<td style="text-align:right;">

403

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_andalasensis

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_andersonii

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_borneensis

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

8

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_brongersmai

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_cantori

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_cardamomensis

</td>

<td style="text-align:right;">

29

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

25

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_erythrurus

</td>

<td style="text-align:right;">

26

</td>

<td style="text-align:right;">

17

</td>

<td style="text-align:right;">

21

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_flavomaculatus

</td>

<td style="text-align:right;">

220

</td>

<td style="text-align:right;">

178

</td>

<td style="text-align:right;">

112

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_gracilis

</td>

<td style="text-align:right;">

38

</td>

<td style="text-align:right;">

37

</td>

<td style="text-align:right;">

38

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_gramineus

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_gumprechti

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_gunaleni

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_hageni

</td>

<td style="text-align:right;">

17

</td>

<td style="text-align:right;">

16

</td>

<td style="text-align:right;">

14

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_honsonensis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_insularis

</td>

<td style="text-align:right;">

126

</td>

<td style="text-align:right;">

108

</td>

<td style="text-align:right;">

83

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_kanburiensis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_labialis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_macrolepis

</td>

<td style="text-align:right;">

16

</td>

<td style="text-align:right;">

13

</td>

<td style="text-align:right;">

14

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_macrops

</td>

<td style="text-align:right;">

64

</td>

<td style="text-align:right;">

46

</td>

<td style="text-align:right;">

42

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_malabaricus

</td>

<td style="text-align:right;">

135

</td>

<td style="text-align:right;">

95

</td>

<td style="text-align:right;">

100

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_malcolmi

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_mcgregori

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

15

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_nebularis

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_phuketensis

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_popeiorum

</td>

<td style="text-align:right;">

47

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

38

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_puniceus

</td>

<td style="text-align:right;">

29

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

29

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_purpureomaculatus

</td>

<td style="text-align:right;">

121

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:right;">

85

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_rubeus

</td>

<td style="text-align:right;">

17

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:right;">

16

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_sabahi

</td>

<td style="text-align:right;">

59

</td>

<td style="text-align:right;">

36

</td>

<td style="text-align:right;">

44

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_salazar

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_schultzei

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_septentrionalis

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_sichuanensis

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_stejnegeri

</td>

<td style="text-align:right;">

1487

</td>

<td style="text-align:right;">

1459

</td>

<td style="text-align:right;">

1278

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_strigatus

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_sumatranus

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_tibetanus

</td>

<td style="text-align:right;">

11

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_trigonocephalus

</td>

<td style="text-align:right;">

120

</td>

<td style="text-align:right;">

76

</td>

<td style="text-align:right;">

87

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_venustus

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_vogeli

</td>

<td style="text-align:right;">

60

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

51

</td>

</tr>

<tr>

<td style="text-align:left;">

Trimeresurus\_yunnanensis

</td>

<td style="text-align:right;">

13

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Tropidolaemus\_huttoni

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Tropidolaemus\_philippensis

</td>

<td style="text-align:right;">

20

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

16

</td>

</tr>

<tr>

<td style="text-align:left;">

Tropidolaemus\_subannulatus

</td>

<td style="text-align:right;">

227

</td>

<td style="text-align:right;">

69

</td>

<td style="text-align:right;">

202

</td>

</tr>

<tr>

<td style="text-align:left;">

Tropidolaemus\_wagleri

</td>

<td style="text-align:right;">

175

</td>

<td style="text-align:right;">

157

</td>

<td style="text-align:right;">

146

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_altaica

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_ammodytes

</td>

<td style="text-align:right;">

717

</td>

<td style="text-align:right;">

689

</td>

<td style="text-align:right;">

635

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_anatolica

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_aspis

</td>

<td style="text-align:right;">

15099

</td>

<td style="text-align:right;">

14518

</td>

<td style="text-align:right;">

5384

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_barani

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_berus

</td>

<td style="text-align:right;">

40482

</td>

<td style="text-align:right;">

35516

</td>

<td style="text-align:right;">

16846

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_darevskii

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_dinniki

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_eriwanensis

</td>

<td style="text-align:right;">

64

</td>

<td style="text-align:right;">

43

</td>

<td style="text-align:right;">

11

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_graeca

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_kaznakovi

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

12

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_latastei

</td>

<td style="text-align:right;">

1866

</td>

<td style="text-align:right;">

1756

</td>

<td style="text-align:right;">

1260

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_lotievi

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

6

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_monticola

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_orlovi

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

3

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_renardi

</td>

<td style="text-align:right;">

682

</td>

<td style="text-align:right;">

622

</td>

<td style="text-align:right;">

596

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_seoanei

</td>

<td style="text-align:right;">

681

</td>

<td style="text-align:right;">

649

</td>

<td style="text-align:right;">

576

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_transcaucasiana

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

12

</td>

<td style="text-align:right;">

18

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_ursinii

</td>

<td style="text-align:right;">

1439

</td>

<td style="text-align:right;">

1405

</td>

<td style="text-align:right;">

178

</td>

</tr>

<tr>

<td style="text-align:left;">

Vipera\_walser

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

10

</td>

</tr>

</tbody>

</table>

### Updating Shiny App Data

``` r
load("/Users/rhettrautsaw/Dropbox/GitHub/ShinyServer/Apps/VenomMaps/data/shiny-data.RData")
combined_distribution<-st_read("~/Dropbox/Projects/2020_MacroCharDisp/VenomMaps/Final/Viper_Distributions.geojson")
occ<-readxl::read_xlsx("~/Dropbox/Projects/2020_MacroCharDisp/VenomMaps/Final/combined_records_v3.xlsx")

occ <- occ %>% group_by(final_species) %>% distinct(latitude,longitude, .keep_all = T)
save.image("~/Dropbox/GitHub/ShinyServer/Apps/VenomMaps/data/shiny-data_2021-09-13.RData")
```

# Summarizing ENM Results

## Averaging Selected Final Models

``` r
summarizeModels<-function(path){
  files<-list.files(Sys.glob(paste0(path,"/*")), pattern="*_avg.asc", full.names = T)
  
  rasters<-stack(files)
  crs(rasters)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  rasters_mean<-mean(rasters, na.rm=T)
  
  new_raster_name<-paste0(path,"/",gsub(".asc", ".tif",gsub(".*\\/","", files[1])))
  
  writeRaster(rasters_mean, new_raster_name, format="GTiff")

}


species<-list.dirs(recursive=F, full.names=F)[-143]

for(i in species){
  print(i)
  summarizeModels(paste0(i,"/Final_Models"))
  summarizeModels(paste0(i,"/Final_Models_bioclim"))
  summarizeModels(paste0(i,"/Final_Models_topvars"))
}
```

## Gathering Final Model Statistics

``` r
# Gathering all Maxent Results
files<-list.files(Sys.glob(paste0("*/Final_Models*/*")), pattern="maxentResults.csv", full.names = T)
results<-read.csv(files[1])
results<-results[nrow(results),]
results$file<-files[1]
for(i in 2:length(files)){
  tmp<-read.csv(files[i])
  tmp<-tmp[nrow(tmp),]
  tmp$file<-files[i]
  results<-merge(results, tmp, all=T)
}
write.csv(results, "Final_Results/allModelResults.csv")


# Gathering Final Model Evaluations
files<-list.files(Sys.glob(paste0("*/Final_Models*evaluation")), pattern="fm_evaluation_results.csv", full.names = T)
results<-read.csv(files[1])
results$file<-files[1]
for(i in 2:length(files)){
  tmp<-read.csv(files[i])
  tmp$file<-files[i]
  results<-merge(results, tmp, all=T)
}
write.csv(results, "Final_Results/fm_evaulation_results.csv")


# Gathering Calibration AIC results
files<-list.files(Sys.glob(paste0("*/Calibration_results*")), pattern="best_candidate_models_OR_AICc.csv", full.names = T)
results<-read.csv(files[1])
results$file<-files[1]
for(i in 2:length(files)){
  tmp<-read.csv(files[i])
  tmp$file<-files[i]
  results<-merge(results, tmp, all=T)
}
write.csv(results, "Final_Results/best_candidate_models_OR_AICc.csv", row.names = F)
```

## Summary Statistics

``` r
# ANALYZING/PLOTTING RESULTS
data<-read_xlsx("allModelResults.xlsx", col_types = c(rep("text",7), rep("numeric", 234)))
# Mean of multiple models
data_mean<-data %>% group_by(Species, testSet) %>% summarize_if(is.numeric, mean, na.rm=T)

# Filtering datasets
data_all<-data %>% filter(testSet=="all") 
data_bioclim<-data %>% filter(testSet=="bioclim") 
data_combo<-data %>% filter(testSet=="topvars")
data_mean_all<-data_mean %>% filter(testSet=="all") 
data_mean_bioclim<-data_mean %>% filter(testSet=="bioclim") 
data_mean_combo<-data_mean %>% filter(testSet=="topvars") 

# Determining which Environmental Set was the best
## Count Number of Times each Environmental Set was Selected
data_all %>% group_by(Species) %>% count(envSet) %>% count(envSet) %>% ungroup() %>% count(envSet)

## Count Mean Number of Models Selected per Environmental Set
data_all %>% group_by(Species) %>% count(envSet) %>% ungroup() %>% group_by(envSet) %>% summarize(mean=mean(n))

## Count Number of Samples per Selected Environmental Set
data_all %>% select(Species, envSet, X.Training.samples) %>% distinct() %>% group_by(envSet) %>% summarize(mean=mean(X.Training.samples))

# AUC Statistics
## All Models
summary(data_mean_all$Training.AUC)
data_all %>% group_by(envSet) %>% summarize(median=median(Training.AUC))

## Bioclim AUC statistics
summary(data_mean_bioclim$Training.AUC)

## Combo AUC statistics
summary(data_mean_combo$Training.AUC)
```

## Plotting AUC and Variable Contributions

``` r
## Summarize AUC
A <- ggviolin(data_all, x="envSet", y="Training.AUC", fill="envSet", palette="npg", trim = T, add="boxplot", 
              ylab="AUC", xlab="Environmental Data Sets", ylim=c(0.5, 1)) + scale_x_discrete(labels=c('Combo', 'Bioclim', 'Topography', 'Landcover')) +
  theme_pubclean() + rremove("legend")

B <- ggviolin(data_mean, x="testSet", y="Training.AUC", fill="testSet", palette=c("gray", get_palette("npg", 2)), trim = T, add="boxplot",
              ylab="AUC", xlab="Comparative Data Sets", ylim=c(0.5, 1), order=c("all", "topvars", "bioclim")) + scale_x_discrete(labels=c('All', 'Combo-only', 'Bioclim-only')) +
  theme_pubclean() + rremove("legend")
  #annotate("text", x=1, y=1.01,  label=paste0(round(range(data_mean[data2$testSet=="all",]$Training.AUC), 3), collapse ="-")) + 
  #annotate("text", x=2, y=1.01,  label=paste0(round(range(data2[data2$testSet=="bioclim",]$Training.AUC), 3), collapse ="—")) + 
  #annotate("text", x=3, y=1.01,  label=paste0(round(range(data2[data2$testSet=="topvars",]$Training.AUC), 3), collapse ="—"))

# Environmental Contributions
# envContrib<-colMeans(data_mean_combo[,59:103], na.rm = T)
# names(envContrib)<-gsub(".contribution", "",colnames(data_mean_combo[,59:103]))
# envContrib<-as.data.frame(envContrib)
# colnames(envContrib)[1]<-"Model Contribution"
# 
# envPermu<-colMeans(data_mean_combo[,14:58], na.rm = T)
# names(envPermu)<-gsub(".permutation.importance", "", colnames(data_mean_combo[,14:58]))
# envPermu<-as.data.frame(envPermu)
# colnames(envPermu)[1]<-"Permutation Importance"

tmp<-as.matrix(data_mean_combo[,59:103])
tmp[is.nan(tmp)] <- 0
tmp<-as.data.frame(tmp)
envContrib<-colMeans(tmp[,1:45], na.rm = T)
names(envContrib)<-gsub(".contribution", "",colnames(data_mean_combo[,59:103]))
envContrib<-as.data.frame(envContrib)
colnames(envContrib)[1]<-"Model Contribution"

tmp<-as.matrix(data_mean_combo[,14:58])
tmp[is.nan(tmp)] <- 0
tmp<-as.data.frame(tmp)
envPermu<-colMeans(tmp[,1:45], na.rm = T)
names(envPermu)<-gsub(".permutation.importance", "", colnames(data_mean_combo[,14:58]))
envPermu<-as.data.frame(envPermu)
colnames(envPermu)[1]<-"Permutation Importance"

envContrib2<-cbind(envContrib, envPermu)
rownames(envContrib2) <- gsub("wc2.0_bio_30s_","bioclim_",gsub("_1KMmd_GMTED", "", gsub("consensus_full_class", "landcover",rownames(envContrib2))))

library(pheatmap)
library(ggplotify)
C<-as.ggplot(pheatmap(t(as.matrix(envContrib2)), cluster_rows = F, cluster_cols = F, angle_col=315))

(A+B)/C + plot_layout(heights = c(2,1)) + plot_annotation(tag_levels = 'A')
```

## Plotting Example Niche Models

``` r
example_plots<-c("Agkistrodon_bilineatus", "Bothriechis_schlegelii", "Bothrops_jararaca", "Cerrophidion_godmani", "Crotalus_cerastes", "Lachesis_muta", "Metlapilcoatlus_nummifer", "Porthidium_nasutum", "Sistrurus_miliarius")
for(i in 1:length(example_plots)){
  tmp_shp<-st_read(paste0("~/Dropbox/GitHub/ShinyServer/Apps/VenomMaps/data/distributions/", example_plots[i],".geojson"), crs=4326)
  tmp<-raster(paste0("Final_Models_topvars/",example_plots[i],"_avg.tif"))
  A<-ggplot() + layer_spatial(tmp) + scale_fill_gradientn("P(Occurrence)", colours = rev(terrain.colors(255)), na.value=NA, limits=c(0,1)) + new_scale_fill() +
    geom_sf(tmp_shp, mapping=aes(), fill=NA) + 
    labs(title=example_plots[i]) + theme_pubclean() +
    annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)
  
  assign(paste0("Fig2.",i),A)
}

(Fig2.1 + Fig2.2 + Fig2.3)/(Fig2.4 + Fig2.5 + Fig2.6)/(Fig2.7 + Fig2.8 + Fig2.9) + plot_layout(widths = c(1,1,1))
```

## Plotting Example Bioclim vs. Combo Model

``` r
bioclim<-raster("Final_Models_bioclim/Agkistrodon_piscivorus_avg.tif")
combo<-raster("Final_Models_topvars/Agkistrodon_piscivorus_avg.tif")
A<-ggplot() + layer_spatial(bioclim) + scale_fill_gradientn("P(Occurrence)", colours = rev(terrain.colors(255)), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  labs(title="Agkistrodon piscivorus bioclim") + theme_pubclean() +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)
B<-ggplot() + layer_spatial(combo) + scale_fill_gradientn("P(Occurrence)", colours = rev(terrain.colors(255)), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  labs(title="Agkistrodon piscivorus combo") + theme_pubclean() +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)

A+B
```

## Summarizing Mean Combo Models (Supplemental Table 1)

``` r
# Create Supplemental Summary Table
data3.1<-data_combo %>% count(Species)

data3.2<-data_combo %>% select(c(1,5,6)) %>% distinct() %>% 
  group_by(Species) %>% summarize_all(~paste(., collapse = ','))

#data3<-data %>% group_by(Species, testSet) %>% select(Species, testSet, envSet) %>% distinct() %>% summarize_if(is.character, paste, collapse="; ")

data4<-merge(data3.1, data3.2)
data4<-merge(data4, data_mean_combo, by="Species")

write.csv(data4, "Final_Model_Results.csv", row.names = F)
```

# References

  - Anderson CG, Greenbaum E. 2012. Phylogeography of Northern
    Populations of the Black-Tailed Rattlesnake (Crotalus molossus Baird
    And Girard, 1853), With the Revalidation of C. ornatus Hallowell,
    1854. Herpetological Monographs 26: 19–57.
  - Blair C, Bryson RW, Linkem CW, Lazcano D, Klicka J, McCormack JE.
    2018. Cryptic diversity in the Mexican highlands: Thousands of UCE
    loci help illuminate phylogenetic relationships, species limits and
    divergence times of montane rattlesnakes (Viperidae: Crotalus ).
    Molecular Ecology Resources 0–2.
  - Bryson RW, Linkem CW, Dorcas ME, Lathrop A, Jones JM, Alvarado-Díaz
    J, Grünwald CI, Murphy RW. 2014. Multilocus species delimitation in
    the Crotalus triseriatus species group (serpentes: Viperidae:
    Crotalinae), with the description of two new species. Zootaxa 3826:
    475–496.
  - Bryson RW, Murphy RW, Lathrop A, Lazcano-Villareal D. 2011.
    Evolutionary drivers of phylogeographical diversity in the highlands
    of Mexico: A case study of the Crotalus triseriatus species group of
    montane rattlesnakes. Journal of Biogeography 38: 697–710.
  - Campbell and Lamar. 2004. The Venomous Reptiles of the Western
    Hemisphere.
  - Carbajal-Márquez RA, Cedeño-Vázquez JR, Martínez-Arce A, Neri-Castro
    E, Machkour-M’Rabet SC. 2020. Accessing cryptic diversity in
    Neotropical rattlesnakes (Serpentes: Viperidae: Crotalus) with the
    description of two new species. Zootaxa 4729: 451–481.
  - Carrasco PA, Grazziotin FG, Cruz Farfán RS, Koch C, Antonio Ochoa J,
    Scrocchi GJ, Leynaud GC, Chaparro JC. 2019. A new species of
    Bothrops (Serpentes: Viperidae: Crotalinae) from Pampas del Heath,
    southeastern Peru, with comments on the systematics of the Bothrops
    neuwiedi species group. Zootaxa 4565: 301–344.
  - Dal Vechio F, Prates I, Grazziotin FG, Zaher H, Rodrigues MT. 2018.
    Phylogeography and historical demography of the arboreal pit viper
    Bothrops bilineatus (Serpentes, Crotalinae) reveal multiple
    connections between Amazonian and Atlantic rain forests. Journal of
    Biogeography 45: 2415–2426.
  - Davis MA, Douglas MR, Collyer ML, Douglas ME. 2016. Deconstructing a
    species-complex: geometric morphometric and molecular analyses
    define species in the Western Rattlesnake (Crotalus viridis). PLOS
    ONE 11: e0146166.
  - Diniz-Sousa R, Moraes J do N, Rodrigues-da-Silva TM, Oliveira CS,
    Caldeira CA da S. 2020. A brief review on the natural history,
    venomics and the medical importance of bushmaster (Lachesis) pit
    viper snakes. Toxicon: X 7: 100053.
  - Douglas ME, Douglas MR, Schuett GW, Porras LW, Thomason BL. 2007.
    Genealogical Concordance between Mitochondrial and Nuclear DNAs
    Supports Species Recognition of the Panamint Rattlesnake (Crotalus
    mitchellii stephensi). Copeia 2007: 920–932.
  - Grünwald CI, Jones JM, Franz-Chávez H, Ahumada-Carrillo IT. 2015. A
    new species of Ophryacus (Serpentes: Viperidae: Crotalinae) from
    eastern Mexico with comments of the taxonomy of related pitvipers.
    Mesoamerican Herpetology 2: 388–416.
  - Heimes. 2016. Snakes of Mexico.
  - Jadin RC, Smith EN, Campbell JA. 2011. Unravelling a tangle of
    Mexican serpents: A systematic revision of highland pitvipers.
    Zoological Journal of the Linnean Society 163: 943–958.
  - Mason AJ, Grazziotin FG, Zaher H, Lemmon AR, Moriarty Lemmon E,
    Parkinson CL. 2019. Reticulate evolution in nuclear Middle America
    causes discordance in the phylogeny of palm‐pitvipers (Viperidae:
    Bothriechis). Journal of Biogeography 46: 833–844.
  - Meik JM, Schaack S, Flores-Villela O, Streicher JW. 2018.
    Integrative taxonomy at the nexus of population divergence and
    speciation in insular speckled rattlesnakes. Journal of Natural
    History 52: 989–1016.
  - Meik JM, Streicher JW, Lawing AM, Flores-Villela O, Fujita MK. 2015.
    Limitations of climatic data for inferring species boundaries:
    Insights from speckled rattlesnakes. PLoS ONE 10: 1–19.
  - Pace
  - Reptile Database. 2020
  - Roll et al. 2017
  - Salazar-Valenzuela CD. 2016. Diversification in the Neotropics:
    Insights from Demographic and Phylogenetic Patterns of Lancehead
    Pitvipers (Bothrops spp.). The Ohio State University.
  - Saldarriaga-Córdoba M, Parkinson CL, Daza JM, Wüster W, Sasa M.
    2017. Phylogeography of the Central American lancehead Bothrops
    asper (SERPENTES: VIPERIDAE). PLoS ONE 12: 1–21.
  - Smith EN, Ferrari-Castro JA. 2008. A new species of jumping pitviper
    of the genus Atropoides (Serpentes: Viperidae: Crotalinae) from the
    Sierra de Botaderos and the Sierra la Muralla, Honduras. Zootaxa 68:
    57–68.
  - Timms J, Chaparro JC, Venegas PJ, Salazar-Valenzuela D, Scrooch G,
    Cuevas J, Leynaud G, Carrasco PA. 2019. A new species of pitviper of
    the genus Bothrops (Serpentes: Viperidae: Crotalinae) from the
    Central Andes of South America. Zootaxa 4656: 99–120.
  - Yañez-Arenas C, Castaño-Quintero S, Rioja-Nieto R, Rodríguez-Medina
    K, Chiappa-Carrara X. 2020. Assessing the Relative Role of
    Environmental Factors That Limit the Distribution of the Yucatan
    Rattlesnake (Crotalus tzabcan). Journal of Herpetology 54: 216.
