This repository contains the code used to

1.  Update occurrence record taxonomy for New World Crotalinae from several sources.
2.  Construct updated digital distributions all 159 species of New World Crotalinae.
3.  Use distributions to further update occurrence record taxonomy.

# Table of Contents

-   [Prep](#prep)
    -   [Synonyms List](#synonyms-list)
-   [Collecting Occurrence Records](#collecting-occurrence-records)
    -   [Filtering Databases](#filtering-databases)
    -   [Combining Databases](#combining-databases)
    -   [Cleaning `subspecies` Column](#cleaning-subspecies-column)
    -   [Updating Taxonomy](#updating-taxonomy)
-   [Constructing Distribution Maps](#constructing-distribution-maps)
    -   [Combining Databases](#combining-databases-1)
    -   [qGIS](#qgis)
    -   [Smooth & Clip Distributions](#smooth--clip-distributions)
    -   [Extract Country Information for Each Species](#extract-country-information-for-each-species)
-   [Spatial Join Correction](#spatial-join-correction)
    -   [Reference Phylogeny (bind.tip2.R)](#reference-phylogeny-bindtip2r)
    -   [Occurence Point Cleaner (occ_cleaner.R)](#occurence-point-cleaner-occ_cleanerr)
    -   [Manual Check](#manual-check)
-   [Final Occurrence Counts](#final-occurrence-counts)
-   [Ecological Niche Modelling](#ecological-niche-modelling)
    -   [Spatial Bias File](#spatial-bias-file)
    -   [spThin criteria creation](#spthin-criteria-creation)
-   [Summarizing ENM Results](#summarizing-enm-results)
    -   [Gathering Final Models & Statistics](#gathering-final-models--statistics)
    -   [Analyzing Results](#analyzing-results)
    -   [Calibration Models](#calibration-models)
    -   [Final Models](#final-models)
    -   [Variable Contributions](#variable-contributions)
    -   [Plotting Example Niche Models](#plotting-example-niche-models)
    -   [Plotting Example Theshold Models](#plotting-example-theshold-models)
    -   [Plotting Example Bioclim vs. Combo Model](#plotting-example-bioclim-vs-combo-model)
    -   [Plotting Example Bioclim vs. Topo vs. Combo Models](#plotting-example-bioclim-vs-topo-vs-combo-models)
-   [References](#references)


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
library(rnaturalearth)
library(RColorBrewer)
library(patchwork)
```

## Synonyms List

I obtained a list of all synonyms for species in “Crotalinae” from the
[Reptile Database (May 2020)](http://www.reptile-database.org/) (Thanks
Peter Uetz!) and combined this with a full list of viper species to
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

`ID, species, subspecies, locality, latitude, longitude, accuracy, source, voucher, issues/notes`

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
combined_data<-Reduce(function(x, y) merge(x, y, all=TRUE), list(gbif2, brazil_atlas2, herpmapper2, bison2, bioweb2, custom2)) #idigbio2
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

combined_data3<-combined_data2 %>% filter(!str_detect(species," NA$")) %>% filter(!str_detect(species," sp.$")) %>% filter(!grepl("fossil",source))
print(paste("Count number of records:", nrow(combined_data3)))
print(paste("Count number of NW records:", nrow(combined_data3 %>% filter(species %in% unique(synonyms[synonyms$Geography=="NewWorld",]$Synonym)))))

# write_csv(combined_data2,"occ_databases/z_combined_records/combined_records_v0.csv")
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
  filter(!is.na(latitude))
  
# write_csv(combined_data5,"occ_databases/z_combined_records/combined_records_v0.1.csv")
```

Taxonomy was then manually checked and further cleaned
(combined_records_v1).

# Constructing Distribution Maps

Preliminary distribution maps were obtained from
[IUCN](https://www.iucnredlist.org/resources/spatial-data-download),
[Heimes 2016 (Snakes of
Mexico)](https://www.chimaira.de/herpetofauna-mexicana-vol-1-snakes-of-mexico.html?___store=english&___from_store=default)
provided by Christoph Grünwald and Jason Jones of
[HerpMX](http://herp.mx/), [Roll et
al. 2017](https://www.nature.com/articles/s41559-017-0332-2), and custom
distribution maps made previously.

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
  
  st_write(distribution_clip, paste0("Final/distributions/",names[i],".geojson"), layer=names[i], driver="GeoJSON")
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

#countries<-st_read("~/Dropbox/GIS_Shapefiles/Basemaps/GADM_Administrative/countries.shp")
countries<-ne_countries(scale=10, returnclass = "sf")

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

## Occurence Point Cleaner (occ_cleaner.R)

Next I developed `occ_cleaner` which will do all the work for me
including overlapping the points with the distributions, calculating
phylogenetic distance, and filtering to minimize phylogenetic and
geographic distance.

To run this function you need:

1.  Range maps in class `sf` in crs: 3857 or other projected CRS
    -   `ranges<-st_read(dsn="file.geojson")`
    -   `ranges<-st_transform(ranges,crs=3857)`
2.  Phylogeny in class `phylo`
    -   `tree<-read.tree("tree.newick")`
3.  Occurrence records in class `sf` in crs:3857 or other projected CRS
    -   `points<-read.csv("file.csv")`
    -   `points<-st_as_sf(points,coords=c("longitude","latitude"), crs=4326)`
    -   `points<-st_transform(points, crs=3857)`

I also took this opportunity to combine all my distribution maps into
one file (`Viper_Distributions.geojson`). *WARNING: occ_cleaner can take
a very long time with many records. To speed this up I ran `occ_cleaner`
on a super computer with 80 cores and 1500 GB of RAM*

``` r
source("occ_cleaner.R")

files<-list.files("./Final/Distributions",pattern=".geojson",full.names = T)
ranges<-st_read(dsn=files[1])
ranges$NAME_0<-NULL
#st_write(tmp, files[1],delete_dsn = T)
for(i in 2:length(files)){
  tmp<-st_read(dsn=files[i])
  tmp$NAME_0<-NULL
  #st_write(tmp, files[i],delete_dsn = T)
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
Author Year, Genus genus species species). Final version –> v4

# Final Occurrence Counts

``` r
occ<-read_xlsx("Final/combined_records_v4.xlsx")
occ_clean <- occ %>% filter(!grepl("POTENTIAL_DUBIOUS_RECORD",flag_detailed))

# print(paste("Final record count:", nrow(occ_clean)))
# print(paste("Final NW record count:", nrow(occ_clean %>% filter(final_species %in% unique(gsub(" ","_",synonyms[synonyms$Geography=="NewWorld",]$Species))))))
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

### Updating Shiny App Data

``` r
load("~/Dropbox/GitHub/ShinyServer/Apps/VenomMaps/shiny_support_material/shiny-data_2021-10-13.RData")
combined_distribution<-st_read("~/Dropbox/Projects/2020_MacroCharDisp/VenomMaps/Final/Viper_Distributions.geojson")
occ<-readxl::read_xlsx("~/Dropbox/Projects/2020_MacroCharDisp/VenomMaps/Final/combined_records_v4.xlsx")

occ <- occ %>% group_by(final_species) %>% distinct(latitude,longitude, .keep_all = T) %>% filter(!str_detect(ID,"^HM"))
save.image("~/Dropbox/GitHub/ShinyServer/Apps/VenomMaps/shiny_support_material/shiny-data_2021-10-13.RData")
```

# Ecological Niche Modelling

## Spatial Bias File

<https://scottrinnan.wordpress.com/2015/08/31/how-to-construct-a-bias-file-with-r-for-use-in-maxent-modeling/>

``` r
library(dismo) # interface with MaxEnt
library(raster) # spatial data manipulation
library(MASS) # for 2D kernel density function
library(magrittr) # for piping functionality, i.e., %>%
library(maptools) # reading shapefiles
library(tidyverse)
library(sf)
occ<-read_csv("autokuenm/combined_records_v4_clean.csv") %>% distinct(final_species, latitude, longitude) #mutate(fakeID="all") %>% dplyr::select(fakeID, everything()) %>%
env<-raster("autokuenm/tifs/wc2.0_bio_30s_01.tif")
occ2<-st_as_sf(occ, coords=c("longitude","latitude"), crs=4326)
occ3<-st_crop(occ2, extent(env))
occ.ras <- rasterize(as_Spatial(occ3), env, 1)
presences <- which(values(occ.ras) == 1)
pres.locs <- coordinates(occ.ras)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occ.ras), ncol(occ.ras)))
dens.ras <- raster(dens)
plot(dens.ras)
writeRaster(dens.ras, "autokuenm/SpatialBias.tif")

dens.ras<-raster("autokuenm/SpatialBias.tif")
new<-projectRaster(dens.ras, env)
new[is.na(new[])] <- 0 
# crs(dens.ras)<-crs(env)
# new<-resample(dens.ras, env)
bias_new<-round(new*1000, digits = 3)
writeRaster(bias_new, "autokuenm/SpatialBias_resamp.tif")
```

## spThin criteria creation

``` r
## First looked at the relationship between the distribution size and the number of occurrence records
suppressMessages(occ <- read_csv("combined_records_v4_clean.csv"))
dist<-sf::st_read("../Viper_Distributions.geojson",crs=4326, quiet=T)
dist$area<-as.numeric(st_area(st_transform(dist, crs=3857)))
tmp3<-st_set_geometry(dist,NULL)
tmp3=tmp3[,c("Species","area")] %>% group_by(Species) %>% summarize(area=sum(area))
tmp2=occ %>% group_by(final_species) %>% summarize(count=n()) %>% rename(Species=final_species)
tmp4=left_join(tmp2,tmp3)
mod<-lm(log(count)~log(area), tmp4)
ggscatter(tmp4, "area", "count", conf.int=T, add="reg.line", xlab="Log10(Area) [m^2]", ylab="Log10(Occurrences)") + xscale("log10", .format = T) + yscale("log10") + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))


## Then used the occurrence records to build a linear model with a maximum spatial thinning parameter of 25 km
tmp2=occ %>% group_by(final_species) %>% summarize(count=n())
tmp<-data.frame(V1=c(1,30000), V2=c(0,25))
mod<-lm(V2~log(V1), data=tmp)
as.numeric(floor(predict(mod, data.frame(V1=c(NO_OCC_RECORDS)))))
```

See information on
[autokuenm](https://github.com/RhettRautsaw/VenomMaps/tree/master/code/autokuenm).

Models were output from MaxEnt in raw format; however, have been
converted to cloglog format and further threshold using a 10th
Percentile Training Presence.

# Summarizing ENM Results

## Gathering Final Models & Statistics

``` bash
mkdir Final_Results
cd Final_Results
mkdir all bioclim topvars
cd ..
for i in `cat list.success`
do echo $i
cp $i/Final_Models/FinalModelResults.csv Final_Results/all/${i}_results.csv
cp $i/Final_Models/FinalModelResults_best.csv Final_Results/all/${i}_results_best.csv
cp $i/Final_Models/FinalModelResults_mean.csv Final_Results/all/${i}_results_mean.csv
cp $i/Final_Models/*.tif Final_Results/all/
cp $i/Final_Models_bioclim/FinalModelResults.csv Final_Results/bioclim/${i}_results.csv
cp $i/Final_Models_bioclim/FinalModelResults_best.csv Final_Results/bioclim/${i}_results_best.csv
cp $i/Final_Models_bioclim/FinalModelResults_mean.csv Final_Results/bioclim/${i}_results_mean.csv
cp $i/Final_Models_bioclim/*.tif Final_Results/bioclim/
cp $i/Final_Models_topvars/FinalModelResults.csv Final_Results/topvars/${i}_results.csv
cp $i/Final_Models_topvars/FinalModelResults_best.csv Final_Results/topvars/${i}_results_best.csv
cp $i/Final_Models_topvars/FinalModelResults_mean.csv Final_Results/topvars/${i}_results_mean.csv
cp $i/Final_Models_topvars/*.tif Final_Results/topvars/
done
```

``` r
files<-list.files(Sys.glob("Final_Results/*"), pattern="*results.csv", full.names = T)
results<-read.csv(files[1])
results$file<-files[1]
for(i in 2:length(files)){
  tmp<-read.csv(files[i])
  tmp$file<-files[i]
  results<-merge(results, tmp, all=T)
}
write.csv(results, "Final_Results/combined_results.csv", quote=F)

files_best<-list.files(Sys.glob("Final_Results/*"), pattern="*results_best.csv", full.names = T)
results<-read.csv(files_best[1])
results$file<-files_best[1]
for(i in 2:length(files_best)){
  tmp<-read.csv(files_best[i])
  tmp$file<-files_best[i]
  results<-merge(results, tmp, all=T)
}
write.csv(results, "Final_Results/combined_results_best.csv", quote=F)

files_mean<-list.files(Sys.glob("Final_Results/*"), pattern="*results_mean.csv", full.names = T)
results<-read.csv(files_mean[1])
results$file<-files_mean[1]
for(i in 2:length(files_mean)){
  tmp<-read.csv(files_mean[i])
  tmp$file<-files_mean[i]
  results<-merge(results, tmp, all=T)
}
write.csv(results, "Final_Results/combined_results_mean.csv", quote=F)
```

## Analyzing Results

``` r
data<-read_xlsx("2022-02_Final_Results/combined_results.xlsx", sheet=2)

# Find best compSet (i.e. THE FINAL MODELS)
data2=data %>% group_by(Species) %>% filter(Omission_rate_at_5. == min(Omission_rate_at_5.)) %>% filter(Mean_AUC_ratio==max(Mean_AUC_ratio))
# added this information to the avg_models sheet of the results xlsx
```

## Calibration Models

``` r
data<-read_xlsx("2022-02_Final_Results/combined_results.xlsx", sheet=2)

# Calculate summary stats for different env sets when comparing all possible models
data_all<-data %>% filter(compSet=="all")
## AUC
summary.stats<-data_all %>% group_by(envSet) %>% get_summary_stats(Training.AUC)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
A1<-ggviolin(data_all, x="envSet", y="Training.AUC", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="AUC", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  theme_pubclean() + rremove("legend")
A2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n","mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## AUC Ratio
summary.stats<-data_all %>% group_by(envSet) %>% get_summary_stats(Mean_AUC_ratio)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
B1<-ggviolin(data_all, x="envSet", y="Mean_AUC_ratio", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="AUC Ratio", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  theme_pubclean() + rremove("legend")
B2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Omission Rate
summary.stats<-data_all %>% group_by(envSet) %>% get_summary_stats(Omission_rate_at_5.)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
C1<-ggviolin(data_all, x="envSet", y="Omission_rate_at_5.", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="Omission Rate at 5%", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  theme_pubclean() + rremove("legend")
C2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Number of Training Records
summary.stats<-data_all %>% group_by(envSet) %>% get_summary_stats(X.Training.samples)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
D1<-ggviolin(data_all, x="envSet", y="X.Training.samples", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="Number of Training Records", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  yscale("log2") + theme_pubclean() + rremove("legend")
D2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

svg("2022-02_Final_Results/ceval_envSet.svg", width=15, height=6)
((A1+B1+C1+D1)+ plot_layout(ncol=4))/((A2+B2+C2+D2) + plot_layout(ncol=4)) + plot_layout(heights = c(5,1))
dev.off()


# Calculate summary stats for different comparative sets
## AUC
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(Training.AUC)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
A1<-ggviolin(data, x="compSet", y="Training.AUC", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="AUC", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  theme_pubclean() + rremove("legend")
A2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## AUC Ratio
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(Mean_AUC_ratio)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
B1<-ggviolin(data, x="compSet", y="Mean_AUC_ratio", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="AUC Ratio", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  theme_pubclean() + rremove("legend")
B2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Omission Rate
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(Omission_rate_at_5.)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
C1<-ggviolin(data, x="compSet", y="Omission_rate_at_5.", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="Omission Rate at 5%", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  theme_pubclean() + rremove("legend")
C2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Number of Training Records
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(X.Training.samples)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
D1<-ggviolin(data, x="compSet", y="X.Training.samples", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="Number of Training Records", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  yscale("log2") + theme_pubclean() + rremove("legend")
D2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

svg("2022-02_Final_Results/ceval_compSet.svg", width=15, height=6)
((A1+B1+C1+D1)+ plot_layout(ncol=4))/((A2+B2+C2+D2) + plot_layout(ncol=4)) + plot_layout(heights = c(5,1))
dev.off()
```

## Final Models

``` r
data<-read_xlsx("2022-02_Final_Results/combined_results.xlsx", sheet=3)

# Calculate summary stats for different env sets when comparing all possible models
data_final<-data %>% filter(finalModel=="TRUE")
## AUC
summary.stats<-data_final %>% group_by(envSet) %>% get_summary_stats(Training.AUC)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
A1<-ggviolin(data_final, x="envSet", y="Training.AUC", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="AUC", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  theme_pubclean() + rremove("legend")
A2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## AUC Ratio
summary.stats<-data_final %>% group_by(envSet) %>% get_summary_stats(Mean_AUC_ratio)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
B1<-ggviolin(data_final, x="envSet", y="Mean_AUC_ratio", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="AUC Ratio", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  theme_pubclean() + rremove("legend")
B2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Omission Rate
summary.stats<-data_final %>% group_by(envSet) %>% get_summary_stats(Omission_rate_at_5.)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
C1<-ggviolin(data_final, x="envSet", y="Omission_rate_at_5.", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="Omission Rate at 5%", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  theme_pubclean() + rremove("legend")
C2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Number of Training Records
summary.stats<-data_final %>% group_by(envSet) %>% get_summary_stats(X.Training.samples)
summary.stats$envSet <- factor(summary.stats$envSet, c("topvars","bioclim","topo"))
D1<-ggviolin(data_final, x="envSet", y="X.Training.samples", add="boxplot", 
             order = c("topvars","bioclim","topo"),fill="envSet", palette="npg", 
             ylab="Number of Training Records", xlab="Environmental Data Sets") +
  scale_x_discrete(labels=c('Combination', 'Bioclim', 'Topography')) +
  yscale("log2") + theme_pubclean() + rremove("legend")
D2<-ggsummarytable(summary.stats, x="envSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

svg("2022-02_Final_Results/feval_envSet.svg", width=15, height=6)
((A1+B1+C1+D1)+ plot_layout(ncol=4))/((A2+B2+C2+D2) + plot_layout(ncol=4)) + plot_layout(heights = c(5,1))
dev.off()


# Calculate summary stats for different comparative sets
## AUC
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(Training.AUC)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
A1<-ggviolin(data, x="compSet", y="Training.AUC", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="AUC", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  theme_pubclean() + rremove("legend")
A2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## AUC Ratio
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(Mean_AUC_ratio)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
B1<-ggviolin(data, x="compSet", y="Mean_AUC_ratio", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="AUC Ratio", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  theme_pubclean() + rremove("legend")
B2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Omission Rate
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(Omission_rate_at_5.)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
C1<-ggviolin(data, x="compSet", y="Omission_rate_at_5.", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="Omission Rate at 5%", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  theme_pubclean() + rremove("legend")
C2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

## Number of Training Records
summary.stats<-data %>% group_by(compSet) %>% get_summary_stats(X.Training.samples)
summary.stats$compSet <- factor(summary.stats$compSet, c("all","bioclim","topvars"))
D1<-ggviolin(data, x="compSet", y="X.Training.samples", add="boxplot", 
             order = c("all","bioclim","topvars"),fill="compSet", palette="npg", 
             ylab="Number of Training Records", xlab="Comparative Sets") +
  scale_x_discrete(labels=c('All', 'Combination Models\nOnly', 'Bioclim Models\nOnly')) +
  yscale("log2") + theme_pubclean() + rremove("legend")
D2<-ggsummarytable(summary.stats, x="compSet", 
                   y=c("n", "mean","sd", "median","iqr"), digits = 2) + 
  clean_table_theme()

svg("2022-02_Final_Results/feval_compSet.svg", width=15, height=6)
((A1+B1+C1+D1)+ plot_layout(ncol=4))/((A2+B2+C2+D2) + plot_layout(ncol=4)) + plot_layout(heights = c(5,1))
dev.off()
```

## Variable Contributions

``` r
data<-read_xlsx("2022-02_Final_Results/combined_results.xlsx", sheet=3) %>% filter(compSet=="topvars")

tmp<-as.matrix(data[,str_detect(colnames(data), "contribution")])
tmp[is.na(tmp)] <- 0
tmp<-as.data.frame(tmp)
envContrib<-colMeans(tmp, na.rm = T)
names(envContrib)<-gsub(".contribution", "",names(envContrib))
envContrib<-as.data.frame(envContrib)
colnames(envContrib)[1]<-"Model Contribution"

tmp<-as.matrix(data[,str_detect(colnames(data), "permutation")])
tmp[is.na(tmp)] <- 0
tmp<-as.data.frame(tmp)
envPermu<-colMeans(tmp, na.rm = T)
names(envPermu)<-gsub(".permutation.importance", "", names(envPermu))
envPermu<-as.data.frame(envPermu)
colnames(envPermu)[1]<-"Permutation Importance"

envContrib2<-merge(envContrib, envPermu, by=0) %>% column_to_rownames("Row.names")
rownames(envContrib2) <- gsub("wc2.0_bio_30s_","bioclim_",gsub("_1KMmd_GMTED", "", gsub("consensus_full_class", "landcover",rownames(envContrib2))))

library(pheatmap)
library(ggplotify)
svg("2022-02_Final_Results/envContrib.svg", height=3, width=15)
C<-as.ggplot(pheatmap(t(as.matrix(envContrib2)), cluster_rows = F, cluster_cols = F, angle_col=315))
dev.off()
```

## Plotting Example Niche Models

``` r
example_plots<-c("Agkistrodon_bilineatus", "Bothriechis_schlegelii", "Bothrops_jararaca", "Cerrophidion_godmani", "Crotalus_cerastes", "Lachesis_muta", "Metlapilcoatlus_nummifer", "Porthidium_nasutum", "Sistrurus_miliarius")

overview<-ne_countries(scale=110, continent = c("North America", "South America"), returnclass = "sf") %>% filter(admin!="Greenland")

overview_map<-ggplot() + geom_sf(overview, mapping=aes(), fill="gray75", color="gray75") + 
  theme_void()

countries<-ne_countries(scale = 10, returnclass = "sf", continent = c("North America", "South America")) 
world<-ne_coastline(scale=10, returnclass = "sf")

rasCol=colorRampPalette(brewer.pal(9,"Greens"))(100)

for(i in 1:length(example_plots)){
  tmp_shp<-st_read(paste0("~/Dropbox/GitHub/ShinyServer/Apps/VenomMaps/data/distributions/", example_plots[i],".geojson"), crs=4326)
  tmp<-raster(paste0("2022-02_Final_Results/final_models/",example_plots[i],"_avg_cloglog.tif"))
  xlim_p=c(st_bbox(tmp_shp)[1]-1, st_bbox(tmp_shp)[3]+1)
  ylim_p=c(st_bbox(tmp_shp)[2]-1, st_bbox(tmp_shp)[4]+1)
  A<-ggplot() + layer_spatial(tmp) + scale_fill_gradientn("P(Occurrence)", colours = rasCol, na.value=NA, limits=c(0,1)) + 
    new_scale_fill() + # rev(terrain.colors(255))
    geom_sf(countries, mapping=aes(), fill=NA, color="gray75") + 
    geom_sf(tmp_shp, mapping=aes(), fill=NA, color="black", ) + coord_sf(xlim=xlim_p, ylim=ylim_p, expand = FALSE) +
    labs(title=example_plots[i]) + theme_void() + theme(panel.border = element_rect(colour = "black", fill=NA, size=5))
    #annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)
  
  #assign(paste0("Fig2.",i),A)
  assign(paste0(example_plots[i]),A)
}

#(Fig2.1 + Fig2.2 + Fig2.3)/(Fig2.4 + Fig2.5 + Fig2.6)/(Fig2.7 + Fig2.8 + Fig2.9) + plot_layout(widths = c(1,1,1))

layout <- '
AADDEE
AADDEE
BBDDFF
BBDDFF
CCDDGG
CCDDGG
'

svg("2022-02_Final_Results/example_enms.svg", width=11, height=8.5)
wrap_plots(A=Crotalus_cerastes, B=Cerrophidion_godmani, C=Bothriechis_schlegelii, D=overview_map, E=Sistrurus_miliarius, F=Metlapilcoatlus_nummifer, G=Bothrops_jararaca, design = layout)
dev.off()
```

## Plotting Example Theshold Models

``` r
example_plots<-c("Agkistrodon_bilineatus", "Bothriechis_schlegelii", "Bothrops_jararaca", "Cerrophidion_godmani", "Crotalus_cerastes", "Lachesis_muta", "Metlapilcoatlus_nummifer", "Porthidium_nasutum", "Sistrurus_miliarius")

overview<-ne_countries(scale=110, continent = c("North America", "South America"), returnclass = "sf") %>% filter(admin!="Greenland")

overview_map<-ggplot() + geom_sf(overview, mapping=aes(), fill="gray75", color="gray75") + 
  theme_void()

countries<-ne_countries(scale = 10, returnclass = "sf", continent = c("North America", "South America")) 
world<-ne_coastline(scale=10, returnclass = "sf")

rasCol=colorRampPalette(brewer.pal(9,"Greens"))(100)

for(i in 1:length(example_plots)){
  tmp_shp<-st_read(paste0("~/Dropbox/GitHub/ShinyServer/Apps/VenomMaps/data/distributions/", example_plots[i],".geojson"), crs=4326)
  tmp<-raster(paste0("2022-02_Final_Results/final_models/",example_plots[i],"_avg_cloglog_p10b.tif"))
  xlim_p=c(st_bbox(tmp_shp)[1]-1, st_bbox(tmp_shp)[3]+1)
  ylim_p=c(st_bbox(tmp_shp)[2]-1, st_bbox(tmp_shp)[4]+1)
  A<-ggplot() + layer_spatial(tmp) + scale_fill_gradientn("P(Occurrence)", colours = rasCol, na.value=NA, limits=c(0,1)) + 
    new_scale_fill() + # rev(terrain.colors(255))
    geom_sf(countries, mapping=aes(), fill=NA, color="gray75") + 
    geom_sf(tmp_shp, mapping=aes(), fill=NA, color="black", ) + coord_sf(xlim=xlim_p, ylim=ylim_p, expand = FALSE) +
    labs(title=example_plots[i]) + theme_void() + theme(panel.border = element_rect(colour = "black", fill=NA, size=5))
    #annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering)
  
  #assign(paste0("Fig2.",i),A)
  assign(paste0(example_plots[i]),A)
}

#(Fig2.1 + Fig2.2 + Fig2.3)/(Fig2.4 + Fig2.5 + Fig2.6)/(Fig2.7 + Fig2.8 + Fig2.9) + plot_layout(widths = c(1,1,1))

layout <- '
AADDEE
AADDEE
BBDDFF
BBDDFF
CCDDGG
CCDDGG
'

svg("2022-02_Final_Results/example_enms_thresh.svg", width=11, height=8.5)
wrap_plots(A=Crotalus_cerastes, B=Cerrophidion_godmani, C=Bothriechis_schlegelii, D=overview_map, E=Sistrurus_miliarius, F=Metlapilcoatlus_nummifer, G=Bothrops_jararaca, design = layout)
dev.off()
```

## Plotting Example Bioclim vs. Combo Model

``` r
bioclim<-raster("2022-02_Final_Results/bioclim/Agkistrodon_piscivorus_avg_cloglog.tif")
combo<-raster("2022-02_Final_Results/topvars/Agkistrodon_piscivorus_avg_cloglog.tif")
rivers<-ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
rivers<-st_crop(rivers,extent(bioclim))

xlim_p=c(st_bbox(bioclim)[1], st_bbox(bioclim)[3])
ylim_p=c(st_bbox(bioclim)[2], st_bbox(bioclim)[4])

A<-ggplot() + layer_spatial(bioclim) + scale_fill_gradientn("P(Occurrence)", colours = colorRampPalette(brewer.pal(9,"Greens"))(100), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  geom_sf(countries, mapping=aes(), fill=NA, color="gray75") +
  coord_sf(xlim=xlim_p, ylim=ylim_p, expand=F) + 
  labs(title="Bioclim Model") + theme_void() + rremove("legend")
B<-ggplot() + layer_spatial(bioclim) + scale_fill_gradientn("P(Occurrence)", colours = colorRampPalette(brewer.pal(9,"Greens"))(100), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  geom_sf(countries, mapping=aes(), fill=NA, color="gray75") + geom_sf(rivers, mapping=aes(), color="black") + 
  coord_sf(xlim=xlim_p, ylim=ylim_p, expand=F)  + 
  labs(title="Bioclim Model + Rivers") + theme_void() + rremove("legend")
C<-ggplot() + layer_spatial(combo) + scale_fill_gradientn("P(Occurrence)", colours = colorRampPalette(brewer.pal(9,"Greens"))(100), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  geom_sf(countries, mapping=aes(), fill=NA, color="gray75") + 
  coord_sf(xlim=xlim_p, ylim=ylim_p, expand=F) + 
  labs(title="Combination Model") + theme_void() 

svg("2022-02_Final_Results/Apisc-envComp.svg", width=11, height=8.5)
overview_map+A+B+C+plot_layout(ncol = 4)
dev.off()
```

## Plotting Example Bioclim vs. Topo vs. Combo Models

Using M. nummifer as an example.

``` r
topo<-raster("2022-02_Final_Results/all/Metlapilcoatlus_nummifer_avg_cloglog.tif")
bioclim<-raster("2022-02_Final_Results/bioclim/Metlapilcoatlus_nummifer_avg_cloglog.tif")
combo<-raster("2022-02_Final_Results/topvars/Metlapilcoatlus_nummifer_avg_cloglog.tif")

xlim_p=c(st_bbox(bioclim)[1]-1, st_bbox(bioclim)[3]+1)
ylim_p=c(st_bbox(bioclim)[2]-1, st_bbox(bioclim)[4]+1)

A<-ggplot() + layer_spatial(topo) + scale_fill_gradientn("P(Occurrence)", colours = colorRampPalette(brewer.pal(9,"Greens"))(100), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  geom_sf(countries, mapping=aes(), fill=NA, color="gray75") +
  coord_sf(xlim=xlim_p, ylim=ylim_p, expand=F) + 
  labs(title="Topographic Model") + theme_void() + rremove("legend")
B<-ggplot() + layer_spatial(bioclim) + scale_fill_gradientn("P(Occurrence)", colours = colorRampPalette(brewer.pal(9,"Greens"))(100), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  geom_sf(countries, mapping=aes(), fill=NA, color="gray75") + 
  coord_sf(xlim=xlim_p, ylim=ylim_p, expand=F)  + 
  labs(title="Bioclim Model") + theme_void() + rremove("legend")
C<-ggplot() + layer_spatial(combo) + scale_fill_gradientn("P(Occurrence)", colours = colorRampPalette(brewer.pal(9,"Greens"))(100), na.value=NA, limits=c(0,1)) + new_scale_fill() +
  geom_sf(countries, mapping=aes(), fill=NA, color="gray75") + 
  coord_sf(xlim=xlim_p, ylim=ylim_p, expand=F) + 
  labs(title="Combination Model") + theme_void() 

svg("2022-02_Final_Results/Mnumm-envComp.svg", width=11, height=8.5)
overview_map+A+B+C+plot_layout(ncol = 4)
dev.off()
```

# References

-   Anderson CG, Greenbaum E. 2012. Phylogeography of Northern
    Populations of the Black-Tailed Rattlesnake (Crotalus molossus Baird
    And Girard, 1853), With the Revalidation of C. ornatus
    Hallowell, 1854. Herpetological Monographs 26: 19–57.
-   Blair C, Bryson RW, Linkem CW, Lazcano D, Klicka J, McCormack
    JE. 2018. Cryptic diversity in the Mexican highlands: Thousands of
    UCE loci help illuminate phylogenetic relationships, species limits
    and divergence times of montane rattlesnakes (Viperidae: Crotalus ).
    Molecular Ecology Resources 0–2.
-   Bryson RW, Linkem CW, Dorcas ME, Lathrop A, Jones JM, Alvarado-Díaz
    J, Grünwald CI, Murphy RW. 2014. Multilocus species delimitation in
    the Crotalus triseriatus species group (serpentes: Viperidae:
    Crotalinae), with the description of two new species. Zootaxa 3826:
    475–496.
-   Bryson RW, Murphy RW, Lathrop A, Lazcano-Villareal D. 2011.
    Evolutionary drivers of phylogeographical diversity in the highlands
    of Mexico: A case study of the Crotalus triseriatus species group of
    montane rattlesnakes. Journal of Biogeography 38: 697–710.
-   Campbell and Lamar. 2004. The Venomous Reptiles of the Western
    Hemisphere.
-   Carbajal-Márquez RA, Cedeño-Vázquez JR, Martínez-Arce A, Neri-Castro
    E, Machkour-M’Rabet SC. 2020. Accessing cryptic diversity in
    Neotropical rattlesnakes (Serpentes: Viperidae: Crotalus) with the
    description of two new species. Zootaxa 4729: 451–481.
-   Carrasco PA, Grazziotin FG, Cruz Farfán RS, Koch C, Antonio Ochoa J,
    Scrocchi GJ, Leynaud GC, Chaparro JC. 2019. A new species of
    Bothrops (Serpentes: Viperidae: Crotalinae) from Pampas del Heath,
    southeastern Peru, with comments on the systematics of the Bothrops
    neuwiedi species group. Zootaxa 4565: 301–344.
-   Dal Vechio F, Prates I, Grazziotin FG, Zaher H, Rodrigues MT. 2018.
    Phylogeography and historical demography of the arboreal pit viper
    Bothrops bilineatus (Serpentes, Crotalinae) reveal multiple
    connections between Amazonian and Atlantic rain forests. Journal of
    Biogeography 45: 2415–2426.
-   Davis MA, Douglas MR, Collyer ML, Douglas ME. 2016. Deconstructing a
    species-complex: geometric morphometric and molecular analyses
    define species in the Western Rattlesnake (Crotalus viridis). PLOS
    ONE 11: e0146166.
-   Diniz-Sousa R, Moraes J do N, Rodrigues-da-Silva TM, Oliveira CS,
    Caldeira CA da S. 2020. A brief review on the natural history,
    venomics and the medical importance of bushmaster (Lachesis) pit
    viper snakes. Toxicon: X 7: 100053.
-   Douglas ME, Douglas MR, Schuett GW, Porras LW, Thomason BL. 2007.
    Genealogical Concordance between Mitochondrial and Nuclear DNAs
    Supports Species Recognition of the Panamint Rattlesnake (Crotalus
    mitchellii stephensi). Copeia 2007: 920–932.
-   Grünwald CI, Jones JM, Franz-Chávez H, Ahumada-Carrillo IT. 2015. A
    new species of Ophryacus (Serpentes: Viperidae: Crotalinae) from
    eastern Mexico with comments of the taxonomy of related pitvipers.
    Mesoamerican Herpetology 2: 388–416.
-   Heimes. 2016. Snakes of Mexico.
-   Jadin RC, Smith EN, Campbell JA. 2011. Unravelling a tangle of
    Mexican serpents: A systematic revision of highland pitvipers.
    Zoological Journal of the Linnean Society 163: 943–958.
-   Mason AJ, Grazziotin FG, Zaher H, Lemmon AR, Moriarty Lemmon E,
    Parkinson CL. 2019. Reticulate evolution in nuclear Middle America
    causes discordance in the phylogeny of palm‐pitvipers (Viperidae:
    Bothriechis). Journal of Biogeography 46: 833–844.
-   Meik JM, Schaack S, Flores-Villela O, Streicher JW. 2018.
    Integrative taxonomy at the nexus of population divergence and
    speciation in insular speckled rattlesnakes. Journal of Natural
    History 52: 989–1016.
-   Meik JM, Streicher JW, Lawing AM, Flores-Villela O, Fujita MK. 2015.
    Limitations of climatic data for inferring species boundaries:
    Insights from speckled rattlesnakes. PLoS ONE 10: 1–19.
-   Pace
-   Reptile Database. 2020
-   Roll et al. 2017
-   Salazar-Valenzuela CD. 2016. Diversification in the Neotropics:
    Insights from Demographic and Phylogenetic Patterns of Lancehead
    Pitvipers (Bothrops spp.). The Ohio State University.
-   Saldarriaga-Córdoba M, Parkinson CL, Daza JM, Wüster W, Sasa M.
    2017. Phylogeography of the Central American lancehead Bothrops
    asper (SERPENTES: VIPERIDAE). PLoS ONE 12: 1–21.
-   Smith EN, Ferrari-Castro JA. 2008. A new species of jumping pitviper
    of the genus Atropoides (Serpentes: Viperidae: Crotalinae) from the
    Sierra de Botaderos and the Sierra la Muralla, Honduras. Zootaxa 68:
    57–68.
-   Timms J, Chaparro JC, Venegas PJ, Salazar-Valenzuela D, Scrooch G,
    Cuevas J, Leynaud G, Carrasco PA. 2019. A new species of pitviper of
    the genus Bothrops (Serpentes: Viperidae: Crotalinae) from the
    Central Andes of South America. Zootaxa 4656: 99–120.
-   Yañez-Arenas C, Castaño-Quintero S, Rioja-Nieto R, Rodríguez-Medina
    K, Chiappa-Carrara X. 2020. Assessing the Relative Role of
    Environmental Factors That Limit the Distribution of the Yucatan
    Rattlesnake (Crotalus tzabcan). Journal of Herpetology 54: 216.
