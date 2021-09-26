library(dplyr)
library(ape)
library(phytools)
library(sf)
library(nngeo)
library(doParallel)
library(foreach)
library(purrr)

# To run this function you need:
# 1. Range maps in class `sf` in crs: 3857 or other projected CRS
#     - ranges<-st_read(dsn="file.geojson")
#     - ranges<-st_transform(ranges,crs=3857)
# 2. Phylogeny in class `phylo`
#     - tree<-read.tree("tree.newick")
# 3. Occurrence records in class `sf` in crs:3857 or other projected CRS
#     - points<-read.csv("file.csv")
#     - points<-st_as_sf(points,coords=c("longitude","latitude"), crs=4326)
#     - points<-st_transform(points, crs=3857)


## Rhett's Setup Example Code
# setwd("~/Dropbox/Projects/2020_MacroCharDisp/Occurence_RangeMaps")
# files<-list.files("./Final/geojson",pattern=".geojson",full.names = T)
# ranges<-st_read(dsn="Final/geojson/Agkistrodon_bilineatus.geojson")
# for(i in 2:159){
#   tmp<-st_read(dsn=files[i])
#   ranges<-rbind(ranges,tmp)
# }
# ranges<-st_transform(ranges,crs=3857)
# 
# tree<-read.tree("NW_Crotalinae.newick")
# 
# occ<-read_csv("occ_databases/z_combined_records/combined_records_cleaned_v1_w_GJV.csv")
# occ_spdf<-occ[complete.cases(occ$latitude),] # Remove NA
# occ_spdf_4326<-st_as_sf(occ_spdf, coords=c("longitude","latitude"), crs=4326)
# # occ_spdf_4326<-occ_spdf_4326[1:50,] # Test set before doing the whole she-bang
# points<-st_transform(occ_spdf_4326,crs=3857)
# maxdist=50000
# k=20
# parallel=6
# maxphydist=5
# pt_id="ID"
# pt_species="new_species"
# range_species="Species"

occ_cleaner<-function(points, ranges, tree, maxdist=50000, k=20, parallel=4, maxphydist=5,
                      pt_id="ID", pt_species="species2", range_species="Species"){
  
  ranges$id<-rownames(ranges)
  
  # Calculate Geographic Overlap/Distance
    print(paste0("Determining nearest ", k," ranges within ", maxdist, " meters"))
    print("Please wait, this could take a while...")
    overlap<-st_nn(points, ranges, k=k, maxdist=maxdist, parallel=parallel, returnDist=T)
    overlap$nn[!sapply(overlap$nn, length)] <- NA   
    overlap$dist[!sapply(overlap$dist, length)] <- NA
    
    points2<-data.frame(points[rep(seq_len(nrow(points)), sapply(overlap$nn,length)), , drop = FALSE], row.names=NULL)
    
    overlap <- overlap %>% set_names("id", "geoDist") %>% map_df(., unlist) 
    overlap<-cbind(as.data.frame(ranges)[overlap$id,],overlap) %>% dplyr::select(-id, -geometry)
    
    overlap<-cbind(points2,overlap)
    rownames(overlap) <- NULL

  # Calculate Phylogenetic Divergence
    print("Calculating Phylogenetic Divergence between recorded and overlapping species")
    registerDoParallel(parallel-2)
    tmp<-foreach(i=1:nrow(overlap), .combine=rbind) %dopar% {
      if(overlap[i,pt_species] %in% tree$tip.label & overlap[i,range_species] %in% tree$tip.label){
        phyDist<-fastDist(tree,overlap[i,pt_species], overlap[i,range_species])/2
      }else{
        phyDist<-NA
      }
    }
    
    overlap<-cbind(overlap,phyDist=tmp)
    
  # Determine if pt_species matches range_species...
  # If not, determine which species represents the closest across geography and phylogeny.
  # Phylogenetic distance is minimized first. If multiple matches/ties, then geographi distance is minimized. 
  # If multiple matches/ties remain then the record is flagged to be manually checked. 
  # Records are also flagged if the phylogenetic distance to the updated species is greater than maxphydist
    print("Filtering to minimize phylogenetic and geographic distance")
    print(paste0("Flagging points with multiple hits or with a phylogenetic divergence greater than ", maxphydist))
    overlap2<-overlap %>% group_by_(pt_id) %>% 
      filter(phyDist==min(phyDist) | is.na(phyDist)) %>% 
      filter(geoDist==min(geoDist) | is.na(geoDist)) %>%
      filter(if_else((!!as.name(pt_species)) %in% (!!as.name(range_species)), (!!as.name(range_species))==(!!as.name(pt_species)), phyDist==min(phyDist))) %>% # need to add an OR NA statement here. Perhaps change to filter_if
      mutate(n=n(),
             phyDist=if_else((!!as.name(range_species))==(!!as.name(pt_species)), 0, phyDist),
             flag=if_else((n>1 | phyDist>maxphydist), "CHECK_THIS_RECORD", ""),
             change=if_else((!!as.name(range_species))!=(!!as.name(pt_species)), "SPECIES_UPDATED", ""),
             flag=if_else((phyDist==0 & geoDist > 0), if_else(flag=="", "POINT_OUTSIDE_DISTRIBUTION", paste(flag,"POINT_OUTSIDE_DISTRIBUTION",sep="; ")), flag))
    # overlap2<-overlap %>% group_by_(pt_id) %>% mutate(geoRank=rank(round(geoDist,2), ties.method = "max"),
    #                                                   phyRank=rank(round(phyDist), ties.method = "max"),
    #                                                   totRank=rank(geoRank+phyRank)) %>%
    #                                            filter(totRank==min(totRank)) %>%
    #                                            #filter(if_else((!!as.name(pt_species)) %in% (!!as.name(range_species)), (!!as.name(range_species))==(!!as.name(pt_species)), totRank==min(totRank))) %>%
    #                                            mutate(n=n(),
    #                                                   flag=if_else((n>1 | phyDist>maxphydist), "CHECK_THIS_RECORD", ""),
    #                                                   change=if_else((!!as.name(range_species))!=(!!as.name(pt_species)), "SPECIES_UPDATED", "")) %>%
    #                                            filter(if_else((!!as.name(pt_species)) %in% (!!as.name(range_species)), (!!as.name(range_species))==(!!as.name(pt_species)), phyDist==min(phyDist))) %>%
    #                                            mutate(n=n(),
    #                                                   flag=if_else((phyDist==0 & geoDist > 0), "POINT_OUTSIDE_DISTRIBUTION", flag)) %>%
    #                                            ungroup() %>% dplyr::select(-geoRank, -phyRank, -totRank)
    # This step has been added into
    # overlap2 <- overlap %>% group_by(ID) %>%
    #   filter(if_else(species2 %in% species3, species3==species2, phyDist==min(phyDist))) %>%
    #   mutate(n=n(),
    #          phyDist=if_else(species3==species2, 0, phyDist),
    #          flag=if_else(species3==species2, gsub("CHECK_THIS_RECORD", "", flag), flag))
    # 
    # overlap3 <- overlap2 %>% group_by(ID) %>%
    #   summarize(overlapping_species = paste(species3, collapse = ','))
    # 
    # overlap4 <- overlap2 %>% group_by(ID) %>%
    #   mutate(species3=if_else(n>1, species2, species3),
    #          subspecies3=if_else(n>1, subspecies2, subspecies3)) %>%
    #   distinct()
    # 
    # overlap5<-merge(overlap4,overlap3)
    # 
    # overlap6<-overlap5 %>% group_by(ID) %>%
    #   mutate(overlapping_species=if_else(n==1, "", overlapping_species),
    #          change=if_else(species3!=species2, "SPECIES_UPDATED", ""),
    #          flag=if_else(n>1, "POTENTIAL_DUBIOUS_RECORD; UNCERTAIN_SPECIES_ID; DISTRIBUTION_ERROR", flag))
  
  return(overlap2)
}

# overlap<-occ_cleaner(points, ranges, tree, maxdist=50000, k=20, parallel=8, maxphydist=5,
#                      pt_id="ID", pt_species="new_species", range_species="Species")

