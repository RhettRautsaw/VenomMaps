library(raster)

#############################
# Convert raw <--> logistic #
#############################
# R = raw output
# L = logistic output
# h = entropy

raw2log<-function(R,h){
  l<-(exp(h)*R)/(1+(exp(h)*R))
  return(l)
}

log2raw<-function(l,h){
  R=-((exp(-h)*l)/(l-1))
  return(R)
}

# raster<-raster("../data/enms/Agkistrodon_bilineatus_avg.tif")
# entropy<-readxl::read_xlsx("../supplemental_material/SupplementalTable1.xlsx", sheet=4) %>%
#   filter(Species=="Agkistrodon_bilineatus") %>% select(Entropy) %>% as.numeric()
# raster2<-log2raw(raster, entropy)

###########################################
# Convert continuous --> threshold binary #
###########################################
# r = raster
# t = threshold value (ranging from 0-1)

threshold<-function(r,t=0.5){
  m<-matrix(c(0,t,0,
              t,1,1), 
            ncol=3, byrow = TRUE)
  result<-reclassify(r,m)
  return(result)
}

# raster<-raster("../data/enms/Agkistrodon_bilineatus_avg.tif")
# raster2<-threshold(raster, 0.5)
