# files <- list.files("data/enms", "*.tif", full.names = T)

stack_diff_extents<-function(files){
  xmin = xmax = ymin = ymax = NA
  
  for(i in files){
    tmp<-raster::raster(i)
    e<-raster::extent(tmp)
    
    if(is.na(xmin) | e[1]<xmin){
      xmin=e[1]
    }
    if(is.na(xmax) | e[2]>xmax){
      xmax=e[2]
    }
    if(is.na(ymin) | e[3]<ymin){
      ymin=e[3]
    }
    if(is.na(ymax) | e[4]>ymax){
      ymax=e[4]
    }
  }
  
  e<-extent(c(xmin,xmax,ymin,ymax))
  
  tmp<-raster::raster(files[1])
  tmp<-extend(tmp,e,value=NA)
  
  for(i in files[-1]){
    tmp2<-raster::raster(i)
    tmp2<-extend(tmp2,e,value=NA)
    tmp<-stack(tmp,tmp2)
  }
  
  return(tmp)
}

# combined_niches<-stack_diff_extents(files)
# 
# names(combined_niches)<-gsub("_avg","",names(combined_niches))
# 
# 
# rm(files)
# rm(stack_diff_extents)













# library(raster)
# library(leaflet)
# library(sf)
# 
# setwd("~/Desktop/")
# 
# tmp<-raster("Crotalus_cerberus_avg.asc")
# 
# crs(tmp)<-CRS("+init=epsg:4326")
# 
# #brPal <- colorRampPalette(c('#008B00FF', '#008B0000', '#008B0000'), alpha=T)
# #brPal <- colorRampPalette(c('#008B00FF', '#008B0000'), alpha=T)
# 
# # Default raster plot = rev(terrain.colors(255))
# 
# brPal <- colorRampPalette(c('#00A600FF', '#61C500BF', '#E6E40280', '#ECB17640', '#F2F2F200'), alpha=T)
# pal <- brPal(255)
# 
# data<-cbind(1:255,rep(1,255))
# barplot(data[,2]~data[,1],col=pal,axes=F, axisnames=F,xlab="",ylab="",space=0,border=NA)
# 
# dist<-read_sf("~/Dropbox/Projects/2020_MacroCharDisp/Occurence_RangeMaps/Final/geojson/Crotalus_cerberus.geojson")
# 
# #leaflet() %>% addTiles() %>% addRasterImage(tmp, colors = rev(pal)) %>% addPolygons(data=dist)
# 
# leaflet() %>% addTiles()  %>% #addWMSTiles('http://ows.mundialis.de/services/service?', layers='TOPO-WMS', group="Topography") %>% 
#     addRasterImage(tmp, colors = rev(pal)) #%>% addPolygons(data=dist)
# 
# plot(tmp)

