# files <- list.files("data/enms", "*.tif", full.names = T)

stack_diff_extents<-function(files){
  rasters <- lapply(files, terra::rast)
  extents <- lapply(rasters, terra::ext)
  combined_extent <- Reduce(terra::union, extents)
  
  extended <- lapply(rasters, terra::extend, y = combined_extent, fill = NA)
  do.call(c, extended)
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
