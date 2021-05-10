files <- list.files("data", "*.asc", full.names = T)

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

combined_niches<-stack_diff_extents(files)

names(combined_niches)<-gsub("_avg","",names(combined_niches))
