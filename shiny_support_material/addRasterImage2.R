epsg4326 <- "+proj=longlat +datum=WGS84 +no_defs"
epsg3857 <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs" # nolint

addRasterImage2 <- function(
  map,
  x,
  colors = if (raster::is.factor(x)) "Set1" else "Spectral",
  opacity = 1,
  attribution = NULL,
  layerId = NULL,
  group = NULL,
  project = FALSE,
  method = c("auto", "bilinear", "ngb"),
  maxBytes = 4 * 1024 * 1024,
  options = tileOptions(),
  data = getMapData(map)
) {
  stopifnot(inherits(x, "RasterLayer"))
  
  options$opacity <- opacity
  options$attribution <- attribution
  
  raster_is_factor <- raster::is.factor(x)
  method <- match.arg(method)
  if (method == "auto") {
    if (raster_is_factor) {
      method <- "ngb"
    } else {
      method <- "bilinear"
    }
  }
  
  if (project) {
    # if we should project the data
    projected <- projectRasterForLeaflet(x, method)
    
    # if data is factor data, make the result factors as well.
    if (raster_is_factor) {
      projected <- raster::as.factor(projected)
    }
  } else {
    # do not project data
    projected <- x
  }
  
  bounds <- raster::extent(
    raster::projectExtent(
      raster::projectExtent(x, crs = sp::CRS(epsg3857)),
      crs = sp::CRS(epsg4326)
    )
  )
  
  if (!is.function(colors)) {
    if (method == "ngb") {
      # 'factors'
      colors <- colorFactor(colors, domain = NULL, na.color = "#00000000", alpha = TRUE)
    } else {
      # 'numeric'
      colors <- colorNumeric(colors, domain = NULL, na.color = "#00000000", alpha = TRUE)
    }
  }
  
  tileData <- raster::values(projected) %>% colors() %>% col2rgb(alpha = TRUE) %>% as.raw()
  dim(tileData) <- c(4, ncol(projected), nrow(projected))
  pngData <- png::writePNG(tileData)
  if (length(pngData) > maxBytes) {
    stop(
      "Raster image too large; ", length(pngData), " bytes is greater than maximum ",
      maxBytes, " bytes"
    )
  }
  encoded <- base64enc::base64encode(pngData)
  uri <- paste0("data:image/png;base64,", encoded)
  
  latlng <- list(
    list(raster::ymax(bounds), raster::xmin(bounds)),
    list(raster::ymin(bounds), raster::xmax(bounds))
  )
  
  invokeMethod(map, data, "addRasterImage", uri, latlng, opacity, attribution, layerId, group, options) %>%
    expandLimits(
      c(raster::ymin(bounds), raster::ymax(bounds)),
      c(raster::xmin(bounds), raster::xmax(bounds))
    )
}



gridOptions <- function(
  tileSize = 256,
  updateWhenIdle = NULL,
  zIndex = 1,
  minZoom = 0,
  maxZoom = NULL,
  ...
) {
  filterNULL(list(
    tileSize = tileSize, updateWhenIdle = updateWhenIdle, zIndex = zIndex,
    minZoom = minZoom, maxZoom = maxZoom,
    ...
  ))
}