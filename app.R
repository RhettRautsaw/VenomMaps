library(shiny)
library(leaflet)
library(leaflet.extras)
library(leaflet.extras2)
library(RColorBrewer)
library(sf)
library(raster)
library(scales)
library(lattice)
library(tidyr)
library(dplyr)
library(readr)
library(ggpubr)
library(ggtext)
library(patchwork)
library(sp)
# library(rgdal)
# library(htmltools)
# library(htmlwidgets)

# Tabbed page example
# https://github.com/rstudio/shiny-examples/tree/master/063-superzip-example

infoGenerator<-function(info, full_info){
    final=list()
    for(i in 1:nrow(info)){
        tmp<-paste0(
            '<br><span style = "font-size:20pt;"><b><i>',gsub("_"," ",info[i,1]),' </b></i></span><span style="font-size:10pt;">',info[i,4],'</span><br>',
            '<span style = "font-size:15pt;">', info[i,8],'</span><br>',
            '<span style = "font-size:15pt;"><strong>Family:</strong> Viperidae (',info[i,7],')</span><br>',
            '<span style = "font-size:15pt;"><strong>Subspecies:</strong> ',info[i,3],'</span><br>',
            '<span style = "font-size:15pt;"><strong>',info[i,10],':</strong> ',info[i,9],' mm</span><br>',
            '<span style = "font-size:15pt;"><strong>Distribution</strong>: ', info[i,16],'</span><br><br>'
        )
        tmp2<-gghistogram(full_info, x="max_length_mm", fill="gray", xlab="Max Length (mm)") + 
            geom_vline(xintercept=as.numeric(info[i,9]), color="red", linetype="dashed") +
            annotate(geom='text', label=as.character(info[i,9]), x=as.numeric(info[i,9]), y=Inf, hjust=-0.5, vjust=5, size=8)+
            labs(title=tmp)+
            theme(plot.title = element_textbox_simple(lineheight = 1.5))
        final[[i]]<-tmp2
    }
    final
}


load("data/shiny-data2.RData")
source("support_scripts/addRasterImage2.R")

# info<-read_csv("data/ViperInfo.csv")
# source("support_scripts/stack_diff_extents.R")
# crs(combined_niches)<-CRS("+init=epsg:4326")
# #brPal <- colorRampPalette(c('#008B00FF', '#008B0000', '#008B0000'), alpha=T)
# brPal <- colorRampPalette(c('#00A600FF', '#61C500BF', '#E6E40280', '#ECB17640', '#F2F2F200'), alpha=T)
# pal <- brPal(255)

species_list<-unique(sort(as.vector(info$species)))
countries_list<-sort(unique(unlist(strsplit(paste(unlist(info$countries), collapse="; "), "; "))))

# UI
ui <- navbarPage("VenomMaps", id="nav",
    
    ###################################
    ######### Interactive Map #########
    ###################################
    
    tabPanel("Interactive map",
        div(class="outer",
            # Include custom styles and javascript
            tags$head(includeCSS("styles.css"), includeScript("gomap.js")),
            # If not using custom CSS, set height of leafletOutput to a number instead of percent
            leafletOutput("map", width="100%", height="100%"),
            
            absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                draggable = TRUE, top = 60, left = "auto", right = 20, bottom = 50,
                width = 330, style = "overflow-y: scroll;",
                img(src = "VenomMaps.png", height = 300, width = 250, style="display: block; margin-left: auto; margin-right: auto"),
                
                p("Use this \n panel to see lists of species by country. Then select your species of interest and explore the different tabs!"),
                selectizeInput(inputId = "countries", label = h4("Country:"), choices = countries_list,
                               multiple = TRUE),
                selectizeInput(inputId = "species", label = h4("Species:"), choices = species_list, 
                               multiple = TRUE),
                
                h5("____________________________________"),
                
                checkboxInput("niche", "Niche Model", FALSE),
                
                checkboxInput("nodist", "Clear Distribution", FALSE),
                
                h5("____________________________________"),
                
                checkboxInput("points", "Occurrence Points", FALSE), # and Heatmap
                p("Gray points are those beyond their distribution"),
                p("Brown points are those with altered species IDs"),
                
                h5("____________________________________"),
                
                actionButton("update", "Update"),
                downloadButton("downloadDist", "Download Shapefile"),
                
                h5("____________________________________"),
                
                h5("Acknowledgements:"),
                p("Distributions of Old World taxa were obtained from ", a(href="https://www.nature.com/articles/s41559-017-0332-2", "Roll et al. 2017")),
                p("Information on max snake lengths was obtained from ", a(href="https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.12398", "Feldman et al. 2015")),
                p("Information on descriptive citations, common names, and subspecies was obtained from the ", a(href="https://reptile-database.reptarium.cz/", "Reptile Database"))
            )
        )
    ),
    
    ###################################
    ####### General Information #######
    ###################################
    
    tabPanel("General Information",
        plotOutput("infoPlot", height = "1000px")
    )

    ###################################
    ######### Phylogeography ##########
    ###################################
    
    # tabPanel("Phylogeography")
    
)


# Server
server<-function(input, output, session) {
    
    ###################################
    ######### Interactive Map #########
    ###################################
    
    output$map<-renderLeaflet({
        leaflet() %>% 
            addMapPane("background", zIndex = 0) %>%        # Level 1: bottom
            addMapPane("polygons", zIndex = 1) %>%          # Level 2: middle
            addMapPane("rasters", zIndex = 100000) %>%      # Level 3: middle
            addMapPane("points", zIndex = 440) %>%          # Level 4: middle
            addMapPane("labels", zIndex = 450) %>%          # Level 5: top
            addWMSTiles('http://ows.mundialis.de/services/service?', layers='TOPO-WMS', group="Topography", options = pathOptions(pane = "background")) %>%
            addProviderTiles(providers$Esri.WorldStreetMap, group="Open Street Map", options = pathOptions(pane = "background")) %>%
            addProviderTiles(providers$Esri.WorldTerrain, group="Terrain", options = pathOptions(pane = "background")) %>%
            addProviderTiles(providers$Esri.WorldImagery, group="Satellite", options = pathOptions(pane = "background")) %>%
            addProviderTiles(providers$Stamen.TonerLines, group="Boundaries", options = pathOptions(pane = "labels")) %>%
            addProviderTiles(providers$Stamen.TonerLabels, group="Labels", options = pathOptions(pane = "labels")) %>%
            addLayersControl(
                baseGroups=c("Topography", "Open Street Map", "Terrain", "Satellite"),
                overlayGroups=c("Boundaries","Labels"),
                options=layersControlOptions(collapsed=F, position="bottomleft")
            ) %>%
            hideGroup("Labels") %>%
            addScaleBar(position = c("topleft"), options = scaleBarOptions()) %>% 
            addEasyprint(options=easyprintOptions(exportOnly = T, sizeModes = list("A4Landscape")))
    })
    
    # Update species list by country selection
    new_species_list<-reactive({
        if(is.null(input$countries)){
            info %>% select(species) %>% arrange(species) %>% distinct() %>% pull(species)
        }else{
            info %>% filter(grepl(paste(input$countries,collapse="|"), countries)) %>%
                     select(species) %>% arrange(species) %>% distinct() %>% pull(species)
        }
    })
    
    observe({
        updateSelectizeInput(session,"species", choices = new_species_list(),
                             selected=head(new_species_list(),1))
    })
    
    # Filter distributions
    distribution<-reactive({
        if(is.null(input$species)){
            as(combined_distribution,"Spatial")
        }else{
            as(combined_distribution %>% filter(Species %in% input$species),"Spatial")
        }
    })
    
    # Filter enms
    niche<-reactive({
        if((input$niche) & any(input$species %in% names(combined_niches))){
            #raster::subset(combined_niches, input$species)
            tmp<-combined_niches[[which(names(combined_niches) %in% input$species)]]
            max(tmp, na.rm = TRUE)
        }else{
            vector(mode="numeric", length=0)
        }
    })
    
    # Download geojson button
    output$downloadDist <- downloadHandler(
        filename = function() {
            paste0("distribution.geojson")
        },
        content = function(con) {
            st_write(as(distribution(),"sf"), dsn = con, layer=con, driver="GeoJSON")
        }
    )
    
    # Update map with distribution/points
    observeEvent(input$update, {
        distribution<-distribution()
        bbox<-st_bbox(as(distribution,"sf")) %>% as.vector()
        sp_factpal<-colorFactor(sp_pal,levels=sort(distribution$Subspecies), ordered=T, reverse=T)
        
        leafletProxy("map", data = distribution) %>%
            clearGroup("distribution") %>% clearGroup("occpoints") %>% clearHeatmap() %>% clearMarkerClusters() %>%  removeControl("distribution") %>% clearImages() %>%
            addPolygons(data=distribution, color="black", weight=3, fillColor = ~sp_factpal(Subspecies), fillOpacity = 0.5, group="distribution", options = pathOptions(pane = "polygons")) %>%
            addLegend(data=distribution, position = "topleft", pal=sp_factpal, values = ~Subspecies, layerId = "distribution") %>%
            fitBounds(bbox[1],bbox[2],bbox[3],bbox[4])
        
        if(input$points){
            # Filter points
            pointsData<-occ %>% filter(Species==input$species) %>%
                mutate(labs=paste0( '<p>', prov, " ", 
                                    ID, '<p></p>',
                                    "recorded: ", new_species,'</p><p>',
                                    "updated: ", Species, '</p><p>',
                                    "accuracy: ", accuracy,' m</p>'))
            labs<-pointsData$labs
            if(nrow(pointsData)!=0){
                leafletProxy("map", data = distribution) %>%
                    addCircleMarkers(data=pointsData,group="occpoints",
                                    lng=~as.numeric(longitude),
                                    lat=~as.numeric(latitude),
                                    radius = 4,
                                    color = ~point_pal(flag),
                                    stroke = FALSE, fillOpacity = 0.8,
                                    options = pathOptions(pane = "points"),
                                    #clusterOptions = markerClusterOptions(),
                                    label = lapply(labs, htmltools::HTML)) #%>%
                    # addHeatmap(data=pointsData,
                    #             lng=~as.numeric(longitude),
                    #             lat=~as.numeric(latitude),
                    #             radius = 30,
                    #             blur=60)
            }
        }
        if(input$niche){
            if(length(niche())!=0){
                if(input$nodist){
                    leafletProxy("map") %>%
                        clearGroup("distribution") %>% removeControl("distribution") %>%
                        addRasterImage2(niche(), colors = rev(pal), options = tileOptions(pane = "rasters"))
                }else{
                    leafletProxy("map") %>%
                        addRasterImage2(niche(), colors = rev(pal), options = tileOptions(pane = "rasters"))   
                }
            }
        }
    })
    
    ###################################
    ####### General Information #######
    ###################################
    
    output$infoPlot <- renderPlot({
        if(is.null(input$species)){
            tmpInfo<-info[1,]
        }else{
            tmpInfo<-info %>% filter(grepl(paste(input$species,collapse="|"), species))
        }
        
        tmp<-infoGenerator(tmpInfo, info)
        
        wrap_plots(tmp, ncol=ifelse(length(tmp)<3,length(tmp),3))
    })

    ###################################
    ######### Phylogeography ##########
    ###################################
    
}

# Run the application 
shinyApp(ui = ui, server = server)
