library(shiny)
library(leaflet)
library(leaflet.extras)
library(leaflet.extras2)
library(RColorBrewer)
library(sf)
library(raster)
library(tidyr)
library(dplyr)
library(readr)
library(ggpubr)
library(ggtext)
library(patchwork)
library(sp)
library(shinybusy)
library(readxl)
# library(scales)
# library(lattice)
# library(rgdal)
# library(htmltools)
# library(htmlwidgets)

# Tabbed page example
# https://github.com/rstudio/shiny-examples/tree/master/063-superzip-example

infoGenerator<-function(info, full_info){
    final=list()
    for(i in 1:nrow(info)){
        tmp<-paste0(
            '<br><span style = "font-size:20pt;"><b><i>',gsub("_"," ",info[i,"species"]),' </b></i></span><span style="font-size:10pt;">',info[i,"description_citation"],'</span><br>',
            '<span style = "font-size:15pt;">', info[i,"common_name"],'</span><br>',
            '<span style = "font-size:15pt;"><strong>Family:</strong> Viperidae (',info[i,"subfamily"],')</span><br>',
            '<span style = "font-size:15pt;"><strong>Subspecies:</strong> ',info[i,"subspecies"],'</span><br>',
            '<span style = "font-size:15pt;"><strong>',info[i,"length_measure"],':</strong> ',info[i,"max_length_mm"],' mm</span><br>',
            '<span style = "font-size:15pt;"><strong>Distribution</strong>: ', info[i,"countries"],'</span><br><br>'
        )
        tmp2<-gghistogram(full_info, x="max_length_mm", fill="gray", xlab="Max Length (mm)") + 
            geom_vline(xintercept=as.numeric(info[i,"max_length_mm"]), color="red", linetype="dashed") +
            annotate(geom='text', label=as.character(info[i,"max_length_mm"]), x=as.numeric(info[i,"max_length_mm"]), y=Inf, hjust=-0.5, vjust=5, size=8)+
            labs(title=tmp)+
            theme(plot.title = element_textbox_simple(lineheight = 1.5))
        final[[i]]<-tmp2
    }
    final
}

# Load Data
sp_pal<-RColorBrewer::brewer.pal(11, "BrBG")
point_pal <- colorFactor(c("red2", "gray75", "darkturquoise", "khaki2"), domain = c("dubious", "geoDist", "noLand", "updated"), na.color = "black")

load("shiny_support_material/shiny-data_2021-10-25.RData")
source("shiny_support_material/addRasterImage2.R")

#info<-read_csv("shiny_support_material/ViperInfo.csv")
info<-read_xlsx("supplemental_material/SupplementalTable1.xlsx", sheet=1)
source("shiny_support_material/stack_diff_extents.R")

#brPal <- colorRampPalette(c('#008B00FF', '#008B0000', '#008B0000'), alpha=T)
brPal <- colorRampPalette(c('#00A600FF', '#61C500BF', '#E6E40280', '#ECB17640', '#F2F2F200'), alpha=T)
pal <- brPal(255)

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

                tags$div(HTML("<br><center><p><a href=\"https://github.com/RhettRautsaw/VenomMaps\"><img src=\"https://img.shields.io/badge/GitHub-RhettRautsaw/VenomMaps-blue\"></a></p></center>")),
                tags$div(HTML("<center><p><a href=\"\"><img src=\"https://img.shields.io/badge/Published%20in-Scientific%20Data-blue\"></a></p></center>")),
                tags$div(HTML("<center><p><a href=\"https://creativecommons.org/licenses/by-nc-sa/4.0/\"><img src=\"https://img.shields.io/badge/License-CC%20BY-blue\"></a></p></center>")),
                tags$div(HTML("<center><p><a href=\"https://zenodo.org/badge/latestdoi/366181795\"><img src=\"https://zenodo.org/badge/366181795.svg\" alt=\"DOI\"></a></p></center>")),
                
                tags$div(HTML("<center><p>Select a species and hit \"Update\"</p><p>You can also filter species by country or add occurrence records/niche models with the check boxes.</p></center>")),
                selectizeInput(inputId = "countries", label = h4("Country:"), choices = countries_list,
                               multiple = TRUE),
                selectizeInput(inputId = "species", label = h4("Species:"), choices = species_list, 
                               multiple = TRUE),
                
                actionButton("update", "Update"),
                
                h5("____________________________________"),
                
                checkboxInput("niche", "Niche Model", FALSE),
                
                checkboxInput("nodist", "Clear Distribution", FALSE),
                
                p("Select \"Update\" again after checking this box"),
                p("Please be patient. These large files can take a while to plot."),
                
                h5("____________________________________"),
                
                checkboxInput("points", "Occurrence Points", FALSE), # and Heatmap
                p("Select \"Update\" again after checking this box"),
                p("Red: points are dubious/questionable"),
                p("Gray: points are beyond their distribution"),
                p("Blue: points are not on land"),
                p("Yellow: species ID has been altered"),
                
                h5("____________________________________"),

                downloadButton("downloadDist", "Download Shapefile"),
                
                h5("____________________________________"),
                
                h5("Acknowledgements:"),
                p("Distributions of Old World Viperidae were obtained from ", a(href="https://www.nature.com/articles/s41559-017-0332-2", "GARD 1.1")),
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
    
    ###################################
    ############# Venom ###############
    ###################################
    
    # tabPanel("Venom")
    
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
    # niche<-reactive({
    #     if((input$niche) & any(input$species %in% names(combined_niches))){
    #         #raster::subset(combined_niches, input$species)
    #         tmp<-combined_niches[[which(names(combined_niches) %in% input$species)]]
    #         max(tmp, na.rm = TRUE)
    #     }else{
    #         vector(mode="numeric", length=0)
    #     }
    # })
    
    # Load ENMs
    niche<-reactive({
        if((input$niche) & length(list.files("data/enms_p", paste0(input$species,"_avg.tif", collapse = "|"), full.names = T))>1){
            files<-list.files("data/enms", paste0(input$species,"_avg.tif", collapse = "|"), full.names = T)
            combined_niches<-stack_diff_extents(files)
            names(combined_niches)<-gsub("_avg","",names(combined_niches))
            # crs(combined_niches)<-CRS("+init=epsg:4326")
            max(combined_niches, na.rm = TRUE)
        }else if((input$niche) & length(list.files("data/enms", paste0(input$species,"_avg.tif", collapse = "|"), full.names = T))==1){
            files<-list.files("data/enms", paste0(input$species,"_avg.tif", collapse = "|"), full.names = T)
            combined_niches<-raster(files)
            names(combined_niches)<-gsub("_avg","",names(combined_niches))
            combined_niches
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
        show_modal_spinner()
        distribution<-distribution()
        bbox<-st_bbox(as(distribution,"sf")) %>% as.vector()
        sp_factpal<-colorFactor(sp_pal,levels=sort(distribution$Subspecies), ordered=T, reverse=T)
        
        leafletProxy("map", data = distribution) %>%
            clearGroup("distribution") %>% clearGroup("occpoints") %>% clearHeatmap() %>% clearMarkerClusters() %>%  removeControl("distribution") %>% clearImages() %>%
            addPolygons(data=distribution, color="black", weight=3, fillColor = ~sp_factpal(Subspecies), fillOpacity = 0.75, group="distribution", options = pathOptions(pane = "polygons")) %>%
            addLegend(data=distribution, position = "topleft", pal=sp_factpal, values = ~Subspecies, layerId = "distribution") %>%
            fitBounds(bbox[1],bbox[2],bbox[3],bbox[4])
        
        if(input$points){
            # Filter points
            pointsData<-occ %>% filter(final_species %in% input$species) %>%
                mutate(labs=paste0( ID, '<p></p>',
                                    "recorded: ", taxonomy_updated_species,'</p><p>',
                                    "updated: ", final_species, '</p><p>',
                                    "accuracy: ", accuracy_m,' m</p><p>',
                                    "flags: ", flag_detailed, '</p>'))
            labs<-pointsData$labs
            if(nrow(pointsData)!=0){
                leafletProxy("map", data = distribution) %>%
                    addCircleMarkers(data=pointsData,group="occpoints",
                                    lng=~as.numeric(longitude),
                                    lat=~as.numeric(latitude),
                                    radius = 4,
                                    fill = TRUE, fillColor = ~point_pal(flag), fillOpacity = 1, weight=1,
                                    stroke = TRUE, color="black", opacity=1,
                                    options = pathOptions(pane = "points"),
                                    #clusterOptions = markerClusterOptions(),
                                    popup = lapply(labs, htmltools::HTML)) #%>%
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
                        addRasterImage2(niche(), colors = rev(pal), options = tileOptions(pane = "rasters"), 
                                        maxBytes = 800*1024*1024, project = T)
                }else{
                    leafletProxy("map") %>%
                        addRasterImage2(niche(), colors = rev(pal), options = tileOptions(pane = "rasters"), 
                                        maxBytes = 800*1024*1024, project = T)   
                }
            }
        }
        remove_modal_spinner()
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
