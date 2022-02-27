#!/usr/bin/env Rscript
# Copyright 2020 Rhett M. Rautsaw
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#=====================#
##### DESCRIPTION #####
#=====================#
#
# autokuenm.R is designed to take occurrence records, distribution maps, and environmental data and prepare data for kuenm.
# 
# Specifically, this script will:
# 1. Subset occurrence records by a species name
# 2. Subset distribution maps by a species name (or create hull from points)
# 3. Overlay the distribution/hull with bioregions to define M areas
# 4. Partition occurrence records for training and testing (ENMeval)
# 5. Create exhaustive combination of all environmental variables or perform heuristic search on sets
#     - Heuristic search involves identifying non-colinear variables with high permutation importance
# 6. Run kuenm
# 7. Summarize results
#
#
#=============#
##### USE #####
#=============#
#
# autokuenm.R -o occ.csv -e env.txt -d ahull -sp Agkistrodon_bilineatus \
#             -spcol1 Species -spcol2 Species -latcol latitude -loncol longitude -b bias.tif -bg 2000 \
#             -bias spatial_bias.tif -g G_variables -c 8 -m 62000
#
#
#======================#
##### REQUIREMENTS #####
#======================#

packages <- c("argparse","readr","dplyr","stringr","parallel","memuse",
              "sf","raster","dismo","rJava","virtualspecies", "concaveman", "spThin",
              "rnaturalearth", "rnaturalearthdata", "usdm","ENMeval","kuenm")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  print("Installing required packages")
  install.packages(packages[!installed_packages], repos="https://cloud.r-project.org")
}


#=========================#
##### SETUP ARGUMENTS #####
#=========================#

suppressPackageStartupMessages(library("argparse"))
mem<-as.numeric(gsub(" GiB","",memuse::Sys.meminfo()[["totalram"]]))*1000
cpu<-parallel::detectCores()

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option
parser$add_argument("-o", default="combined_records_v4_clean.csv", 
                    help="csv of occurrence records. Can contain multiple species and will be filtered by -sp and -spcol flags. [default: \"%(default)s\"]")
parser$add_argument("-e", default="env2.txt", 
                    help="For exhaustive search, provide a folder containing all environmental rasters. For heuristic search, provide a txt file containing paths to all environmental rasters divided into Sets. To skip variable selection and use a pre-defined set of environmental variables, follow heuristic search input (i.e. txt file) and turn on -novarsel. Rasters must be in .tif or .asc format [default: \"%(default)s\"]")
parser$add_argument("-d", default="ahull", 
                    help="Distribution shapefile (e.g., GeoJSON or SHP).  If no shapefile is provided, then you can use a minimum convex polygon or concave alpha hull by specifiying 'mcp' or 'ahull'. Distribution file can contain multiple species and will be filtered by -sp and -spcol flags. [default: \"%(default)s\"]")
parser$add_argument("-sp", default="Agkistrodon_bilineatus", 
                    help="Name of species to model. [default: \"%(default)s\"]")
parser$add_argument("-spcol1", default="final_species", 
                    help="Name of column with species ID in occurrence records csv. [default: \"%(default)s\"]")
parser$add_argument("-spcol2", default="Species", 
                    help="Name of column with species ID in distribution attribute table. [default: \"%(default)s\"]")
parser$add_argument("-latcol", default="latitude", 
                    help="Name of column with latitude in occurrence records csv [default: \"%(default)s\"]")
parser$add_argument("-loncol", default="longitude", 
                    help="Name of column with longitude in occurrence records csv [default: \"%(default)s\"]")
parser$add_argument("-b", default="bioregions/wwf_terr_ecos.shp", 
                    help="Bioregions used to assign M areas [default %(default)s]")
parser$add_argument("-bg", type="integer", default=2000, metavar="number",
                    help="Number of random background points to generate for assigning training/testing partitions (Not used in MaxEnt) [default %(default)s]")
parser$add_argument("-bias", default="SpatialBias_resamp.tif",
                    help="Use a spatial bias file to factor out bias within MaxEnt. If empty, a distance filter is applied to occurrence records.")
parser$add_argument("-novarsel", action="store_true", default=FALSE,
                    help="Turn off variable selection if you want to use all environmental variables in -e list [default: \"%(default)s\"]")
parser$add_argument("-g", default="", 
                    help="Directory containing all environmental rasters for projection in kuenm. Do not include if you do not want to project to larger or forward/backward in time. [default: \"%(default)s\"]")
parser$add_argument("-c", type="integer", default=cpu, metavar="number",
                    help="Number of processor threads to use. Matching this number to the number of cores on your computer speeds up some operations. [default: \"%(default)s\"]")
parser$add_argument("-m", type="integer", default=mem, metavar="number",
                    help="Maximum memory (in megabytes) to be used by maxent while creating the models. [default: \"%(default)s\"]")


#========================#
##### READ ARGUMENTS #####
#========================#

args <- parser$parse_args()
occ <- normalizePath(args$o)
env <- normalizePath(args$e)
bioregions<-normalizePath(args$b)

if(file.exists(env) && !dir.exists(env)){
  search<-"heuristic"
  suppressMessages(env.tab<-readr::read_tsv(env))
  env.list<-normalizePath(env.tab$Var)
}else{
  search<-"exhaustive"
  env.list<-list.files(env, ".tif|.asc", full.names = T)
}

if(args$d=="mcp" | args$d=="ahull"){
  dist<-args$d
}else{
  dist<-normalizePath(args$d) 
}

if(file.exists(args$bias)){
  bias<-normalizePath(args$bias)
  bias_logic<-TRUE
}else{
  bias<-"distance thinning"
  bias_logic<-FALSE
}

if(dir.exists(args$g)){
  gvars<-normalizePath(args$g)
  gvars_logic<-TRUE
}else{
  gvars<-"NA"
  gvars_logic<-FALSE
}

cat("Starting autokuenm...\n")
  cat(paste0("\t Species of Interest:\t", args$sp,"\n"))
  cat(paste0("\t Occurrence Records:\t", occ,"\n"))
  cat(paste0("\t Environmental Data:\t", env,"\n"))
    cat(paste0("\t\t Search Type:\t", search,"\n"))
  cat(paste0("\t Distribution Map:\t", dist,"\n"))
  cat(paste0("\t Bioregions (M areas):\t", bioregions, "\n"))
  cat(paste0("\t Background Points:\t", args$bg,"\n"))
  cat(paste0("\t Spatial Bias:\t\t", bias,"\n"))
  cat(paste0("\t Project:\t\t", gvars_logic,"\n"))
    cat(paste0("\t\t Variables:\t", gvars,"\n"))
  cat(paste0("\t Number of Cores:\t", args$c,"\n"))
  cat(paste0("\t Total Memory:\t\t", args$m, " MB\n"))

#=======================#
##### LOAD PACKAGES #####
#=======================#

cat(paste0("\n", Sys.time(), " ::: Loading packages and preparing directory :::\n"))
suppressWarnings(invisible(suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))))


#=================================#
##### SETUP FUNCTIONS & DISMO #####
#=================================#

dismo_jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

# checking if maxent can be run (normally not part of your script)
if (!file.exists(dismo_jar)){
  invisible(file.copy(from="maxent.jar", to=dismo_jar))
}

source("functions.R")

#================================#
##### PREP WORKING DIRECTORY #####
#================================#

if(!dir.exists(args$sp)){
  dir.create(args$sp)
  invisible(file.copy(from="maxent.jar", to=paste0(args$sp,"/maxent.jar")))
}

setwd(args$sp)


#================================#
##### RECORD/ENV DATA PREP #######
#================================#

if(!file.exists("current_step.log")){
  
  #==========================#
  ##### PREP OCC RECORDS #####
  #==========================#
  
  cat(paste0("\n", Sys.time(), " ::: Reading ",args$o, " :::\n"))
  suppressMessages(occ <- read_csv(occ))
  cat(paste0("\n", Sys.time(), " ::: Filtering ",nrow(occ), " records from ", args$o, " for ", args$sp," :::\n"))
  occ2<-occ %>% filter((!!as.name(args$spcol1))==args$sp)
  cat(paste0("\n", Sys.time(), " ::: ", nrow(occ2), " records remaining :::\n"))
  
  #=====================#
  ##### PREP M AREA #####
  #=====================#
  
  if(dist=="mcp"){
    cat(paste0("\n", Sys.time(), " ::: Creating Minimum Convex Polygon :::\n"))
    dist <- st_as_sf(occ2, coords = c(args$loncol,args$latcol), crs=4326) %>%
      st_union() %>%
      st_convex_hull()
  }else if(dist=="ahull"){
    cat(paste0("\n", Sys.time(), " ::: Creating Concave/Alpha Hull :::\n"))
    occ2_sf <- st_as_sf(occ2, coords = c(args$loncol,args$latcol), crs=4326)
    dist <- ahull(occ2_sf)
    st_write(dist,"ahull.geojson", quiet = T)
  }else{
    cat(paste0("\n", Sys.time(), " ::: Reading ",args$d, " :::\n"))
    dist<-sf::st_read(dist,crs=4326, quiet=T)
    cat(paste0("\n", Sys.time(), " ::: Filtering ",args$d, " for ", args$sp," :::\n"))
    dist<-dist %>% filter((!!as.name(args$spcol2))==args$sp) 
  }
  
  cat(paste0("\n", Sys.time(), " ::: Creating M Area from overlapping biogeographic regions :::\n"))
  bioregions<-sf::st_read(bioregions, crs=4326, quiet=T)
  dist2<-st_sf(Marea(bioregions, dist, cov=20))
  dist2<-fill_holes(dist2, Inf)
  st_write(dist2, "M_area.geojson", quiet=T)
  png("M_area.png")
  plot(dist2)
  dev.off()
  
  #=======================#
  ##### PREP ENV DATA #####
  #=======================#
  
  cat(paste0("\n", Sys.time(), " ::: Reading and stacking ",args$e, " :::\n"))
  envStack <- raster::stack(env.list)
  cat(paste0("\n", Sys.time(), " ::: Cropping environmental stack to M area extent for faster masking :::\n"))
  envStack2<-crop(envStack,extent(dist2))
  cat(paste0("\n", Sys.time(), " ::: Masking environmental stack to M area :::\n"))
  envStack2<-mask(envStack2,dist2) 
  
  #=======================#
  ##### SAMPLING BIAS #####
  #=======================#
  
  cat(paste0("\n", Sys.time(), " ::: Removing records not overlapping with environmental data :::\n"))
  tmp<-raster::extract(envStack2,occ2[,c(args$loncol,args$latcol)])
  occ2<-occ2[complete.cases(tmp),]
  occ2<-occ2[,c(args$spcol1,args$loncol,args$latcol)]
  colnames(occ2)<-c("Species","Longitude","Latitude")
  cat(paste0("\n", Sys.time(), " ::: Removing non-unique records :::\n"))
  occ2<-unique(occ2)
  cat(paste0("\n", Sys.time(), " ::: ", nrow(occ2), " records remaining :::\n"))
  
  if(bias_logic){
    bias<-raster(bias)
    crs(bias)<-crs(envStack2)
    cat(paste0("\n", Sys.time(), " ::: Cropping bias file to M area :::\n"))
    bias2<-crop(bias,extent(dist2))
    cat(paste0("\n", Sys.time(), " ::: Masking bias file to M area :::\n"))
    bias2<-mask(bias2,dist2)
    bias3<-projectRaster(bias2,envStack2)
    png("bias.png")
    plot(bias3)
    dev.off()
    writeRaster(bias3,"bias.asc")
  }else{
    cat(paste0("\n", Sys.time(), " ::: Thinning occ records to remove/reduce sampling bias :::\n"))
    ## Spatial Thin - use linear model to choose distance for thinning based on number of records
    d<-floor(as.numeric(predict(sp_mod, data.frame(V1=nrow(occ2)))))
    cat(paste0("\n", Sys.time(), " ::: Thinning distance: ", d, " km :::\n"))
    d2=spThin::thin(occ2, lat.col="Latitude", long.col ="Longitude", spec.col = "Species", thin.par = d, write.files = F, reps=10, locs.thinned.list.return=T)
    occ2<-cbind.data.frame(Species=args$sp, d2[[which(unlist(lapply(d2,nrow))==max(unlist(lapply(d2, nrow))))[1]]])
    cat(paste0("\n", Sys.time(), " ::: ", nrow(occ2), " records remaining :::\n"))
  }
  
  #==============================#
  ##### INDEPENDENT TEST SET #####
  #==============================#
  
  sp_ind_n<-as.integer(0.05*nrow(occ2))
  if(sp_ind_n >= 5){
    cat(paste0("\n", Sys.time(), " ::: Selecting random 5% of records (n=",sp_ind_n,") as independent dataset :::\n"))
    set.seed(12345)
    sp_ind_rows<-sample(1:nrow(occ2), sp_ind_n)
    sp_ind<-occ2[sp_ind_rows,]
    occ3<-occ2[-sp_ind_rows,]
    write.csv(sp_ind,"Sp_ind.csv", row.names=F, quote=F)
    cat(paste0("\n", Sys.time(), " ::: ", nrow(occ3), " records remaining :::\n")) 
  }else{
    warn<-paste0("\n", Sys.time(), " ::: WARNING ::: Too few records to create independent dataset. Partitioning and using training data for final model creation rather than joint data.\n")
    write(warn,"warnings.txt",append=TRUE)
    cat(warn)
    occ3<-occ2
  }

  #=========================================#
  ##### PARTITION TESTING/TRAINING SETS #####
  #=========================================#
  
  cat(paste0("\n", Sys.time(), " ::: Generating ",args$bg," random background points across cropped environmental grid :::\n"))
  bg <- randomPoints(envStack2[[1]], n = args$bg)
  bg <- as.data.frame(bg)
  
  cat(paste0("\n", Sys.time(), " ::: Partitioning occ records for training and testing :::\n"))
  #block <- get.block(occ3[, c("Longitude", "Latitude")], bg)
  block<-get.checkerboard2(occ = occ3[, c("Longitude", "Latitude")], env = envStack2, bg = bg, aggregation.factor = c(2,2))
  bg<-cbind("bg",bg)
  names(bg)<-c("Species","Longitude","Latitude")
  
  png("blocks_records.png", width = 1080, height=720)
  par(mfrow=c(1,2))
  plot(envStack2[[1]], legend = FALSE, main="occurrence records")
  points(occ3[, c("Longitude", "Latitude")], pch = 21, bg = block[[1]])
  plot(envStack2[[1]], legend = FALSE, main="background points")
  points(bg[, c("Longitude", "Latitude")], pch = 21, bg = block[[2]])
  invisible(dev.off())

  cat(paste0("\n", Sys.time(), " ::: Selecting one partition as test dataset and the rest as training data :::\n"))
  sp_train<-occ3[which(block[[1]]!=2),]
  sp_test<-occ3[which(block[[1]]==2),]
  write.csv(sp_test,"Sp_test.csv", row.names=F, quote=F)
  write.csv(sp_train,"Sp_train.csv", row.names=F, quote=F)
  write.csv(occ3,"Sp_joint.csv", row.names=F, quote=F)
  cat(paste0("\n", Sys.time(), " ::: Using ", nrow(sp_train), " records for training and ", nrow(sp_test), " records for testing :::\n"))
  
  #=================================#
  ##### ENV VARIABLE SELECTION  #####
  #=================================#
  
  if(search=="exhaustive"){
    cat(paste0("\n", Sys.time(), " ::: Writing All Cropped Environmental Variables :::\n"))
    dir.create(paste0("crop_envs/"), recursive = T)
    writeRaster(envStack2, paste0("crop_envs/",names(envStack2),".asc"), bylayer=T, format="GTiff")
    write(names(envStack2),paste0("Variables_all.txt"))
    
    cat(paste0("\n", Sys.time(), " ::: Creating Exhaustive Set of Environmental Variables :::\n"))
    # format<-tools::file_ext(env.list[1])
    # if(format=="tif"){format="GTiff"}else{format="ascii"}
    combinations <- kuenm_varcomb(var.dir = "crop_envs/", out.dir = "M_variables",
                                  min.number = 10, in.format = "GTiff",
                                  out.format = "ascii")
  }
  
  if(search=="heuristic"){
    if (args$novarsel){
      cat(paste0("\n", Sys.time(), " ::: Writing All Cropped Environmental Variables :::\n"))
      dir.create(paste0("M_variables/Set_all"), recursive = T)
      writeRaster(envStack2, paste0("M_variables/Set_all/",names(envStack2),".asc"), bylayer=T, format="ascii")
      write(names(envStack2),paste0("Variables_all.txt"))
    } else{
      cat(paste0("\n", Sys.time(), " ::: Creating Heuristic Set of Environmental Variables :::\n"))
      tmpocc <- cbind.data.frame(Species=args$sp, as.data.frame(gridSample(xy = as.matrix(occ3[,2:3]), r = envStack2[[1]], n = 1)))
      for(i in unique(env.tab$Set)){
        cat(paste0("\n", Sys.time(), " ::: Creating Set ",i," :::\n"))
        set<-which(env.tab$Set==i)
        tmpStack<-envStack2[[set]]
        newStack<-SelectVariables(occ = tmpocc[,2:3], stack = tmpStack)
        dir.create(paste0("M_variables/Set_",i), recursive = T)
        writeRaster(newStack, paste0("M_variables/Set_",i,"/",names(newStack),".asc"), bylayer=T, format="ascii")
        write(names(newStack),paste0("Variables_",i,".txt"))
      }
      cat(paste0("\n", Sys.time(), " ::: Creating Set topvars :::\n"))
      newStack<-SelectVariables(occ = tmpocc[,2:3], stack = envStack2)
      dir.create(paste0("M_variables/Set_topvars"), recursive = T)
      writeRaster(newStack, paste0("M_variables/Set_topvars/",names(newStack),".asc"), bylayer=T, format="ascii")
      write(names(newStack),paste0("Variables_topvars.txt"))
    }
  }
  writeLines("1", "current_step.log")
}
  


#=====================#
##### KUENM SETUP #####
#=====================#
sets=c(unique(env.tab$Set),"topvars")
mem_per_thread<-floor(args$m/args$c)
kuenm_args <- paste("threads=1")
if(bias_logic){
  kuenm_args <- paste(kuenm_args, paste0('biasfile="',normalizePath("bias.asc"),'" biastype=3'))
}

maxent_path <- getwd()


#==================================#
##### KUENM CALIBRATION MODELS #####
#==================================#

if(as.numeric(readLines("current_step.log"))<2){
  unlink("Candidate_Models", recursive=TRUE)
  cat(paste0("\n", Sys.time(), " ::: Creating candidate Maxent models for calibration with kuenm :::\n"))
  reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10) # original
  #reg_mult = c(0.1, seq(0.25, 1, 0.25), seq(2, 6, 2), 8, 10)
  #f_clas = "all"
  f_clas = c("l", "q", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh", "pt", "ph", "th", "lqp", "lqt", "lqh", "lpt", "lph", "lth", "qpt", "qph", "qth", "pth", "lqpt", "lqph", "lqth", "lpth", "qpth", "lqpth")
  #args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
  # note that some arguments are fixed in the function and should not be changed
  kuenm_cal(occ.joint = "Sp_joint.csv",
            occ.tra = "Sp_train.csv",
            M.var.dir = "M_variables",
            batch = "Candidate_models",
            out.dir = "Candidate_Models",
            max.memory = mem_per_thread, #args$m,
            reg.mult = reg_mult,
            f.clas = f_clas,
            args = kuenm_args,
            maxent.path = maxent_path,
            wait = FALSE,
            run = FALSE) #TRUE)
  
  if(length(readLines("Candidate_models.sh"))>args$c){
    jobs<-args$c
  }else{
      jobs<-length(readLines("Candidate_models.sh"))
  }
  
  system(paste0("tail -n +2 Candidate_models.sh | parallel --bar --retries 10 -j ", jobs," 'if eval {} 2>&1 | grep Error; then sleep 10; false; fi'"))
  writeLines("2", "current_step.log")
  cat("\n")
}


#===========================================#
##### KUENM EVALUATE CALIBRATION MODELS #####
#===========================================#

if(as.numeric(readLines("current_step.log"))<3){
  unlink("Calibration_results", recursive=TRUE)
  cat(paste0("\n", Sys.time(), " ::: Evaluating candidate Maxent models and selecting best model by ommission rates and AICc :::\n"))
  # paral_proc <- FALSE # make this true to perform pROC calculations in parallel, recommended
  # only if a powerfull computer is used (see function's help)
  # Note, some of the variables used here as arguments were already created for previous function
  cal_eval <- kuenm_ceval(path = "Candidate_Models",
                          occ.joint = "Sp_joint.csv",
                          occ.tra = "Sp_train.csv",
                          occ.test = "Sp_test.csv",
                          batch = "Candidate_models",
                          out.eval = "Calibration_results",
                          threshold = 5,
                          rand.percent = 50,
                          iterations = 500,
                          kept = TRUE,
                          selection = "OR_AICc",
                          parallel.proc = FALSE, silent=FALSE)
  writeLines("3", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<4){
  for(i in sets){
    out.eval=paste0("Calibration_results_",i)
    unlink(out.eval, recursive=TRUE)
    cat(paste0("\n", Sys.time(), " ::: Evaluating candidate Maxent models (", i,") and selecting best model by ommission rates and AICc :::\n"))
    kuenm_ceval2(out.eval=out.eval, set=i, threshold=5)
  }
  writeLines("5", "current_step.log")
}

# if(as.numeric(readLines("current_step.log"))<4){
#   unlink("Calibration_results_bioclim", recursive=TRUE)
#   dir.create("tmp")
#   system("mv Candidate_Models/*land* tmp")
#   system("mv Candidate_Models/*topo* tmp")
#   system("mv Candidate_Models/*topvars* tmp")
#   system("grep bioc Candidate_models.sh > bioc_models.sh")
#   cat(paste0("\n", Sys.time(), " ::: Evaluating candidate Maxent models (bioclim) and selecting best model by ommission rates and AICc :::\n"))
#   cal_eval <- kuenm_ceval(path = "Candidate_Models",
#                           occ.joint = "Sp_joint.csv",
#                           occ.tra = "Sp_train.csv",
#                           occ.test = "Sp_test.csv",
#                           batch = "bioc_models",
#                           out.eval = "Calibration_results_bioclim",
#                           threshold = 5,
#                           rand.percent = 50,
#                           iterations = 500,
#                           kept = TRUE,
#                           selection = "OR_AICc",
#                           parallel.proc = FALSE, silent=FALSE)
#   writeLines("4", "current_step.log")
# }
# 
# if(as.numeric(readLines("current_step.log"))<5){
#   unlink("Calibration_results_topvars", recursive=TRUE)
#   system("mv Candidate_Models/*bioc* tmp")
#   system("mv tmp/*topvars* Candidate_Models")
#   system("grep topvars Candidate_models.sh > topvars_models.sh")
#   cat(paste0("\n", Sys.time(), " ::: Evaluating candidate Maxent models (topvars) and selecting best model by ommission rates and AICc :::\n"))
#   cal_eval <- kuenm_ceval(path = "Candidate_Models",
#                           occ.joint = "Sp_joint.csv",
#                           occ.tra = "Sp_train.csv",
#                           occ.test = "Sp_test.csv",
#                           batch = "topvars_models",
#                           out.eval = "Calibration_results_topvars",
#                           threshold = 5,
#                           rand.percent = 50,
#                           iterations = 500,
#                           kept = TRUE,
#                           selection = "OR_AICc",
#                           parallel.proc = FALSE, silent=FALSE)
#   system("mv tmp/* Candidate_Models")
#   system("rm -r tmp topvars_models.sh bioc_models.sh")
#   writeLines("5", "current_step.log")
# }


#==================================#
##### KUENM SETUP FINAL MODELS #####
#==================================#

if(as.numeric(readLines("current_step.log"))<6){
  unlink("Final_Models", recursive=TRUE)
  cat(paste0("\n", Sys.time(), " ::: Creating Maxent models with optimal/selected parameters :::\n"))
  if(file.exists("Sp_ind.csv")){
    occ_file="Sp_joint.csv"
  }else{
    warn<-paste0("\n", Sys.time(), " ::: WARNING ::: Using training dataset for Final Model. Evaluating with test dataset rather than 'independent set'. \n")
    write(warn,"warnings.txt",append=TRUE)
    cat(warn)
    occ_file="Sp_train.csv"
  }
  #args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
  # "outputgrids=false" which avoids writing grids of replicated models and only writes the 
  # summary of them (e.g., average, median, etc.) when rep.n > 1
  # note that some arguments are fixed in the function and should not be changed
  # Again, some of the variables used here as arguments were already created for previous functions
  kuenm_mod(occ.joint = occ_file,
            M.var.dir = "M_variables",
            out.eval = "Calibration_results",
            batch = "Final_models",
            rep.n = 10,
            rep.type = "Bootstrap",
            jackknife = TRUE,
            out.dir = "Final_Models",
            max.memory = mem_per_thread, #args$m,
            out.format = "raw",
            project = gvars_logic, # NOT PROJECTING BECAUSE IT TAKES FOREVER
            G.var.dir = gvars,
            ext.type = "all",
            write.mess = FALSE,
            write.clamp = FALSE,
            maxent.path = maxent_path,
            args = kuenm_args,
            wait = FALSE,
            run = FALSE) #TRUE)
  
  if(length(readLines("Final_models.sh"))>args$c){
    jobs<-args$c
  }else{
    jobs<-length(readLines("Final_models.sh"))
  }
  
  system(paste0("tail -n +2 Final_models.sh | parallel --bar --retries 10 -j ", jobs," 'if eval {} 2>&1 | grep Error; then sleep 10; false; fi'"))
  
  writeLines("6", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<7){
  if(file.exists("Sp_ind.csv")){
    occ_file="Sp_joint.csv"
  }else{
    warn<-paste0("\n", Sys.time(), " ::: WARNING ::: Using training dataset for Final Model. Evaluating with test dataset rather than 'independent set'. \n")
    write(warn,"warnings.txt",append=TRUE)
    cat(warn)
    occ_file="Sp_train.csv"
  }
  for(i in sets){
    out.eval=paste0("Calibration_results_",i)
    batch=paste0("Final_models_",i)
    out.dir=paste0("Final_Models_",i)
    unlink(out.dir, recursive=TRUE)
    
    cat(paste0("\n", Sys.time(), " ::: Creating Maxent models (",i,") with optimal/selected parameters :::\n"))
    kuenm_mod(occ.joint = occ_file,
              M.var.dir = "M_variables",
              out.eval = out.eval,
              batch = batch,
              rep.n = 10,
              rep.type = "Bootstrap",
              jackknife = TRUE,
              out.dir = out.dir,
              max.memory = mem_per_thread, #args$m,
              out.format = "raw",
              project = gvars_logic, # NOT PROJECTING BECAUSE IT TAKES FOREVER
              G.var.dir = gvars,
              ext.type = "all",
              write.mess = FALSE,
              write.clamp = FALSE,
              maxent.path = maxent_path,
              args = kuenm_args,
              wait = FALSE,
              run = FALSE) #TRUE)
    
    if(length(readLines(paste0(batch,".sh")))>args$c){
      jobs<-args$c
    }else{
      jobs<-length(readLines(paste0(batch,".sh")))
    }
    
    system(paste0("tail -n +2 ",paste0(batch,".sh")," | parallel --bar --retries 10 -j ", jobs," 'if eval {} 2>&1 | grep Error; then sleep 10; false; fi'"))
    
  }
  writeLines("8", "current_step.log")
}


#=====================================#
##### KUENM EVALUATE FINAL MODELS #####
#=====================================#

if(as.numeric(readLines("current_step.log"))<9){
  unlink("Final_Models_evaluation", recursive=TRUE)
  cat(paste0("\n", Sys.time(), " ::: Running final evaluation of Maxent models with independent data :::\n"))
  if(file.exists("Sp_ind.csv")){
    occ_file1="Sp_joint.csv"
    occ_file2="Sp_ind.csv"
  }else{
    warn<-paste0("\n", Sys.time(), " ::: WARNING ::: Using training dataset for Final Model. Evaluating with test dataset rather than 'independent set'. \n")
    write(warn,"warnings.txt",append=TRUE)
    cat(warn)
    occ_file1="Sp_train.csv"
    occ_file2="Sp_test.csv"
  }
  fin_eval <- kuenm_feval(path = "Final_Models",
                          occ.joint = occ_file1,
                          occ.ind = occ_file2,
                          replicates = TRUE,
                          out.eval = "Final_Models_evaluation",
                          threshold = 5,
                          rand.percent = 50,
                          iterations = 500,
                          parallel.proc = FALSE)
  writeLines("9", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<10){
  if(file.exists("Sp_ind.csv")){
    occ_file1="Sp_joint.csv"
    occ_file2="Sp_ind.csv"
  }else{
    warn<-paste0("\n", Sys.time(), " ::: WARNING ::: Using training dataset for Final Model. Evaluating with test dataset rather than 'independent set'. \n")
    write(warn,"warnings.txt",append=TRUE)
    cat(warn)
    occ_file1="Sp_train.csv"
    occ_file2="Sp_test.csv"
  }
  for(i in sets){
    path=paste0("Final_Models_",i)
    out.eval=paste0("Final_Models_",i,"_evaluation")
    unlink(out.eval, recursive=TRUE)
    cat(paste0("\n", Sys.time(), " ::: Running final evaluation of Maxent models (", i,") with independent data :::\n"))
    fin_eval <- kuenm_feval(path = path,
                            occ.joint = occ_file1,
                            occ.ind = occ_file2,
                            replicates = TRUE,
                            out.eval = out.eval,
                            threshold = 5,
                            rand.percent = 50,
                            iterations = 500,
                            parallel.proc = FALSE)
  }
  writeLines("11", "current_step.log")
}

#=====================================#
########### SUMMARIZE MODELS ##########
#=====================================#
# This function will:
# 1. Calculate the mean from the top models
# 2. Collate statistics
# 3. Convert raw output to cloglog
# 4. Threshold using tenth  percentile  training  presence

cat(paste0("\n", Sys.time(), " ::: Summarizing Results :::\n"))
summarizeModels("Final_Models")

for(i in sets){
  cat(paste0("\n", Sys.time(), " ::: Summarizing Results (",i,") :::\n"))
  summarizeModels(paste0("Final_Models_",i))
}
