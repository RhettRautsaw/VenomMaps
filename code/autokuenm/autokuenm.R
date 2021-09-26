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
# 2. Subset distribution maps by a species name (or create convex hull from points)
# 3. Buffer the resulting distribution and crop environmental data by the buffered extent
# 4. Generate random background points across the landscape and create occurrence record partitions for training and testing (ENMeval)
# 5. Create exhaustive combination of all environmental variables or perform heuristic search on sets
#     - Heuristic search involves identifying non-colinear variables with high permutation importance
# 6. Run kuenm
#
#
#=============#
##### USE #####
#=============#
#
# autokuenm.R -o occ.csv -e env.tab -d combined_distributions.geojson -sp Agkistrodon_bilineatus \
#             -spcol1 Species -spcol2 Species -latcol latitude -loncol longitude -b 100 -bg 2000 \
#             -g G_variables -c 8 -m 62000
#
#
#======================#
##### REQUIREMENTS #####
#======================#

packages <- c("argparse","readr","dplyr","stringr","parallel","memuse",
              "sf","raster","dismo","rJava","virtualspecies",
              "usdm","ENMeval","kuenm")

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
parser$add_argument("-o", default="occ.csv", 
                    help="csv of occurrence records. Can contain multiple species and will be filtered by -sp and -spcol flags. [default: \"%(default)s\"]")
parser$add_argument("-e", default="env.tab", 
                    help="For exhaustive search, provide a folder containing all environmental rasters. For heuristic search, provide a text file containing paths to all environmental rasters divided into Sets (Must be in .tif or .asc format). [default: \"%(default)s\"]")
parser$add_argument("-d", default="distributions.geojson", 
                    help="Distribution shapefile (e.g., GeoJSON or SHP) or 'mcp' if no distribution exists and you want a distribution to . Can contain multiple species and will be filtered by -sp and -spcol flags. [default: \"%(default)s\"]")
parser$add_argument("-sp", default="Agkistrodon_bilineatus", 
                    help="Name of species to model. [default: \"%(default)s\"]")
parser$add_argument("-spcol1", default="Species", 
                    help="Name of column with species ID in occurrence records csv. [default: \"%(default)s\"]")
parser$add_argument("-spcol2", default="Species", 
                    help="Name of column with species ID in distribution attribute table. [default: \"%(default)s\"]")
parser$add_argument("-latcol", default="latitude", 
                    help="Name of column with latitude in occurrence records csv [default: \"%(default)s\"]")
parser$add_argument("-loncol", default="longitude", 
                    help="Name of column with longitude in occurrence records csv [default: \"%(default)s\"]")
parser$add_argument("-b", type="integer", default=100, metavar="number",
                    help="Buffer for distribution in km. [default %(default)s]")
parser$add_argument("-bg", type="integer", default=2000, metavar="number",
                    help="Number of random background points to generate for assigning blocks [default %(default)s]")
parser$add_argument("-g", default="NA", 
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

if(file.exists(env) && !dir.exists(env)){
  search<-"heuristic"
  suppressMessages(env.tab<-readr::read_tsv(env))
  env.list<-normalizePath(env.tab$Var)
}else{
  search<-"exhaustive"
  env.list<-list.files(env, ".tif|.asc", full.names = T)
}

if(args$d=="mcp"){
  dist<-args$d
}else{
  dist<-normalizePath(args$d) 
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
    cat(paste0("\t Buffer Distance:\t", args$b, " km\n"))
  cat(paste0("\t Background Points:\t", args$bg,"\n"))
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

'%notin%' <- Negate('%in%')

## Variable inflation factor based layer selection
## This is based on code by Mike Belitz (mbelitz/Odo_SDM_Rproj) and Shelly Gaynor
SelectVariables<-function(stack, occ){
  print("Running basic MaxEnt to determine variable importance")
  model<-dismo::maxent(stack,occ)
  
  print("Extracting Permutation Importance")
  vIMP <- as.data.frame(as.table(model@results)) %>% 
    dplyr::rename(variables = 1, rem = 2, permutation.importance = 3) %>% 
    dplyr::select(variables, permutation.importance) %>%
    dplyr::filter(stringr::str_detect(variables, '.permutation.importance')) %>% 
    dplyr::mutate(Variables = stringr::word(variables,  sep = fixed(".permutation.importance"))) %>%
    dplyr::select(Variables, permutation.importance) %>%
    dplyr::filter(permutation.importance > 0)
  
  print(paste0(names(stack)[names(stack) %notin% vIMP$Variables], " removed due to 0% permutation importance"))
  
  stack<-stack[[which(names(stack) %in% vIMP$Variables)]]
  
  print("Calculating Variance Inflation Factors")
  pVIF <- usdm::vif(stack)
  
  while(max(pVIF$VIF) >= 10){
    # Left join pVIF with the vIMP from the MaxEnt model
    jdf <- suppressMessages(left_join(pVIF, vIMP))
    # select the variable with the highest VIF and lowest model importance
    lowVar <- jdf %>% 
      dplyr::filter(VIF > 10) %>% 
      dplyr::filter(VIF == sort(VIF, decreasing = TRUE)[1] |
                      VIF == sort(VIF, decreasing = TRUE)[2]) %>% 
      dplyr::filter(permutation.importance == min(permutation.importance))
    # Drop layer from stack
    stack <- raster::dropLayer(stack, as.character(lowVar$Variables))
    print(paste0(as.character(lowVar$Variables), " has been removed"))
    #print("Re-calculating Variance Inflation Factors")
    pVIF <- usdm::vif(stack)
  }
  
  print("Re-running basic MaxEnt")
  model<-dismo::maxent(stack,occ)
  
  print("Re-extracting Permutation Importance")
  vIMP <- as.data.frame(as.table(model@results)) %>% 
    dplyr::rename(variables = 1, rem = 2, permutation.importance = 3) %>% 
    dplyr::select(variables, permutation.importance) %>%
    dplyr::filter(stringr::str_detect(variables, '.permutation.importance')) %>% 
    dplyr::mutate(Variables = stringr::word(variables,  sep = fixed(".permutation.importance"))) %>%
    dplyr::select(Variables, permutation.importance) %>%
    dplyr::filter(permutation.importance > 1)
  
  print(paste0(names(stack)[names(stack) %notin% vIMP$Variables], " removed due to < 1% permutation importance"))
  
  stack<-stack[[which(names(stack) %in% vIMP$Variables)]]
  
  print(names(stack))
  return(stack)
}


#============================================#
##### MAKE NEW DIRECTORY AND COPY MAXENT #####
#============================================#

#unlink(args$sp, recursive=T)
if(!dir.exists(args$sp)){
  dir.create(args$sp)
  invisible(file.copy(from="maxent.jar", to=paste0(args$sp,"/maxent.jar")))
}

setwd(args$sp)


#================================#
##### RECORD/ENV DATA PREP #######
#================================#

if(!file.exists("current_step.log")){
  
  #===================================#
  ##### READ & FILTER OCC RECORDS #####
  #===================================#
  
  cat(paste0("\n", Sys.time(), " ::: Reading ",args$o, " :::\n"))
  suppressMessages(occ <- read_csv(occ))
  cat(paste0("\n", Sys.time(), " ::: Filtering ",nrow(occ), " records from ", args$o, " for ", args$sp," :::\n"))
  occ2<-occ %>% filter((!!as.name(args$spcol1))==args$sp)
  cat(paste0("\n", Sys.time(), " ::: ", nrow(occ2), " records remaining :::\n"))
  
  
  #===============================================#
  ##### READ, FILTER, AND BUFFER DISTRIBUTION #####
  #===============================================#
  
  if(dist=="mcp"){
    dist <- st_as_sf(occ2, coords = c(args$loncol,args$latcol), crs=4326) %>%
      st_union() %>%
      st_convex_hull()
  }else{
    cat(paste0("\n", Sys.time(), " ::: Reading ",args$d, " :::\n"))
    dist<-sf::st_read(dist,crs=4326, quiet=T)
    cat(paste0("\n", Sys.time(), " ::: Filtering ",args$d, " for ", args$sp," :::\n"))
    dist<-dist %>% filter((!!as.name(args$spcol2))==args$sp) 
  }
  cat(paste0("\n", Sys.time(), " ::: Buffering filtered distribution by ", args$b," km :::\n"))
  dist2<-st_transform(st_buffer(st_transform(dist,crs=3857),(args$b*1000)),crs=4326)
  
  
  #=================================================================#
  ##### STACK ENV DATA, CROP, AND MASK BY BUFFERED DISTRIBUTION #####
  #=================================================================#
  
  cat(paste0("\n", Sys.time(), " ::: Reading and stacking ",args$e, " :::\n"))
  envStack <- raster::stack(env.list)
  cat(paste0("\n", Sys.time(), " ::: Cropping environmental stack by buffered distribution extent :::\n"))
  envStack2<-crop(envStack,extent(dist2))
  ## NOT SURE IF I SHOULD MASK OR NOT.
  #envStack2<-mask(envStack2,dist2) # TESTS SUGGEST DO NOT MASK, BUT LEFT CODE HERE ANYWAYS
  
  
  #======================================#
  ##### SUBSAMPLE OCC RECORDS BY ENV #####
  #======================================#
  
  cat(paste0("\n", Sys.time(), " ::: Filtering occ records to 1 per environmental grid cell :::\n"))
  # ISNA<-which(is.na(raster::extract(envStack2,occ2[,c(args$loncol,args$latcol)])[,1]))
  # if(length(ISNA)==0){}else{
  #   occ3<-occ2[-ISNA,]
  # }
  tmp<-raster::extract(envStack2,occ2[,c(args$loncol,args$latcol)])
  occ2<-occ2[complete.cases(tmp),]
  tmpocc <- as.data.frame(gridSample(xy = as.matrix(occ2[,c(args$loncol,args$latcol)]), r = envStack2[[1]], n = 1))
  if(nrow(tmpocc)<20){
    warn<-paste0("\n", Sys.time(), " ::: WARNING ::: Low number of records...proceeding with spatial autocorrelation by allowing more than 1 record per grid cell.\n")
    write(warn,"warnings.txt",append=TRUE)
    cat(warn)
    occ2<-cbind(args$sp,occ2[,c(args$loncol,args$latcol)])
    cat(paste0("\n", Sys.time(), " ::: ", nrow(occ2), " records remaining :::\n"))
  }else{
    occ2<-cbind(args$sp,tmpocc)
    cat(paste0("\n", Sys.time(), " ::: ", nrow(occ2), " records remaining :::\n"))
  }
  colnames(occ2)<-c("Species", "Longitude", "Latitude")
  occ2$Latitude<-as.numeric(occ2$Latitude)
  occ2$Longitude<-as.numeric(occ2$Longitude)
  
  
  #============================================================#
  ##### CREATE SUBSET OF OCC RECORDS FOR INDEPENDENT TEST  #####
  #============================================================#
  
  cat(paste0("\n", Sys.time(), " ::: Selecting random 5% of records as independent dataset :::\n"))
  sp_ind_n<-as.integer(0.05*nrow(occ2))
  if(sp_ind_n >= 6){
    set.seed(12345)
    sp_ind_rows<-sample(1:nrow(occ2), sp_ind_n)
    sp_ind<-occ2[sp_ind_rows,]
    occ3<-occ2[-sp_ind_rows,]
    write.csv(sp_ind,"Sp_ind.csv", row.names=F, quote=F)
  }else{
    warn<-paste0("\n", Sys.time(), " ::: WARNING ::: Too few records. Testing and Independent datasets will actually be identical.\n")
    write(warn,"warnings.txt",append=TRUE)
    cat(warn)
    occ3<-occ2
    # set.seed(12345)
    # sp_ind_rows<-sample(1:nrow(occ2), 6)
    # sp_ind<-occ2[sp_ind_rows,]
    # occ3<-occ2[-sp_ind_rows,]
    # write.csv(sp_ind,"Sp_ind.csv", row.names=F, quote=F)
  }
  
  
  #=========================================================================#
  ##### DEFINE BACKGROUND POINTS AND PARTITION OCC RECORDS INTO BLOCKS  #####
  #=========================================================================#
  
  cat(paste0("\n", Sys.time(), " ::: Generating ",args$bg," random background points across cropped environmental grid :::\n"))
  bg <- randomPoints(envStack2[[1]], n = args$bg)
  bg <- as.data.frame(bg)
  
  cat(paste0("\n", Sys.time(), " ::: Partitioning occ records for training and testing :::\n"))
  #block <- get.block(occ3[, c("Longitude", "Latitude")], bg)
  block<-get.checkerboard2(occ = occ3[, c("Longitude", "Latitude")], env = envStack2, bg.coords = bg, aggregation.factor = c(2,2))
  bg<-cbind("bg",bg)
  names(bg)<-c("Species","Longitude","Latitude")
  
  png("blocks_records.png", width = 1080, height=720)
  par(mfrow=c(1,2))
  plot(envStack2[[1]], legend = FALSE, main="occurrence records")
  points(occ3[, c("Longitude", "Latitude")], pch = 21, bg = block$occ.grp)
  plot(envStack2[[1]], legend = FALSE, main="background points")
  points(bg[, c("Longitude", "Latitude")], pch = 21, bg = block$bg.grp)
  invisible(dev.off())
  
  
  #===============================================#
  ##### SUBSET ONE BLOCK FOR TESTING/TRAINING #####
  #===============================================#
  
  cat(paste0("\n", Sys.time(), " ::: Selecting one partition as test dataset and the rest as training data :::\n"))
  set.seed(12345)
  grp<-sample(c(1,2,3,4), 1)
  sp_train<-occ3[which(block$occ.grp!=grp),]
  sp_test<-occ3[which(block$occ.grp==grp),]
  write.csv(sp_test,"Sp_test.csv", row.names=F, quote=F)
  write.csv(sp_train,"Sp_train.csv", row.names=F, quote=F)
  write.csv(occ3,"Sp_joint.csv", row.names=F, quote=F)
  cat(paste0("\n", Sys.time(), " ::: Using ", nrow(sp_train), " records for training and ", nrow(sp_test), " records for testing :::\n"))
  if(!file.exists("Sp_ind.csv")){
    cat(paste0("\n", Sys.time(), " ::: Independent dataset (Sp_ind.csv) will be identical to testing set due to low records :::\n"))
    write.csv(sp_test,"Sp_ind.csv", row.names=F, quote=F)
  }
  
  
  #========================================#
  ##### EXHAUSTIVE ENVIRONMENTAL SETS  #####
  #========================================#
  
  if(search=="exhaustive"){
    cat(paste0("\n", Sys.time(), " ::: Creating Exhaustive Set of Environmental Variables :::\n"))
    format<-tools::file_ext(env.list[1])
    if(format=="tif"){format="GTiff"}else{format="ascii"}
    combinations <- kuenm_varcomb(var.dir = env, out.dir = "M_variables",
                                  min.number = 2, in.format = format,
                                  out.format = "ascii")
  }
  
  #========================================#
  ##### HEURISTIC ENVIRONMENTAL SETS  #####
  #========================================#
  
  if(search=="heuristic"){
    cat(paste0("\n", Sys.time(), " ::: Creating Heuristic Set of Environmental Variables :::\n"))
    for(i in unique(env.tab$Set)){
      cat(paste0("\n", Sys.time(), " ::: Creating Set ",i," :::\n"))
      set<-which(env.tab$Set==i)
      tmpStack<-envStack2[[set]]
      newStack<-SelectVariables(occ = occ3[,2:3], stack = tmpStack)
      dir.create(paste0("M_variables/Set_",i), recursive = T)
      writeRaster(newStack, paste0("M_variables/Set_",i,"/",names(newStack),".asc"), bylayer=T, format="ascii")
      write(names(newStack),paste0("Variables_",i,".txt"))
    }
    cat(paste0("\n", Sys.time(), " ::: Creating Set topvars :::\n"))
    newStack<-SelectVariables(occ = occ3[,2:3], stack = envStack2)
    dir.create(paste0("M_variables/Set_topvars"), recursive = T)
    writeRaster(newStack, paste0("M_variables/Set_topvars/",names(newStack),".asc"), bylayer=T, format="ascii")
    write(names(newStack),paste0("Variables_topvars.txt"))
  }
  writeLines("1", "current_step.log")
}


#=====================#
##### KUENM SETUP #####
#=====================#
mem_per_thread<-floor(args$m/args$c)
kuenm_args <- paste0("threads=1") #,args$c)
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

# cat(paste0("\n", Sys.time(), " ::: Creating mxe cache files for variables for quicker model testing Maxent models for calibration with kuenm :::\n"))
# system("mkdir tmp")
# cache_script<-c("java -jar ./maxent.jar environmentallayers=./M_variables/SET samplesfile=./Sp_joint.csv outputdirectory=./tmp betamultiplier=0.1 autofeature=false linear=true quadratic=false product=false threshold=false hinge=false threads=1 cache=TRUE extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=false plots=false pictures=false outputformat=raw warnings=false visible=false redoifexists autorun")
# for(i in list.files("M_variables")){
#   s<-gsub("SET", i, cache_script)
#   system(s)
# }
# system("rm -r tmp")

#system(paste0("head -2 Candidate_models.sh | tail -1 | parallel -j 1 '{}'"))
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
                        occ.joint = "Sp_train.csv",
                        occ.tra = "Sp_train.csv",
                        occ.test = "Sp_test.csv",
                        batch = "Candidate_models",
                        out.eval = "Calibration_results",
                        threshold = 5,
                        rand.percent = 50,
                        iterations = 500,
                        kept = TRUE,
                        selection = "OR_AICc",
                        parallel.proc = FALSE)
writeLines("3", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<4){
unlink("Calibration_results_bioclim", recursive=TRUE)
dir.create("tmp")
system("mv Candidate_Models/*land* tmp")
system("mv Candidate_Models/*topo* tmp")
system("mv Candidate_Models/*topvars* tmp")
system("grep bioc Candidate_models.sh > bioc_models.sh")
cat(paste0("\n", Sys.time(), " ::: Evaluating candidate Maxent models (bioclim) and selecting best model by ommission rates and AICc :::\n"))
cal_eval <- kuenm_ceval(path = "Candidate_Models",
                        occ.joint = "Sp_train.csv",
                        occ.tra = "Sp_train.csv",
                        occ.test = "Sp_test.csv",
                        batch = "bioc_models",
                        out.eval = "Calibration_results_bioclim",
                        threshold = 5,
                        rand.percent = 50,
                        iterations = 500,
                        kept = TRUE,
                        selection = "OR_AICc",
                        parallel.proc = FALSE)
writeLines("4", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<5){
unlink("Calibration_results_topvars", recursive=TRUE)
system("mv Candidate_Models/*bioc* tmp")
system("mv tmp/*topvars* Candidate_Models")
system("grep topvars Candidate_models.sh > topvars_models.sh")
cat(paste0("\n", Sys.time(), " ::: Evaluating candidate Maxent models (topvars) and selecting best model by ommission rates and AICc :::\n"))
cal_eval <- kuenm_ceval(path = "Candidate_Models",
                        occ.joint = "Sp_train.csv",
                        occ.tra = "Sp_train.csv",
                        occ.test = "Sp_test.csv",
                        batch = "topvars_models",
                        out.eval = "Calibration_results_topvars",
                        threshold = 5,
                        rand.percent = 50,
                        iterations = 500,
                        kept = TRUE,
                        selection = "OR_AICc",
                        parallel.proc = FALSE)
system("mv tmp/* Candidate_Models")
system("rm -r tmp topvars_models.sh bioc_models.sh")
writeLines("5", "current_step.log")
}


#==================================#
##### KUENM SETUP FINAL MODELS #####
#==================================#

if(as.numeric(readLines("current_step.log"))<6){
unlink("Final_Models", recursive=TRUE)
cat(paste0("\n", Sys.time(), " ::: Creating Maxent models with optimal/selected parameters :::\n"))
#args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions
kuenm_mod(occ.joint = "Sp_joint.csv",
          M.var.dir = "M_variables",
          out.eval = "Calibration_results",
          batch = "Final_models",
          rep.n = 10,
          rep.type = "Bootstrap",
          jackknife = TRUE,
          out.dir = "Final_Models",
          max.memory = mem_per_thread, #args$m,
          out.format = "logistic",
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
unlink("Final_Models_bioclim", recursive=TRUE)
cat(paste0("\n", Sys.time(), " ::: Creating Maxent models (bioclim) with optimal/selected parameters :::\n"))
kuenm_mod(occ.joint = "Sp_joint.csv",
          M.var.dir = "M_variables",
          out.eval = "Calibration_results_bioclim",
          batch = "Final_models_bioclim",
          rep.n = 10,
          rep.type = "Bootstrap",
          jackknife = TRUE,
          out.dir = "Final_Models_bioclim",
          max.memory = mem_per_thread, #args$m,
          out.format = "logistic",
          project = gvars_logic, # NOT PROJECTING BECAUSE IT TAKES FOREVER
          G.var.dir = gvars,
          ext.type = "all",
          write.mess = FALSE,
          write.clamp = FALSE,
          maxent.path = maxent_path,
          args = kuenm_args,
          wait = FALSE,
          run = FALSE) #TRUE)

if(length(readLines("Final_models_bioclim.sh"))>args$c){
  jobs<-args$c
}else{
  jobs<-length(readLines("Final_models_bioclim.sh"))
}

system(paste0("tail -n +2 Final_models_bioclim.sh | parallel --bar --retries 10 -j ", jobs," 'if eval {} 2>&1 | grep Error; then sleep 10; false; fi'"))

writeLines("7", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<8){
unlink("Final_Models_topvars", recursive=TRUE)
cat(paste0("\n", Sys.time(), " ::: Creating Maxent models (topvars) with optimal/selected parameters :::\n"))
kuenm_mod(occ.joint = "Sp_joint.csv",
          M.var.dir = "M_variables",
          out.eval = "Calibration_results_topvars",
          batch = "Final_models_topvars",
          rep.n = 10,
          rep.type = "Bootstrap",
          jackknife = TRUE,
          out.dir = "Final_Models_topvars",
          max.memory = mem_per_thread, #args$m,
          out.format = "logistic",
          project = gvars_logic, # NOT PROJECTING BECAUSE IT TAKES FOREVER
          G.var.dir = gvars,
          ext.type = "all",
          write.mess = FALSE,
          write.clamp = FALSE,
          maxent.path = maxent_path,
          args = kuenm_args,
          wait = FALSE,
          run = FALSE) #TRUE)

if(length(readLines("Final_models_topvars.sh"))>args$c){
  jobs<-args$c
}else{
  jobs<-length(readLines("Final_models_topvars.sh"))
}

system(paste0("tail -n +2 Final_models_topvars.sh | parallel --bar --retries 10 -j ", jobs," 'if eval {} 2>&1 | grep Error; then sleep 10; false; fi'"))

writeLines("8", "current_step.log")
}


#=====================================#
##### KUENM EVALUATE FINAL MODELS #####
#=====================================#

if(as.numeric(readLines("current_step.log"))<9){
unlink("Final_Models_evaluation", recursive=TRUE)
cat(paste0("\n", Sys.time(), " ::: Running final evaluation of Maxent models with independent data :::\n"))
fin_eval <- kuenm_feval(path = "Final_Models",
                        occ.joint = "Sp_joint.csv",
                        occ.ind = "Sp_ind.csv",
                        replicates = TRUE,
                        out.eval = "Final_Models_evaluation",
                        threshold = 5,
                        rand.percent = 50,
                        iterations = 500,
                        parallel.proc = FALSE)
writeLines("9", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<10){
unlink("Final_Models_bioclim_evaluation", recursive=TRUE)
cat(paste0("\n", Sys.time(), " ::: Running final evaluation of Maxent models (bioclim) with independent data :::\n"))
fin_eval <- kuenm_feval(path = "Final_Models_bioclim",
                        occ.joint = "Sp_joint.csv",
                        occ.ind = "Sp_ind.csv",
                        replicates = TRUE,
                        out.eval = "Final_Models_bioclim_evaluation",
                        threshold = 5,
                        rand.percent = 50,
                        iterations = 500,
                        parallel.proc = FALSE)
writeLines("10", "current_step.log")
}

if(as.numeric(readLines("current_step.log"))<11){
unlink("Final_Models_topvars_evaluation", recursive=TRUE)
cat(paste0("\n", Sys.time(), " ::: Running final evaluation of Maxent models (topvars) with independent data :::\n"))
fin_eval <- kuenm_feval(path = "Final_Models_topvars",
                        occ.joint = "Sp_joint.csv",
                        occ.ind = "Sp_ind.csv",
                        replicates = TRUE,
                        out.eval = "Final_Models_topvars_evaluation",
                        threshold = 5,
                        rand.percent = 50,
                        iterations = 500,
                        parallel.proc = FALSE)
writeLines("11", "current_step.log")
}


#======================#
##### ARCHIVE CODE #####
#======================#

# parser$add_argument("-colinear", type="double", default=0.75, metavar="number",
#                     help="Not used for exhaustive search. Multicollinearity cutoff. Variables which have correlation coefficient greater than this number will be sorted and the first value selected [default %(default)s]")
# parser$add_argument("-axes", type="integer"r  , default=2, metavar="number",
#                     help="Not used for exhaustive search. PCA is used to select the most important variables for the final Set. This represents the number of PC axes to use. [default %(default)s]")
# parser$add_argument("-contrib", type="integer", default=5, metavar="number",
#                     help="Not used for exhaustive search. PCA is used to select the most important variables for the final Set. This represents the relative percent contribution a variable must have to be included. [default %(default)s]")
# parser$add_argument("-k", action="store_true", default=FALSE,
#                     help="Keep calibration models [default off]")


# #=====================================================#
# ##### REMOVE COLINEAR ENV VARIABLES FOR EACH SET  #####
# #=====================================================#
# 
# for(i in unique(env.tab$Set)){
#   cat(paste0("\n", Sys.time(), " ::: Creating Set ",i," and removing collinear variables > ",args$colinear," :::\n"))
#   set<-which(env.tab$Set==i)
#   tmpStack<-envStack2[[set]]
#   png(paste0("Colinearity_",i,".png"), width=1080, height=720)
#   vars<-removeCollinearity(tmpStack, plot = TRUE, select.variables = FALSE, multicollinearity.cutoff = args$colinear)
#   dev.off()
#   vars<-lapply(vars, function(v){sort(v)})
#   vars<-unlist(lapply(vars, `[[`, 1))
#   tmpStack<-tmpStack[[which(names(tmpStack) %in% vars)]]
#   dir.create(paste0("M_variables/Set_",i), recursive = T)
#   writeRaster(tmpStack, paste0("M_variables/Set_",i,"/",names(tmpStack),".asc"), bylayer=T, format="ascii")
#   write(names(tmpStack),paste0("Variables_",i,".txt"))
# }
# 
# 
# #=========================================================#
# ##### REMOVE COLINEAR ENV VARIABLES FOR COMBINED SET  #####
# #=========================================================#
# 
# cat(paste0("\n", Sys.time(), " ::: Creating Set topvars and removing collinear variables > ",args$colinear," :::\n"))
# png("Colinearity_combined.png", width=1080, height=720)
# vars<-removeCollinearity(envStack2, plot = TRUE, select.variables = FALSE, multicollinearity.cutoff = args$colinear)
# invisible(dev.off())
# vars<-lapply(vars, function(v){sort(v)})
# vars<-unlist(lapply(vars, `[[`, 1))
# envStack2<-envStack2[[which(names(envStack2) %in% vars)]]
# 
# ## source("../useful_functions.R")
# ## tmp<-raster.pca(envStack,5)
# ## pcaStack<-tmp$rasters
# ## dir.create("M_variables/Set_1", recursive = T)
# ## writeRaster(pcaStack, paste0("M_variables/Set_1/",names(pcaStack),".asc"), bylayer=T, format="ascii")
# 
# 
# #=============================================#
# ##### USE PCA TO SELECT TOP ENV VARIABLES #####
# #=============================================#
# 
# cat(paste0("\n", Sys.time(), " ::: Using PCA to identify the most important (top) variables :::\n"))
# tmp<-raster::extract(envStack2, occ3[,c("Longitude","Latitude")])
# cat(paste0("\n", Sys.time(), " ::: Removing variables with 0 variance :::\n"))
# ISNA<-which(apply(tmp, 2, var)==0)
# if(length(ISNA)==0){}else{
#   tmp<-tmp[,-ISNA]
# }
# tmp<-cbind(occ3,tmp)
# PCA <- prcomp(tmp[,c(4:ncol(tmp))], center=TRUE, scale=TRUE)
# cat(paste0("\n", Sys.time(), " ::: Determining variables in first ", args$axes, " PC axes with greater than ", args$contrib, "% relative contribution :::\n"))
# loadings <- PCA$rotation
# loadings_relative_A <- t(t(abs(loadings))/rowSums(t(abs(loadings))))*100
# topvars<-vector()
# for(i in 1:args$axes){
#   contrib<-loadings_relative_A[,i]
#   vars<-names(which(contrib > args$contrib))
#   topvars<-c(topvars,vars)
# }
# topvars<-sort(unique(topvars))
# tmpStack<-envStack2[[which(names(envStack2) %in% topvars)]]
# dir.create(paste0("M_variables/Set_topvars"), recursive = T)
# writeRaster(tmpStack, paste0("M_variables/Set_topvars/",names(tmpStack),".asc"), bylayer=T, format="ascii")
# write(names(tmpStack),"Variables_topvars.txt")


# sets_var <- "Set1" # a vector of various sets can be used
# out_mop <- "MOP_results"
# percent <- 10
# paral <- FALSE # make this true to perform MOP calculations in parallel, recommended
# # only if a powerfull computer is used (see function's help)
# # Two of the variables used here as arguments were already created for previous functions
# kuenm_mmop(G.var.dir = G_var_dir, M.var.dir = M_var_dir, sets.var = sets_var, out.mop = out_mop,
#            percent = percent, parallel = paral)
