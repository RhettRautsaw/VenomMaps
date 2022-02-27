#### AUTOKUENM SUPPLEMENTAL FUNCTIONS ####

#================#
##### Negate #####
#================#

'%notin%' <- Negate('%in%')

#=========================#
#### Format Converter #####
#=========================#
# R = raw format raster
# l = logistic format raster
# c = cloglog format raster
# h = entropy

raw2log<-function(R,h){
  l<-(exp(h)*R)/(1+(exp(h)*R))
  return(l)
}

raw2clog<-function(R,h){
  c<-1-(exp(-exp(h)*R))
  return(c)
}

log2raw<-function(l,h){
  R=-((exp(-h)*l)/(l-1))
  return(R)
}


#===================#
#### Alpha Hull #####
#===================#
ahull <- function(occ_sf, concave_distance_lim = 5000, buffer_dist=50000){
  occ_sf<-st_transform(occ_sf, crs=3857)
  poly=ne_countries(scale="large", returnclass = "sf") %>% st_transform(crs=3857) %>% st_buffer(0) %>% st_union()
  ahull=concaveman(occ_sf, length_threshold=concave_distance_lim, concavity=2) %>% st_buffer(buffer_dist) %>% 
    st_intersection(poly) %>% st_buffer(0) %>% st_union() %>% st_transform(crs=4326) 
  return(ahull)  
}


#===========================#
##### smoothr fill_holes ####
#===========================#

hole_area <- function(x, crs) {
  sf::st_area(sf::st_sfc(sf::st_polygon(list(x)), crs = crs))
}
fill_holes <- function(x, threshold) {
  UseMethod("fill_holes")
}

#' @export
fill_holes.sfc <- function(x, threshold) {
  stopifnot(inherits(threshold, c("units", "numeric")),
            length(threshold) == 1)
  # check geometry types and get units of feature size
  if (!all(sf::st_is(x, c("POLYGON", "MULTIPOLYGON")))) {
    stop("fill_holes() only works for polygon features.")
  }
  
  # zero threshold returns the input features unchanged
  thresh_nounits <- as.numeric(threshold)
  if (thresh_nounits == 0) {
    return(x)
  } else if (thresh_nounits < 0) {
    stop("threshold cannont be negative")
  }
  
  # convert threshold to crs units
  if (is.na(sf::st_crs(x))) {
    area_units <- units::set_units(1, "m2")
  } else {
    area_units <- units::set_units(1, units(sf::st_area(x[1])),
                                   mode = "standard")
  }
  threshold <- units::set_units(threshold, area_units, mode = "standard")
  
  x_crs <- sf::st_crs(x)
  
  # loop over features, potentially multipart
  for (i in seq_along(x)) {
    # split up multipart geometries
    singles <- sf::st_cast(x[i], "POLYGON")
    # loop over single features
    for (j in seq_along(singles)) {
      # skip features with no holes
      if (length(singles[[j]]) > 1) {
        # area of holes
        a <- vapply(singles[[j]][-1], hole_area, 1, crs = x_crs)
        # assign units
        a <- units::set_units(a, area_units, mode = "standard")
        # remove holes not passing threshold
        singles[[j]] <- sf::st_polygon(singles[[j]][c(TRUE, a >= threshold)])
      }
    }
    # recombine
    if (length(singles) != 1) {
      singles <- sf::st_combine(singles)
    }
    x[[i]] <- singles[[1]]
  }
  # remove empty geometries
  sf::st_sfc(x)
}

#' @export
fill_holes.sf <- function(x, threshold) {
  sf::st_set_geometry(x, fill_holes(sf::st_geometry(x), threshold = threshold))
}


#=============#
#### Marea ####
#=============#
# M Area from Bioregions and polygon overlap
# cov = % theshold for bioregion covered by polygon

Marea<-function(bioregions, dist, cov=25){
  # Calculate area of each bioregion
  bioregions$TMPID = rownames(bioregions)
  bioregions$TMPAREA = as.numeric(st_area(st_buffer(st_transform(bioregions,crs=3857),0)))
  # Intersect bioregions with distribution and calculate area
  i<-st_intersection(st_buffer(st_transform(bioregions,crs=3857),0), st_buffer(st_transform(dist,crs=3857),0))
  i$TMPAREA2<-as.numeric(st_area(i))
  # Calculate % of bioregion covered by distribution and select bioregions over threshold
  m <- i %>% mutate(TMPPERC=(TMPAREA2/TMPAREA)*100) 
  m <- m %>% filter(TMPPERC>=cov) %>% dplyr::select(TMPID)
  # Extract these full bioregions from original bioregions file and union them together
  m<-st_set_geometry(m,NULL)
  m <- as.vector(m[,1])
  m2<-bioregions %>% filter(TMPID %in% m)
  m3<-st_union(m2)
  return(m3)
}


#====================#
#### spThin Model ####
#====================#
# Species with more occurrences get a larger distance for trimming
tmp<-data.frame(V1=c(1,30000),V2=c(0,25))
sp_mod<-lm(V2~log(V1), data=tmp)
rm(tmp)


#==============================#
#### VIF Variable Selection ####
#==============================#
# This is based on code by Mike Belitz (mbelitz/Odo_SDM_Rproj) and Shelly Gaynor
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

#================================#
#### MODIFIED KUENM FUNCTIONS ####
#================================#
#' Evaluation of candidate Maxent models during calibration
#'
#' @description kuenm_ceval evaluates candidate models in terms of statistical
#' significance (partial ROC), prediction ability (omission rates), and model complexity (AICc).
#' After evaluation, this function selects the best models based on user-defined criteria.
#'
#' @param path (character) directory in which folders containig calibration models are being created
#' or were created.
#' @param occ.joint (character) the name of csv file with training and testing occurrences combined;
#' columns must be: species, longitude, latitude.
#' @param occ.tra (character) the name of the csv file with the training occurrences;
#' columns as in occ.joint.
#' @param occ.test (character) the name of the csv file with the evaluation occurrences;
#' columns as in occ.joint.
#' @param batch (character) name of the batch file (bash for Unix) with the code to create all candidate models
#' for calibration.
#' @param out.eval (character) name of the folder where evaluation results will be written.
#' @param threshold (numeric) the percentage of training data omission error allowed (E); default = 5.
#' @param rand.percent (numeric) the percentage of data to be used for the bootstraping process
#' when calculating partial ROCs; default = 50.
#' @param iterations (numeric) the number of times that the bootstrap is going to be repeated;
#' default = 500.
#' @param kept (logical) if FALSE, all candidate models will be erased after evaluation, default = TRUE.
#' @param selection (character) model selection criterion, can be "OR_AICc", "AICc", or "OR";
#' OR = omission rates. Default = "OR_AICc", which means that among models that are statistically significant
#' and that present omission rates below the \code{threshold}, those with delta AICc up to 2 will be
#' selected. See details for other selection criteria.
#' @param parallel.proc (logical) if TRUE, pROC calculations will be performed in parallel using the available
#' cores of the computer. This will demand more RAM and almost full use of the CPU; hence, its use
#' is more recommended in high-performance computers. Using this option will speed up the analyses
#' only if models are large RasterLayers or if \code{iterations} are more than 5000. Default = FALSE.
#' @param silent (logical) if FALSE, report when evaluation is stalled waiting for a model run to
#' complete. Default = TRUE.
#'
#' @return A list with three dataframes containing results from the calibration process and a scatterplot
#' of all models based on the AICc values and omission rates. In addition, a folder, in the
#' working directory, containing a csv file with information about models meeting the user-defined
#' selection criterion, another csv file with a summary of the evaluation and selection process,
#' an extra csv file containing all the statistics of model performance (pROC, AICc, and omission
#' rates) for all candidate models, a png scatterplot of all models based on the AICc values and
#' rates, and an HTML file sumarizing all the information produced after evaluation for helping with
#' further interpretation.
#'
#' @details This function is used after or during the creation of Maxent candidate models for calibration.
#'
#' Other selecton criteria are described below:
#' If "AICc" criterion is chosen, all significant models with delta AICc up to 2 will be selected
#' If "OR" is chosen, the 10 first significant models with the lowest omission rates will be selected.
#'
#' @usage
#' kuenm_ceval(path, occ.joint, occ.tra, occ.test, batch, out.eval,
#'             threshold = 5, rand.percent = 50, iterations = 500,
#'             kept = TRUE, selection = "OR_AICc", parallel.proc = FALSE)
#'
#' @export
#'
#' @examples
#' # To run this function the kuenm_cal function needs te be used first. This previous function will
#' # create the models that kuenm_ceval evaluates.
#'
#' # Variables with information to be used as arguments.
#' occ_joint <- "aame_joint.csv"
#' occ_tra <- "aame_train.csv"
#' batch_cal <- "Candidate_models"
#' out_dir <- "Candidate_Models"
#' occ_test <- "aame_test.csv"
#' out_eval <- "Calibration_results"
#' threshold <- 5
#' rand_percent <- 50
#' iterations <- 100
#' kept <- TRUE
#' selection <- "OR_AICc"
#' paral_proc <- FALSE # make this true to perform pROC calculations in parallel
#'
#' cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
#'                         out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
#'                         kept = kept, selection = selection, parallel.proc = paral_proc)


kuenm_ceval <- function(path, occ.joint, occ.tra, occ.test, batch, out.eval, threshold = 5,
                        rand.percent = 50, iterations = 500, kept = TRUE,
                        selection = "OR_AICc", parallel.proc = FALSE, silent = TRUE) {
  #Checking potential issues
  if (missing(path)) {
    stop(paste("Argument path is not defined, this is necessary for reading the",
               "\ncandidate models created with the kuenm_cal function."))
  }
  if (!dir.exists(path)) {
    stop(paste(path, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (!file.exists(occ.joint)) {
    stop(paste(occ.joint, "does not exist in the working directory, check file name",
               "\nor extension, example: species_joint.csv"))
  }
  if (!file.exists(occ.tra)) {
    stop(paste(occ.tra, "does not exist in the working directory, check file name",
               "\nor extension, example: species_train.csv"))
  }
  if (!file.exists(occ.test)) {
    stop(paste(occ.test, "does not exist in the working directory, check file name",
               "\nor extension, example: species_test.csv"))
  }
  if (missing(out.eval)) {
    stop(paste("Argument out.eval is not defined, this is necessary for creating a folder",
               "\nwith the outputs of this function."))
  }
  if (missing(batch)) {
    stop(paste("Argument batch is not defined, this is necessary for evaluating",
               "\ncandidate models in the order in which they are created."))
  }
  
  
  #####
  #Slash
  if(.Platform$OS.type == "unix") {
    out <- "outputd.\\S*/"
  } else {
    out <- "outputd.\\S*\\\\"
  }
  
  #Data
  ##Source of initial information for model evaluation order
  if(.Platform$OS.type == "unix") {
    bat <- readLines(paste(batch, ".sh", sep = "")) #reading the batch file written to create the calibration models
    
  } else {
    bat <- readLines(paste(batch, ".bat", sep = "")) #reading the batch file written to create the calibration models
    
  }
  
  ###Recognizing the folders names and separating them for different procedures
  fol <- gregexpr("outputd.*\"", bat)
  fold <- regmatches(bat, fol)
  folde <- unlist(fold)
  
  ext <- gregexpr(out, bat)
  extr <- regmatches(bat, ext)
  extra <- unlist(extr)
  extract <- unique(extra)
  
  folder <- gsub(extract, "", folde, fixed = T) #names of all the calibration models folders
  
  folder_a <- gregexpr("M_.*all", folder)
  folder_al <- regmatches(folder, folder_a)
  folder_all <- unlist(folder_al) #folders with the models for calculating AICcs
  
  folder_c <- gregexpr("M_.*cal", folder)
  folder_ca <- regmatches(folder, folder_c)
  folder_cal <- unlist(folder_ca) #folder with the models for calculating pROCs and omission rates
  
  ##Models
  ###For AICc
  dir_names <- as.vector(paste(getwd(), "/", path, "/", folder_all, sep = "")) #vector of the subdirectories with the models
  
  ###For pROC and omission rates
  dir_names1 <- as.vector(paste(getwd(), "/", path, "/", folder_cal, sep = "")) #vector of the subdirectories with the models
  
  ###Names of the models to be evaluated
  mod_nam <- as.vector(gsub("_all", "", folder_all, fixed = TRUE)) #names of the models (taken from the folders names)
  
  ##Complete set and calibration and evaluation occurrences
  oc <- read.csv(occ.joint) #read all occurrences
  oc <- oc[, -1] #erase species name column
  
  occ <- read.csv(occ.tra) #read calibration occurrences
  occ <- occ[, -1] #erase species name column
  
  occ1 <- read.csv(occ.test) #read test occurrences
  occ1 <- occ1[, -1] #erase species name column
  
  #####
  #pROCs, omission rates, and AICcs calculation
  cat("\nPartial ROCs, omission rates, and AICcs calculation, please wait...\n")
  
  aiccs <- list() #empty list of AICc results
  proc_res <- list() #empty list of pROC values
  om_rates <- vector() #empty vector of omision rates
  
  if(.Platform$OS.type == "unix") {
    pb <- txtProgressBar(min = 0, max = length(dir_names), style = 3) #progress bar
  } else {
    pb <- winProgressBar(title = "Progress bar", min = 0, max = length(dir_names),
                         width = 300) #progress bar
  }
  options(list(show.error.messages = FALSE, suppressWarnings = TRUE))
  
  for(i in 1:length(dir_names)) {
    if (!silent) {
      message(sprintf("%d/%d : %s", i, length(dir_names), dir_names[i]))
    }
    Sys.sleep(0.1)
    if(.Platform$OS.type == "unix") {
      setTxtProgressBar(pb, i)
    } else {
      setWinProgressBar(pb, i, title = paste(round(i / length(dir_names) * 100, 2),
                                             "% of the evaluation process has finished"))
    }
    
    
    #AICc calculation
    lbds <- as.vector(list.files(dir_names[i], pattern = ".lambdas",
                                 full.names = TRUE)) #lambdas file
    lambdas <- try(readLines(lbds), silent = TRUE)
    
    par_num <- try(n.par(lambdas), silent = TRUE) #getting the number of parameters for each model
    
    mods <- list.files(dir_names[i], pattern = "*.asc$", full.names = TRUE) #name of ascii model
    mod <- try(raster::raster(mods), silent = TRUE)
    
    aicc <- try(kuenm_aicc(occ = oc, model = mod, npar = par_num), silent = TRUE)
    aiccs[[i]] <- aicc
    
    #If needed, waiting for the model to be created
    aicc_class <- class(aicc)
    
    waiting_for_aicc <- FALSE
    
    while (aicc_class == "try-error") {
      
      if ((!silent) && (!waiting_for_aicc)) {
        waiting_for_aicc <- TRUE
        message(sprintf("Waiting for aicc on %s\n", dir_names[i]))
        Sys.sleep(0.1)
      }
      
      lbds <- as.vector(list.files(dir_names[i], pattern = ".lambdas",
                                   full.names = TRUE)) #lambdas file
      lambdas <- try(readLines(lbds), silent = TRUE)
      
      par_num <- try(n.par(lambdas), silent = TRUE) #getting the number of parameters for each model
      
      mods <- list.files(dir_names[i], pattern = "*.asc$", full.names = TRUE) #name of ascii model
      mod <- try(raster::raster(mods), silent = TRUE)
      
      aicc <- try(kuenm_aicc(occ = oc, model = mod, npar = par_num), silent = TRUE)
      
      aiccs[[i]] <- aicc
      
      aicc_class <- class(aicc)
      if(aicc_class == "data.frame") {
        break()
      }
      
      #For avoiding infinite loops when models cannot be created
      mxlog <- as.vector(list.files(dir_names[i], pattern = ".log",
                                    full.names = TRUE)) #maxent log file
      llin <- try(readLines(mxlog), silent = TRUE)
      
      llin_class <- class(llin)
      
      waiting_for_llin <- FALSE
      
      while (llin_class == "try-error") {
        
        if ((!silent) && (!waiting_for_llin)) {
          waiting_for_llin <- TRUE
          message(sprintf("Waiting for llin on %s\n", dir_names[i]))
          Sys.sleep(0.1)
        }
        
        mxlog <- as.vector(list.files(dir_names[i], pattern = ".log",
                                      full.names = TRUE)) #maxent log file
        llin <- try(readLines(mxlog), silent = TRUE)
        llin_class <- class(llin)
        
        if(llin_class != "try-error") {
          loglin <- tolower(llin)
          
          e <- gregexpr("error", loglin)
          ee <- regmatches(loglin, e)
          eee <- unlist(ee)
          
          if(length(eee) > 0) {
            cat(paste("\nModel in folder", dir_names[i], "was not created because of an error.",
                      "\nCheck your files and software. NA values will be returned.\n"))
            
            ac <- data.frame(NA, NA, NA, NA)
            names(ac) <- c("AICc", "delta.AICc", "w.AIC", "nparam")
            aiccs[[i]] <- ac
          }
          
          break()
        }
      }
    }
    
    #pROCs and omission rates calculation
    mods1 <- list.files(dir_names1[i], pattern = "*.asc$", full.names = TRUE) #name of ascii model
    mod1 <- try(raster::raster(mods1), silent = TRUE)
    
    proc <- try(kuenm_proc(occ.test = occ1, model = mod1, threshold = threshold, # pRoc
                           rand.percent = rand.percent, iterations = iterations,
                           parallel = parallel.proc),
                silent = TRUE)
    
    proc_res[[i]] <- proc[[1]]
    
    
    
    om_rates[i] <- try(kuenm_omrat(model = mod1, threshold = threshold, # omission rates
                                   occ.tra = occ, occ.test = occ1), silent = TRUE)
    
    #If needed, waiting for the model to be created
    proc_class <- class(proc)
    
    waiting_for_proc <- FALSE
    
    while (proc_class == "try-error") {
      
      if ((!silent) && (!waiting_for_proc)) {
        waiting_for_proc <- TRUE
        message(sprintf("Waiting for proc on %s\n", dir_names[i]))
        Sys.sleep(0.1)
      }
      
      
      mods1 <- list.files(dir_names1[i], pattern = "*.asc$", full.names = TRUE) #name of ascii model
      mod1 <- try(raster::raster(mods1), silent = TRUE)
      
      proc <- try(kuenm_proc(occ.test = occ1, model = mod1, threshold = threshold, # pRoc
                             rand.percent = rand.percent, iterations = iterations,
                             parallel = parallel.proc),
                  silent = TRUE)
      
      proc_res[[i]] <- proc[[1]]
      
      om_rates[i] <- try(kuenm_omrat(model = mod1, threshold = threshold, # omission rates
                                     occ.tra = occ, occ.test = occ1), silent = TRUE)
      
      proc_class <- class(proc)
      if(proc_class == "list") {
        break()
      }
      
      #For avoiding infinite loops when models cannot be created
      mxlog <- as.vector(list.files(dir_names1[i], pattern = ".log",
                                    full.names = TRUE)) #maxent log file
      llin <- try(readLines(mxlog), silent = TRUE)
      llin_class <- class(llin)
      
      while (llin_class == "try-error") {
        mxlog <- as.vector(list.files(dir_names1[i], pattern = ".log",
                                      full.names = TRUE)) #maxent log file
        llin <- try(readLines(mxlog), silent = TRUE)
        llin_class <- class(llin)
        
        if(llin_class != "try-error") {
          loglin <- tolower(llin)
          
          e <- gregexpr("error", loglin)
          ee <- regmatches(loglin, e)
          eee <- unlist(ee)
          
          if(length(eee) > 0) {
            cat(paste("\nModel in folder", dir_names1[i], "was not created because of an error.",
                      "\nCheck your files and software. NA values will be returned.\n"))
            
            p_roc <- rep(NA, 2)
            names(p_roc) <- c(paste("Mean_AUC_ratio_at_", threshold, "%", sep = ""), "pval_pROC")
            auc_ratios <- rep(NA, 4)
            names(auc_ratios) <- c("Iteration", paste("AUC_at_", 100 - threshold, "%", sep = ""),
                                   paste("AUC_at_", threshold, "%", sep = ""), "AUC_ratio")
            proc_res[[i]] <- list(p_roc, auc_ratios)
            
            or <- NA
            names(or) <- "om_rate_5%"
            om_rates[i] <- or
            
            break()
          }
        }
      }
    }
    
    #Erasing calibration models after evaluating them if kept = FALSE
    if(kept == FALSE) {
      unlink(dir_names[i], recursive = T)
      unlink(dir_names1[i], recursive = T)
    }
  }
  if(.Platform$OS.type != "unix") {
    suppressMessages(close(pb))
  }
  n.mod <- i
  
  ##Erasing main folder of candidate models if kept = FALSE
  if(kept == FALSE) {
    unlink(path, recursive = T)
    cat("\nAll candidate models were deleted\n")
  }else{
    cat("\nAll candidate models were kept\n")
  }
  
  ##Creating the final tables
  ###From AICc analyses
  aiccs <- do.call(rbind, aiccs) #joining tables
  
  for (i in 1:length(aiccs[, 1])) {
    aiccs[i, 2] <- (aiccs[i, 1] - min(aiccs[, 1], na.rm = TRUE))
    aiccs[i, 3] <- (exp(-0.5 * aiccs[i,2])) / (sum(exp(-0.5 * aiccs[,2]), na.rm = TRUE))
  }
  
  ###From pROC analyses
  proc_res1 <- do.call(rbind, proc_res) #joining tables of the pROC results
  proc_res_m <- data.frame(mod_nam, proc_res1) #adding a new column with the number of AUC ratios interations < 1
  
  #####
  om_rates <- as.numeric(as.character(om_rates))
  #Joining the results
  ku_enm_eval <- data.frame(proc_res_m, om_rates, aiccs)
  colnames(ku_enm_eval) <- c("Model", "Mean_AUC_ratio", "pval_pROC",#changing column names in the final table
                             paste("Omission_rate_at_", threshold, "%", sep = ""), "AICc",
                             "delta_AICc", "W_AICc", "num_parameters")
  
  #####
  #Choosing the best models
  cat("\nSelecting the best candidate models...\n")
  
  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
    if(selection == "OR_AICc") {
      ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
      ku_enm_best <- na.omit(ku_enm_bes[which(ku_enm_bes[, 4] <= threshold / 100), ])
      if(length(ku_enm_best[, 4]) != 0) {
        for (i in 1:length(ku_enm_best[,1])) {
          ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
          ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE))
        }
        ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
        
      }else {
        cat(paste("\nNone of the significant candidate models met the omission rate criterion,",
                  "\nmodels with the smallest omission rate and lowest AICc will be presented\n"))
        
        ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ]
        for (i in 1:length(ku_enm_best[, 1])) {
          ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
          ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[i, 6]), na.rm = TRUE))
        }
        ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
      }
    }
    
    if(selection == "AICc") {
      ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
      ku_enm_best <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]
      if(length(ku_enm_best[, 6]) != 0) {
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
      }else {
        cat(paste("\nNone of the significant candidate models met the AICc criterion,",
                  "\ndelta AICc will be recalculated for significant models\n"))
        
        for (i in 1:length(ku_enm_best[, 6])) {
          ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
          ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE))
        }
        ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
      }
    }
    
    if(selection == "OR") {
      ku_enm_b <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
      ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
      ku_enm_bes1 <- ku_enm_b[ku_enm_b[, 3] <= 0.05, ]
      ku_enm_best <- ku_enm_bes1[ku_enm_bes1[, 4] <= threshold / 100, ]
      if(length(ku_enm_best[, 4]) != 0) {
        if(length(ku_enm_best[, 4]) > 10) {
          ku_enm_best <- ku_enm_best[order(ku_enm_best[, 4]), ][1:10, ]
        }else {
          ku_enm_best <- ku_enm_best[order(ku_enm_best[, 4]), ]
        }
      }else {
        cat(paste("\nNone of the significant candidate models met the omission rate criterion,",
                  "\nmodels with the smallest omission rate will be presented\n"))
        
        ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ][1:10, ]
      }
    }
  }else {
    cat("\nNo valid model selection criterion has been defined,\n
        no file containing best models was created.\n
        Select your best models from the complete list.\n")
  }
  
  #####
  #Statistics of the process
  ##Counting
  ku_enm_sign <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
  ku_enm_sign <- ku_enm_sign[ku_enm_sign[, 3] <= 0.05, ]
  
  ku_enm_or <- ku_enm_eval[ku_enm_eval[, 4] <= threshold / 100, ]
  
  ku_enm_AICc <- ku_enm_eval[!is.na(ku_enm_eval[, 6]), ]
  ku_enm_AICc <- ku_enm_AICc[ku_enm_AICc[, 6] <= 2, ]
  
  ku_enm_best_OR <- ku_enm_sign[ku_enm_sign[, 4] <= threshold / 100, ]
  
  ku_enm_best_AICc <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]
  
  ku_enm_best_OR_AICc <- ku_enm_bes[ku_enm_bes[, 4] <= threshold / 100, ]
  if(length(ku_enm_best_OR_AICc[, 4]) != 0) {
    for (i in 1:length(ku_enm_best_OR_AICc[, 1])) {
      ku_enm_best_OR_AICc[i, 6] <- (ku_enm_best_OR_AICc[i, 5] - min(ku_enm_best_OR_AICc[, 5],
                                                                    na.rm = TRUE))
      ku_enm_best_OR_AICc[i, 7] <- (exp(-0.5 * ku_enm_best_OR_AICc[i, 6])) / (sum(exp(-0.5 * ku_enm_best_OR_AICc[, 6]),
                                                                                  na.rm = TRUE))
    }
    ku_enm_best_OR_AICc <- ku_enm_best_OR_AICc[ku_enm_best_OR_AICc[, 6] <= 2, ]
  }
  
  ##Preparing the table
  r_names <- c("All candidate models", "Statistically significant models", "Models meeting omission rate criteria",
               "Models meeting AICc criteria", "Statistically significant models meeting omission rate criteria",
               "Statistically significant models meeting AICc criteria",
               "Statistically significant models meeting omission rate and AICc criteria")
  statis <- c(length(ku_enm_eval[, 1]),
              length(ku_enm_sign[, 3]),
              length(ku_enm_or[, 4]),
              length(ku_enm_AICc[, 6]),
              length(ku_enm_best_OR[, 4]),
              length(ku_enm_best_AICc[, 6]),
              length(ku_enm_best_OR_AICc[, 2]))
  
  ku_enm_stats <- data.frame(r_names, statis)
  colnames(ku_enm_stats) <- c("Criteria", "Number of models")
  
  #####
  #Writing the results
  ##csv files
  cat("\nWriting kuenm_ceval results...\n")
  dir.create(out.eval)
  
  name <- paste(out.eval, "calibration_results.csv", sep = "/")
  name0 <- paste(out.eval, "calibration_stats.csv", sep = "/")
  name1 <- paste(out.eval, "best_candidate_models_OR_AICc.csv", sep = "/")
  name2 <- paste(out.eval, "best_candidate_models_AICc.csv", sep = "/")
  name3 <- paste(out.eval, "best_candidate_models_OR.csv", sep = "/")
  
  
  write.csv(ku_enm_eval, file = name, eol = "\n", na = "NA", row.names = FALSE)
  write.csv(ku_enm_stats, file = name0, eol = "\n", na = "NA", row.names = FALSE)
  
  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR"){
    if(selection == "OR_AICc"){
      write.csv(ku_enm_best, file = name1, eol = "\n", na = "NA", row.names = FALSE)
    }
    if(selection == "AICc"){
      write.csv(ku_enm_best, file = name2, eol = "\n", na = "NA", row.names = FALSE)
    }
    if(selection == "OR"){
      write.csv(ku_enm_best, file = name3, eol = "\n", na = "NA", row.names = FALSE)
    }
  }
  
  ##Plot
  png(paste(out.eval, "calibration_figure.png", sep = "/"), width = 80, height = 80,
      units = "mm", res = 600)
  
  par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.58)
  dat <- na.omit(ku_enm_eval)
  if (nrow(dat) > 0) {
    plot(na.omit(ku_enm_eval)[, 4]~log(na.omit(ku_enm_eval)[, 5]),
         xlab = "Natural logarithm of AICc", ylab = paste("Omission rates at",
                                                          paste(threshold, "%", sep = ""),
                                                          "threshold value", sep = " "),
         las = 1, col = "grey35")
    
    points(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 4]~log(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 5]),
           col = "red1", pch = 19, cex = 1.1)
    
    if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
      if(selection == "OR_AICc") {
        points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
               col = "dodgerblue1", pch = 17, cex = 1.4)
        legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
               pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("dodgerblue1", "red1", "gray35"), bty = "n",
               inset = c(0.01, 0))
      }
      if(selection == "AICc") {
        points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
               col = "darkorchid1", pch = 17, cex = 1.4)
        legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
               pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("darkorchid1", "red1", "gray35"), bty = "n",
               inset = c(0.01, 0))
      }
      if(selection == "OR") {
        points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
               col = "orange2", pch = 17, cex = 1.4)
        legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
               pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("orange2", "red1", "gray35"), bty = "n",
               inset = c(0.01, 0))
      }
    }
    
    dev.off()
  }
  
  ##html file
  ###Writing the html file
  html_eval(path = out.eval, file.name = "calibration_results")
  
  ##Retuning objects
  ###Dataframes in a list
  list_res <- list(ku_enm_stats, ku_enm_best, ku_enm_eval)
  names(list_res) <- c("Summary", "Best models",
                       "All models")
  
  
  if (nrow(dat) > 0) {
    ###Plot
    ku_enm_plot <- {
      par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.85)
      plot(na.omit(ku_enm_eval)[, 4]~log(na.omit(ku_enm_eval)[, 5]),
           xlab = "Natural logarithm of AICc", ylab = paste("Omission rates at",
                                                            paste(threshold, "%", sep = ""),
                                                            "threshold value", sep = " "),
           las = 1, col = "grey35")
      
      points(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 4]~log(na.omit(ku_enm_eval[!ku_enm_eval[, 1] %in% ku_enm_bes[, 1], ])[, 5]),
             col = "red1", pch = 19, cex = 1.1)
      
      if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
        if(selection == "OR_AICc") {
          points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
                 col = "dodgerblue1", pch = 17, cex = 1.4)
          legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
                 pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("dodgerblue1", "red1", "gray35"), bty = "n",
                 inset = c(0.01, 0))
        }
        if(selection == "AICc") {
          points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
                 col = "darkorchid1", pch = 17, cex = 1.4)
          legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
                 pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("darkorchid1", "red1", "gray35"), bty = "n",
                 inset = c(0.01, 0))
        }
        if(selection == "OR") {
          points(na.omit(ku_enm_best)[, 4]~log(na.omit(ku_enm_best)[, 5]),
                 col = "orange2", pch = 17, cex = 1.4)
          legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
                 pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), col = c("orange2", "red1", "gray35"), bty = "n",
                 inset = c(0.01, 0))
        }
      }
    }
  }
  
  #####
  #Finalizing the function
  cat("\nProcess finished\n")
  cat(paste("A folder containing results of the calibration of", n.mod,
            "\ncandidate models has been written\n", sep = " "))
  
  cat(paste("\nThe folder", out.eval, "contains:\n", sep = " "))
  cat("   -A html file and its dependencies that sum all the results, check\n")
  cat(paste("    ", "calibration_results.html\n", sep = ""))
  
  cat("   -Two csv files with all models' calibration results and stats,\n")
  
  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR") {
    if(selection == "OR_AICc"){
      cat("    and an aditional csv file containing the best models selected by OR and AICc.\n")
    }
    if(selection == "AICc") {
      cat("    and an aditional csv file containing the best models selected by AICc.\n")
    }
    if(selection == "OR") {
      cat("    and an aditional csv file containing the best models selected by OR.\n")
    }
  }
  
  cat(paste("\nCheck your working directory!!!", getwd(), "\n", sep = "    "))
  
  options(list(show.error.messages = TRUE, suppressWarnings = FALSE))
  
  return(list_res)
}

#' Helper function to calculate the AICc values (number of parameters).
#'
#' @param x An object derived from reading the lambdas file created for Maxent.
#' Use \code{\link[base]{readLines}} function to read the file.
#'
#' @export

n.par <- function(x) {
  lambdas <- x[1:(length(x) - 4)]
  countNonZeroParams <- function(x) {
    if (strsplit(x, split = ", ")[[1]][2] != "0.0")
      1
  }
  no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(no.params)
}

kuenm_ceval2<-function(out.eval="Calibration_results_topvars", set="topvars", threshold=5){
  cal_res<-read.csv("Calibration_results/calibration_results.csv")
  ku_enm_eval=cal_res[str_detect(cal_res$Model,set),]
  ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
  ku_enm_best <- na.omit(ku_enm_bes[which(ku_enm_bes[, 4] <= threshold / 100), ])
  if(length(ku_enm_best[, 4]) != 0) {
    for (i in 1:length(ku_enm_best[,1])) {
      ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
      ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE))
    }
    ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
    ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
    
  }else {
    cat(paste("\nNone of the significant candidate models met the omission rate criterion,",
              "\nmodels with the smallest omission rate and lowest AICc will be presented\n"))
    
    ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ]
    for (i in 1:length(ku_enm_best[, 1])) {
      ku_enm_best[i, 6] <- (ku_enm_best[i, 5] - min(ku_enm_best[, 5], na.rm = TRUE))
      ku_enm_best[i, 7] <- (exp(-0.5 * ku_enm_best[i, 6])) / (sum(exp(-0.5 * ku_enm_best[i, 6]), na.rm = TRUE))
    }
    ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
    ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
  }
  ku_enm_sign <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
  ku_enm_sign <- ku_enm_sign[ku_enm_sign[, 3] <= 0.05, ]
  
  ku_enm_or <- ku_enm_eval[ku_enm_eval[, 4] <= threshold / 100, ]
  
  ku_enm_AICc <- ku_enm_eval[!is.na(ku_enm_eval[, 6]), ]
  ku_enm_AICc <- ku_enm_AICc[ku_enm_AICc[, 6] <= 2, ]
  
  ku_enm_best_OR <- ku_enm_sign[ku_enm_sign[, 4] <= threshold / 100, ]
  
  ku_enm_best_AICc <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]
  
  ku_enm_best_OR_AICc <- ku_enm_bes[ku_enm_bes[, 4] <= threshold / 100, ]
  if(length(ku_enm_best_OR_AICc[, 4]) != 0) {
    for (i in 1:length(ku_enm_best_OR_AICc[, 1])) {
      ku_enm_best_OR_AICc[i, 6] <- (ku_enm_best_OR_AICc[i, 5] - min(ku_enm_best_OR_AICc[, 5],
                                                                    na.rm = TRUE))
      ku_enm_best_OR_AICc[i, 7] <- (exp(-0.5 * ku_enm_best_OR_AICc[i, 6])) / (sum(exp(-0.5 * ku_enm_best_OR_AICc[, 6]),
                                                                                  na.rm = TRUE))
    }
    ku_enm_best_OR_AICc <- ku_enm_best_OR_AICc[ku_enm_best_OR_AICc[, 6] <= 2, ]
  }
  
  ##Preparing the table
  r_names <- c("All candidate models", "Statistically significant models", "Models meeting omission rate criteria",
               "Models meeting AICc criteria", "Statistically significant models meeting omission rate criteria",
               "Statistically significant models meeting AICc criteria",
               "Statistically significant models meeting omission rate and AICc criteria")
  statis <- c(length(ku_enm_eval[, 1]),
              length(ku_enm_sign[, 3]),
              length(ku_enm_or[, 4]),
              length(ku_enm_AICc[, 6]),
              length(ku_enm_best_OR[, 4]),
              length(ku_enm_best_AICc[, 6]),
              length(ku_enm_best_OR_AICc[, 2]))
  
  ku_enm_stats <- data.frame(r_names, statis)
  colnames(ku_enm_stats) <- c("Criteria", "Number of models")
  
  #####
  #Writing the results
  ##csv files
  cat("\nWriting kuenm_ceval results...\n")
  dir.create(out.eval)
  name <- paste(out.eval, "calibration_results.csv", sep = "/")
  name0 <- paste(out.eval, "calibration_stats.csv", sep = "/")
  name1 <- paste(out.eval, "best_candidate_models_OR_AICc.csv", sep = "/")
  write.csv(ku_enm_eval, file = name, eol = "\n", na = "NA", row.names = FALSE)
  write.csv(ku_enm_stats, file = name0, eol = "\n", na = "NA", row.names = FALSE)
  write.csv(ku_enm_best, file = name1, eol = "\n", na = "NA", row.names = FALSE)
}

kuenm_aicc <- function (occ, model, npar) {
  if (missing(occ)) {
    stop("Argument occ must be defined, see function's help.")
  }
  if (missing(model)) {
    stop("Argument model must be defined, see function's help.")
  }
  if (!class(model)[1] %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
    stop("model must be a RasterLayer or RasterStack object. See function's help.")
  }
  if (missing(npar)) {
    stop("Argument npar must be defined, see function's help.")
  }
  if (dim(model)[3] != length(npar)) {
    stop("Number of models to evaluate must correspond with length of npar vector.")
  }
  
  AIC.valid <- npar < nrow(occ)
  if (dim(model)[3] == 0) {
    res <- data.frame(cbind(AICc = NA, delta.AICc = NA,
                            w.AIC = NA, parameters = npar))
    warning("Cannot calculate AICc when model = FALSE... returning NA's.")
  } else {
    vals <- raster::extract(model, occ)
    probsum <- sum(raster::values(model), na.rm = TRUE)
    lval<-log(t(t(vals)/probsum))
    lval[which(is.infinite(lval))]<-NA
    LL <- colSums(lval, na.rm = TRUE)
    AICc <- (2 * npar - 2 * LL) + (2 * (npar) * (npar + 1)/(nrow(occ) - npar - 1))
    AICc[AIC.valid == FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if (sum(is.na(AICc)) == length(AICc)) {
      warning("AICc not valid: too many parameters, or likelihood = Inf... returning NA.")
      res <- data.frame(cbind(AICc, delta.AICc = NA, w.AIC = NA,
                              parameters = npar))
    } else {
      delta.AICc <- (AICc - min(AICc, na.rm = TRUE))
      w.AIC <- (exp(-0.5 * delta.AICc))/(sum(exp(-0.5 *
                                                   delta.AICc), na.rm = TRUE))
      res <- data.frame(AICc, delta.AICc, w.AIC, parameters = npar)
      rownames(res) <- NULL
    }
  }
  rownames(res) <- NULL
  return(res)
}

#=========================#
#### Threshold Rasters ####
#=========================#
# This is based on code by Cecina Babich Morrow
# https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/#:~:text=present%20(P10).-,Minimum%20training%20presence,suitability%20value%20for%20the%20species.
raster_threshold <- function(input_raster, points = NULL, type = NULL, threshold = NULL, binary = FALSE) {
  if (!is.null(points)) {
    pointVals <- raster::extract(input_raster, points)
    if (type == "mtp") {
      threshold <- min(na.omit(pointVals))
    } else if (type == "p10") {
      if (length(pointVals) < 10) {
        p10 <- floor(length(pointVals) * 0.9)
      } else {
        p10 <- ceiling(length(pointVals) * 0.9)
      }
      threshold <- rev(sort(pointVals))[p10]
    }
  }
  raster_thresh <- input_raster
  raster_thresh[raster_thresh < threshold] <- NA
  if (binary) {
    raster_thresh[raster_thresh >= threshold] <- 1
  }
  return(raster_thresh)
}

#========================#
#### Summarize Models ####
#========================#
summarizeModels<-function(path, occ=occ2){
  set=str_split(path,"_",simplify=T)[3]
  if(is.na(set)){
    set="all"
    set2=""
  }else{
    set2=paste0("_",set)
  }
  
  # Collate maxentResults
  print("Collating maxentResults")
  files<-list.files(Sys.glob(paste0(path,"/*")), pattern="maxentResults.csv", full.names = T)
  results<-read.csv(files[1])
  results<-results[nrow(results),]
  results$Model<-str_split(files[1],"/",simplify = T)[2]
  results$Set<-set
  if(length(files)>1){
    for(i in 2:length(files)){
      tmp<-read.csv(files[i])
      tmp<-tmp[nrow(tmp),]
      tmp$Model<-str_split(files[i],"/",simplify = T)[2]
      tmp$Set<-set
      results<-merge(results, tmp, all=T)
    } 
  }
  
  # Merge Evaluation Results
  print("Merging with Evaluation Results")
  files<-list.files(Sys.glob(paste0(path,"_evaluation")), pattern="fm_evaluation_results.csv", full.names = T)
  results2<-read.csv(files[1])
  results3<-merge(results2,results,all=T)
  
  # Merge Calibration Results
  print("Merging with Calibration Results")
  files<-list.files(Sys.glob(paste0("Calibration_results",set2)), pattern="best_candidate_models_OR_AICc.csv", full.names = T)
  results4<-read.csv(files[1])
  colnames(results4)[2:ncol(results4)]<-paste0("Cal_",colnames(results4)[2:ncol(results4)])
  results5<-merge(results4,results3,all=T)
  results5<-results5 %>% dplyr::select(Species, Model, Set, Entropy, everything(.))
  write.csv(results5, paste0(path,"/FinalModelResults.csv"), row.names = F, quote=F)
  
  # Find Best Models
  results6=results5 %>% filter(!is.na(Omission_rate_at_5.) & !is.na(Mean_AUC_ratio)) %>% 
    filter(Omission_rate_at_5. == min(Omission_rate_at_5.)) %>% 
    filter(Mean_AUC_ratio==max(Mean_AUC_ratio))
  write.csv(results6, paste0(path,"/FinalModelResults_best.csv"), row.names = F, quote=F)
  
  best_models<-results6$Model
  
  # Calculate mean and export to combined files
  results_mean<-results6 %>% group_by(Species) %>% summarize_if(is.numeric, mean, na.rm=T)
  write.csv(results_mean, paste0(path,"/FinalModelResults_mean.csv"), row.names = F, quote=F)
  
  # Get mean raster
  print("Calculating Mean Raster")
  files<-list.files(Sys.glob(paste0(path,"/",best_models)), pattern="*_avg.asc", full.names = T)
  rasters<-stack(files)
  crs(rasters)<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  rasters_mean<-mean(rasters, na.rm=T)
  new_raster_name<-paste0(path,"/",gsub(".asc", ".tif",gsub(".*\\/","", files[1])))
  writeRaster(rasters_mean, new_raster_name, format="GTiff")
  
  # Convert raw mean raster to cloglog
  print("Converting Raw Output to Cloglog Format")
  h<-as.numeric(results_mean$Entropy)
  c<-raw2clog(rasters_mean, h)
  c<-raw2clog(rasters_mean, h)
  r_prefix<-gsub(".tif", "", new_raster_name)
  writeRaster(c, paste0(r_prefix,"_cloglog.tif"), format="GTiff")
  
  # P10 Threshold
  print("Using 10th Percentile Training Presence to Threshold")
  sdm_p10<-raster_threshold(c, occ[,2:3], type="p10")
  sdm_p10b<-raster_threshold(c, occ[,2:3], type="p10", binary = T)
  writeRaster(sdm_p10, paste0(r_prefix, "_cloglog_p10.tif"), bylayer=T, format="GTiff")
  writeRaster(sdm_p10b, paste0(r_prefix, "_cloglog_p10b.tif"), bylayer=T, format="GTiff")
  
}



