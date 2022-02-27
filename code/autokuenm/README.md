# autokuenm
## Rhett M. Rautsaw

Species Distribution Modeling in one easy script. autokuenm.R is designed to take occurrence records and environmental data and prepare it for kuenm.

Specifically, this script will:
1. Subset occurrence records by a species name
2. Subset distribution maps by a species name (or create hull from points)
3. Overlay the distribution/hull with bioregions to define M areas
4. Partition occurrence records for training and testing (ENMeval)
5. Create exhaustive combination of all environmental variables or perform heuristic search on sets
	- Heuristic search involves identifying non-colinear variables with high permutation importance
6. Run kuenm
7. Summarize results

# Installation

autokuenm requires the following packages in R as well as gnu-parallel. 
```
packages <- c("argparse","readr","dplyr","stringr","tidyr","parallel","memuse",
			"sf","raster","dismo","rJava","virtualspecies", "concaveman", "spThin",
			"rnaturalearth", "rnaturalearthdata", "usdm","ENMeval","kuenm")
```

Many of these packages can be installed via anaconda to create an isolated autokuenm environment.
```
conda create -n r_spatial r-argparse r-readr r-dplyr r-stringr r-tidyr r-memuse r-sf r-raster r-dismo r-devtools r-rgdal r-rjava r-concaveman r-rnaturalearth r-rnaturalearthdata r-geos r-lwgeom parallel
conda activate r_spatial
```

`ENMeval`, `virtualspecies`, `usdm`,`spThin`, and `kuenm` will have to be installed manually as they are not available in anaconda.
```
# In R
install.packages("ENMeval")
install.packages("virtualspecies")
install.packages("usdm")
install.packages("spThin")
devtools::install_github("marlonecobos/kuenm")
devtools::install_github("ropensci/rnaturalearthhires")
```

# Arguments
| Flag		| Default					| Description																																																																|
|-----------|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -o		| occ.csv					| csv of occurrence records. Can contain multiple species and will be filtered by -sp and -spcol1 flags. Latitude and Longitude columns can be set with -latcol and -loncol																									|
| -e		| env.tab					| For exhaustive search, provide a folder containing all environmental rasters. For heuristic search, provide a text file containing paths to all environmental rasters divided into Sets (Must be in .tif or .asc format).													|
| -d		| ahull						| Distribution shapefile (e.g., GeoJSON or SHP).  If no shapefile is provided, then you can use a minimum convex polygon or concave alpha hull by specifiying 'mcp' or 'ahull'. Distribution file can contain multiple species and will be filtered by -sp and -spcol flags.|
| -sp		| Agkistrodon_bilineatus	| Name of species to model.																																																													|
| -spcol1	| Species					| Name of column with species ID in occurrence records csv.																																																					|
| -spcol2	| Species					| Name of column with species ID in distribution attribute table.																																																			|
| -latcol	| latitude					| Name of column with latitude in occurrence records csv																																																					|
| -loncol	| longitude					| Name of column with longitude in occurrence records csv																																																					|
| -b		| wwf_terr_ecos.shp			| Bioregions used to assign M areas.																																																										|
| -bg		| 2000						| Number of random background points to generate for assigning training/testing partitions																																													|
| -bias		| SpatialBias_resamp.tif	| Use a spatial bias file to factor out bias within MaxEnt. If empty, a distance filter is applied to occurrence records.																																					|
| -novarsel	| FALSE						| Turn off variable selection if you want to use all environmental variables in -e list.																																													|
| -g		| NA						| Directory containing all environmental rasters for projection in kuenm. Do not include if you do not want to project to larger or forward/backward in time.																												|
| -c		| `system cores`			| Number of processor threads to use. Matching this number to the number of cores on your computer speeds up some operations.																																				|
| -m		| `system memory`			| Maximum memory (in megabytes) to be used by maxent while creating the models.																																																|


The script is run from command line like so:
```
autokuenm.R -o occ.csv -e env.tab -d ahull 
    -sp Agkistrodon_bilineatus -spcol1 Species -spcol2 Species 
    -latcol latitude -loncol longitude -b 100 -bg 2000 
    -bias bias_file.tif -g G_variables 
    -c 8 -m 62000
```

# Other Functions/Resources
In addition to `autokuenm` several additional functions are provided in `functions.R`. 

- Format Converters: A series of functions to convert raw MaxEnt output into logistic and cloglog formats.
	- `raw2log`
	- `raw2clog`
	- `log2raw`
- `ahull`: A function to use [`concaveman`]() to calculate concave alpha hulls
- `Marea`: A function to create M-areas from a distribution/hull and biogeographic regions
- `SelectVariables`: A function to select the top contributing variables and remove collinearity through iterative examination of MaxEnt permutation importance and VIF.
- `raster_threshold`: A function to perform raster thresholding based on a given value or via standard thesholding methodologies including minimum training presence and 10th percentile training presence.
	- Function by [Cecina Babich Morrow](https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/#:~:text=present%20(P10).-,Minimum%20training%20presence,suitability%20value%20for%20the%20species)
- `summarizeModels`: A function to summarize final kuenm results and identify the model with the lowest omission rate and highest AUC ratio. 

# Troubleshooting

## Common Error Messages
### Couldn't get file lock
```
Jun 28, 2021 5:04:42 PM java.util.prefs.FileSystemPreferences syncWorld
WARNING: Couldn't flush user prefs: java.util.prefs.BackingStoreException: Couldn't get file lock.
Jun 28, 2021 5:05:12 PM java.util.prefs.FileSystemPreferences syncWorld
WARNING: Couldn't flush user prefs: java.util.prefs.BackingStoreException: Couldn't get file lock.
```
This error message sometimes displays itself while maxent is running. Usually over and over again. Don't worry, maxent is still running properly. You can make this error message go away be removing the lock file in `~/.java/.userPrefs/.user.lock.*`

### Error reading file *.mxe
```
Error: Error reading file /PATH/TO/M_variables/Set_*/maxent.cache/*.mxe
```
This error message sometimes displays itself when running the calibration kuenm models in parallel. Maxent creates cache files (.mxe, .info) for faster access of the environmental layers during subsequent model runs. When running in parallel, multiple maxent runs will try to create these cache files simultaneously, causing one of the runs to fail. But have no fear, because I have set parallel to retry failed models 3 times. It should work the second time around since the cache file will have already been made by another model run. 