# autokuenm
## Rhett M. Rautsaw

Ecological Niche Models in one easy script. autokuenm.R is designed to take occurrence records, distribution maps, and environmental data and prepare it for kuenm.

Specifically, this script will:
1. Subset occurrence records by a species name
2. Subset distribution maps by a species name (or create convex hull from points)
3. Buffer the resulting distribution and crop environmental data by the buffered extent
4. Thin occurrence records to one point per environmental grid cell
5. Generate random background points and partition occurrence records for training and testing (`ENMeval`)
6. Create exhaustive combination of all environmental variables or perform heuristic search on sets
    - Heuristic search involves identifying non-colinear variables with high `maxent` permutation importance
7. Run `kuenm`

# Installation

autokuenm requires the following packages in R as well as gnu-parallel. 
```
packages <- c("argparse","readr","dplyr","stringr","parallel","memuse",
              "sf","raster","dismo","rJava","virtualspecies",
              "usdm","ENMeval","kuenm")
```

Many of these packages can be installed via anaconda to create an isolated autokuenm environment.
```
conda create -n r_spatial r-argparse r-readr r-dplyr r-stringr r-memuse r-sf r-raster r-dismo r-devtools r-rgdal r-rjava parallel
conda activate r_spatial
```

`ENMeval`, `virtualspecies`, `usdm`, and `kuenm` will have to be installed manually as they are not available in anaconda.
```
# In R
install.packages("ENMeval")
install.packages("virtualspecies")
install.packages("usdm")
devtools::install_github("marlonecobos/kuenm")
```

# Arguments
| Flag      | Default                | Description                                                                                                                                                                                                                                                                                      |
|-----------|------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -o        | occ.csv                | csv of occurrence records. Can contain multiple species and will be filtered by -sp and -spcol1 flags. Latitude and Longitude columns can be set with -latcol and -loncol                                                                                                                        |
| -e        | env.tab                | For exhaustive search, provide a folder containing all environmental rasters. For heuristic search, provide a text file containing paths to all environmental rasters divided into Sets (Must be in .tif or .asc format).                                                                        |
| -d        | distributions.geojson  | Distribution shapefile (e.g., GeoJSON or SHP). Alternatively you can specify "mcp" here if you do not have a distribution for your species. In this case, a convex hull will be created from the occurrence records. Can contain multiple species and will be filtered by -sp and -spcol2 flags. |
| -sp       | Agkistrodon_bilineatus | Name of species to model.                                                                                                                                                                                                                                                                        |
| -spcol1   | Species                | Name of column with species ID in occurrence records csv.                                                                                                                                                                                                                                        |
| -spcol2   | Species                | Name of column with species ID in distribution attribute table.                                                                                                                                                                                                                                  |
| -latcol   | latitude               | Name of column with latitude in occurrence records csv                                                                                                                                                                                                                                           |
| -loncol   | longitude              | Name of column with longitude in occurrence records csv                                                                                                                                                                                                                                          |
| -b        | 100                    | Buffer for distribution in km.                                                                                                                                                                                                                                                                   |
| -bg       | 2000                   | Number of random background points to generate for assigning blocks                                                                                                                                                                                                                              |
| -g        | NA                     | Directory containing all environmental rasters for projection in kuenm. Do not include if you do not want to project to larger or forward/backward in time.                                                                                                                                      |
| -c        | `system cores`         | Number of processor threads to use. Matching this number to the number of cores on your computer speeds up some operations.                                                                                                                                                                      |
| -m        | `system memory`        | Maximum memory (in megabytes) to be used by maxent while creating the models.                                                                                                                                                                                                                    |


The script is run from command line like so:
```
autokuenm.R -o occ.csv -e env.tab -d combined_distributions.geojson 
    -sp Agkistrodon_bilineatus -spcol1 Species -spcol2 Species 
    -latcol latitude -loncol longitude -b 100 -bg 2000 -g G_variables 
    -c 8 -m 62000
```

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