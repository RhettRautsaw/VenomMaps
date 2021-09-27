# VenomMaps: Updated Distribution Maps and Niche Models for New World Pitvipers (Crotalinae)
<img align="right" src="www/VenomMaps.png" width=150>

### Rhett M. Rautsaw, Gustavo Jiménez, Christoph Grünwald, Erich P. Hofmann, Laura Alencar, Marcio Martins, Paola Carrasco, Tiffany Doan, & Christopher L. Parkinson

# Table of Contents

- Cleaned occurrence records for Viperidae from [GBIF](https://www.gbif.org/), [HerpMapper](https://www.herpmapper.org/), [Bison](https://bison.usgs.gov/), [Brazil Snake Atlas](https://bioone.org/journals/South-American-Journal-of-Herpetology/volume-14/issue-sp1/SAJH-D-19-00120.1/Atlas-of-Brazilian-Snakes--Verified-Point-Locality-Maps-to/10.2994/SAJH-D-19-00120.1.short), [BioWeb Ecuador](https://bioweb.bio/), and custom databases/georeferencing. 
  - [`data/occurrence`](https://github.com/RhettRautsaw/VenomMaps/tree/master/data/occurrence)
- Updated distribution maps in `geojson` format for all New World Pitvipers as well as distribution maps for Old World Vipers from [Roll et al. 2017](https://www.nature.com/articles/s41559-017-0332-2).
  - [`data/distributions`](https://github.com/RhettRautsaw/VenomMaps/tree/master/data/distributions)
- Final niche models for all New World pitvipers
  - [`data/enms`](https://github.com/RhettRautsaw/VenomMaps/tree/master/data/enms)
- Code used to clean occurrence records, construct distribution maps, and summarize enm results.
  - [`code/*`](https://github.com/RhettRautsaw/VenomMaps/tree/master/code)
- Code used to construct ecological niche models
  - [`code/autokuenm/*`](https://github.com/RhettRautsaw/VenomMaps/tree/master/code/autokuenm)
- Shiny App to visualize results
  - [`app.R`](https://RhettRautsaw.app/shiny/VenomMaps/)


# Shiny App

You can view the VenomMaps Shiny App for free on my Digital Ocean Droplet: 
[RhettRautsaw.app/shiny/VenomMaps](https://RhettRautsaw.app/shiny/VenomMaps/)

## Running the Application Locally

This app can also be run through R:

```R
library(shiny)

# Easiest way is to use runGitHub
runGitHub("VenomMaps", "RhettRautsaw")

# Run a tar or zip file directly
runUrl("https://github.com/RhettRautsaw/VenomMaps/archive/master.tar.gz")
runUrl("https://github.com/RhettRautsaw/VenomMaps/archive/master.zip")
```

To run a Shiny app from a subdirectory in the repo or zip file, you can use the `subdir` argument. This repository happens to contain another copy of the app in `inst/shinyapp/`.

```R
runGitHub("VenomMaps", "RhettRautsaw", subdir = "inst/shinyapp/")

runUrl("https://github.com/RhettRautsaw/VenomMaps/archive/master.tar.gz",
  subdir = "inst/shinyapp/")
```

## Download the Data

You can clone this GitHub repository to run the Shiny App or to download all the `geojson` files.

```R
# First clone the repository with git. If you have cloned it into
# ~/VenomMaps, first go to that directory, then use runApp().
setwd("~/VenomMaps")
runApp()
```

# Todo List (Additions)

- Update distribution maps for Old World Vipers and other venomous snake species.
- Add phylogeographic tracing/playback
- Venom Information
- Diet Information