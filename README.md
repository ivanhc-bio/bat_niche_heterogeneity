# Elevational range as a candidate pathway linking topographic heterogeneity to climatic niche breadth in Neotropical bats

Data and code to reproduce the analyses of the manuscript submitted to *Ecography* (53 Neotropical bat species; phylogenetic generalized least squares and phylogenetic path analysis across 500 posterior trees).

* Manuscript: \[citation to be added upon publication]
* Archived version: \[Zenodo DOI]
* Occurrence data source: GBIF.org (17 September 2026) GBIF Occurrence Download, https://doi.org/10.15468/dl.rpjhx4

## Repository structure

|Folder|Content|
|-|-|
|`data/`|Curated inputs (see `data/DATA\_DICTIONARY.md`)|
|`scripts/`|R scripts, in the order of the workflow below|
|`results/`|Output tables of the analyses reported in the manuscript and Supporting Information|
|`trees/`|`output.nex`: 500 posterior trees (pruned to the 53 species)|

## Workflow (run in this order)

Scripts have a `base\_dir` variable at the top (currently `D:/capitulo2/reanalyses\_17sep`): **edit it to your local copy** before running. The TRI raster path (`tri\_path`) must also point to your local copy of `tri\_1KMmd\_GMTEDmd.tif`. Script comments are in Spanish.

|Step|Script|What it does|
|-|-|-|
|0|`00\_restaurar\_registros\_desde\_gbif.R`|Restores the files with coordinates from the GBIF download DOI (needed by scripts that use coordinates)|
|1|`gbif\_download\_clean\_thin.R`|Downloads GBIF records, cleans, filters by elevation and spatially thins (10 km)|
|2|`revision\_visual\_manual.R`|Manual visual review support (QGIS review was done outside R)|
|3|`auditar\_registros\_cautivos.R`, `excluir\_registro\_cautiverio.R`|Flags and excludes captive-origin records after checking source collection records|
|4|`extraccion\_bioclim\_heterogeneidad.R`|Extracts CHELSA v2.1 bioclimatic variables, TRI and local roughness per record|
|5|`metricas\_elipsoide.R`, `diagnostico\_outliers.R`, `actualizar\_tras\_revision\_outliers.R`|Minimum-volume ellipsoids (bio5, bio6, bio16, bio17), centroid distances, environmental outlier review|
|6|`agregacion\_elevacion.R`, `integrar\_gremio\_rango.R`, `variables\_geograficas\_dieta.R`, `emparejamiento\_arbol.R`|Species-level table (`data/species\_table\_final.csv`), range area, guild, tree-tip matching|
|7|`pgls\_v7\_principal.R`, `phylopath\_v7.R`|Main PGLS and phylogenetic path analysis (500 trees)|
|8|`colinealidad\_v7.R`|Collinearity checks: VIF, commonality, suppression, eight path candidates, leave-one-out, directed subsampling|
|9|`robustez\_v7.R`, `reconstruir\_ANTES\_outliers.R`|Robustness: n ≥ 30, 16 focal species, outliers reinstated, range area, null-model SES|
|10|`sensibilidad\_calidad\_v7.R`|Sensitivity to coordinate uncertainty > 10 km and to U.S. records of *Dasypterus ega*|
|11|`robustez\_null\_model\_geometrico\_v3.R`, `sensibilidad\_thinning.R`|Terrain-availability null model; thinning-distance sensitivity|
|12|`generar\_tableS1\_pipeline\_curacion.R`, `generar\_figS3\_thinning\_sensitivity.R`|Supporting tables/figures derived from the pipeline|

## Phylogenies

`trees/output.nex` is a posterior sample of 500 trees from the mammal supertree of Upham, Esselstyn \& Jetz (2019), subsampled and pruned from VertLife.org on 19 September 2026. Tip-label matching is in `data/species\_tip\_label\_matching.csv`. Scripts expect it at `<base\_dir>/trees/output.nex`.

## External inputs not included

* **CHELSA v2.1** bioclimatic rasters (Karger et al. 2017) and the **GMTED2010-based TRI** layer (Amatulli et al. 2018); both are public (the scripts read them directly).

## Software

R with the packages `rgbif`, `spThin`, `terra`, `ntbox`, `ape`, `caper`, `phylopath`, `dplyr`, `readr`, `tidyr`, `tibble`. Random seed 2026 is fixed wherever minimum-volume ellipsoids are estimated; ellipsoids are stochastic estimators, so repeated estimation can change volumes slightly (median 1.6%).

## Data licensing and attribution

GBIF records are **not redistributed** here: individual records carry different licences (CC0, CC BY 4.0 and CC BY-NC 4.0), so the repository contains only the `gbifID` of each record and derived values (climate and terrain values per record, species-level metrics). The records themselves, including coordinates, can be retrieved from the GBIF download (https://doi.org/10.15468/dl.rpjhx4) and restored with `scripts/00\_restaurar\_registros\_desde\_gbif.R`. Please cite the download DOI and the original datasets when reusing them. Code is released under the MIT licence (`LICENSE`); derived tables and results under CC BY 4.0. The posterior trees (`trees/output.nex`) derive from Upham et al. (2019) and VertLife.org; please cite them.

## Contact

Iván Hernández-Chávez
Laboratorio de Mastozoología Evolutiva y Colecciones Científicas, Facultad de Ciencias, UNAM
ivanhc@ciencias.unam.mx
Iván Alejandro Hernández Chávez

Laboratorio de Mastozoología Evolutiva y Colecciones Científicas, Facultad de Ciencias, Universidad Nacional Autónoma de México

ivanhc@ciencias.unam.mx

