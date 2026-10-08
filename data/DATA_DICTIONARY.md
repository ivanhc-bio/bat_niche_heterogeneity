# Data dictionary

| File | Rows | Description |
|---|---|---|
| `record_ids_final.csv` | 11,516 | Final curated records (53 species) after thinning, manual review and exclusions, identified by `gbifID`. No coordinates (see README). Columns: `species`, `gbifID`, `countryCode`, `coordinateUncertaintyInMeters`, `year`, `license` (licence of the original record), `elevation_extracted` (m, from the 1-km GMTED-derived layer) |
| `environmental_values_per_record.csv` | 11,516 | One row per `gbifID`, same order as `record_ids_final.csv`. CHELSA v2.1 `bio1`–`bio19` (1981–2010), `tri` (Terrain Ruggedness Index, 1 km), `bio1_rugosidad_local` and `bio12_rugosidad_local` (local SD in a 3×3 cell window) |
| `species_table_final.csv` | 53 | Species-level data used in all models: `n`, `niche_volume`, `dist_centroide_media`, `dist_centroide_sd` (Mahalanobis distance to ellipsoid centroid), `tri_media`, `bio12_rugosidad_media`, `elev_mediana`, `elev_rango_robusto` (2.5–97.5th percentile spread), `family`, `guild`, `area_rango_km2` (convex hull) |
| `species_tip_label_matching.csv` | 53 | Matching of species names to tip labels of the Upham et al. (2019) phylogeny |
| `curation_pipeline_per_species.csv` | 53 | Records per species at each curation step (Table S1) |
| `thinning_sensitivity.csv` | 53 | Retained records at 1, 2, 5 and 10 km thinning (Fig. S3) |
| `null_model_SES_elevrange.csv` | 53 | Standardized effect size of elevational range relative to the terrain-availability null model |
| `captive_candidates_review.csv` | 12 | `gbifID`s flagged as possible captive-origin specimens, with the review decision (all retained after checking source collection records) |

Column names in scripts and result files are in Spanish (e.g. `bio12_rugosidad_media` = mean local precipitation heterogeneity; `elev_rango_robusto` = robust elevational range).
