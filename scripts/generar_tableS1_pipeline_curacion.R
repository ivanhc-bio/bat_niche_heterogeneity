#####################################################################
# TABLE S1: pipeline de curación de datos por especie
#
# Combina los conteos por etapa que gbif_download_clean_thin.R ya
# guardó (n_raw, n_after_coord_clean, n_after_elev_filter, n_after_thin)
# con el conteo FINAL real (post revisión visual manual + remoción de
# outliers ambientales), tomado de tu archivo final de análisis.
#####################################################################

pkgs <- c("dplyr", "readr")
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

# Conteos por etapa (ya generados por el script principal)
conteos <- read_csv(file.path(base_dir, "base_spp_actualizado.csv"), show_col_types = FALSE) %>%
  distinct(species, n_raw, n_after_coord_clean, n_after_elev_filter, n_after_thin)

# Conteo final real por especie, desde el archivo maestro de ocurrencias
# (ya refleja: revisión visual manual + remoción de los 2 outliers
# ambientales finales, es decir, es posterior a n_after_thin)
ocurrencias_finales <- read_csv(
  file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv"),
  show_col_types = FALSE)

n_final <- ocurrencias_finales %>%
  count(species, name = "n_final")

tabla_s1 <- conteos %>%
  left_join(n_final, by = "species") %>%
  mutate(n_final = ifelse(is.na(n_final), 0, n_final),
         removidos_revision_visual_y_outliers = n_after_thin - n_final) %>%
  rename(
    `Raw GBIF records`                    = n_raw,
    `After coordinate cleaning`           = n_after_coord_clean,
    `After elevation filter`              = n_after_elev_filter,
    `After 10-km spatial thinning`        = n_after_thin,
    `Removed (manual visual review + environmental outliers)` = removidos_revision_visual_y_outliers,
    `Final n (used in analyses)`          = n_final
  ) %>%
  arrange(species)

write_csv(tabla_s1, file.path(out_dir, "TableS1_curation_pipeline.csv"))

# --- Chequeos de sanidad --------------------------------------------
message(sprintf("Especies en la tabla: %d", nrow(tabla_s1)))
message(sprintf("Suma de 'Final n': %d  (debe ser 11,517 si coincide con lo ya reportado)",
                 sum(tabla_s1$`Final n (used in analyses)`)))
negativos <- tabla_s1 %>% filter(`Removed (manual visual review + environmental outliers)` < 0)
if (nrow(negativos) > 0) {
  message("⚠ ATENCIÓN: hay especies donde el n final es MAYOR que el n post-thinning -- ")
  message("  revisa si 'TODAS_LAS_ESPECIES_final_manual.csv' es realmente posterior a thinning,")
  message("  o si hay un desfase de nombres entre archivos:")
  print(negativos %>% select(species, `After 10-km spatial thinning`, `Final n (used in analyses)`))
}
