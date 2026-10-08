#####################################################################
# 00. RESTAURAR LOS ARCHIVOS CON COORDENADAS A PARTIR DE LA DESCARGA DE GBIF
#
# Este repositorio NO redistribuye los registros de GBIF (licencias
# individuales, incl. CC BY-NC). Contiene solo los gbifID y valores
# derivados. Para correr los scripts que necesitan coordenadas:
#   1. Descarga el ZIP de GBIF Occurrence Download (formato SIMPLE_CSV):
#      https://doi.org/10.15468/dl.rpjhx4
#   2. Ajusta las rutas de abajo y corre este script.
# Salidas (mismo formato que usan los demas scripts):
#   <res_dir>/TODAS_LAS_ESPECIES_final_manual.csv   (registros finales + coordenadas)
#   <out_dir>/valores_ambientales_por_registro.csv  (species, lon, lat, bio1-19, tri, rugosidades)
#####################################################################
pkgs <- c("dplyr", "readr")
invisible(lapply(pkgs, library, character.only = TRUE))

repo_dir <- "."                                   # carpeta de este repositorio
gbif_zip <- "ruta/a/la_descarga_de_GBIF.zip"      # <- ajusta
base_dir <- "D:/capitulo2/reanalyses_17sep"       # <- ajusta (igual que en los demas scripts)
res_dir  <- file.path(base_dir, "resultados_por_especie")
out_dir  <- file.path(base_dir, "resultadosfinales")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ids <- read_csv(file.path(repo_dir, "data", "record_ids_final.csv"),
                col_types = cols(.default = "c"), show_col_types = FALSE)
env <- read_csv(file.path(repo_dir, "data", "environmental_values_per_record.csv"),
                col_types = cols(.default = "c"), show_col_types = FALSE)
stopifnot(nrow(ids) == 11516, nrow(env) == 11516, all(ids$gbifID == env$gbifID))

tmp <- tempfile(); dir.create(tmp)
unzip(gbif_zip, exdir = tmp)
f <- list.files(tmp, pattern = "\\.csv$", full.names = TRUE)[1]
gb <- read_tsv(f, col_types = cols(.default = "c"), quote = "", show_col_types = FALSE) %>%
  select(gbifID, decimalLatitude, decimalLongitude) %>% distinct(gbifID, .keep_all = TRUE)

occ <- ids %>% left_join(gb, by = "gbifID")
stopifnot(!anyNA(occ$decimalLatitude), !anyNA(occ$decimalLongitude))   # todos los gbifID deben estar en la descarga

valores <- env %>%
  left_join(occ %>% select(gbifID, decimalLongitude, decimalLatitude), by = "gbifID") %>%
  select(species, decimalLongitude, decimalLatitude, everything(), -gbifID)

write_csv(occ,     file.path(res_dir, "TODAS_LAS_ESPECIES_final_manual.csv"))
write_csv(valores, file.path(out_dir, "valores_ambientales_por_registro.csv"))
message("Listo: 11,516 registros restaurados con coordenadas.")
message("Nota: se restauran solo las columnas necesarias para los analisis; las demas columnas del registro de GBIF pueden volver a obtenerse por gbifID.")
