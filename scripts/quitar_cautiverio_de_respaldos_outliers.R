#####################################################################
# Quita el registro de cautiverio (Carollia perspicillata, US) de los
# RESPALDOS "ANTES_de_revision_outliers", para que la robustez
# "con outliers" (robustez_con_outliers_ambientales.R) difiera del
# análisis principal SOLO en los 2 outliers ambientales y no también
# en el registro del zoológico. Hace copia de seguridad de cada archivo.
#####################################################################
library(dplyr); library(readr)

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

occ_antes_csv <- file.path(base_dir, "resultados_por_especie",
                           "TODAS_LAS_ESPECIES_final_manual_ANTES_de_revision_outliers.csv")
val_antes_csv <- file.path(out_dir, "valores_ambientales_por_registro_ANTES_de_revision_outliers.csv")

occ_antes <- read_csv(occ_antes_csv, show_col_types = FALSE)
val_antes <- read_csv(val_antes_csv, show_col_types = FALSE)

cautiverio <- occ_antes %>% filter(species == "Carollia perspicillata", countryCode == "US")
stopifnot(nrow(cautiverio) == 1)   # si falla, avísame: no es el mismo registro
print(cautiverio %>% select(species, stateProvince, decimalLatitude, decimalLongitude, locality))

file.copy(occ_antes_csv, sub("\\.csv$", "_CON_cautiverio.csv", occ_antes_csv), overwrite = FALSE)
file.copy(val_antes_csv, sub("\\.csv$", "_CON_cautiverio.csv", val_antes_csv), overwrite = FALSE)

k <- cautiverio %>% select(species, decimalLongitude, decimalLatitude)
occ_nuevo <- anti_join(occ_antes, k, by = c("species", "decimalLongitude", "decimalLatitude"))
val_nuevo <- anti_join(val_antes, k, by = c("species", "decimalLongitude", "decimalLatitude"))

message(sprintf("occ ANTES: %d -> %d | valores ANTES: %d -> %d",
                nrow(occ_antes), nrow(occ_nuevo), nrow(val_antes), nrow(val_nuevo)))
stopifnot(nrow(occ_antes) - nrow(occ_nuevo) == 1, nrow(val_antes) - nrow(val_nuevo) == 1)

write_csv(occ_nuevo, occ_antes_csv)
write_csv(val_nuevo, val_antes_csv)
message("Listo.")
