#####################################################################
# Excluye ÚNICAMENTE el registro con problema de calidad de dato
# (espécimen en cautiverio), y refresca heterogeneidad + métricas del
# elipsoide + elevación para la especie afectada.
#
# Diagnóstico (verificado sobre TODAS_LAS_ESPECIES_final_manual.csv):
# de los 8 registros con countryCode == "US" (fuera del Neotrópico
# geográfico), solo 1 se excluye:
#   - Carollia perspicillata: 1 registro -- "Point Defiance Zoo and
#     Aquarium", Washington. Es un espécimen EN CAUTIVERIO, no una
#     ocurrencia silvestre -- se excluye por calidad de dato, no por
#     biogeografía.
# Los otros 7 (Dasypterus ega, sur de Texas) SE CONSERVAN: son
# ocurrencias silvestres genuinas y, aunque caen en la región Neártica,
# truncar la especie a solo su porción neotropical sesgaría a la baja
# su volumen de nicho estimado (ver nota de Limitaciones en el
# manuscrito sobre especies cuya distribución cruza al Neártico, como
# Choeronycteris mexicana; Arroyo-Cabrales et al. 1987).
#
# Ninguna otra especie tiene registros con countryCode fuera del
# Neotrópico (se revisaron las 53; todos los demás códigos de país
# corresponden a México, Centro y Sudamérica, o a las Antillas, que
# sí forman parte del área de estudio).
#
# Fix mínimo: Carollia perspicillata pasa de 1078 a 1077 registros
# (de 11,517 a 11,516 en el total del dataset).
#####################################################################

library(dplyr)
library(readr)
library(ntbox)

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

especie_afectada <- "Carollia perspicillata"

# ------------------------------------------------------------------
# 1. QUITAR EL REGISTRO EN CAUTIVERIO (con respaldo)
# ------------------------------------------------------------------
occ_csv <- file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv")
occ_backup <- file.path(base_dir, "resultados_por_especie",
                         "TODAS_LAS_ESPECIES_final_manual_ANTES_de_excluir_cautiverio.csv")

occ <- read_csv(occ_csv, show_col_types = FALSE)
file.copy(occ_csv, occ_backup, overwrite = TRUE)

registro_cautiverio <- occ %>%
  filter(species == especie_afectada, countryCode == "US")

message("Registro identificado como espécimen en cautiverio:")
print(registro_cautiverio %>% select(species, stateProvince, decimalLatitude, decimalLongitude, locality))
stopifnot(nrow(registro_cautiverio) == 1)

occ_actualizado <- occ %>%
  anti_join(registro_cautiverio %>% select(species, decimalLongitude, decimalLatitude),
            by = c("species", "decimalLongitude", "decimalLatitude"))

message(sprintf("\nOcurrencias: %d -> %d (se quitó %d registro)",
                 nrow(occ), nrow(occ_actualizado), nrow(occ) - nrow(occ_actualizado)))

write_csv(occ_actualizado, occ_csv)
message("Respaldo del archivo original guardado en: ", occ_backup)

# ------------------------------------------------------------------
# 2. QUITARLO DE LA TABLA DE VALORES AMBIENTALES (con respaldo)
# ------------------------------------------------------------------
valores_csv <- file.path(out_dir, "valores_ambientales_por_registro.csv")
valores_backup <- file.path(out_dir, "valores_ambientales_por_registro_ANTES_de_excluir_cautiverio.csv")

valores <- read_csv(valores_csv, show_col_types = FALSE)
file.copy(valores_csv, valores_backup, overwrite = TRUE)

valores_actualizado <- valores %>%
  anti_join(registro_cautiverio %>% select(species, decimalLongitude, decimalLatitude),
            by = c("species", "decimalLongitude", "decimalLatitude"))

message(sprintf("Valores ambientales: %d -> %d (se quitó %d registro)",
                 nrow(valores), nrow(valores_actualizado), nrow(valores) - nrow(valores_actualizado)))

write_csv(valores_actualizado, valores_csv)

# ------------------------------------------------------------------
# 3. REFRESCAR heterogeneidad_por_especie.csv
# ------------------------------------------------------------------
heterogeneidad_por_especie <- valores_actualizado %>%
  group_by(species) %>%
  summarise(
    n = n(),
    tri_media = mean(tri, na.rm = TRUE),
    tri_sd    = sd(tri, na.rm = TRUE),
    bio1_rugosidad_media  = mean(bio1_rugosidad_local, na.rm = TRUE),
    bio12_rugosidad_media = mean(bio12_rugosidad_local, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(heterogeneidad_por_especie, file.path(out_dir, "heterogeneidad_por_especie.csv"))
message("heterogeneidad_por_especie.csv refrescado")

# ------------------------------------------------------------------
# 4. REFRESCAR métricas del elipsoide + distancias por registro
#    (solo para Carollia perspicillata)
# ------------------------------------------------------------------
vars_finales <- c("bio5", "bio6", "bio16", "bio17")
nivel_elipsoide <- 0.95

metricas_previas <- read_csv(file.path(out_dir, "metricas_elipsoide_por_especie.csv"), show_col_types = FALSE)
distancias_previas <- read_csv(file.path(out_dir, "distancias_centroide_por_registro.csv"), show_col_types = FALSE)

set.seed(2026)

datos_sp <- valores_actualizado %>%
  filter(species == especie_afectada) %>%
  select(all_of(vars_finales)) %>%
  na.omit()

n <- nrow(datos_sp)

cc <- ntbox::cov_center(data = datos_sp, mve = TRUE, level = nivel_elipsoide, vars = vars_finales)
in_elip <- ntbox::inEllipsoid(centroid = cc$centroid, eShape = cc$covariance,
                               env_data = datos_sp, level = nivel_elipsoide)
dists <- in_elip$mh_dist

distancias_sp <- valores_actualizado %>%
  filter(species == especie_afectada) %>%
  select(species, decimalLongitude, decimalLatitude, all_of(vars_finales)) %>%
  na.omit() %>%
  mutate(mh_dist = dists, dentro_del_elipsoide = in_elip$in_Ellipsoid == 1)

metrica_sp <- tibble(species = especie_afectada, n = n, niche_volume = cc$niche_volume,
                      dist_centroide_media = mean(dists), dist_centroide_sd = sd(dists),
                      metodo = "mve")

message(sprintf("%-28s n=%-4d vol=%-15.1f dist_media=%-8.2f dist_sd=%-8.2f",
                 especie_afectada, n, cc$niche_volume, mean(dists), sd(dists)))

metricas_elipsoide <- metricas_previas %>%
  filter(species != especie_afectada) %>%
  bind_rows(metrica_sp) %>%
  arrange(species)

distancias_centroide_por_registro <- distancias_previas %>%
  filter(species != especie_afectada) %>%
  bind_rows(distancias_sp) %>%
  arrange(species)

write_csv(metricas_elipsoide, file.path(out_dir, "metricas_elipsoide_por_especie.csv"))
write_csv(distancias_centroide_por_registro, file.path(out_dir, "distancias_centroide_por_registro.csv"))

# ------------------------------------------------------------------
# 5. tabla_final_analisis.csv no necesita refrescar elevación aquí:
#    el registro excluido es de Carollia perspicillata, cuya elevación
#    reportada no se ve afectada por un único registro de >1000; si
#    quieres paridad exacta, re-ejecuta el bloque de elevación de
#    agregacion_elevacion.R solo para esta especie.
# ------------------------------------------------------------------
tabla_final <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)

tabla_final_actualizada <- tabla_final %>%
  filter(species != especie_afectada) %>%
  bind_rows(
    metricas_elipsoide %>% filter(species == especie_afectada) %>%
      left_join(heterogeneidad_por_especie %>% select(-n), by = "species") %>%
      left_join(tabla_final %>% filter(species == especie_afectada) %>%
                  select(species, elev_mediana, elev_rango_robusto, elev_p2_5, elev_p97_5),
                by = "species")
  ) %>%
  arrange(species)

write_csv(tabla_final_actualizada, file.path(out_dir, "tabla_final_analisis.csv"))

message("\nListo. Especie afectada: ", especie_afectada)
message("Comparación antes/después (debería ser prácticamente idéntico, es 1 registro de 1078):")
print(metricas_previas %>% filter(species == especie_afectada) %>% select(species, n, niche_volume))
print(metricas_elipsoide %>% filter(species == especie_afectada) %>% select(species, n, niche_volume))
message("\nSi todo luce bien, vuelve a correr el PGLS/phylopath de siempre con tabla_final_analisis.csv")
message("(usa las 500 filogenias ya guardadas -- no hace falta re-muestrear el árbol).")
message("\nDataset total: 11,517 -> 11,516 registros. Actualiza ese número en Resultados del manuscrito.")
