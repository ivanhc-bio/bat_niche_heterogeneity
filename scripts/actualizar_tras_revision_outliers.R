#####################################################################
# Aplica las decisiones de la revisión de outliers ambientales:
#   - Artibeus phaeotis: se conservan todos los puntos (sin cambios)
#   - Balantiopteryx plicata: se retira el punto 2 (el de Belice,
#     bio5 inusualmente bajo)
#   - Diphylla ecaudata: se retira el punto 5 (el de Chiapas, muy
#     alejado geográfica y climáticamente del resto)
#
# Identifica los puntos exactos a partir de diagnostico_outliers_
# ambientales.csv (ya calculado y ordenado por distancia al centroide,
# así que no hay riesgo de re-identificarlos mal por redondeo), los
# quita del archivo maestro y de la tabla de valores ambientales, y
# refresca heterogeneidad + métricas del elipsoide + distancias por
# registro. NO vuelve a tocar GBIF ni CHELSA -- todo esto es rápido.
#####################################################################

library(dplyr)
library(readr)
library(ntbox)

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

# ------------------------------------------------------------------
# 1. IDENTIFICAR LOS PUNTOS EXACTOS A ELIMINAR
# ------------------------------------------------------------------
diagnostico <- read_csv(file.path(out_dir, "diagnostico_outliers_ambientales.csv"),
                         show_col_types = FALSE)

punto_bp <- diagnostico %>%
  filter(species == "Balantiopteryx plicata") %>%
  arrange(desc(mh_dist)) %>%
  slice(2)

punto_de <- diagnostico %>%
  filter(species == "Diphylla ecaudata") %>%
  arrange(desc(mh_dist)) %>%
  slice(5)

puntos_a_eliminar <- bind_rows(punto_bp, punto_de) %>%
  select(species, decimalLongitude, decimalLatitude)

message("Puntos identificados para eliminar:")
print(puntos_a_eliminar)

# ------------------------------------------------------------------
# 2. QUITARLOS DEL ARCHIVO MAESTRO DE OCURRENCIAS (con respaldo)
# ------------------------------------------------------------------
occ_csv <- file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv")
occ_backup <- file.path(base_dir, "resultados_por_especie",
                         "TODAS_LAS_ESPECIES_final_manual_ANTES_de_revision_outliers.csv")

occ <- read_csv(occ_csv, show_col_types = FALSE)
file.copy(occ_csv, occ_backup, overwrite = TRUE)  # respaldo antes de sobreescribir

occ_actualizado <- occ %>%
  anti_join(puntos_a_eliminar, by = c("species", "decimalLongitude", "decimalLatitude"))

message(sprintf("Ocurrencias: %d -> %d (se quitaron %d registros)",
                 nrow(occ), nrow(occ_actualizado), nrow(occ) - nrow(occ_actualizado)))

write_csv(occ_actualizado, occ_csv)
message("Respaldo del archivo original guardado en: ", occ_backup)

# ------------------------------------------------------------------
# 3. QUITARLOS DE LA TABLA DE VALORES AMBIENTALES (con respaldo)
# ------------------------------------------------------------------
valores_csv <- file.path(out_dir, "valores_ambientales_por_registro.csv")
valores_backup <- file.path(out_dir, "valores_ambientales_por_registro_ANTES_de_revision_outliers.csv")

valores <- read_csv(valores_csv, show_col_types = FALSE)
file.copy(valores_csv, valores_backup, overwrite = TRUE)

valores_actualizado <- valores %>%
  anti_join(puntos_a_eliminar, by = c("species", "decimalLongitude", "decimalLatitude"))

message(sprintf("Valores ambientales: %d -> %d (se quitaron %d registros)",
                 nrow(valores), nrow(valores_actualizado), nrow(valores) - nrow(valores_actualizado)))

write_csv(valores_actualizado, valores_csv)

# ------------------------------------------------------------------
# 4. REFRESCAR heterogeneidad_por_especie.csv (sin volver a tocar CHELSA)
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
message("✅ heterogeneidad_por_especie.csv refrescado")

# ------------------------------------------------------------------
# 5. REFRESCAR métricas del elipsoide + distancias por registro
# ------------------------------------------------------------------
vars_finales <- c("bio5", "bio6", "bio16", "bio17")
nivel_elipsoide <- 0.95

especies <- sort(unique(valores_actualizado$species))
resultados <- list()
distancias_registro <- list()

set.seed(2026)  # fija el submuestreo aleatorio del método MVE para que sea reproducible

for (sp in especies) {

  datos_sp <- valores_actualizado %>%
    filter(species == sp) %>%
    select(all_of(vars_finales)) %>%
    na.omit()

  n <- nrow(datos_sp)

  if (n <= length(vars_finales)) {
    resultados[[sp]] <- tibble(species = sp, n = n, niche_volume = NA,
                                dist_centroide_media = NA, dist_centroide_sd = NA,
                                metodo = "insuficiente")
    next
  }

  cc <- tryCatch({
    list(res = ntbox::cov_center(data = datos_sp, mve = TRUE, level = nivel_elipsoide,
                                  vars = vars_finales),
         metodo = "mve")
  }, error = function(e) {
    tryCatch(
      list(res = ntbox::cov_center(data = datos_sp, mve = FALSE, level = nivel_elipsoide,
                                    vars = vars_finales),
           metodo = "covarianza_clasica"),
      error = function(e2) NULL
    )
  })

  if (is.null(cc)) {
    resultados[[sp]] <- tibble(species = sp, n = n, niche_volume = NA,
                                dist_centroide_media = NA, dist_centroide_sd = NA,
                                metodo = "error")
    next
  }

  res <- cc$res

  in_elip <- ntbox::inEllipsoid(
    centroid = res$centroid,
    eShape   = res$covariance,
    env_data = datos_sp,
    level    = nivel_elipsoide
  )

  dists <- in_elip$mh_dist

  distancias_registro[[sp]] <- valores_actualizado %>%
    filter(species == sp) %>%
    select(species, decimalLongitude, decimalLatitude, all_of(vars_finales)) %>%
    na.omit() %>%
    mutate(mh_dist = dists, dentro_del_elipsoide = in_elip$in_Ellipsoid == 1)

  resultados[[sp]] <- tibble(
    species = sp, n = n, niche_volume = res$niche_volume,
    dist_centroide_media = mean(dists), dist_centroide_sd = sd(dists),
    metodo = cc$metodo
  )

  message(sprintf("%-28s n=%-4d vol=%-15.1f dist_media=%-8.2f dist_sd=%-8.2f (%s)",
                   sp, n, res$niche_volume, mean(dists), sd(dists), cc$metodo))
}

metricas_elipsoide <- bind_rows(resultados)
write_csv(metricas_elipsoide, file.path(out_dir, "metricas_elipsoide_por_especie.csv"))

distancias_por_registro_completo <- bind_rows(distancias_registro)
write_csv(distancias_por_registro_completo, file.path(out_dir, "distancias_centroide_por_registro.csv"))

tabla_final <- metricas_elipsoide %>%
  left_join(heterogeneidad_por_especie %>% select(-n), by = "species")
write_csv(tabla_final, file.path(out_dir, "tabla_final_analisis.csv"))

message("\n✅ Todo refrescado tras la revisión de outliers:")
message("   - metricas_elipsoide_por_especie.csv")
message("   - distancias_centroide_por_registro.csv")
message("   - heterogeneidad_por_especie.csv")
message("   - tabla_final_analisis.csv")
