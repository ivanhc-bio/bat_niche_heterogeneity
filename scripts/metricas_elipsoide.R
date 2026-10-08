#####################################################################
# Métricas del elipsoide de nicho por especie: volumen, y media/SD de
# las distancias de Mahalanobis de cada registro al centroide.
# Variables: bio5, bio6, bio16, bio17 (elegidas por baja colinealidad
# y relevancia fisiológica/de recursos para murciélagos)
#####################################################################

library(ntbox)
library(dplyr)
library(readr)
library(purrr)

base_dir <- "D:/capitulo2/reanalyses_17sep"
valores  <- read_csv(file.path(base_dir, "resultadosfinales", "valores_ambientales_por_registro.csv"),
                      show_col_types = FALSE)

vars_finales <- c("bio5", "bio6", "bio16", "bio17")
nivel_elipsoide <- 0.95  # proporción de puntos incluidos en el elipsoide

especies <- sort(unique(valores$species))

resultados <- list()
distancias_registro <- list()

set.seed(2026)  # fija el submuestreo aleatorio del método MVE para que sea reproducible

for (sp in especies) {

  datos_sp <- valores %>%
    filter(species == sp) %>%
    select(all_of(vars_finales)) %>%
    na.omit()

  n <- nrow(datos_sp)

  if (n <= length(vars_finales)) {
    message(sprintf("%-28s n=%-4d -> insuficiente (<= %d variables), se omite",
                     sp, n, length(vars_finales)))
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
    message(sprintf("%-28s n=%-4d -> cov_center falló con ambos métodos", sp, n))
    resultados[[sp]] <- tibble(species = sp, n = n, niche_volume = NA,
                                dist_centroide_media = NA, dist_centroide_sd = NA,
                                metodo = "error")
    next
  }

  res <- cc$res

  # Distancias de Mahalanobis de cada registro al centroide del elipsoide,
  # vía la función oficial de ntbox (en vez de una réplica manual)
  in_elip <- ntbox::inEllipsoid(
    centroid = res$centroid,
    eShape   = res$covariance,
    env_data = datos_sp,
    level    = nivel_elipsoide
  )

  dists <- in_elip$mh_dist

  # Guardamos también el detalle por registro (coordenadas + distancia +
  # dentro/fuera del elipsoide al nivel elegido), para uso posterior
  distancias_registro[[sp]] <- valores %>%
    filter(species == sp) %>%
    select(species, decimalLongitude, decimalLatitude, all_of(vars_finales)) %>%
    na.omit() %>%
    mutate(
      mh_dist = dists,
      dentro_del_elipsoide = in_elip$in_Ellipsoid == 1
    )

  resultados[[sp]] <- tibble(
    species = sp,
    n = n,
    niche_volume = res$niche_volume,
    dist_centroide_media = mean(dists),
    dist_centroide_sd = sd(dists),
    metodo = cc$metodo
  )

  message(sprintf("%-28s n=%-4d vol=%-15.1f dist_media=%-8.2f dist_sd=%-8.2f (%s)",
                   sp, n, res$niche_volume, mean(dists), sd(dists), cc$metodo))
}

metricas_elipsoide <- bind_rows(resultados)

write_csv(metricas_elipsoide,
          file.path(base_dir, "resultadosfinales", "metricas_elipsoide_por_especie.csv"))

message("\n✅ Listo: ", file.path(base_dir, "resultadosfinales", "metricas_elipsoide_por_especie.csv"))

# Distancias al centroide, registro por registro, todas las especies
distancias_por_registro_completo <- bind_rows(distancias_registro)
write_csv(distancias_por_registro_completo,
          file.path(base_dir, "resultadosfinales", "distancias_centroide_por_registro.csv"))
message("✅ Distancias por registro (todas las especies): ",
        file.path(base_dir, "resultadosfinales", "distancias_centroide_por_registro.csv"))

# ------------------------------------------------------------------
# Unir de una vez con la tabla de heterogeneidad (TRI + rugosidad
# climática) que ya tienes, para dejar una sola tabla lista para análisis
# ------------------------------------------------------------------
heterogeneidad <- read_csv(file.path(base_dir, "resultadosfinales", "heterogeneidad_por_especie.csv"),
                            show_col_types = FALSE)

tabla_final <- metricas_elipsoide %>%
  left_join(heterogeneidad %>% select(-n), by = "species")

write_csv(tabla_final, file.path(base_dir, "resultadosfinales", "tabla_final_analisis.csv"))
message("✅ Tabla combinada (elipsoide + heterogeneidad) lista en: ",
        file.path(base_dir, "resultadosfinales", "tabla_final_analisis.csv"))
