#####################################################################
# Diagnóstico de outliers AMBIENTALES (no geográficos) para las
# especies con SD de distancia al centroide desproporcionadamente
# alta respecto a su media: Artibeus phaeotis, Balantiopteryx plicata,
# Diphylla ecaudata.
#
# Identifica qué registros puntuales están "jalando" la varianza,
# usando la misma clasificación dentro/fuera del elipsoide que ya
# usa ntbox::inEllipsoid() (nivel 0.95), y los ordena por qué tan
# lejos están del centroide.
#####################################################################

library(ntbox)
library(dplyr)
library(readr)

base_dir <- "D:/capitulo2/reanalyses_17sep"
valores  <- read_csv(file.path(base_dir, "resultadosfinales", "valores_ambientales_por_registro.csv"),
                      show_col_types = FALSE)

vars_finales <- c("bio5", "bio6", "bio16", "bio17")
nivel_elipsoide <- 0.95

especies_revisar <- c("Artibeus phaeotis", "Balantiopteryx plicata", "Diphylla ecaudata")

diagnostico <- list()

for (sp in especies_revisar) {

  datos_sp <- valores %>%
    filter(species == sp) %>%
    select(species, decimalLongitude, decimalLatitude, all_of(vars_finales)) %>%
    na.omit()

  cc <- ntbox::cov_center(data = datos_sp %>% select(all_of(vars_finales)),
                           mve = TRUE, level = nivel_elipsoide, vars = vars_finales)

  in_elip <- ntbox::inEllipsoid(
    centroid = cc$centroid,
    eShape   = cc$covariance,
    env_data = datos_sp %>% select(all_of(vars_finales)),
    level    = nivel_elipsoide
  )

  datos_sp_dx <- datos_sp %>%
    mutate(
      mh_dist = in_elip$mh_dist,
      dentro_del_elipsoide_95 = in_elip$in_Ellipsoid == 1
    ) %>%
    arrange(desc(mh_dist))

  n_fuera <- sum(!datos_sp_dx$dentro_del_elipsoide_95)

  message(sprintf("\n%s: %d de %d registros caen FUERA del elipsoide al 95%%",
                   sp, n_fuera, nrow(datos_sp_dx)))
  message("Los 5 registros más alejados del centroide:")
  print(datos_sp_dx %>% select(decimalLongitude, decimalLatitude, all_of(vars_finales), mh_dist) %>% head(5))

  diagnostico[[sp]] <- datos_sp_dx
}

diagnostico_completo <- bind_rows(diagnostico)

out_path <- file.path(base_dir, "resultadosfinales", "diagnostico_outliers_ambientales.csv")
write_csv(diagnostico_completo, out_path)

message("\n✅ Diagnóstico completo guardado en: ", out_path)
message("   (ordenado por especie y por distancia al centroide, de mayor a menor)")
