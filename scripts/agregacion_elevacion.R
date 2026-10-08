#####################################################################
# Agregación por especie de la elevación ya extraída (capa Amatulli
# et al. 2018, 1 km) -- mediana (moderador de afinidad elevacional) y
# rango elevacional robusto (percentil 2.5-97.5, en vez del rango
# crudo min-max, sensible a un solo punto atípico).
#
# NO vuelve a extraer nada -- usa la columna 'elevation_extracted' ya
# presente en tu archivo maestro (post-thinning, post-revisión de
# outliers), y la une a tabla_final_analisis.csv.
#####################################################################

library(dplyr)
library(readr)

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

occ <- read_csv(file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv"),
                 show_col_types = FALSE)

stopifnot("elevation_extracted" %in% names(occ))

# ------------------------------------------------------------------
# AGREGACIÓN POR ESPECIE
# ------------------------------------------------------------------
elevacion_por_especie <- occ %>%
  filter(!is.na(elevation_extracted)) %>%
  group_by(species) %>%
  summarise(
    n_elev = n(),
    elev_mediana = median(elevation_extracted),
    elev_p2_5  = quantile(elevation_extracted, 0.025),
    elev_p97_5 = quantile(elevation_extracted, 0.975),
    elev_rango_robusto = elev_p97_5 - elev_p2_5,
    elev_min_crudo = min(elevation_extracted),   # de referencia, no se usa como respuesta
    elev_max_crudo = max(elevation_extracted),
    .groups = "drop"
  )

write_csv(elevacion_por_especie, file.path(out_dir, "elevacion_por_especie.csv"))
message("✅ elevacion_por_especie.csv guardado (", nrow(elevacion_por_especie), " especies)")

# ------------------------------------------------------------------
# CHEQUEO RÁPIDO: posible "regla elevacional de Rapoport"
# (especies de mayor elevación mediana con rango más angosto, por
# restricción geométrica de área disponible en altura, no necesariamente
# por heterogeneidad ambiental)
# ------------------------------------------------------------------
rapoport_check <- cor.test(elevacion_por_especie$elev_mediana,
                            elevacion_por_especie$elev_rango_robusto)
message("\nChequeo Rapoport (elev. mediana vs. rango robusto):")
print(rapoport_check)

# ------------------------------------------------------------------
# UNIR A LA TABLA FINAL DE ANÁLISIS
# ------------------------------------------------------------------
tabla_final <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)

tabla_final_actualizada <- tabla_final %>%
  left_join(elevacion_por_especie %>%
              select(species, elev_mediana, elev_rango_robusto, elev_p2_5, elev_p97_5),
            by = "species")

write_csv(tabla_final_actualizada, file.path(out_dir, "tabla_final_analisis.csv"))

message("\n✅ tabla_final_analisis.csv actualizada con elev_mediana y elev_rango_robusto")
message("   (elev_mediana: moderador continuo de afinidad elevacional para el PGLS)")
message("   (elev_rango_robusto: nueva variable de respuesta -- amplitud de nicho vertical)")
