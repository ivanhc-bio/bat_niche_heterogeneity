#####################################################################
# Integra a tabla_final_analisis.csv:
#   - Gremio trófico y familia, desde base_spp.csv (con la corrección
#     de Uroderma magnirostrum -> Frugivorous y Trachops cirrhosus ->
#     Carnivorous, que estaban intercambiadas)
#   - Área del rango geográfico (convex hull, km2), ya calculada en
#     rango_geografico_por_especie.csv
#####################################################################

library(dplyr)
library(readr)

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

# ------------------------------------------------------------------
# 1. GREMIO TRÓFICO + FAMILIA (deduplicado, con la corrección)
# ------------------------------------------------------------------
base_spp <- read_csv(file.path(base_dir, "base_spp.csv"), show_col_types = FALSE)

gremio_familia <- base_spp %>%
  distinct(species, family, guild) %>%
  mutate(
    guild = case_when(
      species == "Uroderma magnirostrum" ~ "Frugivorous",
      species == "Trachops cirrhosus"    ~ "Carnivorous",
      TRUE ~ guild
    )
  )

stopifnot(nrow(gremio_familia) == n_distinct(gremio_familia$species))  # 1 fila por especie

write_csv(gremio_familia, file.path(out_dir, "gremio_familia_por_especie.csv"))
message("✅ gremio_familia_por_especie.csv guardado (con la corrección Uroderma/Trachops)")
print(gremio_familia %>% count(guild, sort = TRUE))

# ------------------------------------------------------------------
# 2. UNIR A LA TABLA FINAL (gremio, familia, y rango geográfico si ya existe)
# ------------------------------------------------------------------
tabla_final <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)

tabla_final_actualizada <- tabla_final %>%
  dplyr::select(-any_of(c("family", "guild", "area_rango_km2"))) %>%  # por si ya existían de una corrida anterior
  left_join(gremio_familia, by = "species")

rango_geo_path <- file.path(out_dir, "rango_geografico_por_especie.csv")
if (file.exists(rango_geo_path)) {
  rango_geografico <- read_csv(rango_geo_path, show_col_types = FALSE)
  tabla_final_actualizada <- tabla_final_actualizada %>%
    left_join(rango_geografico %>% dplyr::select(species, area_rango_km2), by = "species")
  message("✅ Área de rango geográfico también integrada")
} else {
  message("ℹ No encontré rango_geografico_por_especie.csv todavía -- corre primero")
  message("  variables_geograficas_dieta.R si quieres integrarlo en este mismo paso.")
}

write_csv(tabla_final_actualizada, file.path(out_dir, "tabla_final_analisis.csv"))

message("\n✅ tabla_final_analisis.csv actualizada.")
message("   Verificación rápida -- Uroderma magnirostrum / Trachops cirrhosus:")
print(tabla_final_actualizada %>%
        filter(species %in% c("Uroderma magnirostrum", "Trachops cirrhosus")) %>%
        dplyr::select(species, family, guild))
