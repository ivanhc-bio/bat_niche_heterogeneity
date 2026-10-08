#####################################################################
# 1) Tamaño de rango geográfico por especie (área del polígono convexo
#    de los puntos de ocurrencia finales, en km²) -- covariable de
#    control para el efecto de heterogeneidad sobre el nicho.
# 2) Plantilla de gremio trófico/dieta por especie -- clasificación
#    preliminar a NIVEL DE GÉNERO basada en conocimiento general de
#    quirópteros neotropicales. *** DEBES VERIFICARLA CONTRA TU LIBRO
#    DE REFERENCIA *** antes de usarla -- especialmente para especies
#    con dietas mixtas/oportunistas dentro de un género.
#####################################################################

pkgs <- c("sf", "dplyr", "readr")
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed], dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

occ <- read_csv(file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv"),
                 show_col_types = FALSE)

# ------------------------------------------------------------------
# 1. ÁREA DEL POLÍGONO CONVEXO (km²), en una proyección de área
#    equivalente centrada en el continente americano (no en grados,
#    que distorsiona área)
# ------------------------------------------------------------------
proj_area_equivalente <- "+proj=cea +lat_ts=0 +lon_0=-75 +datum=WGS84 +units=km"

calcular_area_hull <- function(df) {
  if (nrow(df) < 3) return(NA_real_)  # un polígono necesita >=3 puntos no colineales
  pts <- sf::st_as_sf(df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
    sf::st_transform(proj_area_equivalente)
  hull <- sf::st_convex_hull(sf::st_union(pts))
  as.numeric(sf::st_area(hull))  # km^2, porque units=km en la proyección
}

rango_geografico <- occ %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  group_by(species) %>%
  group_modify(~ tibble(area_rango_km2 = calcular_area_hull(.x), n_puntos_hull = nrow(.x))) %>%
  ungroup()

n_insuficientes <- sum(is.na(rango_geografico$area_rango_km2))
if (n_insuficientes > 0) {
  message(sprintf("⚠ %d especie(s) con <3 puntos, no se pudo calcular el hull:", n_insuficientes))
  print(rango_geografico %>% filter(is.na(area_rango_km2)))
}

write_csv(rango_geografico, file.path(out_dir, "rango_geografico_por_especie.csv"))
message("✅ rango_geografico_por_especie.csv guardado")

# ------------------------------------------------------------------
# 2. PLANTILLA DE GREMIO TRÓFICO (clasificación preliminar a revisar)
# ------------------------------------------------------------------
especies <- sort(unique(occ$species))

genero <- function(sp) strsplit(sp, " ")[[1]][1]

clasificacion_generos <- c(
  # Frugívoros
  "Artibeus" = "Frugivoro", "Carollia" = "Frugivoro", "Sturnira" = "Frugivoro",
  "Platyrrhinus" = "Frugivoro", "Uroderma" = "Frugivoro", "Chiroderma" = "Frugivoro",
  "Enchisthenes" = "Frugivoro", "Dermanura" = "Frugivoro",
  # Nectarivoros/polinivoros
  "Glossophaga" = "Nectarivoro", "Anoura" = "Nectarivoro",
  # Hematofagos
  "Desmodus" = "Hematofago", "Diphylla" = "Hematofago",
  # Carnivoros / insectivoros de emboscada (gleaning)
  "Trachops" = "Carnivoro", "Micronycteris" = "Insectivoro_gleaning",
  # Piscivoro/insectivoro
  "Noctilio" = "Piscivoro_insectivoro",
  # Insectivoros aereos (Molossidae, Vespertilionidae, Emballonuridae, Mormoopidae)
  "Molossus" = "Insectivoro_aereo", "Eumops" = "Insectivoro_aereo",
  "Nyctinomops" = "Insectivoro_aereo", "Promops" = "Insectivoro_aereo",
  "Cynomops" = "Insectivoro_aereo", "Molossops" = "Insectivoro_aereo",
  "Myotis" = "Insectivoro_aereo", "Eptesicus" = "Insectivoro_aereo",
  "Lasiurus" = "Insectivoro_aereo", "Dasypterus" = "Insectivoro_aereo",
  "Rhogeessa" = "Insectivoro_aereo",
  "Saccopteryx" = "Insectivoro_aereo", "Rhynchonycteris" = "Insectivoro_aereo",
  "Balantiopteryx" = "Insectivoro_aereo", "Peropteryx" = "Insectivoro_aereo",
  "Diclidurus" = "Insectivoro_aereo",
  "Pteronotus" = "Insectivoro_aereo"
)

plantilla_dieta <- tibble(
  species = especies,
  genero = sapply(especies, genero),
  gremio_trofico_PRELIMINAR = clasificacion_generos[sapply(especies, genero)],
  verificado = FALSE
)

sin_clasificar <- plantilla_dieta %>% filter(is.na(gremio_trofico_PRELIMINAR))
if (nrow(sin_clasificar) > 0) {
  message("\n⚠ Géneros sin clasificación preliminar (revisar manualmente):")
  print(sin_clasificar)
}

write_csv(plantilla_dieta, file.path(out_dir, "plantilla_gremio_trofico.csv"))
message("\n✅ plantilla_gremio_trofico.csv guardada -- REVÍSALA contra tu bibliografía,")
message("   especialmente Trachops (a veces clasificado más ampliamente como carnívoro/omnívoro)")
message("   y cualquier especie con dieta documentada como mixta/oportunista.")
message("   Cuando la confirmes, cambia 'verificado' a TRUE y la unimos a tabla_final_analisis.csv.")
