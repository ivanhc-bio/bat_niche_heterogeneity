#####################################################################
# Análisis de sensibilidad: ¿qué tanto cambia el n final por especie
# según la distancia de thinning espacial?
#
# Corre el thinning a varias distancias candidatas (por defecto
# 1, 2, 5 y 10 km) sobre los datos YA limpios y filtrados por
# elevación (los mismos que usó el pipeline principal, antes del
# thinning de 10 km), y arma una tabla comparativa por especie.
#
# Requiere que ya hayas corrido el script principal al menos una vez
# (usa los CSVs "_clean.csv" guardados en gbif_clean/, y reconstruye
# el filtro de elevación con el mismo DEM). Si tu sesión de R sigue
# abierta desde entonces, reutiliza dem_terra directamente; si no,
# lo reconstruye igual que la primera vez (tarda ~1 min).
#####################################################################

pkgs <- c("dplyr", "readr", "spThin", "terra", "elevatr", "sf", "purrr", "tidyr")
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed], dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------------------------------------------------
# CONFIGURACIÓN
# ------------------------------------------------------------------
base_dir   <- "D:/capitulo2/reanalyses_17sep"
clean_dir  <- file.path(base_dir, "gbif_clean")
input_csv  <- file.path(base_dir, "base_spp.csv")
out_csv    <- file.path(base_dir, "sensibilidad_thinning_distancias.csv")

# Distancias candidatas a evaluar (km). Agrega/quita las que quieras.
thin_distances <- c(1, 2, 5, 10)

# Umbral mínimo de referencia para marcar especies "en riesgo" a 10 km
# (ajústalo cuando definas cuántas variables usarás en el elipsoide;
# una regla común es apuntar a varias veces el número de variables).
min_n_target <- 30

elev_ref <- read_csv(input_csv, show_col_types = FALSE) %>%
  distinct(species, elev_min_ref, elev_max_ref)

# ------------------------------------------------------------------
# DEM: reutiliza el de tu sesión si ya existe; si no, lo reconstruye
# ------------------------------------------------------------------
if (!exists("dem_terra")) {
  message("No encontré 'dem_terra' en la sesión, reconstruyendo el DEM (una sola vez)...")

  clean_files <- list.files(clean_dir, pattern = "_clean\\.csv$", full.names = TRUE)
  all_coords <- map_dfr(clean_files, ~ read_csv(.x, show_col_types = FALSE) %>%
                           select(decimalLongitude, decimalLatitude)) %>%
    distinct() %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

  dem_raster <- get_elev_raster(all_coords, z = 5, src = "aws", clip = "bbox")
  dem_terra  <- terra::rast(dem_raster)
} else {
  message("Usando 'dem_terra' ya cargado en la sesión.")
}

# ------------------------------------------------------------------
# LOOP: elevación + thinning a varias distancias, por especie
# ------------------------------------------------------------------
clean_files <- list.files(clean_dir, pattern = "_clean\\.csv$", full.names = TRUE)
species_names <- clean_files %>%
  basename() %>%
  gsub("_clean\\.csv$", "", .) %>%
  gsub("_", " ", .)

resultados <- list()

for (i in seq_along(clean_files)) {

  sp <- species_names[i]
  sp_clean <- read_csv(clean_files[i], show_col_types = FALSE)
  n_clean <- nrow(sp_clean)

  # --- Filtro de elevación (idéntico al del script principal) ------
  ref_row <- elev_ref %>% filter(species == sp)
  emin <- ref_row$elev_min_ref[1]
  emax <- ref_row$elev_max_ref[1]

  if (is.na(emin) || is.na(emax)) {
    sp_elev <- sp_clean
  } else {
    elev_extracted <- terra::extract(
      dem_terra, cbind(sp_clean$decimalLongitude, sp_clean$decimalLatitude)
    )
    sp_elev <- sp_clean %>%
      mutate(elevation_dem = elev_extracted[[ncol(elev_extracted)]]) %>%
      filter(!is.na(elevation_dem), elevation_dem >= emin, elevation_dem <= emax)
  }

  n_elev <- nrow(sp_elev)

  fila <- tibble(species = sp, n_clean = n_clean, n_after_elev_filter = n_elev)

  if (n_elev == 0) {
    for (d in thin_distances) fila[[paste0("n_thin_", d, "km")]] <- 0
    resultados[[sp]] <- fila
    next
  }

  # --- Thinning a cada distancia candidata --------------------------
  for (d in thin_distances) {

    thin_input <- sp_elev %>%
      mutate(species_col = sp) %>%
      select(species_col, decimalLongitude, decimalLatitude) %>%
      as.data.frame()

    thinned <- thin(
      loc.data = thin_input,
      lat.col = "decimalLatitude", long.col = "decimalLongitude",
      spec.col = "species_col", thin.par = d, reps = 100,
      locs.thinned.list.return = TRUE,
      write.files = FALSE, write.log.file = FALSE, verbose = FALSE
    )

    n_thin <- max(sapply(thinned, nrow))
    fila[[paste0("n_thin_", d, "km")]] <- n_thin
  }

  resultados[[sp]] <- fila
  message(sprintf("%-30s clean=%-6d elev=%-6d  ->  %s",
                   sp, n_clean, n_elev,
                   paste(sprintf("%dkm=%d", thin_distances,
                                 unlist(fila[paste0("n_thin_", thin_distances, "km")])),
                         collapse = "  ")))
}

tabla_sensibilidad <- bind_rows(resultados)

# Marca las especies que, incluso a la distancia MÁS FINA evaluada,
# no alcanzan tu umbral mínimo de referencia (esas necesitan atención
# aparte sin importar qué distancia elijas)
col_mas_fina <- paste0("n_thin_", min(thin_distances), "km")
tabla_sensibilidad <- tabla_sensibilidad %>%
  mutate(en_riesgo_incluso_fino = .data[[col_mas_fina]] < min_n_target)

write_csv(tabla_sensibilidad, out_csv)

message("\n✅ Listo. Tabla comparativa guardada en: ", out_csv)
message(sprintf("   %d especie(s) no alcanzan n=%d ni con la distancia más fina evaluada (%d km).",
                 sum(tabla_sensibilidad$en_riesgo_incluso_fino), min_n_target, min(thin_distances)))
