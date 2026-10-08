#####################################################################
# Descarga, limpieza y thinning espacial de registros GBIF
# Autor: generado con Claude
# Input : base_spp.csv (columna "species" con nombres científicos;
#          puede tener especies repetidas para distintos subsets)
# Output: - un .csv limpio/thinned por especie
#         - base_spp.csv actualizado con los conteos por etapa
#####################################################################

# ------------------------------------------------------------------
# 0. PAQUETES
# ------------------------------------------------------------------
pkgs <- c("rgbif", "dplyr", "readr", "purrr", "tidyr", "spThin",
          "CoordinateCleaner", "stringr", "janitor", "elevatr", "terra", "sf")

installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed], dependencies = TRUE)

invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------------------------------------------------
# 1. CONFIGURACIÓN DE RUTAS
# ------------------------------------------------------------------
base_dir     <- "D:/capitulo2/reanalyses_17sep"
input_csv    <- file.path(base_dir, "base_spp.csv")
raw_dir      <- file.path(base_dir, "gbif_raw")        # descarga cruda por especie
clean_dir    <- file.path(base_dir, "gbif_clean")      # sin duplicados
thin_dir     <- file.path(base_dir, "gbif_thinned")    # después de thinning 10 km
flagged_csv  <- file.path(base_dir, "especies_sin_match_o_revisar.csv")
output_csv   <- file.path(base_dir, "base_spp_actualizado.csv")

dir.create(raw_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(clean_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(thin_dir,  showWarnings = FALSE, recursive = TRUE)

thin_km <- 10   # distancia mínima de thinning en km

# --- Filtros adicionales solicitados por el usuario ----------------
# a) Solo especímenes preservados (voucher físico), igual que en QGIS.
#    Si para alguna especie esto deja muy pocos registros, considera
#    ampliar a c("PRESERVED_SPECIMEN", "MATERIAL_SAMPLE").
basis_of_record <- "PRESERVED_SPECIMEN"

# b) Filtro temporal. CHELSA v2.1 (bio1-bio19) representa el periodo de
#    referencia climatológica 1981-2010 (Karger et al. 2017, 2021), así
#    que por defecto restringimos las ocurrencias a esa misma ventana
#    para evitar desajuste temporal entre presencia y clima.
#    Si te quedas con muestras muy pequeñas para alguna especie, una
#    alternativa común en la literatura es ampliar el rango (ej. 1970-2020).
year_start <- 1970
year_end   <- 2026

# ------------------------------------------------------------------
# 2. CREDENCIALES GBIF
# ------------------------------------------------------------------
# occ_download() requiere una cuenta gratuita en https://www.gbif.org/user/profile
# NO pongas tus credenciales directamente en el script. Configúralas una
# sola vez en tu .Renviron (usethis::edit_r_environ()) agregando:
#
#   GBIF_USER=tu_usuario
#   GBIF_PWD=tu_contraseña
#   GBIF_EMAIL=tu_correo@ejemplo.com
#
# y reinicia R. rgbif las detecta automáticamente.

stopifnot(
  "Falta GBIF_USER/GBIF_PWD/GBIF_EMAIL en el entorno. Ver instrucciones arriba." =
    all(nzchar(c(Sys.getenv("GBIF_USER"), Sys.getenv("GBIF_PWD"), Sys.getenv("GBIF_EMAIL"))))
)

# ------------------------------------------------------------------
# 3. LEER LISTA DE ESPECIES (únicas, sin importar subset duplicado)
# ------------------------------------------------------------------
base_spp <- read_csv(input_csv, show_col_types = FALSE)

sp_unique <- base_spp %>%
  distinct(species) %>%
  pull(species) %>%
  sort()

message(sprintf("Especies únicas a descargar: %d", length(sp_unique)))

# ------------------------------------------------------------------
# 4. RESOLVER NOMBRES / SINONIMIAS CONTRA EL BACKBONE DE GBIF
# ------------------------------------------------------------------
# name_backbone_checklist() devuelve, para cada nombre (sea el aceptado
# o un sinónimo), el "usageKey" del taxón ACEPTADO por el backbone.
# Usar ese usageKey como taxonKey en la descarga asegura que se
# incluyan los registros reportados bajo sinónimos, porque GBIF indexa
# cada ocurrencia con el taxonKey aceptado independientemente del
# nombre científico original con el que fue subida.

matches <- name_backbone_checklist(sp_unique) %>%
  as_tibble() %>%
  mutate(original_name = sp_unique)

# Revisa manualmente: sin match, match difuso (fuzzy) o rango != SPECIES
to_review <- matches %>%
  filter(matchType %in% c("NONE", "FUZZY") | rank != "SPECIES" | is.na(usageKey))

if (nrow(to_review) > 0) {
  write_csv(to_review, flagged_csv)
  message(sprintf(
    "⚠ %d nombre(s) requieren revisión manual (sin match exacto / posible sinonimia dudosa). Ver: %s",
    nrow(to_review), flagged_csv
  ))
}

# Nos quedamos con los que sí tienen un usageKey de especie válido
matches_ok <- matches %>%
  filter(!is.na(usageKey), rank == "SPECIES")

taxon_keys <- as.integer(unique(matches_ok$usageKey))
message(sprintf("Taxon keys válidos para descarga: %d", length(taxon_keys)))

# Tabla de referencia nombre original -> taxonKey usado -> nombre aceptado
key_lookup <- matches_ok %>%
  select(original_name, usageKey, accepted_name = species, status) %>%
  mutate(usageKey = as.integer(usageKey))
write_csv(key_lookup, file.path(base_dir, "tabla_equivalencia_nombres_gbif.csv"))

# ------------------------------------------------------------------
# 5. DESCARGA MASIVA (una sola solicitud asíncrona para todas las especies)
# ------------------------------------------------------------------
gbif_dl <- occ_download(
  pred_in("taxonKey", taxon_keys),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred("occurrenceStatus", "PRESENT"),
  pred("basisOfRecord", basis_of_record),
  pred_gte("year", year_start),
  pred_lte("year", year_end),
  format = "SIMPLE_CSV"
)

message("Solicitud de descarga enviada a GBIF. Key: ", gbif_dl)
occ_download_wait(gbif_dl)                      # espera a que GBIF prepare el archivo
dl_path <- occ_download_get(gbif_dl, path = raw_dir, overwrite = TRUE)
occ_raw <- occ_download_import(dl_path)

message(sprintf("Registros crudos descargados en total: %d", nrow(occ_raw)))
write_csv(occ_raw, file.path(raw_dir, "gbif_raw_all_species.csv"))

# ------------------------------------------------------------------
# 6. LIMPIEZA Y THINNING POR ESPECIE
# ------------------------------------------------------------------
# Unimos el taxonKey descargado con el nombre original (puede haber
# más de un nombre original -> mismo usageKey si eran sinónimos)

occ_raw <- occ_raw %>%
  left_join(key_lookup, by = c("taxonKey" = "usageKey"))

# Rango de elevación de referencia por especie (verificado: idéntico
# entre subsets 16 y 48 para una misma especie en tu base_spp.csv)
elev_ref <- base_spp %>%
  distinct(species, elev_min_ref, elev_max_ref)

# --- DEM único para TODA la extensión (evita descargar un ráster ---
# --- gigante por cada una de las 53 especies) -----------------------
# z = 5 (~4.8 km de resolución) es suficiente para un filtro de rango
# altitudinal amplio; si necesitas más precisión sube z, pero ten en
# cuenta que el tamaño de descarga crece muy rápido con el zoom.
dem_z <- 7

message(sprintf("Descargando DEM una sola vez (zoom=%d) para toda la extensión de tus registros...", dem_z))

occ_valid_coords <- occ_raw %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  filter(!(decimalLatitude == 0 & decimalLongitude == 0)) %>%
  distinct(decimalLongitude, decimalLatitude) %>%
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

dem_raster <- get_elev_raster(occ_valid_coords, z = dem_z,
                               src = "aws", clip = "bbox")
dem_terra  <- terra::rast(dem_raster)   # asegura un SpatRaster de terra

message("DEM listo. Comenzando limpieza/filtrado por especie...")

summary_list <- list()

for (sp in sp_unique) {

  # Todos los nombres originales que mapean al mismo taxón aceptado
  keys_for_sp <- key_lookup %>% filter(original_name == sp) %>% pull(usageKey)

  if (length(keys_for_sp) == 0) {
    summary_list[[sp]] <- tibble(species = sp, n_raw = 0,
                                  n_after_coord_clean = 0,
                                  n_after_elev_filter = 0, n_after_thin = 0,
                                  n_pres = 0)
    next
  }

  sp_occ <- occ_raw %>% filter(taxonKey %in% keys_for_sp)
  n_raw <- nrow(sp_occ)

  # --- 6a. Coordenadas válidas y sin duplicados --------------------
  sp_clean <- sp_occ %>%
    filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
    filter(!(decimalLatitude == 0 & decimalLongitude == 0)) %>%
    distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

  # --- (Opcional recomendado, desactivado por defecto) --------------
  # Filtros adicionales de calidad con CoordinateCleaner (coordenadas
  # en capitales, centroides de país/provincia, en instituciones, mar, etc.)
  # Descomenta si los quieres aplicar:
  #
  # flags <- clean_coordinates(
  #   x = sp_clean, lon = "decimalLongitude", lat = "decimalLatitude",
  #   species = "species",
  #   tests = c("capitals", "centroids", "equal", "gbif", "institutions",
  #             "seas", "zeros")
  # )
  # sp_clean <- sp_clean[flags$.summary, ]

  n_after_coord_clean <- nrow(sp_clean)

  if (n_after_coord_clean == 0) {
    summary_list[[sp]] <- tibble(species = sp, n_raw = n_raw,
                                  n_after_coord_clean = 0,
                                  n_after_elev_filter = 0, n_after_thin = 0,
                                  n_pres = 0)
    next
  }

  write_csv(sp_clean, file.path(clean_dir, paste0(gsub(" ", "_", sp), "_clean.csv")))

  # --- 6b. Filtro de elevación (DEM real, no el campo crudo de GBIF) ---
  # El campo "elevation" de GBIF suele venir vacío/poco confiable porque
  # depende de si el colector lo reportó. En su lugar, extraemos la
  # elevación real del DEM descargado una sola vez arriba (dem_terra).
  ref_row <- elev_ref %>% filter(species == sp)
  emin <- ref_row$elev_min_ref[1]
  emax <- ref_row$elev_max_ref[1]

  if (is.na(emin) || is.na(emax)) {
    # Sin rango de referencia bibliográfico: no se filtra por elevación
    sp_elev <- sp_clean
  } else {
    elev_extracted <- terra::extract(
      dem_terra,
      cbind(sp_clean$decimalLongitude, sp_clean$decimalLatitude)
    )

    sp_elev <- sp_clean %>%
      # última columna = valor de elevación, sea que terra incluya o no
      # una columna "ID" al inicio (varía según versión del paquete)
      mutate(elevation_dem = elev_extracted[[ncol(elev_extracted)]]) %>%
      filter(!is.na(elevation_dem),
             elevation_dem >= emin, elevation_dem <= emax)
  }

  n_after_elev_filter <- nrow(sp_elev)

  if (n_after_elev_filter == 0) {
    summary_list[[sp]] <- tibble(species = sp, n_raw = n_raw,
                                  n_after_coord_clean = n_after_coord_clean,
                                  n_after_elev_filter = 0, n_after_thin = 0,
                                  n_pres = 0)
    next
  }

  # --- 6c. Thinning espacial a thin_km (10 km) ----------------------
  sp_clean_thin_input <- sp_elev %>%
    mutate(species_col = sp) %>%
    select(species_col, decimalLongitude, decimalLatitude) %>%
    as.data.frame()

  thinned <- thin(
    loc.data       = sp_clean_thin_input,
    lat.col        = "decimalLatitude",
    long.col       = "decimalLongitude",
    spec.col       = "species_col",
    thin.par       = thin_km,
    reps           = 100,
    locs.thinned.list.return = TRUE,
    write.files    = FALSE,
    write.log.file = FALSE,
    verbose        = FALSE
  )

  # spThin genera `reps` soluciones; nos quedamos con la que retiene más puntos
  best_rep <- thinned[[which.max(sapply(thinned, nrow))]]

  sp_thinned_full <- sp_elev %>%
    semi_join(best_rep, by = c("decimalLongitude" = "Longitude",
                                "decimalLatitude"  = "Latitude"))

  n_after_thin <- nrow(sp_thinned_full)

  write_csv(sp_thinned_full,
            file.path(thin_dir, paste0(gsub(" ", "_", sp), "_thinned10km.csv")))

  summary_list[[sp]] <- tibble(
    species              = sp,
    n_raw                = n_raw,
    n_after_coord_clean  = n_after_coord_clean,
    n_after_elev_filter  = n_after_elev_filter,
    n_after_thin         = n_after_thin,
    n_pres               = n_after_thin   # registros de presencia finales para el ENM
  )

  message(sprintf("%-30s raw=%-6d clean=%-6d elev=%-6d thinned(%dkm)=%-6d",
                   sp, n_raw, n_after_coord_clean, n_after_elev_filter, thin_km, n_after_thin))
}

summary_df <- bind_rows(summary_list)

# ------------------------------------------------------------------
# 7. ACTUALIZAR base_spp.csv CON LOS CONTEOS (respeta filas duplicadas)
# ------------------------------------------------------------------
base_spp_updated <- base_spp %>%
  select(-any_of(c("n_raw", "n_after_coord_clean", "n_after_elev_filter",
                    "n_after_thin", "n_pres"))) %>%
  left_join(summary_df, by = "species")

write_csv(base_spp_updated, output_csv)

message("\n✅ Listo.")
message("  - CSV crudo combinado:   ", file.path(raw_dir, "gbif_raw_all_species.csv"))
message("  - CSVs limpios por spp:  ", clean_dir)
message("  - CSVs thinned por spp:  ", thin_dir)
message("  - Tabla de equivalencias de nombres (sinonimias): ",
        file.path(base_dir, "tabla_equivalencia_nombres_gbif.csv"))
message("  - base_spp.csv actualizado con conteos: ", output_csv)
if (nrow(to_review) > 0) {
  message("  - ⚠ Revisar nombres sin match/fuzzy: ", flagged_csv)
}
