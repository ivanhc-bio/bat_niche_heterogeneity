#####################################################################
# Extracción de variables bioclimáticas (CHELSA v2.1) + heterogeneidad
# topográfica (TRI, Amatulli et al. 2018) y climática (rugosidad local
# de bio01/bio12), a partir de tu base final de ocurrencias.
#
# CHELSA se lee vía streaming (GDAL /vsicurl/) y se recorta al extent
# de tus registros + un margen de seguridad -- NO se descargan las
# capas globales completas (cada bio# global pesa cientos de MB).
#
# La capa TRI de Amatulli et al. (2018) sí debes descargarla tú
# manualmente (es un solo archivo, no tan pesado como para necesitar
# streaming) desde:
#   https://www.earthenv.org/topography
#   (o el DOI alternativo: https://doi.pangaea.de/10.1594/PANGAEA.867115)
# Busca la variable "Terrain Ruggedness Index (TRI)" a 1 km de grano
# (normalmente el archivo se llama algo como tri_1KMmd_GMTEDmd.tif).
# Pon la ruta donde la guardes en `tri_path` abajo.
#####################################################################

pkgs <- c("terra", "dplyr", "readr", "purrr", "stringr", "tidyr")
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed], dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------------------------------------------------
# CONFIGURACIÓN
# ------------------------------------------------------------------
base_dir     <- "D:/capitulo2/reanalyses_17sep"
occ_csv      <- file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv")
tri_path     <- file.path(base_dir, "variables", "tri_1KMmd_GMTEDmd.tif")  # <- ajusta si tu archivo se llama distinto
out_dir      <- file.path(base_dir, "resultadosfinales")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Ventana móvil para el índice de rugosidad climática (en celdas de
# CHELSA, ~1 km/celda). 3 = ventana de 3x3 celdas (~3x3 km), análoga
# a la ventana usada por Amatulli et al. (2018) para TRI.
window_size <- 3

# Margen alrededor del extent de tus puntos (en grados) para que el
# recorte de CHELSA tenga suficiente borde y el focal() no genere NA
# justo en tus registros más periféricos.
buffer_deg <- 0.25

# ------------------------------------------------------------------
# CONFIGURAR GDAL PARA LECTURA EFICIENTE VÍA HTTP (streaming, sin
# descargar los .tif completos)
# ------------------------------------------------------------------
terra::setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
terra::setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
terra::setGDALconfig("GDAL_HTTP_MULTIPLEX", "YES")
terra::setGDALconfig("VSI_CACHE", "TRUE")

# ------------------------------------------------------------------
# 1. LEER OCURRENCIAS FINALES
# ------------------------------------------------------------------
occ <- read_csv(occ_csv, show_col_types = FALSE)

cols_necesarias <- c("species", "decimalLongitude", "decimalLatitude")
faltantes <- setdiff(cols_necesarias, names(occ))
if (length(faltantes) > 0) {
  stop(sprintf(
    "Faltan columnas esperadas (%s) en tu CSV. Columnas disponibles: %s",
    paste(faltantes, collapse = ", "), paste(names(occ), collapse = ", ")
  ))
}

occ <- occ %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude))
message(sprintf("Registros leídos: %d  |  Especies: %d",
                 nrow(occ), n_distinct(occ$species)))

# Extent de trabajo: bbox de todos tus puntos + margen de seguridad
ext_puntos <- terra::ext(
  min(occ$decimalLongitude) - buffer_deg, max(occ$decimalLongitude) + buffer_deg,
  min(occ$decimalLatitude)  - buffer_deg, max(occ$decimalLatitude)  + buffer_deg
)

# ------------------------------------------------------------------
# 2. CHELSA: leer vía streaming y recortar (sin descargar el global)
# ------------------------------------------------------------------
chelsa_base_url <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"

message("\nDescargando (streaming, solo el recorte necesario) las 19 bioclimáticas de CHELSA...")

bio_layers <- list()
for (i in 1:19) {
  var_name <- paste0("bio", i)
  url <- sprintf("%s/CHELSA_bio%d_1981-2010_V.2.1.tif", chelsa_base_url, i)

  r <- tryCatch({
    full <- terra::rast(paste0("/vsicurl/", url))
    terra::crop(full, ext_puntos)
  }, error = function(e) {
    message(sprintf("  ⚠ No pude leer %s (%s). Revisa la URL/tu conexión.", var_name, url))
    NULL
  })

  if (!is.null(r)) {
    names(r) <- var_name
    bio_layers[[var_name]] <- r
    message(sprintf("  ✓ %s recortado (%d x %d celdas)", var_name, nrow(r), ncol(r)))
  }
}

if (length(bio_layers) == 0) {
  stop(paste(
    "No se pudo leer ninguna capa de CHELSA. Verifica en",
    "https://chelsa-climate.org/downloads/ que la URL base siga siendo",
    "la misma, y actualiza 'chelsa_base_url' si cambió."
  ))
}

chelsa_stack <- terra::rast(bio_layers)
message(sprintf("\nStack de CHELSA listo: %d variable(s) recortada(s).", terra::nlyr(chelsa_stack)))

# ------------------------------------------------------------------
# 3. TRI (Amatulli et al. 2018) -- archivo local, se recorta también
# ------------------------------------------------------------------
if (!file.exists(tri_path)) {
  stop(sprintf(
    "No encontré el archivo TRI en: %s\nDescárgalo de https://www.earthenv.org/topography y ajusta 'tri_path'.",
    tri_path
  ))
}

tri_full <- terra::rast(tri_path)
tri_crop <- terra::crop(tri_full, ext_puntos)
names(tri_crop) <- "tri"
message(sprintf("TRI recortado (%d x %d celdas).", nrow(tri_crop), ncol(tri_crop)))

# ------------------------------------------------------------------
# 4. ÍNDICE DE RUGOSIDAD CLIMÁTICA LOCAL (bio01 y bio12)
#    Analogía directa con TRI: SD en una ventana móvil de 3x3 celdas
# ------------------------------------------------------------------
message(sprintf("\nCalculando rugosidad climática local (ventana %dx%d) para bio01 y bio12...",
                 window_size, window_size))

bio1_rough  <- terra::focal(chelsa_stack[["bio1"]],  w = window_size, fun = "sd", na.rm = TRUE)
bio12_rough <- terra::focal(chelsa_stack[["bio12"]], w = window_size, fun = "sd", na.rm = TRUE)
names(bio1_rough)  <- "bio1_rugosidad_local"
names(bio12_rough) <- "bio12_rugosidad_local"

# ------------------------------------------------------------------
# 5. EXTRACCIÓN POR REGISTRO (todas las bio + TRI + rugosidad clim.)
# ------------------------------------------------------------------
message("\nExtrayendo valores en cada registro...")

coords <- cbind(occ$decimalLongitude, occ$decimalLatitude)

extraer_col <- function(raster_layer, coords, colname) {
  r <- terra::extract(raster_layer, coords)
  as.numeric(unlist(r[[colname]]))
}

ex_tri    <- extraer_col(tri_crop, coords, "tri")
ex_rough1 <- extraer_col(bio1_rough, coords, "bio1_rugosidad_local")
ex_rough2 <- extraer_col(bio12_rough, coords, "bio12_rugosidad_local")

ex_bio_df <- as.data.frame(terra::extract(chelsa_stack, coords))
ex_bio_df$ID <- NULL
ex_bio_df[] <- lapply(ex_bio_df, function(col) as.numeric(unlist(col)))

valores_por_registro <- occ %>%
  select(species, decimalLongitude, decimalLatitude) %>%
  bind_cols(as_tibble(ex_bio_df)) %>%
  mutate(
    tri = ex_tri,
    bio1_rugosidad_local = ex_rough1,
    bio12_rugosidad_local = ex_rough2
  )

write_csv(valores_por_registro, file.path(out_dir, "valores_ambientales_por_registro.csv"))

write_csv(valores_por_registro, file.path(out_dir, "valores_ambientales_por_registro.csv"))

# ------------------------------------------------------------------
# 6. RESUMEN POR ESPECIE (para unir con volumen/distancias del elipsoide)
# ------------------------------------------------------------------
heterogeneidad_por_especie <- valores_por_registro %>%
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

message("\n✅ Listo.")
message("   Valores por registro (para selección de variables del elipsoide): ",
        file.path(out_dir, "valores_ambientales_por_registro.csv"))
message("   Heterogeneidad por especie (TRI + rugosidad climática):           ",
        file.path(out_dir, "heterogeneidad_por_especie.csv"))
