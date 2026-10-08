#####################################################################
# ROBUSTECER HALLAZGO 1, PARTE B -- VERSIÓN v3 (corrige el bug real)
#
# Diagnóstico confirmó: terra::extract(dem, puntos, buffer=...) NO
# aplica un buffer real cuando el raster está en coordenadas
# geográficas (WGS84, grados) -- silenciosamente devuelve solo la
# celda del punto (1 celda por punto, siempre). La ruta que sí
# funciona (confirmada: ~366 celdas por punto con buffer de 10 km,
# consistente con pi*10^2 ~ 314 km^2 a 1km de resolución) es:
# proyectar los puntos a un CRS métrico -> construir el buffer ahí
# -> reproyectar a WGS84 -> extraer.
#####################################################################

pkgs <- c("terra", "dplyr", "readr", "purrr", "tibble")
invisible(lapply(pkgs, library, character.only = TRUE))
select <- dplyr::select

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

set.seed(2026)

BUFFER_KM <- 10
N_NULL    <- 999
CRS_METRICO <- "+proj=cea +lat_ts=0 +lon_0=-75"  # mismo CRS equal-area ya usado para area_rango_km2

dem <- terra::rast(file.path(base_dir, "variables", "elevation.tif"))

ocurrencias <- read_csv(
  file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv"),
  show_col_types = FALSE)

especies <- unique(ocurrencias$species)
message(sprintf("Especies a procesar: %d", length(especies)))

calcular_ses <- function(sp) {
  pts_sp <- ocurrencias %>% filter(species == sp)
  n_sp <- nrow(pts_sp)
  if (n_sp < 5) return(NULL)

  pts_vect <- terra::vect(pts_sp, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

  # Buffer real: proyectar -> construir buffer en metros -> reproyectar
  pts_proj <- terra::project(pts_vect, CRS_METRICO)
  buf_proj <- terra::buffer(pts_proj, width = BUFFER_KM * 1000)
  buf_geo  <- terra::project(buf_proj, "EPSG:4326")

  extraidos <- tryCatch(terra::extract(dem, buf_geo), error = function(e) NULL)
  if (is.null(extraidos)) return(NULL)
  elev_disponible <- extraidos[[2]]
  elev_disponible <- elev_disponible[!is.na(elev_disponible)]
  if (length(elev_disponible) < 10) return(NULL)

  elev_obs <- pts_sp$elevation_extracted
  elev_obs <- elev_obs[!is.na(elev_obs)]
  if (length(elev_obs) < 5) return(NULL)
  rango_obs <- as.numeric(diff(quantile(elev_obs, c(0.025, 0.975))))

  rangos_nulos <- map_dbl(1:N_NULL, function(i) {
    muestra <- sample(elev_disponible, size = n_sp, replace = TRUE)
    as.numeric(diff(quantile(muestra, c(0.025, 0.975))))
  })

  media_nula <- mean(rangos_nulos)
  sd_nula <- sd(rangos_nulos)
  ses <- if (sd_nula > 0) (rango_obs - media_nula) / sd_nula else NA_real_
  p_nula <- mean(rangos_nulos >= rango_obs)

  tibble(species = sp, n_ocurrencias = n_sp, n_celdas_disponibles = length(elev_disponible),
         elev_rango_observado = rango_obs, elev_rango_nulo_media = media_nula,
         elev_rango_nulo_sd = sd_nula, SES_elev_range = ses, p_nulo = p_nula)
}

resultados <- map(especies, function(sp) {
  message(sprintf("  %s", sp))
  tryCatch(calcular_ses(sp), error = function(e) {message("   -> error: ", conditionMessage(e)); NULL})
}) %>% compact() %>% bind_rows()

write_csv(resultados, file.path(out_dir, "null_model_SES_elevrange_v3.csv"))
message(sprintf("\nCompletado: %d / %d especies", nrow(resultados), length(especies)))
message(sprintf("Mediana de celdas disponibles por especie (debe ser MUCHO mayor que n_ocurrencias ahora): %s vs n_ocurrencias mediana %s",
                 median(resultados$n_celdas_disponibles), median(resultados$n_ocurrencias)))
message(sprintf("Especies con rango observado fuera del 95%% nulo (p_nulo<0.05 o >0.95): %d",
                 sum(resultados$p_nulo < 0.05 | resultados$p_nulo > 0.95, na.rm = TRUE)))
