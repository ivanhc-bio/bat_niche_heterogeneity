#####################################################################
# Reconstruye los respaldos "ANTES_de_revision_outliers" para que SÍ
# contengan los 2 registros excluidos como outliers ambientales
# (Diphylla ecaudata: gbifID 2630373568; Balantiopteryx plicata:
# gbifID 1211822249), SIN el registro de cautiverio.
#
# Motivo: los respaldos actuales tienen 11,516 filas (= dataset final),
# así que la prueba "con outliers" de robustez_v7.R no reincorporaba nada.
#
# Uso: 1) corre este script;  2) en robustez_v7.R pon
#      CORRER <- c(n30=FALSE, focales16=FALSE, outliers=TRUE, rango=FALSE, ses=FALSE)
#      (ajusta a los nombres que tenga tu vector) y córrelo.
#####################################################################
pkgs <- c("terra", "dplyr", "readr")
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "D:/capitulo2/reanalyses_17sep"
res_dir  <- file.path(base_dir, "resultados_por_especie")
out_dir  <- file.path(base_dir, "resultadosfinales")
tri_path <- file.path(base_dir, "topografia", "tri_1KMmd_GMTEDmd.tif")
# Archivos "final_manual" por especie (post-QGIS, ANTES de excluir los outliers)
f_dip <- file.path(res_dir, "Diphylla_ecaudata_final_manual.csv")
f_bal <- file.path(res_dir, "Balantiopteryx_plicata_final_manual.csv")

occ_csv <- file.path(res_dir, "TODAS_LAS_ESPECIES_final_manual.csv")          # final, 11,516 filas
val_csv <- file.path(out_dir, "valores_ambientales_por_registro.csv")          # final
occ_antes_csv <- file.path(res_dir, "TODAS_LAS_ESPECIES_final_manual_ANTES_de_revision_outliers.csv")
val_antes_csv <- file.path(out_dir, "valores_ambientales_por_registro_ANTES_de_revision_outliers.csv")

lee <- function(f) readr::read_csv(f, col_types = readr::cols(.default = "c"), show_col_types = FALSE)
occ <- lee(occ_csv); val <- lee(val_csv)
stopifnot(nrow(occ) == 11516, nrow(val) == 11516)

# --- 1. Los 2 registros extra -------------------------------------
# Los CSV por especie no traen la columna `species`: se asigna desde el archivo
extras <- bind_rows(lee(f_dip) %>% mutate(species = "Diphylla ecaudata"),
                    lee(f_bal) %>% mutate(species = "Balantiopteryx plicata")) %>%
  filter(!gbifID %in% occ$gbifID)
print(extras %>% select(gbifID, species, decimalLatitude, decimalLongitude, locality, elevation_extracted))
stopifnot(nrow(extras) == 2, !any(grepl("Zoo", extras$locality)))
extras <- extras %>% select(any_of(names(occ)))

# --- 2. Variables ambientales en esos 2 puntos ----------------------
terra::setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
terra::setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
lon <- as.numeric(extras$decimalLongitude); lat <- as.numeric(extras$decimalLatitude)
ext_p <- terra::ext(min(lon) - 0.25, max(lon) + 0.25, min(lat) - 0.25, max(lat) + 0.25)
base_url <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"
capas <- lapply(1:19, function(i) {
  r <- terra::crop(terra::rast(sprintf("/vsicurl/%s/CHELSA_bio%d_1981-2010_V.2.1.tif", base_url, i)), ext_p)
  names(r) <- paste0("bio", i); r })
chelsa <- terra::rast(capas)
tri <- terra::crop(terra::rast(tri_path), ext_p); names(tri) <- "tri"
r1  <- terra::focal(chelsa[["bio1"]],  w = 3, fun = "sd", na.rm = TRUE)
r12 <- terra::focal(chelsa[["bio12"]], w = 3, fun = "sd", na.rm = TRUE)
co <- cbind(lon, lat)
ex <- as.data.frame(terra::extract(chelsa, co)); ex$ID <- NULL
ex[] <- lapply(ex, function(x) as.numeric(unlist(x)))
ex$tri <- as.numeric(unlist(terra::extract(tri, co)[["tri"]]))
ex$bio1_rugosidad_local  <- as.numeric(unlist(terra::extract(r1,  co)[[1]]))
ex$bio12_rugosidad_local <- as.numeric(unlist(terra::extract(r12, co)[[1]]))
val_extra <- bind_cols(tibble::tibble(species = extras$species, decimalLongitude = extras$decimalLongitude,
                                      decimalLatitude = extras$decimalLatitude), tibble::as_tibble(ex)) %>%
  mutate(across(everything(), as.character))
print(val_extra %>% select(species, bio5, bio6, bio16, bio17, tri, bio12_rugosidad_local))

# --- 3. Respaldos ANTES = final + 2 registros -----------------------
for (f in c(occ_antes_csv, val_antes_csv))
  if (file.exists(f)) file.copy(f, sub("\\.csv$", "_SIN2outliers.csv", f), overwrite = FALSE)
occ_antes <- bind_rows(occ, extras)
val_antes <- bind_rows(val, val_extra[, names(val)])
stopifnot(nrow(occ_antes) == 11518, nrow(val_antes) == 11518)
readr::write_csv(occ_antes, occ_antes_csv)
readr::write_csv(val_antes, val_antes_csv)
message("Listo: respaldos ANTES con 11,518 filas (sin cautiverio, con los 2 outliers).")
