#####################################################################
# Revisión visual interactiva de registros GBIF thinned, especie por
# especie, contra tu conocimiento de la distribución (bibliografía).
#
# Flujo: para cada especie se abre un mapa interactivo (leaflet).
# Haces clic sobre los puntos que quieras ELIMINAR (fuera del rango
# conocido), das clic en "Done", y el script filtra y guarda la base
# final + un log de cuántos puntos se quitaron por especie.
#####################################################################

pkgs <- c("mapedit", "mapview", "leaflet", "sf", "dplyr", "readr", "shiny")
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed], dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# ------------------------------------------------------------------
# CONFIGURACIÓN
# ------------------------------------------------------------------
base_dir  <- "D:/capitulo2/reanalyses_17sep"
thin_dir  <- file.path(base_dir, "gbif_thinned")        # input: CSVs ya thinned
final_dir <- file.path(base_dir, "gbif_final_revisado") # output: CSVs revisados
log_file  <- file.path(base_dir, "log_revision_visual.csv")

dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

# Especies a revisar (por defecto, todas las que tengan CSV thinned).
# Si quieres revisar solo algunas primero, reemplaza por un vector manual,
# ej: species_to_review <- c("Artibeus aztecus", "Dasypterus ega")
species_to_review <- list.files(thin_dir, pattern = "_thinned10km\\.csv$") %>%
  gsub("_thinned10km\\.csv$", "", .) %>%
  gsub("_", " ", .) %>%
  sort()

# Si ya revisaste algunas especies en una sesión anterior, no las repitas
ya_revisadas <- list.files(final_dir, pattern = "_final\\.csv$") %>%
  gsub("_final\\.csv$", "", .) %>%
  gsub("_", " ", .)
species_to_review <- setdiff(species_to_review, ya_revisadas)

message(sprintf("Especies pendientes de revisión: %d", length(species_to_review)))

log_list <- list()

# ------------------------------------------------------------------
# LOOP DE REVISIÓN (uno a la vez, es inherentemente manual)
# ------------------------------------------------------------------
for (sp in species_to_review) {

  f <- file.path(thin_dir, paste0(gsub(" ", "_", sp), "_thinned10km.csv"))
  if (!file.exists(f)) next

  dat <- read_csv(f, show_col_types = FALSE) %>%
    mutate(.orig_id = row_number())

  pts <- st_as_sf(dat, coords = c("decimalLongitude", "decimalLatitude"),
                   crs = 4326, remove = FALSE)

  cat("\n==================================================\n")
  cat(sprintf("Especie: %s  (%d registros)\n", sp, nrow(pts)))
  cat("En el visor: haz clic en los puntos a ELIMINAR (fuera del\n")
  cat("rango conocido). Puedes cambiar el mapa base (esquina sup.\n")
  cat("derecha) a OpenTopoMap para ver relieve/contexto. Cuando\n")
  cat("termines, da clic en 'Done'.\n")

  seleccion <- mapedit::selectFeatures(
    pts,
    viewer = shiny::paneViewer()
  )

  ids_eliminar <- if (nrow(seleccion) > 0) seleccion$.orig_id else integer(0)

  dat_final <- dat %>%
    filter(!(.orig_id %in% ids_eliminar)) %>%
    select(-.orig_id)

  write_csv(dat_final, file.path(final_dir, paste0(gsub(" ", "_", sp), "_final.csv")))

  log_list[[sp]] <- tibble(
    species        = sp,
    n_antes_revision = nrow(dat),
    n_eliminados     = length(ids_eliminar),
    n_final          = nrow(dat_final)
  )

  cat(sprintf(">> %s: eliminados %d de %d -> quedan %d\n",
              sp, length(ids_eliminar), nrow(dat), nrow(dat_final)))

  # Guarda el log incrementalmente por si quieres pausar a mitad de sesión
  write_csv(bind_rows(log_list), log_file)
}

message("\n✅ Revisión visual completa (o pausada). Log en: ", log_file)
message("   CSVs finales en: ", final_dir)
