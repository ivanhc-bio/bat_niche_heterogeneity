#####################################################################
# Auditoría de registros potencialmente cautivos / fuera de contexto
# silvestre en la base FINAL (todas las especies).
#
# Por qué: basisOfRecord == "PRESERVED_SPECIMEN" NO excluye animales que
# vivieron en cautiverio (un espécimen de museo de un animal de zoológico
# sigue siendo PRESERVED_SPECIMEN). Este script marca, sin borrar nada,
# los registros que requieren decisión humana, y opcionalmente los quita.
#
# Salidas (en resultadosfinales/):
#   auditoria_registros_marcados.csv   <- revisar a mano (columna `decision`)
#   resumen_auditoria_por_especie.csv
# Paso 2 (opcional, tras revisar): poner EXCLUIR <- TRUE para generar
#   TODAS_LAS_ESPECIES_final_manual_SIN_marcados.csv
#####################################################################
library(dplyr); library(readr); library(stringr)

base_dir <- "D:/capitulo2/reanalyses_17sep"
occ_csv  <- file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv")
out_dir  <- file.path(base_dir, "resultadosfinales")
EXCLUIR  <- FALSE          # TRUE tras revisar auditoria_registros_marcados.csv
umbral_incertidumbre_m <- 10000   # = distancia de thinning (10 km)

occ <- read_csv(occ_csv, show_col_types = FALSE)

# Texto libre donde aparece el contexto del registro (SIMPLE_CSV no trae occurrenceRemarks)
txt <- tibble::tibble(gbifID = occ$gbifID)
txt$texto <- str_c(coalesce(occ$locality, ""), " | ", coalesce(occ$recordedBy, ""), " | ",
                   coalesce(occ$typeStatus, ""))

pat_cautiverio <- regex(paste(
  "\\bzoo\\b", "zool[oó]gico", "zoological (garden|park)", "parque zool",
  "aquarium", "acuario", "captiv", "cautiv", "\\bcage", "jaula", "vivari",
  "mascota", "\\bpet\\b", "rescue", "rescate", "rehabilit", "criadero",
  "zoocriadero", "bioterio", "laboratory colony", "born in", "bred in",
  sep = "|"), ignore_case = TRUE)

# Excepciones conocidas de falsos positivos ("Museu Zoologia", "Pampa de Hospital", etc.)
pat_falsos <- regex("Museu Zoologia|Pampa de Hospital|Hospital de Satipu|Santuario del Manat|Santuario de Fauna y Flora|Santuario Nacional|Wildlife Sanctuary|Ecoparque", ignore_case = TRUE)

marcas <- occ %>%
  mutate(
    flag_texto_cautiverio = str_detect(txt$texto, pat_cautiverio) & !str_detect(txt$texto, pat_falsos),
    flag_establecimiento  = !is.na(establishmentMeans) & !tolower(establishmentMeans) %in% c("native", "native; reintroduced"),
    flag_pais_no_neotropical = countryCode %in% c("US", "CA") | decimalLatitude > 33,  # ver nota Nearctic
    flag_incertidumbre    = !is.na(coordinateUncertaintyInMeters) & coordinateUncertaintyInMeters > umbral_incertidumbre_m
  ) %>%
  mutate(motivo = str_c(
    if_else(flag_texto_cautiverio, "texto_cautiverio;", ""),
    if_else(flag_establecimiento, "establishmentMeans;", ""),
    if_else(flag_pais_no_neotropical, "fuera_Neotropico(US/CA);", ""),
    if_else(flag_incertidumbre, "incertidumbre>10km;", ""))) %>%
  filter(motivo != "") %>%
  select(gbifID, species, countryCode, decimalLatitude, decimalLongitude, locality, recordedBy,
         institutionCode, catalogNumber, establishmentMeans, coordinateUncertaintyInMeters, year, motivo) %>%
  mutate(decision = NA_character_)   # llena: "excluir" / "mantener"

write_csv(marcas, file.path(out_dir, "auditoria_registros_marcados.csv"))
resumen <- marcas %>% mutate(m = str_detect(motivo, "texto_cautiverio")) %>%
  group_by(species) %>%
  summarise(n_marcados = n(), n_texto_cautiverio = sum(m),
            n_fuera_neotropico = sum(str_detect(motivo, "fuera_Neo")),
            n_incert_10km = sum(str_detect(motivo, "incertidumbre")), .groups = "drop") %>%
  arrange(desc(n_marcados))
write_csv(resumen, file.path(out_dir, "resumen_auditoria_por_especie.csv"))
message(sprintf("Marcados: %d registros de %d (texto cautiverio: %d; fuera Neotrópico: %d; incertidumbre > %d m: %d)",
        nrow(marcas), nrow(occ), sum(str_detect(marcas$motivo, "texto_cautiverio")),
        sum(str_detect(marcas$motivo, "fuera_Neo")), umbral_incertidumbre_m,
        sum(str_detect(marcas$motivo, "incertidumbre"))))

if (EXCLUIR) {
  rev <- read_csv(file.path(out_dir, "auditoria_registros_marcados.csv"), show_col_types = FALSE)
  quitar <- rev %>% filter(tolower(decision) == "excluir") %>% pull(gbifID)
  nuevo <- occ %>% filter(!gbifID %in% quitar)
  write_csv(nuevo, sub("\\.csv$", "_SIN_marcados.csv", occ_csv))
  message(sprintf("Excluidos %d; quedan %d. Después repite: valores ambientales -> elipsoides -> PGLS.", length(quitar), nrow(nuevo)))
}
