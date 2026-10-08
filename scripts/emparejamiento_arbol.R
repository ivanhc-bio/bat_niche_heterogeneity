#####################################################################
# Emparejamiento de nombres: tus 53 especies (tabla_final_analisis.csv)
# vs. los tip labels de la muestra de 500 árboles (VertLife/Upham et al.
# 2019, formato Genus_species).
#
# La taxonomía de VertLife (basada en MSW3) puede diferir de la de
# GBIF/backbone en algunos géneros -- p. ej. ya sabemos que usa
# "Lasiurus_ega" en vez de "Dasypterus_ega", y "Dermanura_aztecus" en
# vez de "Artibeus_aztecus" (justo lo opuesto al backbone de GBIF).
#####################################################################

pkgs <- c("ape", "dplyr", "readr", "stringr")
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) install.packages(pkgs[!installed], dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

# ------------------------------------------------------------------
# 1. LEER LOS ÁRBOLES Y SACAR LOS TIP LABELS (se asume el mismo
#    conjunto de puntas en los 500 árboles -- solo cambia topología/
#    longitudes de rama entre muestras posteriores)
# ------------------------------------------------------------------
arboles <- ape::read.nexus(file.path(base_dir, "trees", "output.nex"))
message(sprintf("Árboles leídos: %d", length(arboles)))

tips_arbol <- sort(unique(arboles[[1]]$tip.label))
message(sprintf("Puntas (especies) en el árbol: %d", length(tips_arbol)))

tips_arbol_espacio <- str_replace(tips_arbol, "_", " ")

# ------------------------------------------------------------------
# 2. TUS 53 ESPECIES
# ------------------------------------------------------------------
datos <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)
mis_especies <- sort(unique(datos$species))
message(sprintf("Tus especies: %d", length(mis_especies)))

# ------------------------------------------------------------------
# 3. EQUIVALENCIAS YA CONOCIDAS (de la resolución de sinonimias en
#    GBIF -- ojo, el árbol a veces usa la convención CONTRARIA)
# ------------------------------------------------------------------
equivalencias_conocidas <- tibble::tribble(
  ~species,            ~tip_label,
  "Dasypterus ega",     "Lasiurus_ega",
  "Artibeus aztecus",   "Dermanura_aztecus"
)

# ------------------------------------------------------------------
# 4. EMPAREJAMIENTO AUTOMÁTICO (match exacto directo)
# ------------------------------------------------------------------
match_directo <- tibble(species = mis_especies) %>%
  mutate(tip_label_directo = str_replace(species, " ", "_")) %>%
  mutate(en_arbol = tip_label_directo %in% tips_arbol)

resueltas_directo <- match_directo %>% filter(en_arbol) %>%
  transmute(species, tip_label = tip_label_directo, metodo = "match_directo")

pendientes <- match_directo %>% filter(!en_arbol) %>% pull(species)

# Aplicamos las equivalencias ya conocidas a las pendientes
resueltas_conocidas <- equivalencias_conocidas %>%
  filter(species %in% pendientes) %>%
  mutate(metodo = "equivalencia_conocida")

pendientes <- setdiff(pendientes, resueltas_conocidas$species)

# ------------------------------------------------------------------
# 5. PISTAS PARA LAS QUE SIGUEN PENDIENTES: buscamos tips del árbol
#    que compartan género con la especie pendiente, para que sea
#    fácil decidir manualmente
# ------------------------------------------------------------------
if (length(pendientes) > 0) {
  message(sprintf("\n⚠ %d especie(s) sin resolver automáticamente:\n", length(pendientes)))

  pistas <- map_dfr(pendientes, function(sp) {
    genero <- str_split(sp, " ")[[1]][1]
    candidatos <- tips_arbol[str_starts(tips_arbol, paste0(genero, "_"))]
    # también busca por el epíteto específico en cualquier género del árbol,
    # por si el género cambió pero el epíteto se conservó
    epiteto <- str_split(sp, " ")[[1]][2]
    candidatos_epiteto <- tips_arbol[str_ends(tips_arbol, paste0("_", epiteto))]
    todos <- union(candidatos, candidatos_epiteto)
    tibble(species = sp,
           candidatos_en_arbol = if (length(todos) > 0) paste(todos, collapse = " | ") else "NINGUNO ENCONTRADO")
  })

  print(pistas, n = Inf)
  write_csv(pistas, file.path(out_dir, "pistas_emparejamiento_pendiente.csv"))
  message("\nPistas guardadas en: ", file.path(out_dir, "pistas_emparejamiento_pendiente.csv"))
} else {
  pistas <- tibble(species = character(), candidatos_en_arbol = character())
}

# ------------------------------------------------------------------
# 6. TABLA DE EQUIVALENCIAS (para revisar y completar a mano las
#    filas sin tip_label, usando la columna de pistas)
# ------------------------------------------------------------------
tabla_equivalencias <- bind_rows(resueltas_directo, resueltas_conocidas) %>%
  bind_rows(tibble(species = pendientes, tip_label = NA_character_, metodo = "PENDIENTE_REVISAR")) %>%
  arrange(species)

write_csv(tabla_equivalencias, file.path(out_dir, "equivalencias_especies_arbol.csv"))

message(sprintf("\n✅ Resueltas automáticamente: %d / %d",
                 sum(tabla_equivalencias$metodo != "PENDIENTE_REVISAR"), length(mis_especies)))
message("   Tabla guardada en: ", file.path(out_dir, "equivalencias_especies_arbol.csv"))
if (length(pendientes) > 0) {
  message("   Completa manualmente la columna 'tip_label' para las filas PENDIENTE_REVISAR,")
  message("   usando las pistas de 'pistas_emparejamiento_pendiente.csv', y mándamelas.")
}


library(dplyr)
library(readr)

out_dir <- "D:/capitulo2/reanalyses_17sep/resultadosfinales"

equivalencias <- read_csv(file.path(out_dir, "equivalencias_especies_arbol.csv"), show_col_types = FALSE)

confirmaciones <- tibble::tribble(
  ~species,            ~tip_label,
  "Artibeus phaeotis",  "Dermanura_phaeotis",
  "Artibeus watsoni",   "Dermanura_watsoni",
  "Glossophaga mutica", "Glossophaga_soricina"
)

equivalencias_actualizada <- equivalencias %>%
  rows_update(confirmaciones %>% mutate(metodo = "confirmado_manual"), by = "species")

write_csv(equivalencias_actualizada, file.path(out_dir, "equivalencias_especies_arbol.csv"))

# Verificación: ¿queda alguna especie sin tip_label?
pendientes_finales <- equivalencias_actualizada %>% filter(is.na(tip_label))
print(pendientes_finales)
message(sprintf("Especies resueltas: %d / %d",
                sum(!is.na(equivalencias_actualizada$tip_label)), nrow(equivalencias_actualizada)))


equivalencias_final <- equivalencias_actualizada %>%
  rows_update(
    tibble(species = "Anoura peruana", tip_label = "Anoura_geoffroyi", metodo = "confirmado_manual"),
    by = "species"
  )

write_csv(equivalencias_final, file.path(out_dir, "equivalencias_especies_arbol.csv"))

stopifnot(all(!is.na(equivalencias_final$tip_label)))
message("✅ 53/53 especies emparejadas con el árbol.")
print(equivalencias_final %>% count(metodo))
