#####################################################################
# ROBUSTEZ, ESPECIFICACIÓN v7 (todas las pruebas en un solo script)
#
# Reproduce las pruebas de robustez previas con EXACTAMENTE la misma
# especificación que los análisis principales v7:
#   log_niche_volume   ~ tri_c + bio12_rugosidad_media + log_n
#   log_dist_media     ~ tri_c + bio12_rugosidad_media + log_n
#   log_dist_sd        ~ tri_c + bio12_rugosidad_media + log_n
#   elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n
#   phylopath: mismos 4 candidatos que phylopath_v7.R
#
# IMPORTANTE: TRI_CENTRO es SIEMPRE la media de TRI de las 53 especies
# de tabla_final_analisis.csv (aunque el subconjunto sea menor), para que
# los coeficientes y el pico de TRI sean comparables con los principales.
#
# Pruebas (activa/desactiva con los interruptores de abajo):
#   A  n30        : excluye especies con n < 30
#   B  focales16  : 16 especies de estrategia elevacional contrastante
#   C  outliers   : devuelve los 2 registros ambientales marcados como outliers
#                   (usa los respaldos YA sin el registro de cautiverio)
#   D  rango      : añade log(area_rango_km2) como covariable
#   E  ses        : TRI -> exceso de rango elevacional sobre el modelo nulo (SES)
#
# Salidas (resultadosfinales/), una pareja por prueba:
#   ROBUSTEZ_v7_<prueba>_pgls_resumen.csv   (estimate, IC entre árboles, SE, % sig, lambda, R2)
#   ROBUSTEZ_v7_<prueba>_phylopath_seleccion.csv (cuando aplica)
#   ROBUSTEZ_v7_<prueba>_pico_TRI.csv       (pico de TRI del modelo de rango elevacional)
# Y al final:
#   ROBUSTEZ_v7_comparacion.csv   principal vs. cada prueba (estimate y % sig)
#
# Tiempo aproximado: ~ lo que tardó pgls_v7_principal.R + phylopath_v7.R
# por cada prueba (A-C son las más lentas). Puedes correr una a la vez.
#####################################################################
pkgs <- c("ape", "caper", "phylopath", "dplyr", "readr", "tibble", "tidyr")
faltantes <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(faltantes) > 0) install.packages(faltantes, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))
select <- dplyr::select

# ---------------- INTERRUPTORES ----------------
CORRER <- c(n30 = TRUE, focales16 = TRUE, outliers = TRUE, rango = TRUE, ses = TRUE)
# -----------------------------------------------

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

equivalencias <- read_csv(file.path(out_dir, "equivalencias_especies_arbol.csv"), show_col_types = FALSE)
datos <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)
TRI_CENTRO <- mean(datos$tri_media)
message(sprintf("TRI_CENTRO = %.4f (debe coincidir con el de pgls_v7_principal.R)", TRI_CENTRO))

arboles <- ape::read.nexus(file.path(base_dir, "trees", "output.nex"))

# ------------------------------------------------------------------
# Utilidades
# ------------------------------------------------------------------
preparar <- function(df) {
  df %>%
    left_join(equivalencias %>% select(species, tip_label), by = "species") %>%
    mutate(log_niche_volume = log(niche_volume),
           log_dist_media   = log(dist_centroide_media),
           log_dist_sd      = log(dist_centroide_sd),
           log_n            = log(n),
           tri_c            = tri_media - TRI_CENTRO,
           tri_c2           = tri_c^2) %>%
    rename(tiplabel = tip_label) %>%
    filter(!is.na(tiplabel), !is.na(niche_volume), !is.na(elev_rango_robusto)) %>%
    as.data.frame()
}

MODELOS_V7 <- list(
  log_niche_volume   = log_niche_volume   ~ tri_c + bio12_rugosidad_media + log_n,
  log_dist_media     = log_dist_media     ~ tri_c + bio12_rugosidad_media + log_n,
  log_dist_sd        = log_dist_sd        ~ tri_c + bio12_rugosidad_media + log_n,
  elev_rango_robusto = elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n
)

# PGLS sobre las 500 filogenias; devuelve resumen y (si existe el modelo) pico de TRI
correr_pgls <- function(dat, modelos, etiqueta, pico_modelo = "elev_rango_robusto") {
  message(sprintf("\n=== PGLS [%s]: %d especies ===", etiqueta, nrow(dat)))
  res <- list(); picos <- list()
  for (i in seq_along(arboles)) {
    arbol_i <- arboles[[i]]; arbol_i$node.label <- NULL
    cd <- tryCatch(caper::comparative.data(phy = arbol_i, data = dat, names.col = "tiplabel",
                                            vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE),
                   error = function(e) NULL)
    if (is.null(cd)) next
    for (nm in names(modelos)) {
      fit <- tryCatch(caper::pgls(modelos[[nm]], data = cd, lambda = "ML"), error = function(e) NULL)
      if (is.null(fit)) next
      s <- summary(fit); co <- as.data.frame(s$coefficients); co$termino <- rownames(co)
      res[[length(res) + 1]] <- co %>%
        transmute(arbol = i, modelo = nm, termino, estimate = Estimate, se = `Std. Error`,
                  p_valor = `Pr(>|t|)`, lambda = fit$param["lambda"], r2 = s$r.squared)
      if (nm == pico_modelo && all(c("tri_c", "tri_c2") %in% names(coef(fit)))) {
        b1 <- coef(fit)["tri_c"]; b2 <- coef(fit)["tri_c2"]
        picos[[length(picos) + 1]] <- tibble(arbol = i, b1 = b1, b2 = b2,
                                             TRI_pico = ifelse(b2 < 0, TRI_CENTRO - b1 / (2 * b2), NA_real_))
      }
    }
    if (i %% 100 == 0) message(sprintf("  PGLS %s: %d/%d", etiqueta, i, length(arboles)))
  }
  resumen <- bind_rows(res) %>% filter(termino != "(Intercept)") %>%
    group_by(modelo, termino) %>%
    summarise(n_arboles = n(),
              estimate_mediana = median(estimate),
              estimate_p2_5 = quantile(estimate, 0.025), estimate_p97_5 = quantile(estimate, 0.975),
              se_mediana = median(se), p_mediana = median(p_valor),
              pct_significativos = mean(p_valor < 0.05) * 100,
              lambda_mediana = median(lambda), r2_mediano = median(r2), .groups = "drop") %>%
    arrange(modelo, termino)
  write_csv(resumen, file.path(out_dir, sprintf("ROBUSTEZ_v7_%s_pgls_resumen.csv", etiqueta)))
  if (length(picos) > 0) {
    pk <- bind_rows(picos)
    write_csv(pk, file.path(out_dir, sprintf("ROBUSTEZ_v7_%s_pico_TRI.csv", etiqueta)))
    message(sprintf("  Pico de TRI [%s]: mediana %.2f (IC entre árboles %.2f - %.2f); %d árboles con b2<0",
                    etiqueta, median(pk$TRI_pico, na.rm = TRUE),
                    quantile(pk$TRI_pico, 0.025, na.rm = TRUE), quantile(pk$TRI_pico, 0.975, na.rm = TRUE),
                    sum(pk$b2 < 0)))
  }
  print(resumen %>% select(modelo, termino, estimate_mediana, se_mediana, pct_significativos,
                            lambda_mediana, r2_mediano), n = Inf)
  resumen
}

# Phylopath v7 (mismos 4 candidatos que phylopath_v7.R)
correr_phylopath <- function(dat, etiqueta, extra_vol = NULL) {
  message(sprintf("\n=== PHYLOPATH [%s]: %d especies ===", etiqueta, nrow(dat)))
  rownames(dat) <- dat$tiplabel
  comunes <- list(elev_rango_robusto ~ tri_c + tri_c2 + log_n,
                  bio12_rugosidad_media ~ tri_c,
                  tri_c2 ~ tri_c)
  mods <- define_model_set(
    mediacion_completa   = c(comunes, log_niche_volume ~ elev_rango_robusto + bio12_rugosidad_media + log_n),
    sin_mediacion        = c(comunes, log_niche_volume ~ tri_c + tri_c2 + bio12_rugosidad_media + log_n),
    mediacion_parcial    = c(comunes, log_niche_volume ~ elev_rango_robusto + tri_c + tri_c2 + bio12_rugosidad_media + log_n),
    rutas_independientes = c(comunes, log_niche_volume ~ bio12_rugosidad_media + log_n)
  )
  ganador <- character(length(arboles))
  for (i in seq_along(arboles)) {
    arbol_i <- arboles[[i]]; arbol_i$node.label <- NULL
    aj <- tryCatch(phylo_path(mods, data = dat, tree = arbol_i, model = "lambda"), error = function(e) NULL)
    if (is.null(aj)) next
    s <- tryCatch(summary(aj), error = function(e) NULL)
    if (is.null(s)) next
    ganador[i] <- as.character(s$model[1])
    if (i %% 100 == 0) message(sprintf("  phylopath %s: %d/%d", etiqueta, i, length(arboles)))
  }
  sel <- tibble(modelo_ganador = ganador) %>% filter(modelo_ganador != "") %>%
    count(modelo_ganador, name = "n_arboles") %>%
    mutate(pct = round(100 * n_arboles / sum(n_arboles), 1)) %>% arrange(desc(n_arboles))
  write_csv(sel, file.path(out_dir, sprintf("ROBUSTEZ_v7_%s_phylopath_seleccion.csv", etiqueta)))
  print(sel)
  sel
}

# ------------------------------------------------------------------
# A. n >= 30
# ------------------------------------------------------------------
if (CORRER["n30"]) {
  d <- preparar(datos) %>% filter(n >= 30)
  message(sprintf("n>=30: %d especies", nrow(d)))
  correr_pgls(d, MODELOS_V7, "n30")
  correr_phylopath(d, "n30")
}

# ------------------------------------------------------------------
# B. 16 especies focales
#    (con n = 16 el modelo de rango elevacional tiene 6 coeficientes:
#     es una prueba de CONSISTENCIA de dirección, no de potencia)
# ------------------------------------------------------------------
if (CORRER["focales16"]) {
  especies_16 <- tibble::tribble(
    ~species,                  ~estrategia,
    "Noctilio leporinus",      "Lowland-narrow",
    "Pteronotus gymnonotus",   "Lowland-narrow",
    "Peropteryx macrotis",     "Lowland-narrow",
    "Artibeus phaeotis",       "Lowland-narrow",
    "Glossophaga mutica",      "Lowland-wide",
    "Myotis riparius",         "Lowland-wide",
    "Desmodus rotundus",       "Lowland-wide",
    "Sturnira parvidens",      "Lowland-wide",
    "Eptesicus brasiliensis",  "Montane-wide",
    "Lasiurus blossevillii",   "Montane-wide",
    "Enchisthenes hartii",     "Montane-wide",
    "Sturnira hondurensis",    "Montane-wide",
    "Artibeus aztecus",        "Montane-narrow",
    "Anoura peruana",          "Montane-narrow",
    "Sturnira erythromos",     "Montane-narrow",
    "Myotis oxyotus",          "Montane-narrow"
  )
  d <- preparar(datos) %>% inner_join(especies_16, by = "species")
  message(sprintf("Focales: %d / 16", nrow(d)))
  if (nrow(d) < 16) message("⚠ Faltan: ", paste(setdiff(especies_16$species, d$species), collapse = ", "))
  correr_pgls(d, MODELOS_V7, "focales16")
  correr_phylopath(d, "focales16")
}

# ------------------------------------------------------------------
# C. Outliers ambientales de vuelta (2 registros) -- SIN el registro de cautiverio
#    Requiere los respaldos *_ANTES_de_revision_outliers.csv ya limpiados con
#    quitar_cautiverio_de_respaldos_outliers.R
# ------------------------------------------------------------------
if (CORRER["outliers"]) {
  suppressPackageStartupMessages(library(ntbox))
  valores_antes <- read_csv(file.path(out_dir, "valores_ambientales_por_registro_ANTES_de_revision_outliers.csv"),
                            show_col_types = FALSE)
  occ_antes <- read_csv(file.path(base_dir, "resultados_por_especie",
                                  "TODAS_LAS_ESPECIES_final_manual_ANTES_de_revision_outliers.csv"),
                        show_col_types = FALSE)
  heterogeneidad_antes <- valores_antes %>% group_by(species) %>%
    summarise(n = n(), tri_media = mean(tri, na.rm = TRUE),
              bio1_rugosidad_media = mean(bio1_rugosidad_local, na.rm = TRUE),
              bio12_rugosidad_media = mean(bio12_rugosidad_local, na.rm = TRUE), .groups = "drop")
  elevacion_antes <- occ_antes %>% filter(!is.na(elevation_extracted)) %>% group_by(species) %>%
    summarise(elev_mediana = median(elevation_extracted),
              elev_rango_robusto = quantile(elevation_extracted, 0.975) - quantile(elevation_extracted, 0.025),
              .groups = "drop")
  vars_finales <- c("bio5", "bio6", "bio16", "bio17")
  set.seed(2026)
  re <- list()
  for (sp in sort(unique(valores_antes$species))) {
    dsp <- valores_antes %>% filter(species == sp) %>% select(all_of(vars_finales)) %>% na.omit()
    if (nrow(dsp) <= length(vars_finales)) next
    r <- tryCatch(ntbox::cov_center(data = dsp, mve = TRUE, level = 0.95, vars = vars_finales), error = function(e) NULL)
    if (is.null(r)) next
    dd <- ntbox::inEllipsoid(centroid = r$centroid, eShape = r$covariance, env_data = dsp, level = 0.95)$mh_dist
    re[[sp]] <- tibble(species = sp, niche_volume = r$niche_volume,
                       dist_centroide_media = mean(dd), dist_centroide_sd = sd(dd))
  }
  tabla_antes <- bind_rows(re) %>% left_join(heterogeneidad_antes, by = "species") %>%
    left_join(elevacion_antes, by = "species")
  write_csv(tabla_antes, file.path(out_dir, "ROBUSTEZ_v7_outliers_tabla.csv"))
  d <- preparar(tabla_antes)
  message(sprintf("Con outliers: %d especies; n total = %d (reportado: %d)",
                  nrow(d), sum(d$n), sum(datos$n)))
  correr_pgls(d, MODELOS_V7, "outliers")
  correr_phylopath(d, "outliers")
}

# ------------------------------------------------------------------
# D. Tamaño del rango geográfico como covariable
# ------------------------------------------------------------------
if (CORRER["rango"]) {
  d <- preparar(datos) %>% filter(!is.na(area_rango_km2)) %>% mutate(log_area_rango = log(area_rango_km2))
  mods_rango <- list(
    elev_con_rango = elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n + log_area_rango,
    vol_con_rango  = log_niche_volume   ~ tri_c + bio12_rugosidad_media + log_n + log_area_rango
  )
  correr_pgls(d, mods_rango, "rango", pico_modelo = "elev_con_rango")
}

# ------------------------------------------------------------------
# E. SES del rango elevacional (modelo nulo geométrico)
#    Requiere null_model_SES_elevrange_v3.csv (de robustez_null_model_geometrico_v3.R)
# ------------------------------------------------------------------
if (CORRER["ses"]) {
  ses <- read_csv(file.path(out_dir, "null_model_SES_elevrange_v3.csv"), show_col_types = FALSE)
  d <- preparar(datos) %>%
    left_join(ses %>% select(species, SES_elev_range, p_nulo), by = "species") %>%
    filter(!is.na(SES_elev_range))
  message(sprintf("Especies con SES: %d", nrow(d)))
  mods_ses <- list(
    SES_lineal    = SES_elev_range ~ tri_c + bio12_rugosidad_media + log_n,
    SES_cuadratico = SES_elev_range ~ tri_c + tri_c2 + bio12_rugosidad_media + log_n
  )
  correr_pgls(d, mods_ses, "ses", pico_modelo = "SES_cuadratico")
  message(sprintf("Especies con rango observado > nulo (p<0.05): %d / %d",
                  sum(d$p_nulo < 0.05, na.rm = TRUE), nrow(d)))
}

# ------------------------------------------------------------------
# Comparación final: principal vs. cada prueba
# ------------------------------------------------------------------
principal <- read_csv(file.path(out_dir, "pgls_v7_resumen.csv"), show_col_types = FALSE) %>%
  filter(!grepl("^(D_|E_)", modelo)) %>%
  select(modelo, termino, estimate_principal = estimate_mediana, se_principal = se_mediana,
         pct_sig_principal = pct_significativos, lambda_principal = lambda_mediana)
comp <- principal
for (et in c("n30", "focales16", "outliers", "rango")) {   # incluye pruebas corridas antes (si existe el archivo)
  f <- file.path(out_dir, sprintf("ROBUSTEZ_v7_%s_pgls_resumen.csv", et))
  if (!file.exists(f)) next
  r <- read_csv(f, show_col_types = FALSE) %>%
    select(modelo, termino, estimate_mediana, se_mediana, pct_significativos, lambda_mediana)
  names(r)[3:6] <- paste0(c("estimate_", "se_", "pct_sig_", "lambda_"), et)
  comp <- comp %>% left_join(r, by = c("modelo", "termino"))
}
write_csv(comp, file.path(out_dir, "ROBUSTEZ_v7_comparacion.csv"))
message("\n✅ Listo. Archivos ROBUSTEZ_v7_* en resultadosfinales/")
print(comp, n = Inf, width = Inf)
