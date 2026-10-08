#####################################################################
# SENSIBILIDAD A CALIDAD DE COORDENADAS Y LÍMITE GEOGRÁFICO (v7)
#
# Misma especificación que robustez_v7.R (mismos modelos, TRI_CENTRO,
# 500 árboles, set.seed(2026) en los elipsoides). Rehace TODO desde los
# registros: heterogeneidad, elevación y elipsoides MVE por especie.
#
#   incert10km : sin registros con coordinateUncertaintyInMeters > 10 000
#                (registros SIN dato de incertidumbre se conservan)
#   sinTexas   : sin los registros de Dasypterus ega fuera del Neotrópico
#                (countryCode US)
#   ambos      : las dos exclusiones a la vez
#
# Salidas (resultadosfinales/): ROBUSTEZ_v7_<variante>_pgls_resumen.csv,
#   _phylopath_seleccion.csv, _pico_TRI.csv, _tabla.csv  y
#   ROBUSTEZ_v7_sensibilidad_comparacion.csv (principal vs. variantes).
# Interruptores abajo. Orden de duración: igual que cada prueba de robustez_v7.R.
#####################################################################
pkgs <- c("ape", "caper", "phylopath", "dplyr", "readr", "tibble", "tidyr")
faltantes <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(faltantes) > 0) install.packages(faltantes, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))
select <- dplyr::select


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
# Datos por registro (alineados con valores_ambientales_por_registro.csv)
# ------------------------------------------------------------------
suppressPackageStartupMessages(library(ntbox))
CORRER_CAL <- c(incert10km = TRUE, sinTexas = TRUE, ambos = TRUE)
UMBRAL_M <- 10000

occ <- read_csv(file.path(base_dir, "resultados_por_especie", "TODAS_LAS_ESPECIES_final_manual.csv"),
                show_col_types = FALSE)
val <- read_csv(file.path(out_dir, "valores_ambientales_por_registro.csv"), show_col_types = FALSE)
stopifnot(nrow(occ) == nrow(val), all(occ$species == val$species))
stopifnot(isTRUE(all.equal(as.numeric(occ$decimalLongitude), as.numeric(val$decimalLongitude), tolerance = 1e-6)))
val$gbifID <- occ$gbifID
val$elevation_extracted <- occ$elevation_extracted
val$flag_incert <- !is.na(occ$coordinateUncertaintyInMeters) & occ$coordinateUncertaintyInMeters > UMBRAL_M
val$flag_texas  <- occ$species == "Dasypterus ega" & occ$countryCode %in% "US"
message(sprintf("Registros: %d | incertidumbre > %d m: %d | D. ega fuera del Neotrópico (US): %d",
                nrow(val), UMBRAL_M, sum(val$flag_incert), sum(val$flag_texas)))

tabla_desde_registros <- function(v) {
  het <- v %>% group_by(species) %>%
    summarise(n = n(), tri_media = mean(tri, na.rm = TRUE),
              bio1_rugosidad_media = mean(bio1_rugosidad_local, na.rm = TRUE),
              bio12_rugosidad_media = mean(bio12_rugosidad_local, na.rm = TRUE), .groups = "drop")
  elev <- v %>% filter(!is.na(elevation_extracted)) %>% group_by(species) %>%
    summarise(elev_mediana = median(elevation_extracted),
              elev_rango_robusto = quantile(elevation_extracted, 0.975) - quantile(elevation_extracted, 0.025),
              .groups = "drop")
  vars_finales <- c("bio5", "bio6", "bio16", "bio17")
  set.seed(2026)
  re <- list()
  for (sp in sort(unique(v$species))) {
    dsp <- v %>% filter(species == sp) %>% select(all_of(vars_finales)) %>% na.omit()
    if (nrow(dsp) <= length(vars_finales)) next
    r <- tryCatch(ntbox::cov_center(data = dsp, mve = TRUE, level = 0.95, vars = vars_finales), error = function(e) NULL)
    if (is.null(r)) next
    dd <- ntbox::inEllipsoid(centroid = r$centroid, eShape = r$covariance, env_data = dsp, level = 0.95)$mh_dist
    re[[sp]] <- tibble(species = sp, niche_volume = r$niche_volume,
                       dist_centroide_media = mean(dd), dist_centroide_sd = sd(dd))
  }
  bind_rows(re) %>% left_join(het, by = "species") %>% left_join(elev, by = "species")
}

variantes <- list(
  incert10km = function(v) v %>% filter(!flag_incert),
  sinTexas   = function(v) v %>% filter(!flag_texas),
  ambos      = function(v) v %>% filter(!flag_incert, !flag_texas)
)
for (et in names(variantes)) {
  if (!CORRER_CAL[et]) next
  v <- variantes[[et]](val)
  message(sprintf("\n### %s: %d registros (se quitan %d)", et, nrow(v), nrow(val) - nrow(v)))
  tb <- tabla_desde_registros(v)
  write_csv(tb, file.path(out_dir, sprintf("ROBUSTEZ_v7_%s_tabla.csv", et)))
  d <- preparar(tb)
  message(sprintf("%s: %d especies; n mínimo = %d", et, nrow(d), min(d$n)))
  correr_pgls(d, MODELOS_V7, et)
  correr_phylopath(d, et)
}

# ------------------------------------------------------------------
# Comparación con el análisis principal
# ------------------------------------------------------------------
comp <- read_csv(file.path(out_dir, "pgls_v7_resumen.csv"), show_col_types = FALSE) %>%
  filter(!grepl("^(D_|E_)", modelo)) %>%
  select(modelo, termino, estimate_principal = estimate_mediana, se_principal = se_mediana,
         pct_sig_principal = pct_significativos, lambda_principal = lambda_mediana)
for (et in names(variantes)) {
  f <- file.path(out_dir, sprintf("ROBUSTEZ_v7_%s_pgls_resumen.csv", et))
  if (!file.exists(f)) next
  r <- read_csv(f, show_col_types = FALSE) %>%
    select(modelo, termino, estimate_mediana, se_mediana, pct_significativos, lambda_mediana)
  names(r)[3:6] <- paste0(c("estimate_", "se_", "pct_sig_", "lambda_"), et)
  comp <- comp %>% left_join(r, by = c("modelo", "termino"))
}
write_csv(comp, file.path(out_dir, "ROBUSTEZ_v7_sensibilidad_comparacion.csv"))
message("\n✅ Listo. Archivos ROBUSTEZ_v7_{incert10km,sinTexas,ambos}_* y ROBUSTEZ_v7_sensibilidad_comparacion.csv")
print(comp, n = Inf, width = Inf)
