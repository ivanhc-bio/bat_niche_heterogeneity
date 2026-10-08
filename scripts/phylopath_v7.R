#####################################################################
# PATH ANALYSIS FILOGENÉTICO, ESPECIFICACIÓN v7
#
# Igual que phylopath_mediacion.R pero con la especificación nueva:
#   - TRI entra como par (tri_c, tri_c2) para capturar la curvatura sobre
#     el rango elevacional; tri_c2 ~ tri_c se declara en TODOS los modelos
#     (es una función determinista de tri_c, no una variable independiente).
#   - bio12_rugosidad_media ~ tri_c se declara en TODOS los modelos
#     (la covarianza TRI-precipitación ya mostró ganar 100 % de los árboles).
#   - log_n (esfuerzo de muestreo) es predictor de rango elevacional y
#     de volumen en TODOS los modelos.
# Los 4 candidatos solo difieren en la estructura de log_niche_volume:
#   mediacion_completa, sin_mediacion, mediacion_parcial, rutas_independientes
#
# Salidas (resultadosfinales/):
#   phylopath_v7_seleccion_modelo.csv     % de árboles que gana cada modelo
#   phylopath_v7_cicc_por_arbol.csv       CICc por modelo y árbol (diferencias)
#   phylopath_v7_coeficientes_estandarizados.csv  mediana/IC de cada ruta
#                                         (coeficientes del modelo promediado)
#   phylopath_v7_tabla_rutas.csv          rutas del modelo ganador, escala original (PGLS)
#   phylopath_v7_modelos_candidatos.png   diagrama de los 4 candidatos
#####################################################################
pkgs <- c("phylopath", "ape", "caper", "dplyr", "readr", "tibble")
faltantes <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(faltantes) > 0) install.packages(faltantes, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

equivalencias <- read_csv(file.path(out_dir, "equivalencias_especies_arbol.csv"), show_col_types = FALSE)
datos <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)
TRI_CENTRO <- mean(datos$tri_media)

datos_pp <- datos %>%
  left_join(equivalencias %>% dplyr::select(species, tip_label), by = "species") %>%
  mutate(log_niche_volume = log(niche_volume), log_n = log(n),
         tri_c = tri_media - TRI_CENTRO, tri_c2 = tri_c^2) %>%
  rename(tiplabel = tip_label) %>%
  filter(!is.na(tiplabel), !is.na(niche_volume), !is.na(elev_rango_robusto)) %>%
  as.data.frame()
rownames(datos_pp) <- datos_pp$tiplabel
stopifnot(nrow(datos_pp) == 53)

arboles <- ape::read.nexus(file.path(base_dir, "trees", "output.nex"))

comunes <- list(
  elev_rango_robusto ~ tri_c + tri_c2 + log_n,
  bio12_rugosidad_media ~ tri_c,
  tri_c2 ~ tri_c
)
modelos_causales <- define_model_set(
  mediacion_completa   = c(comunes, log_niche_volume ~ elev_rango_robusto + bio12_rugosidad_media + log_n),
  sin_mediacion        = c(comunes, log_niche_volume ~ tri_c + tri_c2 + bio12_rugosidad_media + log_n),
  mediacion_parcial    = c(comunes, log_niche_volume ~ elev_rango_robusto + tri_c + tri_c2 + bio12_rugosidad_media + log_n),
  rutas_independientes = c(comunes, log_niche_volume ~ bio12_rugosidad_media + log_n)
)

plot_modelos <- tryCatch(plot_model_set(modelos_causales), error = function(e) NULL)
if (!is.null(plot_modelos))
  ggplot2::ggsave(file.path(out_dir, "phylopath_v7_modelos_candidatos.png"), plot_modelos,
                   width = 10, height = 8, dpi = 300)

ganador <- character(length(arboles)); cicc <- list(); coefs_largo <- list(); fallos <- 0
t0 <- Sys.time()
for (i in seq_along(arboles)) {
  arbol_i <- arboles[[i]]; arbol_i$node.label <- NULL
  ajuste <- tryCatch(phylo_path(modelos_causales, data = datos_pp, tree = arbol_i, model = "lambda"),
                     error = function(e) NULL)
  if (is.null(ajuste)) { fallos <- fallos + 1; next }
  s <- tryCatch(summary(ajuste), error = function(e) NULL)
  if (is.null(s)) { fallos <- fallos + 1; next }
  ganador[i] <- as.character(s$model[1])
  cicc[[length(cicc) + 1]] <- tibble(arbol = i, modelo = as.character(s$model), CICc = s$CICc, delta_CICc = s$delta_CICc)

  prom <- tryCatch(average(ajuste), error = function(e) NULL)
  if (!is.null(prom)) {
    m <- as.matrix(prom$coef)
    for (a in rownames(m)) for (b in colnames(m))
      if (!is.na(m[a, b]) && m[a, b] != 0)
        coefs_largo[[length(coefs_largo) + 1]] <- tibble(arbol = i, desde = a, hacia = b, coef = m[a, b])
  }
  if (i %% 50 == 0) message(sprintf("  ...%d/%d (%.1f min)", i, length(arboles),
                                     as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (fallos > 0) message(sprintf("⚠ %d árbol(es) fallaron", fallos))

sel <- tibble(arbol = seq_along(arboles), modelo_ganador = ganador) %>% filter(modelo_ganador != "") %>%
  count(modelo_ganador, name = "n_arboles") %>% mutate(pct = round(100 * n_arboles / sum(n_arboles), 1)) %>%
  arrange(desc(n_arboles))
write_csv(sel, file.path(out_dir, "phylopath_v7_seleccion_modelo.csv"))
write_csv(bind_rows(cicc), file.path(out_dir, "phylopath_v7_cicc_por_arbol.csv"))
message("\n✅ Modelo ganador por árbol:"); print(sel)

resumen_coef <- bind_rows(coefs_largo) %>% group_by(desde, hacia) %>%
  summarise(n_arboles = n(), coef_mediana = median(coef),
            coef_p2_5 = quantile(coef, 0.025), coef_p97_5 = quantile(coef, 0.975), .groups = "drop")
write_csv(resumen_coef, file.path(out_dir, "phylopath_v7_coeficientes_estandarizados.csv"))
message("\n✅ Coeficientes estandarizados (modelo promediado):"); print(resumen_coef, n = Inf)

# Rutas en escala original (PGLS por ecuación) para el modelo ganador más frecuente
modelo_top <- sel$modelo_ganador[1]
message(sprintf("\nTabla de rutas en escala original para el modelo '%s'", modelo_top))
ecuaciones <- list(
  "elev_rango_robusto <- tri_c + tri_c2 + log_n" = elev_rango_robusto ~ tri_c + tri_c2 + log_n,
  "log_niche_volume <- elev + tri_c + tri_c2 + bio12 + log_n" =
    log_niche_volume ~ elev_rango_robusto + tri_c + tri_c2 + bio12_rugosidad_media + log_n,
  "bio12_rugosidad_media <- tri_c" = bio12_rugosidad_media ~ tri_c
)
if (modelo_top != "mediacion_parcial")
  message("⚠ El ganador no es mediacion_parcial: ajusta 'ecuaciones' a la estructura ganadora antes de reportar.")
res <- list()
for (i in seq_along(arboles)) {
  arbol_i <- arboles[[i]]; arbol_i$node.label <- NULL
  cd <- tryCatch(caper::comparative.data(phy = arbol_i, data = datos_pp, names.col = "tiplabel",
                                          vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE), error = function(e) NULL)
  if (is.null(cd)) next
  for (nm in names(ecuaciones)) {
    f <- tryCatch(caper::pgls(ecuaciones[[nm]], data = cd, lambda = "ML"), error = function(e) NULL)
    if (is.null(f)) next
    s <- summary(f); co <- as.data.frame(s$coefficients); co$termino <- rownames(co)
    res[[length(res) + 1]] <- co %>% transmute(arbol = i, ecuacion = nm, termino, estimate = Estimate,
                                               se = `Std. Error`, p_valor = `Pr(>|t|)`, lambda = f$param["lambda"])
  }
}
tabla_rutas <- bind_rows(res) %>% filter(termino != "(Intercept)") %>% group_by(ecuacion, termino) %>%
  summarise(n_arboles = n(), estimate_mediana = median(estimate), estimate_p2_5 = quantile(estimate, 0.025),
            estimate_p97_5 = quantile(estimate, 0.975), se_mediana = median(se),
            pct_significativos = mean(p_valor < 0.05) * 100, lambda_mediana = median(lambda), .groups = "drop")
write_csv(tabla_rutas, file.path(out_dir, "phylopath_v7_tabla_rutas.csv"))
print(tabla_rutas, n = Inf)
