#####################################################################
# PGLS PRINCIPALES, ESPECIFICACIÓN v7 (500 árboles, lambda por ML)
#
# Cambios respecto a pgls_500_arboles.R (respuesta a la revisión del
# manuscrito y a los diagnósticos de no linealidad / muestreo):
#   1. Las distancias al centroide se analizan en escala LOG (distribución
#      muy asimétrica a la derecha en escala original).
#   2. Todos los modelos incluyen log(n) = log(número de registros) para
#      controlar el esfuerzo de muestreo.
#   3. El rango elevacional incluye TRI centrada y TRI^2 (curvatura).
#      TRI se centra en su media de las 53 especies (TRI_CENTRO).
#
# Modelos principales:
#   log_niche_volume ~ tri_c + bio12_rugosidad_media + log_n
#   log_dist_media   ~ tri_c + bio12_rugosidad_media + log_n
#   log_dist_sd      ~ tri_c + bio12_rugosidad_media + log_n
#   elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n
#
# Diagnósticos (se guardan aparte, no son los modelos reportados):
#   D_*: los 3 primeros modelos + tri_c2, para documentar que NO hay
#        curvatura en esas respuestas y por eso solo se modela en el rango.
#   E_elev_sin_curvatura: rango elevacional lineal + log_n (referencia
#        para delta-AIC y para ver cuánto cambia el efecto de TRI).
#
# Salidas (resultadosfinales/):
#   pgls_v7_completo.csv        una fila por árbol x modelo x término
#   pgls_v7_resumen.csv         mediana, IC 2.5-97.5 entre árboles, SE mediano,
#                               % significativos, lambda, R2
#   pgls_v7_AIC.csv             AIC por árbol (principal vs. diagnóstico)
#   pgls_v7_pico_TRI.csv        TRI en que el rango elevacional alcanza su
#                               máximo (por árbol) y su resumen
#####################################################################
pkgs <- c("ape", "caper", "dplyr", "readr", "tibble", "tidyr")
faltantes <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(faltantes) > 0) install.packages(faltantes, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

equivalencias <- read_csv(file.path(out_dir, "equivalencias_especies_arbol.csv"), show_col_types = FALSE)
datos <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)

TRI_CENTRO <- mean(datos$tri_media)       # constante de centrado (se reutiliza en toda la v7)
message(sprintf("TRI_CENTRO = %.4f", TRI_CENTRO))

datos_pgls <- datos %>%
  left_join(equivalencias %>% dplyr::select(species, tip_label), by = "species") %>%
  mutate(log_niche_volume = log(niche_volume),
         log_dist_media   = log(dist_centroide_media),
         log_dist_sd      = log(dist_centroide_sd),
         log_n            = log(n),
         tri_c            = tri_media - TRI_CENTRO,
         tri_c2           = tri_c^2) %>%
  rename(tiplabel = tip_label) %>%
  filter(!is.na(tiplabel), !is.na(niche_volume)) %>%
  as.data.frame()
stopifnot(nrow(datos_pgls) == 53)

arboles <- ape::read.nexus(file.path(base_dir, "trees", "output.nex"))
stopifnot(all(datos_pgls$tiplabel %in% arboles[[1]]$tip.label))

modelos <- list(
  # ---- PRINCIPALES ----
  log_niche_volume   = log_niche_volume   ~ tri_c + bio12_rugosidad_media + log_n,
  log_dist_media     = log_dist_media     ~ tri_c + bio12_rugosidad_media + log_n,
  log_dist_sd        = log_dist_sd        ~ tri_c + bio12_rugosidad_media + log_n,
  elev_rango_robusto = elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n,
  # ---- DIAGNÓSTICOS ----
  D_log_niche_volume = log_niche_volume   ~ tri_c + tri_c2 + bio12_rugosidad_media + log_n,
  D_log_dist_media   = log_dist_media     ~ tri_c + tri_c2 + bio12_rugosidad_media + log_n,
  D_log_dist_sd      = log_dist_sd        ~ tri_c + tri_c2 + bio12_rugosidad_media + log_n,
  E_elev_sin_curvatura = elev_rango_robusto ~ tri_c + bio12_rugosidad_media + elev_mediana + log_n
)

resultados <- list(); aics <- list(); picos <- list(); fallos <- 0
t0 <- Sys.time()
for (i in seq_along(arboles)) {
  arbol_i <- arboles[[i]]; arbol_i$node.label <- NULL
  comp_data <- tryCatch(
    caper::comparative.data(phy = arbol_i, data = datos_pgls, names.col = "tiplabel",
                             vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE),
    error = function(e) NULL)
  if (is.null(comp_data)) { fallos <- fallos + 1; next }

  for (nm in names(modelos)) {
    fit <- tryCatch(caper::pgls(modelos[[nm]], data = comp_data, lambda = "ML"), error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit); co <- as.data.frame(s$coefficients); co$termino <- rownames(co)
    resultados[[length(resultados) + 1]] <- co %>%
      transmute(arbol = i, modelo = nm, termino, estimate = Estimate, se = `Std. Error`,
                p_valor = `Pr(>|t|)`, lambda = fit$param["lambda"], r2 = s$r.squared)
    ll <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
    k  <- length(coef(fit)) + 1
    aics[[length(aics) + 1]] <- tibble(arbol = i, modelo = nm, AIC = -2 * ll + 2 * k)
    if (nm == "elev_rango_robusto") {
      b1 <- coef(fit)["tri_c"]; b2 <- coef(fit)["tri_c2"]
      picos[[length(picos) + 1]] <- tibble(arbol = i, b1 = b1, b2 = b2,
                                           TRI_pico = ifelse(b2 < 0, TRI_CENTRO - b1 / (2 * b2), NA_real_))
    }
  }
  if (i %% 50 == 0) message(sprintf("  ...%d/%d (%.1f min)", i, length(arboles),
                                     as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (fallos > 0) message(sprintf("⚠ %d árbol(es) fallaron", fallos))

completo <- bind_rows(resultados)
write_csv(completo, file.path(out_dir, "pgls_v7_completo.csv"))

resumen <- completo %>% filter(termino != "(Intercept)") %>%
  group_by(modelo, termino) %>%
  summarise(n_arboles = n(),
            estimate_mediana = median(estimate),
            estimate_p2_5 = quantile(estimate, 0.025), estimate_p97_5 = quantile(estimate, 0.975),
            se_mediana = median(se), p_mediana = median(p_valor),
            pct_significativos = mean(p_valor < 0.05) * 100,
            lambda_mediana = median(lambda), r2_mediano = median(r2), .groups = "drop") %>%
  arrange(modelo, termino)
write_csv(resumen, file.path(out_dir, "pgls_v7_resumen.csv"))

aic_tab <- bind_rows(aics) %>% tidyr::pivot_wider(names_from = modelo, values_from = AIC) %>%
  mutate(dAIC_vol_cuad  = D_log_niche_volume - log_niche_volume,
         dAIC_dmedia_cuad = D_log_dist_media - log_dist_media,
         dAIC_dsd_cuad  = D_log_dist_sd - log_dist_sd,
         dAIC_elev_cuad_vs_lineal = elev_rango_robusto - E_elev_sin_curvatura)
write_csv(aic_tab, file.path(out_dir, "pgls_v7_AIC.csv"))

pico_tab <- bind_rows(picos)
write_csv(pico_tab, file.path(out_dir, "pgls_v7_pico_TRI.csv"))

message("\n✅ Resumen de PGLS v7:"); print(resumen, n = Inf)
message("\nMediana de delta-AIC (negativo = el modelo con TRI^2 ajusta mejor):")
print(aic_tab %>% summarise(across(starts_with("dAIC"), \(x) median(x, na.rm = TRUE))))
message("\nTRI en el pico del rango elevacional (mediana, p2.5, p97.5):")
print(quantile(pico_tab$TRI_pico, c(0.5, 0.025, 0.975), na.rm = TRUE))
message(sprintf("(constante de centrado TRI_CENTRO = %.4f)", TRI_CENTRO))
