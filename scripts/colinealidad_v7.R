#####################################################################
# COLINEALIDAD TRI vs. rugosidad de precipitación, ESPECIFICACIÓN v7
# (reemplaza colinealidad_1, _2, _3 y _6; Tablas S9-S13 y Fig. S5)
#
# Cambios respecto a la versión anterior:
#   - TRI entra centrada (tri_c) y, en el rango elevacional, junto con
#     tri_c2. TRI y TRI^2 se tratan como UN bloque ("TRI") en la partición
#     de varianza y en la prueba de supresión.
#   - Todos los modelos llevan log(n) (y elev_mediana en el rango elev.).
#   - Distancias al centroide en escala log.
#   - Los 8 candidatos de phylopath son los 4 de phylopath_v7.R con y sin
#     la arista bio12_rugosidad_media ~ tri_c (en v7 esa arista va en todos
#     los modelos principales; aquí se prueba que no sea una imposición).
#
# Partes (activa/desactiva con CORRER):
#   1 vif        : VIF/GVIF (instantáneo)
#   2 particion  : partición de varianza + supresión (500 árboles x 4 resp. x 4 spec.)
#   3 phylopath  : 8 candidatos con/sin la arista (500 árboles)
#   4 curva      : leave-one-out, remoción voraz y curva de robustez
#                  (la más lenta; piso = 20 como antes)
#
# Salidas (resultadosfinales/):
#   COLIN_V7_1_VIF.csv
#   COLIN_V7_2_completo.csv, COLIN_V7_2_particion.csv, COLIN_V7_2_supresion.csv
#   COLIN_V7_3_seleccion_8modelos.csv, COLIN_V7_3_con_vs_sin_arista.csv
#   COLIN_V7_4_leave_one_out.csv, COLIN_V7_4_trayectoria.csv,
#   COLIN_V7_4_curva_completo.csv, COLIN_V7_4_curva_resumen.csv
#####################################################################
CORRER <- c(vif = TRUE, particion = TRUE, phylopath = TRUE, curva = TRUE)
piso <- 20
curva_tamanos <- c(53, 45, 35, 25, piso)

pkgs <- c("ape", "caper", "phylopath", "car", "dplyr", "readr", "tibble", "tidyr", "purrr")
faltantes <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(faltantes) > 0) install.packages(faltantes, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))
select <- dplyr::select

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")

equivalencias <- read_csv(file.path(out_dir, "equivalencias_especies_arbol.csv"), show_col_types = FALSE)
datos <- read_csv(file.path(out_dir, "tabla_final_analisis.csv"), show_col_types = FALSE)
TRI_CENTRO <- mean(datos$tri_media)
message(sprintf("TRI_CENTRO = %.4f", TRI_CENTRO))

datos_pgls <- datos %>%
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
stopifnot(nrow(datos_pgls) == 53)
rownames(datos_pgls) <- datos_pgls$tiplabel

arboles <- ape::read.nexus(file.path(base_dir, "trees", "output.nex"))
stopifnot(all(datos_pgls$tiplabel %in% arboles[[1]]$tip.label))

r_global <- cor(datos_pgls$tri_media, datos_pgls$bio12_rugosidad_media)
r_tri_tri2 <- cor(datos_pgls$tri_c, datos_pgls$tri_c2)
message(sprintf("r(TRI, rugosidad precip.) = %.3f ;  r(tri_c, tri_c2) = %.3f", r_global, r_tri_tri2))

# ===================================================================
# 1. VIF / GVIF
#    (calculados con lm solo como vehículo de la matriz de diseño)
# ===================================================================
if (CORRER["vif"]) {
  d <- datos_pgls
  d$tri_cruda  <- d$tri_media
  d$tri_cruda2 <- d$tri_media^2
  v <- function(f) { x <- car::vif(lm(f, data = d)); if (is.matrix(x)) x[, 1] else x }

  tabla_vif <- bind_rows(
    tibble(conjunto = "volumen y distancias: tri_c + bio12 + log_n",
           predictor = names(v(log_niche_volume ~ tri_c + bio12_rugosidad_media + log_n)),
           vif = as.numeric(v(log_niche_volume ~ tri_c + bio12_rugosidad_media + log_n))),
    tibble(conjunto = "rango elev.: tri_c + tri_c2 + bio12 + elev_mediana + log_n",
           predictor = names(v(elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n)),
           vif = as.numeric(v(elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n))),
    tibble(conjunto = "rango elev. con TRI SIN centrar (para mostrar el efecto de centrar)",
           predictor = names(v(elev_rango_robusto ~ tri_cruda + tri_cruda2 + bio12_rugosidad_media + elev_mediana + log_n)),
           vif = as.numeric(v(elev_rango_robusto ~ tri_cruda + tri_cruda2 + bio12_rugosidad_media + elev_mediana + log_n)))
  )
  # GVIF del bloque TRI (lineal + cuadrático) tratado como un solo término
  g <- car::vif(lm(elev_rango_robusto ~ poly(tri_c, 2) + bio12_rugosidad_media + elev_mediana + log_n, data = d))
  tabla_gvif <- tibble(conjunto = "rango elev., GVIF con bloque TRI = poly(tri_c, 2)",
                       predictor = rownames(g), vif = g[, "GVIF"], gl = g[, "Df"],
                       gvif_ajustado = g[, "GVIF^(1/(2*Df))"])
  write_csv(bind_rows(tabla_vif, tabla_gvif), file.path(out_dir, "COLIN_V7_1_VIF.csv"))
  message("\n✅ VIF:"); print(tabla_vif, n = Inf); print(tabla_gvif)
}

# ===================================================================
# 2. PARTICIÓN DE VARIANZA + SUPRESIÓN (TRI como bloque)
#    Especificaciones por respuesta: completo, solo_tri, solo_precip, base
#    ('base' = solo los controles log_n [+ elev_mediana]).
# ===================================================================
if (CORRER["particion"]) {
  spec_simple <- function(y) list(
    completo    = as.formula(paste(y, "~ tri_c + bio12_rugosidad_media + log_n")),
    solo_tri    = as.formula(paste(y, "~ tri_c + log_n")),
    solo_precip = as.formula(paste(y, "~ bio12_rugosidad_media + log_n")),
    base        = as.formula(paste(y, "~ log_n")))
  especificaciones <- list(
    log_niche_volume = spec_simple("log_niche_volume"),
    log_dist_media   = spec_simple("log_dist_media"),
    log_dist_sd      = spec_simple("log_dist_sd"),
    elev_rango_robusto = list(
      completo    = elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n,
      solo_tri    = elev_rango_robusto ~ tri_c + tri_c2 + elev_mediana + log_n,
      solo_precip = elev_rango_robusto ~ bio12_rugosidad_media + elev_mediana + log_n,
      base        = elev_rango_robusto ~ elev_mediana + log_n))

  res <- list(); t0 <- Sys.time()
  for (i in seq_along(arboles)) {
    arbol_i <- arboles[[i]]; arbol_i$node.label <- NULL
    cd <- tryCatch(caper::comparative.data(phy = arbol_i, data = datos_pgls, names.col = "tiplabel",
                                            vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE), error = function(e) NULL)
    if (is.null(cd)) next
    for (resp in names(especificaciones)) for (sp in names(especificaciones[[resp]])) {
      fit <- tryCatch(caper::pgls(especificaciones[[resp]][[sp]], data = cd, lambda = "ML"), error = function(e) NULL)
      if (is.null(fit)) next
      s <- summary(fit); co <- as.data.frame(s$coefficients); co$termino <- rownames(co)
      res[[length(res) + 1]] <- co %>% transmute(arbol = i, respuesta = resp, especificacion = sp, termino,
                                                 estimate = Estimate, se = `Std. Error`, p_valor = `Pr(>|t|)`,
                                                 r2 = s$r.squared)
    }
    if (i %% 50 == 0) message(sprintf("  partición: %d/%d (%.1f min)", i, length(arboles),
                                       as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  tc <- bind_rows(res)
  write_csv(tc, file.path(out_dir, "COLIN_V7_2_completo.csv"))

  part <- tc %>% distinct(arbol, respuesta, especificacion, r2) %>%
    pivot_wider(names_from = especificacion, values_from = r2, names_prefix = "r2_") %>%
    mutate(incremento_total = r2_completo - r2_base,
           unico_tri    = r2_completo - r2_solo_precip,
           unico_precip = r2_completo - r2_solo_tri,
           compartido   = incremento_total - unico_tri - unico_precip)
  resumen_part <- part %>% group_by(respuesta) %>%
    summarise(n_arboles = n(), r2_base = median(r2_base), r2_completo = median(r2_completo),
              incremento_total = median(incremento_total), unico_tri = median(unico_tri),
              unico_precip = median(unico_precip), compartido = median(compartido), .groups = "drop")
  write_csv(resumen_part, file.path(out_dir, "COLIN_V7_2_particion.csv"))
  message("\n✅ Partición de varianza (mediana entre árboles; únicos y compartido son sobre el incremento tras los controles):")
  print(resumen_part)

  # Supresión: coeficiente en el modelo 'solo' vs. en el 'completo'
  sup <- function(term, solo) {
    tc %>% filter(termino == term, especificacion %in% c("completo", solo)) %>%
      select(arbol, respuesta, especificacion, estimate) %>%
      pivot_wider(names_from = especificacion, values_from = estimate) %>%
      filter(!is.na(.data[[solo]]), !is.na(completo)) %>%
      group_by(respuesta) %>%
      summarise(n_arboles = n(), coef_solo = median(.data[[solo]]), coef_completo = median(completo),
                cambio_signo_pct = mean(sign(.data[[solo]]) != sign(completo)) * 100, .groups = "drop") %>%
      mutate(variable = term)
  }
  resumen_sup <- bind_rows(sup("tri_c", "solo_tri"), sup("tri_c2", "solo_tri"),
                           sup("bio12_rugosidad_media", "solo_precip")) %>%
    select(respuesta, variable, everything())
  write_csv(resumen_sup, file.path(out_dir, "COLIN_V7_2_supresion.csv"))
  message("\n✅ Supresión (coeficiente solo vs. completo):"); print(resumen_sup, n = Inf)
}

# ===================================================================
# 3. PHYLOPATH: 8 candidatos (4 estructuras de volumen x con/sin arista TRI -> precipitación)
#    'tri_c2 ~ tri_c' va en TODOS (es determinista).
# ===================================================================
if (CORRER["phylopath"]) {
  base_e  <- list(elev_rango_robusto ~ tri_c + tri_c2 + log_n, tri_c2 ~ tri_c)
  con_a   <- c(base_e, list(bio12_rugosidad_media ~ tri_c))
  estr <- list(
    mediacion_completa   = log_niche_volume ~ elev_rango_robusto + bio12_rugosidad_media + log_n,
    sin_mediacion        = log_niche_volume ~ tri_c + tri_c2 + bio12_rugosidad_media + log_n,
    mediacion_parcial    = log_niche_volume ~ elev_rango_robusto + tri_c + tri_c2 + bio12_rugosidad_media + log_n,
    rutas_independientes = log_niche_volume ~ bio12_rugosidad_media + log_n)
  mods <- c(
    setNames(lapply(estr, function(f) c(con_a, f)), paste0(names(estr), "_con_arista")),
    setNames(lapply(estr, function(f) c(base_e, f)), paste0(names(estr), "_sin_arista")))
  modelos_causales <- do.call(define_model_set, mods)
  ganador <- character(length(arboles)); t0 <- Sys.time()
  for (i in seq_along(arboles)) {
    arbol_i <- arboles[[i]]; arbol_i$node.label <- NULL
    aj <- tryCatch(phylo_path(modelos_causales, data = datos_pgls, tree = arbol_i, model = "lambda"), error = function(e) NULL)
    if (is.null(aj)) next
    s <- tryCatch(summary(aj), error = function(e) NULL)
    if (is.null(s)) next
    ganador[i] <- as.character(s$model[1])
    if (i %% 50 == 0) message(sprintf("  phylopath 8 mod.: %d/%d (%.1f min)", i, length(arboles),
                                       as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  tg <- tibble(modelo_ganador = ganador) %>% filter(modelo_ganador != "") %>%
    mutate(con_arista = grepl("_con_arista$", modelo_ganador))
  sel8 <- tg %>% count(modelo_ganador, con_arista, name = "n_arboles") %>%
    mutate(pct = round(100 * n_arboles / sum(n_arboles), 1)) %>% arrange(desc(n_arboles))
  cvs <- tg %>% count(con_arista, name = "n_arboles") %>% mutate(pct = round(100 * n_arboles / sum(n_arboles), 1))
  write_csv(sel8, file.path(out_dir, "COLIN_V7_3_seleccion_8modelos.csv"))
  write_csv(cvs,  file.path(out_dir, "COLIN_V7_3_con_vs_sin_arista.csv"))
  message("\n✅ Ocho candidatos:"); print(sel8, n = Inf); print(cvs)
}

# ===================================================================
# 4. LEAVE-ONE-OUT + REMOCIÓN VORAZ + CURVA DE ROBUSTEZ
#    Criterio de colinealidad: r(tri_media, bio12_rugosidad_media) (como antes).
#    Con n pequeño el modelo de rango elevacional tiene 6 coeficientes:
#    interpretar los tamaños 25 y 20 como consistencia de dirección.
# ===================================================================
if (CORRER["curva"]) {
  loo <- purrr::map_dfr(datos_pgls$species, function(sp) {
    sub <- datos_pgls %>% filter(species != sp)
    r_sin <- cor(sub$tri_media, sub$bio12_rugosidad_media)
    tibble(species = sp, r_sin_esta_especie = r_sin, delta_r = r_global - r_sin)
  }) %>% arrange(desc(abs(delta_r)))
  write_csv(loo, file.path(out_dir, "COLIN_V7_4_leave_one_out.csv"))

  restantes <- datos_pgls$species
  trayectoria <- tibble(paso = 0, n_especies = length(restantes), r_tri_precip = r_global, especie_removida = NA_character_)
  paso <- 0
  while (length(restantes) > piso) {
    paso <- paso + 1
    sub_actual <- datos_pgls %>% filter(species %in% restantes)
    cand <- purrr::map_dfr(restantes, function(sp) {
      sub <- sub_actual %>% filter(species != sp)
      tibble(species = sp, r_si_sale = cor(sub$tri_media, sub$bio12_rugosidad_media))
    })
    peor <- cand %>% slice_min(abs(r_si_sale), n = 1, with_ties = FALSE)
    restantes <- setdiff(restantes, peor$species)
    trayectoria <- bind_rows(trayectoria, tibble(paso = paso, n_especies = length(restantes),
                                                 r_tri_precip = peor$r_si_sale, especie_removida = peor$species))
  }
  write_csv(trayectoria, file.path(out_dir, "COLIN_V7_4_trayectoria.csv"))
  message(sprintf("\nr bajó de %.3f (53 sp.) a %.3f (%d sp.)", r_global, tail(trayectoria$r_tri_precip, 1), piso))

  especies_en_tamano <- function(n) {
    k <- 53 - n
    if (k == 0) return(datos_pgls$species)
    setdiff(datos_pgls$species, trayectoria %>% filter(paso >= 1, paso <= k) %>% pull(especie_removida))
  }
  curva_tamanos <- sort(unique(curva_tamanos[curva_tamanos >= piso & curva_tamanos <= 53]), decreasing = TRUE)
  modelos <- list(
    log_niche_volume   = log_niche_volume   ~ tri_c + bio12_rugosidad_media + log_n,
    log_dist_media     = log_dist_media     ~ tri_c + bio12_rugosidad_media + log_n,
    log_dist_sd        = log_dist_sd        ~ tri_c + bio12_rugosidad_media + log_n,
    elev_rango_robusto = elev_rango_robusto ~ tri_c + tri_c2 + bio12_rugosidad_media + elev_mediana + log_n)

  rc <- list(); t0 <- Sys.time()
  for (n_sub in curva_tamanos) {
    ds <- datos_pgls %>% filter(species %in% especies_en_tamano(n_sub))
    r_sub <- cor(ds$tri_media, ds$bio12_rugosidad_media)
    message(sprintf("\n--- PGLS en subset de %d especies (r = %.3f) ---", n_sub, r_sub))
    for (i in seq_along(arboles)) {
      arbol_i <- ape::drop.tip(arboles[[i]], setdiff(arboles[[i]]$tip.label, ds$tiplabel)); arbol_i$node.label <- NULL
      cd <- tryCatch(caper::comparative.data(phy = arbol_i, data = ds, names.col = "tiplabel",
                                              vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE), error = function(e) NULL)
      if (is.null(cd)) next
      for (nm in names(modelos)) {
        fit <- tryCatch(caper::pgls(modelos[[nm]], data = cd, lambda = "ML"), error = function(e) NULL)
        if (is.null(fit)) next
        s <- summary(fit); co <- as.data.frame(s$coefficients); co$termino <- rownames(co)
        rc[[length(rc) + 1]] <- co %>% transmute(n_especies = n_sub, r_subset = r_sub, arbol = i, modelo = nm, termino,
                                                 estimate = Estimate, se = `Std. Error`, p_valor = `Pr(>|t|)`)
      }
      if (i %% 100 == 0) message(sprintf("  %d/%d árboles (%.1f min en total)", i, length(arboles),
                                          as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }
  tcurva <- bind_rows(rc)
  write_csv(tcurva, file.path(out_dir, "COLIN_V7_4_curva_completo.csv"))
  resumen_curva <- tcurva %>% filter(termino %in% c("tri_c", "tri_c2", "bio12_rugosidad_media", "log_n")) %>%
    group_by(n_especies, r_subset, modelo, termino) %>%
    summarise(n_arboles = n(), estimate_mediana = median(estimate),
              estimate_p2_5 = quantile(estimate, 0.025), estimate_p97_5 = quantile(estimate, 0.975),
              se_mediana = median(se), pct_significativos = mean(p_valor < 0.05) * 100, .groups = "drop") %>%
    arrange(modelo, termino, desc(n_especies))
  write_csv(resumen_curva, file.path(out_dir, "COLIN_V7_4_curva_resumen.csv"))
  message("\n✅ Curva de robustez:"); print(resumen_curva, n = Inf)
}
message("\n✅ Listo. Archivos COLIN_V7_* en resultadosfinales/")
