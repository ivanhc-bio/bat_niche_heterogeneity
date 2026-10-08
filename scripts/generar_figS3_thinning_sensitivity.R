#####################################################################
# FIGURE S3: sensibilidad del n retenido por especie según la
# distancia de thinning espacial (1, 2, 5, 10 km)
#
# Usa directamente sensibilidad_thinning_distancias.csv, ya generado
# por sensibilidad_thinning.R -- solo falta esta pieza de graficado.
#####################################################################

pkgs <- c("dplyr", "readr", "tidyr", "ggplot2")
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "D:/capitulo2/reanalyses_17sep"
out_dir  <- file.path(base_dir, "resultadosfinales")
fig_SI   <- file.path(base_dir, "manuscript_EN", "figures_SI")

tabla <- read_csv(file.path(base_dir, "sensibilidad_thinning_distancias.csv"), show_col_types = FALSE)

# Formato largo: una fila por especie x distancia, con el % de
# registros retenidos respecto al n disponible tras el filtro de
# elevación (n_after_elev_filter), que es el insumo común a las 4
# distancias evaluadas
datos_largos <- tabla %>%
  pivot_longer(cols = starts_with("n_thin_"),
               names_to = "distancia_km", values_to = "n_retenido") %>%
  mutate(distancia_km = as.numeric(gsub("n_thin_|km", "", distancia_km)),
         pct_retenido = 100 * n_retenido / n_after_elev_filter,
         muestra_pequena_10km = n_after_elev_filter > 0 &
           n_retenido[distancia_km == 10] < 30) %>%
  filter(n_after_elev_filter > 0)

# Marca, por especie, si a 10 km queda con n < 30 (para resaltarlas)
n_en_10km <- datos_largos %>% filter(distancia_km == 10) %>%
  transmute(species, pequena_10km = n_retenido < 30)
datos_largos <- datos_largos %>% select(-muestra_pequena_10km) %>%
  left_join(n_en_10km, by = "species")

figS3 <- ggplot(datos_largos, aes(x = distancia_km, y = pct_retenido, group = species,
                                    color = pequena_10km)) +
  geom_line(alpha = 0.5, linewidth = 0.4) +
  geom_point(alpha = 0.6, size = 1.3) +
  scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "#B71C1C"),
                      labels = c("n >= 30 at 10 km", "n < 30 at 10 km"),
                      name = NULL) +
  scale_x_continuous(breaks = c(1, 2, 5, 10)) +
  labs(title = "Spatial thinning sensitivity",
       x = "Thinning distance (km)",
       y = "Records retained (% of post-elevation-filter n)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(fig_SI, "FigS3_thinning_sensitivity.png"), figS3, width = 7, height = 5.5, dpi = 400)
ggsave(file.path(fig_SI, "FigS3_thinning_sensitivity.pdf"), figS3, width = 7, height = 5.5)

message("Listo: FigS3_thinning_sensitivity.png / .pdf")
message(sprintf("Especies con n < 30 ya a 10 km: %d / %d",
                 sum(n_en_10km$pequena_10km), nrow(n_en_10km)))
