library(ggplot2)
library(cowplot)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
output_dir <- file.path(resultsdir, "plots_for_a_thesis", "synthetic_correlations")

# ==============================================================================
# PALETA KOLORÓW I ETYKIETY
# ==============================================================================

color_map <- c(
  "10" = "#99FF99", 
  "5"  = "#FFFF99", 
  "2"  = "#FFCC99", 
  "1"  = "#FFB3B3", 
  "0"  = "#99CCFF", 
  "NA" = "#D3D3D3"
)

legend_labels <- c(
  "10" = "Excellent (min EST/SE >= 10)",
  "5"  = "Good (min EST/SE >= 5)",
  "2"  = "Acceptable (min EST/SE >= 2)",
  "1"  = "Poor (min EST/SE >= 1)",
  "0"  = "Unstable (min EST/SE < 1)",
  "NA" = "Non-available"
)

# ==============================================================================
# GENEROWANIE SZTUCZNYCH DANYCH
# ==============================================================================
set.seed(42)

# Tworzymy wektor grup dla wszystkich 6 kategorii (w tym "NA") -> 6 * 100 = 600 obserwacji
groups <- rep(names(color_map), each = 100)

dummy_dt <- data.frame(
  spearman_rho = rnorm(length(groups), mean = 0, sd = 0.3),
  color_group = factor(
    groups, 
    levels = names(color_map) # POPRAWKA: Używamy wszystkich kluczy z color_map, łącznie z "NA"
  )
)

# ==============================================================================
# WYKRES HISTOGRAMU Z LEGENDĄ (POZIOMA)
# ==============================================================================
p <- ggplot(dummy_dt, aes(x = spearman_rho, fill = color_group)) +
  geom_histogram(binwidth = 0.08, color = "white", linewidth = 0.2, position = "stack") +
  scale_fill_manual(
    values = color_map,
    labels = legend_labels,
    name = "Estimation quality:",
    drop = FALSE
  ) +
  guides(fill = guide_legend(
    title.position = "top",  # Tytuł nad elementami
    title.hjust = 0.5,       # Wyśrodkowanie tytułu
    nrow = 1,                # Wymuszenie 1 wiersza (wszystkie 6 elementów w poziomie)
    byrow = TRUE
  )) +
  labs(
    title = "Sztuczny histogram do wycięcia legendy",
    x = "Estimated Spearman's correlation",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 8.5, color = "black"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(6, 8, 6, 8)
  )

# Zapis do PDF — szerokość 15 cali zapewnia wystarczająco dużo miejsca dla 6 elementów
output_file <- file.path(output_dir, "histogram_with_legend_with_NA_horizontal.pdf")
pdf(output_file, width = 15, height = 5)
print(p)
dev.off()

cat("Zapisano wykres ze sztucznym histogramem i legendą do:", output_file, "\n")
