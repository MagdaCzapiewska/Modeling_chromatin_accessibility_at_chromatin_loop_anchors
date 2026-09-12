library(ggtern)
library(ggplot2)
library(viridis)

# 1. Funkcja do analitycznego wyznaczania gęstości Dirichleta
ddirichlet_3d <- function(p1, p2, p3, alpha) {
  eps <- 1e-7
  p1 <- pmax(p1, eps)
  p2 <- pmax(p2, eps)
  p3 <- pmax(p3, eps)
  
  # Logarytm funkcji Beta dla wielu zmiennych (stała normalizująca)
  log_B <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  
  # Logarytm gęstości
  log_pdf <- (alpha[1] - 1) * log(p1) + 
             (alpha[2] - 1) * log(p2) + 
             (alpha[3] - 1) * log(p3) - log_B
  
  exp(log_pdf)
}

# 2. Funkcja do obliczania teoretycznych korelacji w rozkładzie Dirichleta
dirichlet_correlations <- function(alpha) {
  a0 <- sum(alpha)
  
  corr_12 <- -sqrt((alpha[1] * alpha[2]) / ((a0 - alpha[1]) * (a0 - alpha[2])))
  corr_13 <- -sqrt((alpha[1] * alpha[3]) / ((a0 - alpha[1]) * (a0 - alpha[3])))
  corr_23 <- -sqrt((alpha[2] * alpha[3]) / ((a0 - alpha[2]) * (a0 - alpha[3])))
  
  return(c("p1 vs p2" = corr_12, "p1 vs p3" = corr_13, "p2 vs p3" = corr_23))
}

# 3. Obsługa argumentów z linii poleceń
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 4) {
  pdf_file <- args[1]
  alpha <- as.numeric(args[2:4])
} else {
  pdf_file <- "dirichlet_density.pdf"
  alpha <- c(2, 5, 2)
  cat("Lack of complete set of arguments. Using default:\n")
  cat("  Output file:", pdf_file, "\n")
  cat("  Alpha:", paste(alpha, collapse = ", "), "\n\n")
}

# 4. Obliczenie korelacji i sformatowanie napisu pod trójkątem
correlations <- dirichlet_correlations(alpha)

corr_text <- sprintf(
  "Correlations:  r(p1, p2) = %.2f   |   r(p1, p3) = %.2f   |   r(p2, p3) = %.2f",
  correlations["p1 vs p2"], 
  correlations["p1 vs p3"], 
  correlations["p2 vs p3"]
)

# Wypisanie w konsoli
cat("Theoretical correlations between components:\n")
for (pair in names(correlations)) {
  cat(sprintf("  %s: %.4f\n", pair, correlations[pair]))
}
cat("\n")

# 5. Generowanie gęstej siatki punktów w sympleksie (p1 + p2 + p3 = 1)
grid_res <- 150
vals <- seq(0.001, 0.999, length.out = grid_res)
grid <- expand.grid(p1 = vals, p2 = vals)
grid$p3 <- 1 - grid$p1 - grid$p2

# Odrzucenie punktów leżących poza sympleksem
grid <- grid[grid$p3 > 0, ]

# 6. Obliczenie analitycznej gęstości prawdopodobieństwa
grid$density <- ddirichlet_3d(grid$p1, grid$p2, grid$p3, alpha)

# 7. Wykres ggtern z czarnym podpisem korelacji w kroju czcionki tytułu
p <- ggtern(grid, aes(x = p1, y = p2, z = p3, color = density)) +
  geom_point(size = 1.5, stroke = 0) +
  scale_color_viridis_c(
    option = "D", 
    name = "Density"
  ) +
  guides(color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = unit(6, "cm"),
    barheight = unit(0.4, "cm")
  )) +
  labs(
    title = "Probability density function (Dirichlet)",
    subtitle = sprintf("alpha = (%g, %g, %g)", alpha[1], alpha[2], alpha[3]),
    caption = corr_text,
    x = "p1",
    y = "p2",
    z = "p3"
  ) +
  theme_bw() +
  theme_showarrows() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, color = "black"),
    plot.caption = element_text(size = 10.5, face = "bold", hjust = 0.5, color = "black", margin = margin(t = 8, b = 6))
  )

# 8. Zapis do pliku PDF
pdf(pdf_file, width = 6.5, height = 7)
print(p)
dev.off()

cat("Successfully generated a file:", pdf_file, "\n")
