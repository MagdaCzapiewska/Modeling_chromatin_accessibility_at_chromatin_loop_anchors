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

# 2. Obsługa argumentów z linii poleceń
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

# 3. Generowanie gęstej siatki punktów w sympleksie (p1 + p2 + p3 = 1)
grid_res <- 250
vals <- seq(0.001, 0.999, length.out = grid_res)
grid <- expand.grid(p1 = vals, p2 = vals)
grid$p3 <- 1 - grid$p1 - grid$p2

# Odrzucenie punktów leżących poza sympleksem
grid <- grid[grid$p3 > 0, ]

# 4. Obliczenie analitycznej gęstości prawdopodobieństwa
grid$density <- ddirichlet_3d(grid$p1, grid$p2, grid$p3, alpha)

# 5. Wykres ggtern – rysowanie gęstej chmury punktów dającej efekt ciągłej mapy
p <- ggtern(grid, aes(x = p1, y = p2, z = p3, color = density)) +
  geom_point(size = 0.5, stroke = 0) +
  scale_color_viridis_c(
    option = "D", 
    name = "Density"
  ) +
  labs(
    title = "Probability density function (Dirichlet)",
    subtitle = sprintf("alpha = (%g, %g, %g)", alpha[1], alpha[2], alpha[3]),
    x = "p1",
    y = "p2",
    z = "p3"
  ) +
  theme_bw() +
  theme_showarrows()

# 6. Zapis do pliku PDF
pdf(pdf_file, width = 6, height = 5)
print(p)
dev.off()

cat("Successfully generated a file:", pdf_file, "\n")
