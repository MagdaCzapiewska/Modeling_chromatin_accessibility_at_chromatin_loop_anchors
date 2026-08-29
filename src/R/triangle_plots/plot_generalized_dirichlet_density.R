library(ggtern)
library(ggplot2)
library(viridis)

# 1. Funkcja obliczająca analityczną gęstość Uogólnionego Rozkładu Dirichleta (GDM3)
dgdm_3d <- function(p1, p2, p3, alpha, beta) {
  eps <- 1e-7
  p1 <- pmax(p1, eps)
  p2 <- pmax(p2, eps)
  p3 <- pmax(p3, eps)
  one_minus_p1 <- pmax(1 - p1, eps)
  
  # Logarytmy funkcji Beta (stałe normalizujące)
  log_B1 <- lbeta(alpha[1], beta[1])
  log_B2 <- lbeta(alpha[2], beta[2])
  
  # Logarytm gęstości GDM (3 składniki)
  log_pdf <- (alpha[1] - 1) * log(p1) + 
             (alpha[2] - 1) * log(p2) + 
             (beta[2] - 1) * log(p3) + 
             (beta[1] - alpha[2] - beta[2]) * log(one_minus_p1) - 
             log_B1 - log_B2
  
  exp(log_pdf)
}

# 2. Obsługa argumentów z linii poleceń
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 5) {
  pdf_file <- args[1]
  alpha <- as.numeric(args[2:3])
  beta  <- as.numeric(args[4:5])
} else {
  pdf_file <- "gdm_density.pdf"
  alpha <- c(2, 3)
  beta  <- c(8, 5)
  cat("Lack of complete arguments. Using default:\n")
  cat("  Output file:", pdf_file, "\n")
  cat("  Alpha:", paste(alpha, collapse = ", "), "\n")
  cat("  Beta: ", paste(beta, collapse = ", "), "\n\n")
}

# 3. Generowanie gęstej siatki punktów na sympleksie (p1 + p2 + p3 = 1)
grid_res <- 250
vals <- seq(0.001, 0.999, length.out = grid_res)
grid <- expand.grid(p1 = vals, p2 = vals)
grid$p3 <- 1 - grid$p1 - grid$p2

# Odrzucenie punktów leżących poza sympleksem
grid <- grid[grid$p3 > 0, ]

# 4. Wyznaczenie gęstości dla każdego punktu
grid$density <- dgdm_3d(grid$p1, grid$p2, grid$p3, alpha, beta)

# 5. Generowanie wykresu ggtern
p <- ggtern(grid, aes(x = p1, y = p2, z = p3, color = density)) +
  geom_point(size = 0.6, stroke = 0) +
  scale_color_viridis_c(
    option = "D", 
    name = "Density"
  ) +
  labs(
    title = "Probability density function (generalized Dirichlet)",
    subtitle = sprintf("alpha = (%g, %g), beta = (%g, %g)", 
                       alpha[1], alpha[2], beta[1], beta[2]),
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

cat("Pomyślnie wygenerowano plik PDF:", pdf_file, "\n")
