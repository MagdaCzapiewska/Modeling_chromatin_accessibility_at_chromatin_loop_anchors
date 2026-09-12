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

# 2. Funkcja do obliczania teoretycznych korelacji w rozkładzie GDM
gdm_correlations <- function(alpha, beta) {
  a1 <- alpha[1]; a2 <- alpha[2]
  b1 <- beta[1];  b2 <- beta[2]
  
  Ez1 <- a1 / (a1 + b1)
  E1_z1 <- b1 / (a1 + b1)
  Ez1_sq <- a1 * (a1 + 1) / ((a1 + b1) * (a1 + b1 + 1))
  E1_z1_sq <- b1 * (b1 + 1) / ((a1 + b1) * (a1 + b1 + 1))
  Ez1_1_z1 <- a1 * b1 / ((a1 + b1) * (a1 + b1 + 1))
  
  Ez2 <- a2 / (a2 + b2)
  E1_z2 <- b2 / (a2 + b2)
  Ez2_sq <- a2 * (a2 + 1) / ((a2 + b2) * (a2 + b2 + 1))
  E1_z2_sq <- b2 * (b2 + 1) / ((a2 + b2) * (a2 + b2 + 1))
  Ez2_1_z2 <- a2 * b2 / ((a2 + b2) * (a2 + b2 + 1))
  
  Ep1 <- Ez1
  Ep2 <- E1_z1 * Ez2
  Ep3 <- E1_z1 * E1_z2
  
  Var_p1 <- Ez1_sq - Ez1^2
  Var_p2 <- E1_z1_sq * Ez2_sq - Ep2^2
  Var_p3 <- E1_z1_sq * E1_z2_sq - Ep3^2
  
  Cov_p1_p2 <- Ez1_1_z1 * Ez2 - Ep1 * Ep2
  Cov_p1_p3 <- Ez1_1_z1 * E1_z2 - Ep1 * Ep3
  Cov_p2_p3 <- E1_z1_sq * Ez2_1_z2 - Ep2 * Ep3
  
  r12 <- Cov_p1_p2 / sqrt(Var_p1 * Var_p2)
  r13 <- Cov_p1_p3 / sqrt(Var_p1 * Var_p3)
  r23 <- Cov_p2_p3 / sqrt(Var_p2 * Var_p3)
  
  return(c("p1 vs p2" = r12, "p1 vs p3" = r13, "p2 vs p3" = r23))
}

# 3. Obsługa argumentów z linii poleceń
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 5) {
  pdf_file <- args[1]
  alpha <- as.numeric(args[2:3])
  beta  <- as.numeric(args[4:5])
} else {
  pdf_file <- "gdm_density.pdf"
  alpha <- c(0.5, 3)
  beta  <- c(0.5, 3)
  cat("Lack of complete arguments. Using default:\n")
  cat("  Output file:", pdf_file, "\n")
  cat("  Alpha:", paste(alpha, collapse = ", "), "\n")
  cat("  Beta: ", paste(beta, collapse = ", "), "\n\n")
}

# 4. Obliczenie korelacji i sformatowanie podpisu
correlations <- gdm_correlations(alpha, beta)

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

# 5. Generowanie gęstej siatki punktów na sympleksie (p1 + p2 + p3 = 1)
grid_res <- 150
vals <- seq(0.001, 0.999, length.out = grid_res)
grid <- expand.grid(p1 = vals, p2 = vals)
grid$p3 <- 1 - grid$p1 - grid$p2

# Odrzucenie punktów leżących poza sympleksem
grid <- grid[grid$p3 > 0, ]

# 6. Wyznaczenie gęstości dla każdego punktu
grid$density <- dgdm_3d(grid$p1, grid$p2, grid$p3, alpha, beta)

# 7. Generowanie wykresu ggtern (p2 po lewej (x), p1 u góry (y), p3 po prawej (z))
p <- ggtern(grid, aes(x = p2, y = p1, z = p3, color = density)) +
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
    title = "Probability density function (generalized Dirichlet)",
    subtitle = sprintf("alpha = (%g, %g), beta = (%g, %g)", 
                       alpha[1], alpha[2], beta[1], beta[2]),
    caption = corr_text,
    x = "p2",
    y = "p1",
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

cat("Pomyślnie wygenerowano plik PDF:", pdf_file, "\n")
