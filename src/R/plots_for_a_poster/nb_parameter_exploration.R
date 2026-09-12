library(ggplot2)
library(patchwork)

MU <- 5000
X_MAX <- 12000
x_seq <- seq(0, X_MAX, by = 10)


base_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "#f1f2f6"),
    axis.title.x = element_blank(), 
    
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    
    plot.background = element_rect(fill = "#fbc53105", color = NA),
    plot.margin = margin(12, 6, 12, 6)
  )


p_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "ATAC-seq\nreads\nper cell", 
           size = 5, fontface = "bold", color = "#2c3e50", hjust = 0.5) +
  labs(title = "Parameter") + 
  theme_void() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 14)),
    plot.background = element_rect(fill = "#f8f9fa", color = NA),
    plot.margin = margin(12, 6, 12, 6)
  )

df_05 <- data.frame(x = x_seq, y = dnbinom(x_seq, mu = MU, size = 0.5))
p1 <- ggplot(df_05, aes(x = x, y = y)) +
  geom_area(fill = "#9b59b6", alpha = 0.25) +
  geom_line(color = "#8e44ad", size = 1.2) + 
  geom_vline(xintercept = MU, color = "#e74c3c", linetype = "dashed", size = 0.9) +
  labs(title = "size = 0.5") + 
  scale_x_continuous(breaks = c(0, MU, 10000), labels = c("0", "5000", "10,000")) +
  base_theme

df_10 <- data.frame(x = x_seq, y = dnbinom(x_seq, mu = MU, size = 1.0))
p2 <- ggplot(df_10, aes(x = x, y = y)) +
  geom_area(fill = "#3498db", alpha = 0.25) +
  geom_line(color = "#2980b9", size = 1.2) +
  geom_vline(xintercept = MU, color = "#e74c3c", linetype = "dashed", size = 0.9) +
  labs(title = "size = 1.0") +
  scale_x_continuous(breaks = c(0, MU, 10000), labels = c("0", "5000", "10,000")) +
  base_theme

df_20 <- data.frame(x = x_seq, y = dnbinom(x_seq, mu = MU, size = 2.0))
p3 <- ggplot(df_20, aes(x = x, y = y)) +
  geom_area(fill = "#1abc9c", alpha = 0.25) +
  geom_line(color = "#16a085", size = 1.2) +
  geom_vline(xintercept = MU, color = "#e74c3c", linetype = "dashed", size = 0.9) +
  labs(title = "size = 2.0") +
  scale_x_continuous(breaks = c(0, MU, 10000), labels = c("0", "5000", "10,000")) +
  base_theme

df_inf <- data.frame(x = x_seq, y = dpois(x_seq, lambda = MU))
p4 <- ggplot(df_inf, aes(x = x, y = y)) +
  geom_area(fill = "#e67e22", alpha = 0.25) +
  geom_line(color = "#d35400", size = 1.2) +
  geom_vline(xintercept = MU, color = "#e74c3c", linetype = "dashed", size = 0.9) +
  labs(title = "size = infinity") + 
  scale_x_continuous(breaks = c(0, MU, 10000), labels = c("0", "5000", "10,000")) +
  base_theme

p5 <- ggplot() +
  geom_segment(aes(x = MU, xend = MU, y = 0, yend = 1), color = "#2c3e50", size = 1.8) +
  geom_point(aes(x = MU, y = 1), color = "#2c3e50", size = 4) +
  geom_vline(xintercept = MU, color = "#e74c3c", linetype = "dashed", size = 0.9) +
  labs(title = "Fixed number of reads") + 
  scale_x_continuous(breaks = c(0, MU, 10000), labels = c("0", "5000", "10,000"), limits = c(0, X_MAX)) +
  scale_y_continuous(limits = c(0, 1.05)) +
  base_theme

layout <- (
  p_label | p1 | p2 | p3 | p4 | p5
) + 
  plot_layout(ncol = 6, widths = c(0.6, 1, 1, 1, 1, 1)) +
  plot_annotation(
    title = expression(bold("Total ATAC-seq read counts (negative binomial dist. for " * mu * " = 5000)")),
    theme = theme(
      plot.title.position = "plot", 
      plot.title = element_text(size = 22, hjust = 0.5, color = "#2c3e50", margin = margin(t = 12, b = 4))
    )
  )

pdf_output_path = "nb_parameter_exploration_R.pdf"
pdf(file = pdf_output_path, width = 12.5, height = 2.5)
print(layout)
dev.off()
