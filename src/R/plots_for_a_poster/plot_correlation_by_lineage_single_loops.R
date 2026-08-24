library(data.table)
library(ggplot2)
library(patchwork)

# 1. Konfiguracja ścieżek
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

output_file <- file.path(".", "correlation_by_loop_and_lineage.pdf")
correlation_dir = file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

# 2. Definicja okien czasowych
time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

# 3. Wczytanie i połączenie wszystkich plików korelacji
cor_list <- list()
for (tw in time_windows) {
  file_path <- file.path(correlation_dir, paste0("cor_", tw, ".tsv.gz"))
  if (file.exists(file_path)) {
    dt <- fread(file_path)
    cor_list[[tw]] <- dt
  } else {
    warning(paste("Plik nie istnieje:", file_path))
  }
}
all_cor <- rbindlist(cor_list)

# Filtrowanie pętli na podstawie pliku referencyjnego
all_cor <- all_cor[loop_id %in% loops_data$loop_id]

# Tworzenie czystej nazwy tkanki oraz identyfikatora okno_klaster
all_cor[, tissue := gsub("^[0-9]+_", "", population)]
all_cor[, time_window := factor(time_window, levels = time_windows)]
all_cor[, tw_pop := paste0(time_window, "_", population)]

# 4. Definicja wspólnego korzenia (root)
root_pops <- c(
  "00-02_0_Unknown", "00-02_1_Blastoderm", "00-02_2_Blastoderm", 
  "00-02_3_Blastoderm", "00-02_4_Blastoderm", "00-02_5_Unknown",
  "02-04_6_Blastoderm"
)

# 5. Poprawione mapowanie klastrów do 5 głównych linii rozwojowych (bez spacji i kropek)
lineages <- list(
  Ectodermal = c(
    root_pops,
    "02-04_0_Ectoderm_anlage", "02-04_1_Ectoderm_anlage",
    "04-06_0_Ectoderm_anlage", "04-06_1_Ectoderm_anlage", "04-06_2_Ectoderm_anlage", "04-06_3_Head_ectoderm", "04-06_5_Ectoderm_anlage", "04-06_10_Amnioserosa_anlage",
    "06-08_0_Epidermis_prim", "06-08_2_Ectoderm_anlage", "06-08_6_Head_ectoderm_prim", "06-08_9_Amnioserosa_anlage", "06-08_15_Ectoderm_anlage",
    "08-10_0_Epidermis_prim", "08-10_7_Tracheal_system_prim", "08-10_9_Epidermis_prim", "08-10_11_Amnioserosa",
    "10-12_0_Epidermis", "10-12_9_Tracheal_system", "10-12_14_Amnioserosa", "10-12_20_Salivary_gland",
    "12-14_0_Epidermis", "12-14_7_Tracheal_system", "12-14_14_Amnioserosa", "12-14_20_Salivary_gland", "12-14_25_Epidermis",
    "14-16_3_Epidermis", "14-16_5_Epidermis", "14-16_9_Tracheal_system", "14-16_16_Amnioserosa", "14-16_19_Head_ectoderm", "14-16_20_Salivary_gland", "14-16_24_Epidermis",
    "16-18_3_Epidermis", "16-18_8_Epidermis", "16-18_10_Tracheal_system", "16-18_17_Head_ectoderm", "16-18_18_Head_ectoderm", "16-18_19_Salivary_gland", "16-18_21_Amnioserosa",
    "18-20_2_Epidermis", "18-20_5_Head_ectoderm", "18-20_9_Tracheal_system"
  ),
  Neuroectodermal = c(
    root_pops,
    "04-06_12_Ventral_nerve_cord_prim", "04-06_15_Mesectoderm_anlage",
    "06-08_5_Ventral_nerve_cord_prim", "06-08_12_Mesectoderm_anlage",
    "08-10_1_Ventral_nerve_cord_prim", "08-10_10_Brain_prim", "08-10_13_PNS_and_sense", "08-10_15_Ventral_midline", "08-10_16_Glia",
    "10-12_1_Ventral_nerve_cord", "10-12_3_Ventral_nerve_cord", "10-12_4_Brain", "10-12_6_PNS_and_sense", "10-12_16_Glia", "10-12_18_Ventral_midline",
    "12-14_1_Ventral_nerve_cord", "12-14_4_Ventral_nerve_cord", "12-14_6_Brain", "12-14_8_PNS_and_sense", "12-14_10_Neural", "12-14_17_Glia", "12-14_19_Ventral_midline",
    "14-16_0_Ventral_nerve_cord", "14-16_4_Neural", "14-16_7_Ventral_nerve_cord", "14-16_11_PNS_and_sense", "14-16_12_Brain", "14-16_14_Glia", "14-16_22_Glia",
    "16-18_1_Ventral_nerve_cord", "16-18_6_Neural", "16-18_11_PNS_and_sense", "16-18_14_Ventral_nerve_cord", "16-18_15_Glia", "16-18_20_Neural", "16-18_25_Glia",
    "18-20_6_Neural", "18-20_7_Ventral_nerve_cord"
  ),
  Mesodermal = c(
    root_pops,
    "02-04_2_Mesoderm_anlage",
    "04-06_4_Muscle_prim", "04-06_7_Muscle_prim", "04-06_9_Muscle_prim", "04-06_14_Plasmatocytes", "04-06_17_Muscle_prim",
    "06-08_1_Muscle_prim", "06-08_8_Muscle_prim", "06-08_11_Muscle_prim", "06-08_13_Plasmatocytes", "06-08_14_Muscle_prim",
    "08-10_2_Muscle_prim", "08-10_5_Muscle_prim", "08-10_8_Fat_body", "08-10_14_Plasmatocytes",
    "10-12_5_Somatic_muscle", "10-12_8_Somatic_muscle", "10-12_11_Fat_body", "10-12_12_Plasmatocytes", "10-12_13_Visceral_muscle",
    "12-14_2_Somatic_muscle", "12-14_9_Fat_body", "12-14_12_Visceral_muscle", "12-14_13_Plasmatocytes", "12-14_18_Somatic_muscle", "12-14_22_Plasmatocytes", "12-14_24_Fat_body",
    "14-16_1_Somatic_muscle", "14-16_10_Visceral_muscle", "14-16_13_Fat_body", "14-16_17_Plasmatocytes", "14-16_25_Fat_body", "14-16_26_Somatic_muscle",
    "16-18_2_Somatic_muscle", "16-18_5_Somatic_muscle", "16-18_9_Visceral_muscle", "16-18_13_Fat_body",
    "18-20_0_Somatic_muscle"
  ),
  Endodermal_Digestive = c(
    root_pops,
    "02-04_4_Endoderm_anlage",
    "04-06_6_Midgut_prim", "04-06_8_Hindgut_prim",
    "06-08_3_Midgut_prim", "06-08_4_Foregut_prim", "06-08_10_Hindgut_prim",
    "08-10_3_Midgut_prim", "08-10_4_Foregut_prim", "08-10_12_Hindgut_prim",
    "10-12_2_Midgut", "10-12_7_Pharnyx", "10-12_15_Hindgut", "10-12_17_Malpighian_tubule",
    "12-14_3_Midgut", "12-14_5_Pharnyx", "12-14_15_Hindgut", "12-14_16_Malpighian_tubule", "12-14_23_Midgut",
    "14-16_2_Midgut", "14-16_6_Pharnyx", "14-16_15_Hindgut", "14-16_18_Malpighian_tubule", "14-16_23_Midgut",
    "16-18_0_Midgut", "16-18_7_Pharnyx", "16-18_12_Malpighian_tubule", "16-18_16_Hindgut", "16-18_22_Proventriculus",
    "18-20_1_Pharnyx", "18-20_3_Midgut", "18-20_8_Malpighian_tubule"
  ),
  EXE_Germline = c(
    root_pops,
    "02-04_3_Germ_cell",
    "04-06_11_Yolk", "04-06_18_Germ_cell",
    "06-08_7_Yolk", "06-08_16_Germ_cell",
    "08-10_6_Yolk", "08-10_17_Germ_cell", "08-10_18_Yolk",
    "10-12_10_Yolk", "10-12_19_Germ_cell", "10-12_21_Yolk",
    "12-14_11_Yolk", "12-14_21_Germ_cell",
    "14-16_8_Yolk", "14-16_21_Germ_cell",
    "16-18_4_Yolk", "16-18_23_Germ_cell",
    "18-20_4_Yolk"
  )
)

# 6. Pobranie unikalnych identyfikatorów pętli chromatynowych
unique_loops <- unique(all_cor$loop_id)
total_loops <- length(unique_loops)

# 7. Otwarcie urządzenia PDF (jedna strona = 11x17 cali)
pdf(output_file, width = 11, height = 17)

counter <- 0

# Główna pętla generująca osobną stronę dla każdego loop_id
for (current_loop in unique_loops) {

  counter <- counter + 1
  
  # WYPISYWANIE POSTĘPU W KONSOLI
  # \r sprawia, że tekst nadpisuje się w tej samej linijce, nie śmiecąc w konsoli
  cat(sprintf("\rPrzetwarzanie pętli: %d / %d (ID: %s) ... ", counter, total_loops, current_loop))
  flush.console() # Wymusza natychmiastowe odświeżenie tekstu w konsoli R
  
  # Filtrowanie danych tylko dla bieżącej pętli chromatynowej
  loop_subset <- all_cor[loop_id == current_loop]
  
  # Lista do przechowywania 5 wykresów (linii) dla tej konkretnej pętli
  loop_plots <- list()
  
  for (lineage_name in names(lineages)) {
    # Filtrowanie danych dla danej linii komórkowej
    lineage_data <- loop_subset[tw_pop %in% lineages[[lineage_name]]]
    
    # Rysowanie wykresu punktowego z linią trendu dla tkanki
    p <- ggplot(lineage_data, aes(x = time_window, y = spearman_rho, color = tissue, group = tissue)) +
      #geom_line(linewidth = 0.6, alpha = 0.5) +
      geom_point(size = 2, alpha = 0.9) +
      labs(
        title = paste0("Loop: ", current_loop, " | ", lineage_name, " Lineage"),
        x = "Developmental Time Window",
        y = "Spearman's Rho",
        color = "Tissue"
      ) +
      scale_y_continuous(limits = c(-1, 1)) + # Stały układ osi Y dla łatwiejszego porównywania
      theme_minimal(base_size = 9) +
      theme(
        panel.border = element_rect(fill = NA, color = "grey85"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank()
      )
    
    loop_plots[[lineage_name]] <- p
  }
  
  # Złożenie 5 wykresów w jeden pionowy panel (kolumnę) dla bieżącej pętli
  panel_for_loop <- wrap_plots(loop_plots, ncol = 1)
  
  # Każde wywołanie print() wewnątrz otwartego pdf() tworzy nową stronę
  print(panel_for_loop)
}

# Zamknięcie pliku PDF
dev.off()

message(paste("Gotowe! Wygenerowano wielostronicowy PDF:", output_file))
