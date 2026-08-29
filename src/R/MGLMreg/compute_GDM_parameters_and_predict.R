library(data.table)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

# Rozważamy fitowanie MGLMreg w zależności od czasu (bez tkanki)

# Argumentem przyjmowanym jest okno czasowe, które moze być 00-02, 02-04, 04-06, 06-08, 08-10, 10-12, 12-14, 14-16, 16-18
# Jak dla danej pętli i dla komórek z danego okna czasowego mamy kolumny x_A1, x_A2, x_out, to MGLMfit jest wywołane tak, żeby estymować od najliczniejszych kategorii (najpierw zawsze x_out, a potem x_A1 lub x_A2)
# MGLMreg już wewnątrz sortuje tak, że przewiduje od najliczniejszych

# Dla każdego okna czasowego mamy taki plik

file.path(results, "MGLMfit_GDM_cor", "real_data", "all", "init_1e-6", "cor_{tw}.tsv.gz")
# kolumny: time_window	loop_id	min_est_over_se	spearman_rho	spearman_pvalue	pearson_r	pearson_pvalue	alpha_x_A1_est	alpha_x_A1_SE	alpha_x_A2_est	alpha_x_A2_SE	alpha_x_out_est	alpha_x_out_SE	beta_x_A1_est	beta_x_A1_SE	beta_x_A2_est	beta_x_A2_SE	beta_x_out_est	beta_x_out_SE

# Chcemy wziąć z niego time_window, loop_id, min_est_over_se, spearman_rho, spearman_pvalue, alpha_x_A1_est, alpha_x_A2_est, alpha_x_out_est, beta_x_A1_wst, beta_x_A2_est, beta_x_out_est
# Wiadomo, ze czasami będą NA (przynajmniej dla jednej kategorii w każdej pętli, a czasami nawet i więcej jak coś się nie dofitowało)

# Chcę, żebyś wziął każdą pętlę i policzył estymacje parametrów GDM z tego MGLMreg (biorąc intercept = 1 i time środkowy dla danego okna czasowego)
# Pewnie MGLMreg w tych parametrach ma wskazane, jaki predyktor i jaki parametr, utworzysz parametry dla dokładnie tych kategorii, dla których było fitowanie MGLMfit
# Zapisz to do pliku w kolejnych kolumnach
# Policz tez korelację spearmana pomiędzy A1 i A2 i ją tam zapisz
# zrób też predict

# 

