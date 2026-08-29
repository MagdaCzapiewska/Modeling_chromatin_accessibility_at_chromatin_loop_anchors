library(data.table)

safe_extract <- function(x, name) {
  if (is.null(x)) return(NA_real_)
  if (!(name %in% names(x))) return(NA_real_)
  as.numeric(x[name])
}

rGDM3 <- function(n, fit, sorted_names) {
  alpha <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("alpha_", nm)])
  beta  <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("beta_", nm)])

  p1 <- rbeta(n, alpha[1], beta[1])
  p2 <- rbeta(n, alpha[2], beta[2]) * (1 - p1)
  p3 <- pmax(0, 1 - p1 - p2)

  res <- cbind(p1, p2, p3)
  colnames(res) <- sorted_names
  res
}

compute_gdm_correlations_from_fit <- function(fit, n_sim = 10000, metadata_list = list()) {

  base_dt <- data.table(
    min_est_over_se = NA_real_,
    spearman_rho = NA_real_, spearman_pvalue = NA_real_,
    pearson_r = NA_real_, pearson_pvalue = NA_real_,
    alpha_x_A1_est = NA_real_, alpha_x_A1_SE = NA_real_,
    alpha_x_A2_est = NA_real_, alpha_x_A2_SE = NA_real_,
    alpha_x_out_est = NA_real_, alpha_x_out_SE = NA_real_,
    beta_x_A1_est = NA_real_, beta_x_A1_SE = NA_real_,
    beta_x_A2_est = NA_real_, beta_x_A2_SE = NA_real_,
    beta_x_out_est = NA_real_, beta_x_out_SE = NA_real_
  )
  
  if (length(metadata_list) > 0) {
    base_dt <- cbind(as.data.table(metadata_list), base_dt)
  }
  
  if (is.null(fit) || is.null(fit@estimate)) {
    return(base_dt)
  }
  
  est <- fit@estimate
  se <- fit@SE

  alpha_pars <- names(est)[grep("^alpha_", names(est))]
  fitted_names <- sub("^alpha_", "", alpha_pars)

  if (length(fitted_names) < 2 || !("x_out" %in% fitted_names)) {
    return(base_dt)
  }

  other_fitted <- setdiff(fitted_names, "x_out")[1]

  remaining_name <- setdiff(c("x_A1", "x_A2"), other_fitted)

  sorted_names <- c("x_out", other_fitted, remaining_name)
  
  alpha_est <- setNames(rep(NA_real_, 3), sorted_names)
  beta_est  <- setNames(rep(NA_real_, 3), sorted_names)
  alpha_se  <- setNames(rep(NA_real_, 3), sorted_names)
  beta_se   <- setNames(rep(NA_real_, 3), sorted_names)

  for (nm in sorted_names[1:2]) {
    alpha_est[nm] <- safe_extract(est, paste0("alpha_", nm))
    beta_est[nm]  <- safe_extract(est, paste0("beta_", nm))
  }
  
  alpha_se[sorted_names[1]] <- se[[1]]
  alpha_se[sorted_names[2]] <- se[[2]]
  beta_se[sorted_names[1]]  <- se[[3]]
  beta_se[sorted_names[2]]  <- se[[4]]
  
  est_vec <- c(alpha_est, beta_est)
  se_vec  <- c(alpha_se, beta_se)
  all_params_vec <- c(est_vec, se_vec)

  param_names <- c("alpha_x_A1_est", "alpha_x_A1_SE", 
                   "alpha_x_A2_est", "alpha_x_A2_SE", 
                   "alpha_x_out_est", "alpha_x_out_SE",
                   "beta_x_A1_est",  "beta_x_A1_SE",  
                   "beta_x_A2_est",  "beta_x_A2_SE",  
                   "beta_x_out_est", "beta_x_out_SE")
  
  vals_to_assign <- list(
    alpha_est["x_A1"], alpha_se["x_A1"],
    alpha_est["x_A2"], alpha_se["x_A2"],
    alpha_est["x_out"], alpha_se["x_out"],
    beta_est["x_A1"],  beta_se["x_A1"],
    beta_est["x_A2"],  beta_se["x_A2"],
    beta_est["x_out"], beta_se["x_out"]
  )
  
  base_dt[, (param_names) := vals_to_assign]

  if (sum(!is.na(all_params_vec)) < 8) {
    return(base_dt)
  }
  
  y_sim <- rGDM3(n_sim, fit, sorted_names)
  dt_sim <- as.data.table(y_sim)
  
  ct_s <- suppressWarnings(cor.test(dt_sim$x_A1, dt_sim$x_A2, method = "spearman"))
  ct_p <- suppressWarnings(cor.test(dt_sim$x_A1, dt_sim$x_A2, method = "pearson"))
  
  ratio <- abs(est_vec) / se_vec
  ratio <- ratio[is.finite(ratio)]
  min_ratio_val <- ifelse(length(ratio) == 0, NA_real_, min(ratio))
  
  base_dt[, `:=`(
    spearman_rho = as.numeric(ct_s$estimate),
    spearman_pvalue = ct_s$p.value,
    pearson_r = as.numeric(ct_p$estimate),
    pearson_pvalue = ct_p$p.value,
    min_est_over_se = min_ratio_val
  )]
  
  return(base_dt)
}

compute_gdm_correlations <- function(fit, reads_file, mode = "all", pop_number = "NA", n_sim = 10000, metadata_list = list()) {

  base_dt <- data.table(
    min_est_over_se = NA_real_,
    spearman_rho = NA_real_, spearman_pvalue = NA_real_,
    pearson_r = NA_real_, pearson_pvalue = NA_real_,
    alpha_x_A1_est = NA_real_, alpha_x_A1_SE = NA_real_,
    alpha_x_A2_est = NA_real_, alpha_x_A2_SE = NA_real_,
    alpha_x_out_est = NA_real_, alpha_x_out_SE = NA_real_,
    beta_x_A1_est = NA_real_, beta_x_A1_SE = NA_real_,
    beta_x_A2_est = NA_real_, beta_x_A2_SE = NA_real_,
    beta_x_out_est = NA_real_, beta_x_out_SE = NA_real_
  )
  
  if (length(metadata_list) > 0) {
    base_dt <- cbind(as.data.table(metadata_list), base_dt)
  }
  
  if (is.null(fit) || is.null(fit@estimate)) {
    return(base_dt)
  }
  
  reads <- fread(reads_file, select = c("x_A1", "x_A2", "total_reads", "population"))
  
  if (mode == "pop") {
    reads[, pop_id := sub("^([0-9]+)_.*", "\\1", population)]
    reads <- reads[pop_id == pop_number]
    reads[, pop_id := NULL]
  }
  if ("population" %in% names(reads)) {
    reads[, population := NULL]
  }
  
  reads[, x_out := total_reads - x_A1 - x_A2]
  reads[, total_reads := NULL]
  
  y <- as.matrix(reads[, .(x_A1, x_A2, x_out)])
  col_sums <- colSums(y)
  ord <- order(col_sums, decreasing = TRUE)
  sorted_names <- names(col_sums)[ord]
  rm(reads, y)
  
  est <- fit@estimate
  se <- fit@SE
  
  alpha_est <- setNames(rep(NA_real_, 3), sorted_names)
  beta_est  <- setNames(rep(NA_real_, 3), sorted_names)
  alpha_se  <- setNames(rep(NA_real_, 3), sorted_names)
  beta_se   <- setNames(rep(NA_real_, 3), sorted_names)

  for (nm in sorted_names[1:2]) {
    alpha_est[nm] <- safe_extract(est, paste0("alpha_", nm))
    beta_est[nm]  <- safe_extract(est, paste0("beta_", nm))
  }
  
  alpha_se[sorted_names[1]] <- se[[1]]
  alpha_se[sorted_names[2]] <- se[[2]]
  beta_se[sorted_names[1]]  <- se[[3]]
  beta_se[sorted_names[2]]  <- se[[4]]
  
  est_vec <- c(alpha_est, beta_est)
  se_vec  <- c(alpha_se, beta_se)
  all_params_vec <- c(est_vec, se_vec)

  param_names <- c("alpha_x_A1_est", "alpha_x_A1_SE", 
                   "alpha_x_A2_est", "alpha_x_A2_SE", 
                   "alpha_x_out_est", "alpha_x_out_SE",
                   "beta_x_A1_est",  "beta_x_A1_SE",  
                   "beta_x_A2_est",  "beta_x_A2_SE",  
                   "beta_x_out_est", "beta_x_out_SE")
  
  vals_to_assign <- list(
    alpha_est["x_A1"], alpha_se["x_A1"],
    alpha_est["x_A2"], alpha_se["x_A2"],
    alpha_est["x_out"], alpha_se["x_out"],
    beta_est["x_A1"],  beta_se["x_A1"],
    beta_est["x_A2"],  beta_se["x_A2"],
    beta_est["x_out"], beta_se["x_out"]
  )
  
  base_dt[, (param_names) := vals_to_assign]

  if (sum(!is.na(all_params_vec)) < 8) {
    return(base_dt)
  }
  
  y_sim <- rGDM3(n_sim, fit, sorted_names)
  dt_sim <- as.data.table(y_sim)
  
  ct_s <- suppressWarnings(cor.test(dt_sim$x_A1, dt_sim$x_A2, method = "spearman"))
  ct_p <- suppressWarnings(cor.test(dt_sim$x_A1, dt_sim$x_A2, method = "pearson"))
  
  ratio <- abs(est_vec) / se_vec
  ratio <- ratio[is.finite(ratio)]
  min_ratio_val <- ifelse(length(ratio) == 0, NA_real_, min(ratio))
  
  base_dt[, `:=`(
    spearman_rho = as.numeric(ct_s$estimate),
    spearman_pvalue = ct_s$p.value,
    pearson_r = as.numeric(ct_p$estimate),
    pearson_pvalue = ct_p$p.value,
    min_est_over_se = min_ratio_val
  )]
  
  return(base_dt)
}
