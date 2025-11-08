#!/usr/bin/env Rscript
# HW1 solutions: SCF-1 (part a)
# Compute mean household wealth (networth) and standard errors "by hand"
# i) unweighted classic SE
# ii) weighted classic SE
# iii) Rubins rule for implicates (use sampling weights for within-implicate variance)
# iv) replicate weights only within first implicate
# v) combine implicates + replicate weights

library(haven)
library(dplyr)
library(stringr)


# Define a function to download and import each stata file:

scf_dta_import <-
  function(this_url) {
    this_tf <- tempfile()
    download.file(this_url , this_tf , mode = 'wb')
    this_tbl <- read_dta(this_tf)
    this_df <- data.frame(this_tbl)
    file.remove(this_tf)
    names(this_df) <- tolower(names(this_df))
    this_df
  }

# URLs used in SCF Introduction.R
scf_summary_2022 <- scf_dta_import("https://www.federalreserve.gov/econres/files/scfp2022s.zip")
scf_rw_2022_raw  <- scf_dta_import("https://www.federalreserve.gov/econres/files/scf2022rw1s.zip")

cat("Data downloaded. Preparing variables...\n")

if(!all(c('networth','wgt','y1','yy1') %in% names(scf_summary_2022))) {
  stop('Expected variables networth, wgt, y1, yy1 not found in summary file')
}

# Create working data frame and implicate id
scf <- scf_summary_2022 %>%
  mutate(id_imp = as.numeric(str_sub(y1, -1, -1))) %>%
  arrange(yy1)

# create wgt5 (weights that sum to population when using each implicate)
scf <- scf %>% mutate(wgt5 = wgt * 5,
                      counts = 5)
### Part (i): Unweighted mean and classic SE
cat('\nPART (i) Unweighted mean and classic SE\n')
net <- scf$networth
valid <- !is.na(net)
n_obs <- sum(valid)
mean_unw <- mean(net, na.rm = TRUE)
se_unw <- sd(net, na.rm = TRUE) / sqrt(n_obs)
cat(sprintf('Unweighted: n = %d, mean = %g, se = %g\n', n_obs, mean_unw, se_unw))

### Part (ii): Weighted mean and "classic" weighted SE (naive formula using sampling weights)
cat('\nPART (ii) Weighted mean and classic weighted SE\n')
w <- scf$wgt
# weighted mean
mean_wtd <- sum(scf$networth * w, na.rm = TRUE) / sum(w[!is.na(scf$networth)])
# naive weighted SE: sqrt( sum( (w/sum(w))^2 * (x - mean_wtd)^2 ) )
den <- sum(w[!is.na(scf$networth)])
se_wtd <- sqrt( sum( (w[!is.na(scf$networth)]/den)^2 * ( (scf$networth[!is.na(scf$networth)] - mean_wtd)^2 ) ) )
cat(sprintf('Weighted (naive): mean = %g, se = %g\n', mean_wtd, se_wtd))

### Part (iii): Account for implicates (Rubin's rule) using sampling weights (no replicate weights yet)
cat('\nPART (iii) Rubin\'s rule across implicates (use sampling weights)\n')
M <- length(unique(scf$id_imp))
if(is.na(M) || M < 2) stop('Less than 2 implicates found; cannot apply Rubin\'s rule')

# compute theta (weighted mean) and within variance for each implicate
imp_stats <- scf %>%
  group_by(id_imp) %>%
  summarise(
    theta = sum(networth * w, na.rm = TRUE) / sum(w[!is.na(networth)]),
    # within-imputation variance (naive, using sampling weights like in the assignment)
    var_within = {
      idx <- !is.na(networth)
      wloc <- w[idx]
      xloc <- networth[idx]
      denom <- sum(wloc)
      # compute weighted variance in the "probability weights" sense
      # variance of weighted mean approx: sum( (w/denom)^2 * (x - theta)^2 )
      theta_loc <- sum(xloc * wloc) / denom
      sum((wloc/denom)^2 * (xloc - theta_loc)^2)
    },
    .groups = 'drop'
  )

theta_bar <- mean(imp_stats$theta)
W <- mean(imp_stats$var_within)             # average within-imputation variance
B <- (1/(M-1)) * sum( (imp_stats$theta - theta_bar)^2 )  # between-imputation variance
Tvar <- W + (1 + 1/M) * B
se_rubin_sampling <- sqrt(Tvar)

cat(sprintf('Implicates M = %d\n', M))
cat(sprintf('Theta per implicate: %s\n', paste(round(imp_stats$theta,3), collapse = ', ')))
cat(sprintf('Rubin: theta_bar = %g, SE = %g\n', theta_bar, se_rubin_sampling))

### Part (iv): Replicate weights only within first implicate
cat('\nPART (iv) Replicate weights (only for first implicate)\n')

# Prepare replicate weights table: NA -> 0
rw <- scf_rw_2022_raw %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

# detect replicate patterns: wt1b* and mm*
wt_cols <- grep('^wt1b', names(rw), value = TRUE)
mm_cols <- grep('^mm', names(rw), value = TRUE)
R <- min(length(wt_cols), length(mm_cols))
if(R < 1) stop('No replicate weight columns found')

cat(sprintf('Detected %d replicate columns (using first %d)\n', R, R))

# build combined replicate weights wgt1...wgtR and keep yy1
rw_comb <- rw %>% arrange(yy1)
for(k in seq_len(R)){
  wcol <- paste0('wgt', k)
  rw_comb[[wcol]] <- rw_comb[[wt_cols[k]]] * rw_comb[[mm_cols[k]]]
}
rep_cols <- paste0('wgt', seq_len(R))
rw_comb <- rw_comb %>% select(yy1, all_of(rep_cols))

# merge replicate weights to summary file by yy1
scf_full <- left_join(scf, rw_comb, by = 'yy1')

# use first implicate
first_imp <- min(scf_full$id_imp, na.rm = TRUE)
scf_imp1 <- scf_full %>% filter(id_imp == first_imp)
if(nrow(scf_imp1) == 0) stop('No rows for first implicate')

# point estimate for implicate 1 (using main sampling weight wgt)
theta_imp1 <- sum(scf_imp1$networth * scf_imp1$wgt, na.rm = TRUE) / sum(scf_imp1$wgt[!is.na(scf_imp1$networth)])

# compute replicate estimates theta_r for r in 1:R
theta_r <- numeric(R)
for(k in seq_len(R)){
  rwk <- scf_imp1[[paste0('wgt',k)]]
  denom <- sum(rwk[!is.na(scf_imp1$networth)])
  if(denom == 0) {
    theta_r[k] <- NA
  } else {
    theta_r[k] <- sum(scf_imp1$networth * rwk, na.rm = TRUE) / denom
  }
}
valid_theta_r <- !is.na(theta_r)
var_rep_imp1 <- mean( (theta_r[valid_theta_r] - theta_imp1)^2 )
se_rep_imp1 <- sqrt(var_rep_imp1)
cat(sprintf('Implicate %d: theta = %g, replicate-based SE = %g (R=%d used)\n', first_imp, theta_imp1, se_rep_imp1, sum(valid_theta_r)))

### Part (v): Combine implicates and replicate weights (Rubin-style using replicate variances as within variances)
cat('\nPART (v) Combine implicates + replicate weights\n')

# For each implicate: compute theta_m and replicate-based variance
imp_theta <- scf_full %>% group_by(id_imp) %>%
  summarise(theta_m = sum(networth * wgt, na.rm = TRUE) / sum(wgt[!is.na(networth)]), .groups='drop')

replicate_variance_for_implicate <- function(df_imp, R, rep_prefix = 'wgt'){
  theta_hat <- sum(df_imp$networth * df_imp$wgt, na.rm = TRUE) / sum(df_imp$wgt[!is.na(df_imp$networth)])
  thetas_r <- numeric(R)
  for(k in seq_len(R)){
    rwk <- df_imp[[paste0(rep_prefix,k)]]
    denom <- sum(rwk[!is.na(df_imp$networth)])
    if(denom == 0) {
      thetas_r[k] <- NA
    } else {
      thetas_r[k] <- sum(df_imp$networth * rwk, na.rm = TRUE) / denom
    }
  }
  thetas_r <- thetas_r[!is.na(thetas_r)]
  if(length(thetas_r) == 0) return(NA_real_)
  mean((thetas_r - theta_hat)^2)
}

var_m_vec <- numeric(M)
for(m in seq_len(M)){
  dfm <- scf_full %>% filter(id_imp == m)
  var_m_vec[m] <- replicate_variance_for_implicate(dfm, R)
}

W_rep <- mean(var_m_vec, na.rm = TRUE)   # average within-imputation variance from replicates
B_imp <- (1/(M-1)) * sum( (imp_theta$theta_m - mean(imp_theta$theta_m))^2 )
Tvar_combined <- W_rep + (1 + 1/M) * B_imp
se_combined <- sqrt(Tvar_combined)

cat(sprintf('Combined result: theta_bar = %g\n', mean(imp_theta$theta_m)))
cat(sprintf('Within (replicate) variance W = %g\n', W_rep))
cat(sprintf('Between variance B = %g\n', B_imp))
cat(sprintf('Total variance T = %g, SE = %g\n', Tvar_combined, se_combined))

cat('\nDone. Results summary:\n')
res_tab <- tibble::tibble(
  method = c('Unweighted (i)', 'Weighted naive (ii)', 'Rubin sampling (iii)', 'Replicates imp1 (iv)', 'Rubin + replicates (v)'),
  estimate = c(mean_unw, mean_wtd, theta_bar, theta_imp1, mean(imp_theta$theta_m)),
  se = c(se_unw, se_wtd, se_rubin_sampling, se_rep_imp1, se_combined)
)
print(res_tab)

# ---------------------------
# Part (c): use survey + mitools for correct mean & SE
# ---------------------------
cat('\nPART (c) Using survey + mitools to compute mean and SE\n')
library(survey)
library(mitools)

# Build a list of implicates (each implicate is a dataframe with replicate weights merged)
scf_list_full <- split(scf_full, factor(scf_full$id_imp))
scf_list_full <- lapply(scf_list_full, function(df) df %>% arrange(yy1))

# Prepare replicate weights matrix (rows by unique yy1) as in SCF Introduction.R
# Use the combined replicate table we already built: rw_comb (yy1 + wgt1..wgtR)
repweights_mat <- rw_comb[ , -1]

# Create svrepdesign using the provided replicate weights and the imputation list
scf_design <- svrepdesign(weights = ~ wgt,
                          repweights = repweights_mat,
                          data = imputationList(scf_list_full),
                          scale = 1,
                          rscales = rep(1 / (R-1), R),
                          mse = FALSE,
                          type = "other",
                          combined.weights = TRUE)

# compute mean(networth) with survey functions across implicates and combine with mitools
svy_mean_list <- with(scf_design, svymean(~ networth))
mi_mean <- MIcombine(svy_mean_list)

cat('survey+mitools point estimate and SE:\n')
cat(sprintf('Estimate = %g, SE = %g\n', coef(mi_mean), SE(mi_mean)))

# add survey result to results table and print
res_tab <- dplyr::bind_rows(res_tab,
                            tibble::tibble(method = 'survey+mitools (c)',
                                           estimate = as.numeric(coef(mi_mean)),
                                           se = as.numeric(SE(mi_mean))))
print(res_tab)

# ---------------------------
# Problem 2(a): mean(networth) and SE by group (hhsex and racecl5), and weighted N by group
# Using survey + mitools (recommended)
# ---------------------------
cat('\nProblem 2(a): mean(networth) and SE by `hhsex` and `racecl5` plus weighted counts\n')

# 1) By hhsex
if('hhsex' %in% names(scf_full)){
  mi_hhsex_mean <- MIcombine(with(scf_design, svyby(~ networth, ~ hhsex, svymean)))
  cat('\nMean(networth) by hhsex (survey+mitools):\n')
  print(coef(mi_hhsex_mean))
  print(SE(mi_hhsex_mean))
  # weighted population counts by hhsex (use counts=5 per row so totals reflect population counts)
  mi_hhsex_pop <- MIcombine(with(scf_design, svyby(~ counts, ~ hhsex, svytotal)))
  cat('\nWeighted population counts (by hhsex):\n')
  print(coef(mi_hhsex_pop))
} else {
  warning('hhsex not found in data; skipping hhsex group calculations')
}

# 2) By racecl5
if('racecl5' %in% names(scf_full)){
  mi_race_mean <- MIcombine(with(scf_design, svyby(~ networth, ~ racecl5, svymean)))
  cat('\nMean(networth) by racecl5 (survey+mitools):\n')
  print(coef(mi_race_mean))
  print(SE(mi_race_mean))
  mi_race_pop <- MIcombine(with(scf_design, svyby(~ counts, ~ racecl5, svytotal)))
  cat('\nWeighted population counts (by racecl5):\n')
  print(coef(mi_race_pop))
} else {
  warning('racecl5 not found in data; skipping racecl5 group calculations')
}
 
# Save group results to CSV files for convenience
if(exists('mi_hhsex_mean')){
  hh_coef <- coef(mi_hhsex_mean)
  hh_se <- SE(mi_hhsex_mean)
  hh_pop <- as.numeric(coef(mi_hhsex_pop))
  hh_df <- data.frame(group = names(hh_coef), estimate = as.numeric(hh_coef), se = as.numeric(hh_se), pop = hh_pop)
  write.csv(hh_df, file = 'c:/Users/takao/Documents/wu-seminar-distribution/hw1/hhsex_networth_by_group.csv', row.names = FALSE)
}
if(exists('mi_race_mean')){
  race_coef <- coef(mi_race_mean)
  race_se <- SE(mi_race_mean)
  race_pop <- as.numeric(coef(mi_race_pop))
  race_df <- data.frame(group = names(race_coef), estimate = as.numeric(race_coef), se = as.numeric(race_se), pop = race_pop)
  write.csv(race_df, file = 'c:/Users/takao/Documents/wu-seminar-distribution/hw1/racecl5_networth_by_group.csv', row.names = FALSE)
}

# Save results table to CSV
results_csv <- "c:/Users/takao/Documents/wu-seminar-distribution/hw1/hw1_results.csv"
write.csv(res_tab, results_csv, row.names = FALSE)
cat(sprintf('\nSaved summary results to %s\n', results_csv))

# End of script
