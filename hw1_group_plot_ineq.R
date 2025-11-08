#!/usr/bin/env Rscript
# Compute mean(log(networth)) by homeownership (own) using survey + mitools across implicates
# Saves CSV and PNG

library(haven)
library(dplyr)
library(stringr)
library(survey)
library(mitools)
library(ggplot2)

scf_dta_import <- function(this_url) {
  this_tf <- tempfile()
  download.file(this_url, this_tf, mode = 'wb')
  this_tbl <- haven::read_dta(this_tf)
  this_df <- data.frame(this_tbl)
  file.remove(this_tf)
  names(this_df) <- tolower(names(this_df))
  this_df
}

cat('Downloading SCF summary and replicate weights...\n')
scf_summary_2022 <- scf_dta_import('https://www.federalreserve.gov/econres/files/scfp2022s.zip')
scf_rw_2022_raw  <- scf_dta_import('https://www.federalreserve.gov/econres/files/scf2022rw1s.zip')

names(scf_summary_2022) <- tolower(names(scf_summary_2022))
names(scf_rw_2022_raw) <- tolower(names(scf_rw_2022_raw))

scf <- scf_summary_2022 %>% mutate(id_imp = as.numeric(str_sub(y1, -1, -1))) %>% arrange(yy1)

# Make sure variable 'own' exists
if(!('own' %in% names(scf))) stop('Variable own not found in summary file')

# We'll compute log(networth) after shifting to positive values (to handle negatives)
min_nw <- min(scf$networth, na.rm = TRUE)
shift <- 0
if(min_nw <= 0) shift <- abs(min_nw) + 1
cat(sprintf('Shifting networth by %g to compute log(networth + shift)\n', shift))
scf <- scf %>% mutate(log_net = log(networth + shift))

# Prepare replicate weights
rw <- scf_rw_2022_raw %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
wt_cols <- grep('^wt1b', names(rw), value = TRUE)
mm_cols <- grep('^mm', names(rw), value = TRUE)
R <- min(length(wt_cols), length(mm_cols))
if(R < 1) stop('No replicate weight columns found')

rw_comb <- rw %>% arrange(yy1)
for(k in seq_len(R)){
  wcol <- paste0('wgt', k)
  rw_comb[[wcol]] <- rw_comb[[wt_cols[k]]] * rw_comb[[mm_cols[k]]]
}
rep_cols <- paste0('wgt', seq_len(R))
rw_comb <- rw_comb %>% select(yy1, all_of(rep_cols))

scf_full <- left_join(scf, rw_comb, by = 'yy1')

# Build imputation list
scf_list_full <- split(scf_full, factor(scf_full$id_imp))
scf_list_full <- lapply(scf_list_full, function(df) df %>% arrange(yy1))
repweights_mat <- rw_comb[ , -1]

cat('Building svrepdesign...\n')
scf_design <- svrepdesign(weights = ~ wgt,
                          repweights = repweights_mat,
                          data = imputationList(scf_list_full),
                          scale = 1,
                          rscales = rep(1 / (R-1), R),
                          mse = FALSE,
                          type = 'other',
                          combined.weights = TRUE)

# compute mean(log_net) by own
mi_log_mean <- MIcombine(with(scf_design, svyby(~ log_net, ~ own, svymean)))
est <- coef(mi_log_mean)
se <- SE(mi_log_mean)

groups <- names(est)
df_out <- data.frame(group = groups, estimate = as.numeric(est), se = as.numeric(se), stringsAsFactors = FALSE)
df_out$lower95 <- df_out$estimate - 1.96 * df_out$se
df_out$upper95 <- df_out$estimate + 1.96 * df_out$se

out_csv <- 'c:/Users/takao/Documents/wu-seminar-distribution/hw1/group_by_own_log_networth.csv'
write.csv(df_out, out_csv, row.names = FALSE)
cat(sprintf('Saved group log-networth results to %s\n', out_csv))

# For plotting, convert back to levels by exponentiating and subtracting shift: E[log(X+shift)] isn't E[X], but plotting mean log is an inequality indicator.
plot_file <- 'c:/Users/takao/Documents/wu-seminar-distribution/hw1/group_by_own_log_networth.png'
# Plot mean log values (not exponentiated)
gg <- ggplot(df_out, aes(x = factor(group), y = estimate)) +
  geom_col(fill = '#d6604d') +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2) +
  labs(title = 'Mean log(networth + shift) by homeownership', x = 'own', y = 'Mean log(networth + shift)') +
  theme_minimal()

ggsave(plot_file, gg, width = 7, height = 4)
cat(sprintf('Saved inequality plot to %s\n', plot_file))

cat('Done.\n')
