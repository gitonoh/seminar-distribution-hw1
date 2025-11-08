#!/usr/bin/env Rscript
# Compute mean(networth) and 95% CIs by a chosen household characteristic
# The script picks the first available variable from a candidate list (excluding hhsex and racecl5)

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

# normalize names
names(scf_summary_2022) <- tolower(names(scf_summary_2022))
names(scf_rw_2022_raw) <- tolower(names(scf_rw_2022_raw))

scf <- scf_summary_2022 %>% mutate(id_imp = as.numeric(str_sub(y1, -1, -1))) %>% arrange(yy1)
scf <- scf %>% mutate(wgt5 = wgt * 5, counts = 5)

# Candidate grouping variables (try in this order)
candidates <- c('ownhome','homeown','own_primary','own', 'x30022', 'region', 'age', 'agecl5', 'agecl', 'agegrp', 'educ', 'educyrs', 'edcl', 'hhage')
# exclude variables already used
exclude <- c('hhsex','racecl5')

found_var <- NULL
for(var in candidates){
  if(var %in% names(scf) && !(var %in% exclude)){
    found_var <- var
    break
  }
}

if(is.null(found_var)){
  # fallback: pick the first non-NA categorical-like variable not excluded
  others <- setdiff(names(scf), c('networth','wgt','y1','yy1','id_imp','wgt5','counts', exclude))
  # prefer factor/character variables
  for(var in others){
    if(is.character(scf[[var]]) || is.factor(scf[[var]])){
      found_var <- var; break
    }
  }
}

if(is.null(found_var)) stop('No suitable grouping variable found in dataset')

cat(sprintf('Selected grouping variable: %s\n', found_var))

# Prepare replicate weights like in hw1_solutions.R
rw <- scf_rw_2022_raw %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
wt_cols <- grep('^wt1b', names(rw), value = TRUE)
mm_cols <- grep('^mm', names(rw), value = TRUE)
R <- min(length(wt_cols), length(mm_cols))
if(R < 1) stop('No replicate weight columns found in replicate file')

rw_comb <- rw %>% arrange(yy1)
for(k in seq_len(R)){
  wcol <- paste0('wgt', k)
  rw_comb[[wcol]] <- rw_comb[[wt_cols[k]]] * rw_comb[[mm_cols[k]]]
}
rep_cols <- paste0('wgt', seq_len(R))
rw_comb <- rw_comb %>% select(yy1, all_of(rep_cols))

scf_full <- left_join(scf, rw_comb, by = 'yy1')

# If chosen var is numeric with many unique values, bin into quartiles
var_data <- scf_full[[found_var]]
group_var_name <- found_var
if(is.numeric(var_data) && length(unique(var_data[!is.na(var_data)])) > 10){
  # create 4 quantile bins
  scf_full[[paste0(found_var, '_grp')]] <- cut(var_data, breaks = unique(quantile(var_data, probs = seq(0,1,0.25), na.rm = TRUE)), include.lowest = TRUE)
  group_var_name <- paste0(found_var, '_grp')
}

# build imputation list and replicate matrix
scf_list_full <- split(scf_full, factor(scf_full$id_imp))
scf_list_full <- lapply(scf_list_full, function(df) df %>% arrange(yy1))
repweights_mat <- rw_comb[ , -1]

cat('Building svrepdesign (this may take a moment)...\n')
scf_design <- svrepdesign(weights = ~ wgt,
                          repweights = repweights_mat,
                          data = imputationList(scf_list_full),
                          scale = 1,
                          rscales = rep(1 / (R-1), R),
                          mse = FALSE,
                          type = 'other',
                          combined.weights = TRUE)

# compute means by group and combine
formula_group <- as.formula(paste('~', group_var_name))
mi_group_mean <- MIcombine(with(scf_design, svyby(~ networth, formula_group, svymean)))
est <- coef(mi_group_mean)
se  <- SE(mi_group_mean)

# build output data.frame
groups <- names(est)
df_out <- data.frame(group = groups, estimate = as.numeric(est), se = as.numeric(se), stringsAsFactors = FALSE)
df_out$lower95 <- df_out$estimate - 1.96 * df_out$se
df_out$upper95 <- df_out$estimate + 1.96 * df_out$se

out_csv <- paste0('c:/Users/takao/Documents/wu-seminar-distribution/hw1/group_by_', found_var, '_networth.csv')
write.csv(df_out, out_csv, row.names = FALSE)
cat(sprintf('Saved group results to %s\n', out_csv))

# Plot
plot_file <- paste0('c:/Users/takao/Documents/wu-seminar-distribution/hw1/group_by_', found_var, '_networth.png')
gg <- ggplot(df_out, aes(x = reorder(group, estimate), y = estimate/1e6)) +
  geom_col(fill = '#2b83ba') +
  geom_errorbar(aes(ymin = lower95/1e6, ymax = upper95/1e6), width = 0.2) +
  coord_flip() +
  labs(title = paste('Mean net worth by', found_var), x = found_var, y = 'Mean net worth (million USD)') +
  theme_minimal()
ggsave(plot_file, gg, width = 8, height = 4)
cat(sprintf('Saved bar chart to %s\n', plot_file))

cat('\nDone.\n')
