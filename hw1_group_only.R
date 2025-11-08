library(haven); library(dplyr); library(stringr); library(survey); library(mitools)

scf_dta_import <- function(this_url) {
  this_tf <- tempfile()
  download.file(this_url, this_tf, mode = 'wb')
  this_tbl <- read_dta(this_tf)
  this_df <- data.frame(this_tbl)
  file.remove(this_tf)
  names(this_df) <- tolower(names(this_df))
  this_df
}
scf_summary_2022 <- scf_dta_import("https://www.federalreserve.gov/econres/files/scfp2022s.zip")
scf_rw_2022_raw  <- scf_dta_import("https://www.federalreserve.gov/econres/files/scf2022rw1s.zip")

scf <- scf_summary_2022 %>% mutate(id_imp = as.numeric(str_sub(y1, -1, -1))) %>% arrange(yy1) %>% mutate(wgt5 = wgt*5, counts = 5)

rw <- scf_rw_2022_raw %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
wt_cols <- grep('^wt1b', names(rw), value = TRUE)
mm_cols <- grep('^mm', names(rw), value = TRUE)
R <- min(length(wt_cols), length(mm_cols))
rw_comb <- rw %>% arrange(yy1)
for(k in seq_len(R)){
  rw_comb[[paste0('wgt',k)]] <- rw_comb[[wt_cols[k]]] * rw_comb[[mm_cols[k]]]
}
rw_comb <- rw_comb %>% select(yy1, starts_with('wgt'))
scf_full <- left_join(scf, rw_comb, by='yy1')

scf_list_full <- split(scf_full, factor(scf_full$id_imp))
scf_list_full <- lapply(scf_list_full, function(df) df %>% arrange(yy1))
repweights_mat <- rw_comb[ , -1]

scf_design <- svrepdesign(weights = ~ wgt,
                          repweights = repweights_mat,
                          data = imputationList(scf_list_full),
                          scale = 1,
                          rscales = rep(1 / (R-1), R),
                          mse = FALSE,
                          type = "other",
                          combined.weights = TRUE)

# hhsex
cat('hhsex mean and SE:\n')
mi_hh <- MIcombine(with(scf_design, svyby(~ networth, ~ hhsex, svymean)))
print(coef(mi_hh))
print(SE(mi_hh))
cat('hhsex pop counts:\n')
mi_hh_pop <- MIcombine(with(scf_design, svyby(~ counts, ~ hhsex, svytotal)))
print(coef(mi_hh_pop))

# racecl5
if('racecl5' %in% names(scf_full)){
  cat('racecl5 mean and SE:\n')
  mi_r <- MIcombine(with(scf_design, svyby(~ networth, ~ racecl5, svymean)))
  print(coef(mi_r))
  print(SE(mi_r))
  cat('racecl5 pop counts:\n')
  mi_r_pop <- MIcombine(with(scf_design, svyby(~ counts, ~ racecl5, svytotal)))
  print(coef(mi_r_pop))
} else {
  cat('racecl5 not present\n')
}
