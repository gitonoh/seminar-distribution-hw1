# HW1 (SCF-1) — Part (b) report

This short report summarizes why we use sampling weights, implicates, and replicate weights in the SCF and explains how point estimates and standard errors change across the stepwise calculations you performed in part (a). It also reproduces the numeric results from the `hw1_solutions.R` run and contains quick instructions for reproducing the results locally.

## Where to find the code
- `hw1_solutions.R` — script that performs parts (i)–(v) step-by-step and prints a small results table. Path: `c:/Users/takao/Documents/wu-seminar-distribution/hw1/hw1_solutions.R`

## Quick reproduction (PowerShell / Windows)
1. Install required R packages (one-off):

```powershell
Rscript -e "install.packages(c('haven','dplyr','stringr','tibble'), repos='https://cloud.r-project.org')"
```

2. Run the script:

```powershell
Rscript "c:/Users/takao/Documents/wu-seminar-distribution/hw1/hw1_solutions.R"
```

The script downloads the SCF summary and replicate-weight files from the Federal Reserve website (internet required) and prints the results to the console.

## Numeric results (run performed in this workspace)

The run produced the following table (estimate units: US dollars):

| Method | Estimate (mean networth) | SE |
|---|---:|---:|
| (i) Unweighted classic | 19,956,403 | 726,839 |
| (ii) Weighted naive | 1,059,457 | 15,011.7 |
| (iii) Rubin (sampling weights for within) | 19,547,139 | 932,432 |
| (iv) Replicate weights (implicate 1) | 1,041,815 | 78,574.7 |
| (v) Rubin + replicates (combined) | 1,059,457 | 97,073 |

Notes:
- The unweighted estimate (i) treats the sample as if it were an SRS and gives a sample-level average; because very wealthy households are present and potentially over-sampled, the unweighted average is very large.
- The weighted estimate (ii) uses `wgt` to produce a population-representative average; this is the usual estimand for SCF population descriptions and is much smaller than the unweighted sample mean.
- (iii) shows Rubin's combination across M = 5 implicates using sampling-weight-based within variances; the example illustrates how imputation uncertainty increases variance when imputations differ.
- (iv) shows replicate-based SE for implicate 1 computed from the SCF replicate weights (R = 999 used). This captures design-based sampling variability and is larger than the naive weighted SE.
- (v) shows the recommended combination: compute replicate-based within-imputation variances for each implicate, then use Rubin's formula to combine (adds between-imputation variance). The result is the point estimate on the population-weighted scale and a design+imputation-aware SE.

## Conceptual explanation (condensed)
Below is a concise, practical explanation of the three building blocks used in parts (a)–(c) and how they affect point estimates and standard errors.

- Sampling weights (`wgt`): adjust for unequal selection probabilities, nonresponse and post-stratification so that estimates generalize to the target population. Using weights changes the estimand from a sample-average to a population-average and typically reduces the influence of over-sampled subgroups (e.g., very wealthy households in SCF).

- Implicates (multiple imputation): the SCF provides M completed datasets to reflect uncertainty about imputed values (especially for wealth components). Compute your statistic in each implicate and combine with Rubin's rules: the overall point estimate is the average of implicate estimates; the total variance equals the average within-imputation variance plus (1 + 1/M) times the between-imputation variance. Ignoring implicates underestimates uncertainty.

- Replicate weights: capture sampling variance due to the complex survey design (strata, clustering, multi-stage selection). For each set of replicate weights recompute the estimator; the variability across replicate estimates (with the replication-specific scaling) gives a design-correct estimate of variance. Replicate-based SEs are generally larger than naive weighted SEs because they account for design effects.

How estimates and SEs change when you add each component:

- Unweighted → Weighted: the point estimate moves to reflect the population (often smaller for SCF networth); SE computed naively for the sample is not appropriate for the population estimand.
- Add implicates (Rubin): point estimate becomes the implicate average; SE increases because between-imputation variability is added to within-imputation variance.
- Add replicate weights: SE becomes design-consistent and typically increases relative to naive formulas; combining replicate-based within-imputation variances with Rubin's rule yields the full design+imputation SE reported by `survey`+`mitools`.

Recommended practice: compute population-weighted estimates, use all implicates and the provided replicate weights, and prefer the `survey`+`mitools` implementation (or the `convey` workflow for inequality measures) to obtain design- and imputation-correct standard errors.

## Note on discrepancy between manual and `survey` results

During part (v) you saw a hand-coded combined SE of about 97,073 while the `survey` + `mitools` calculation in part (c) returned a smaller SE of about 46,960. The main reasons for this difference are:

- Scaling and replication rules: the `survey` package applies specific replication scaling (rscales/scale) and internal variance formulas depending on the replication type (bootstrap, BRR, jackknife). My manual part (v) used a straightforward mean-of-squared-differences across replicate estimates as the within-imputation variance; that does not necessarily match the exact scaling `survey` applies for SCF replicate weights.
- Combined.weights / multiplication factors: SCF replicate weights in the `scf_rw_2022` table include multiplication factors. When building the `svrepdesign` the `combined.weights = TRUE` argument and any rscales influence how the replicates get normalized. Small differences in how you compute the denominators produce noticeable differences in SE for high-variance quantities.
- Implementation details: `mitools::MIcombine()` and `survey::svrepdesign()` perform the canonical combination of imputation and replication variances; they are the recommended approach when you want results consistent with Fed publications.

## Group results — race (variable: `racecl5`)

The SCF public recode `racecl5` groups households into five race/ethnicity categories (public codebook recodes). Below I show the point estimate (mean net worth), the standard error (design+imputation combined, from `survey`+`mitools`), and the weighted population total for each group as produced by the analysis scripts. The label mapping used here follows the SCF public codebook conventions.

| racecl5 | Label | Mean net worth (USD) | SE | Weighted population |
|---:|---|---:|---:|---:|
| 1 | White (non-Hispanic) | 1,361,785.11 | 63,203.01 | 87,725,120.04 |
| 2 | Black / African-American | 211,511.74 | 39,885.20 | 15,118,564.71 |
| 3 | Hispanic / Latino | 227,523.38 | 26,517.11 | 12,284,803.88 |
| 4 | Asian | 1,810,102.26 | 319,343.65 | 5,181,423.80 |
| 5 | Other (AI/AN, NH/PI, Other) | 389,473.83 | 51,448.57 | 10,996,476.95 |

Interpretation: White and Asian households show much higher mean net worth than Black, Hispanic, or Other households in these SCF estimates. Asian households have the highest point estimate but also the largest SE (reflecting smaller sample size and more variability). Black and Hispanic households show substantially lower mean net worth with relatively smaller absolute SEs but large relative uncertainty compared with their point estimates. These differences reflect persistent wealth disparities and also larger uncertainty for smaller subgroups in the SCF public sample.


## Group results — homeownership (variable: `own`)

We selected a different household characteristic available in the public SCF summary file: `own` (homeownership indicator). Below are the exact MIcombined point estimates and 95% confidence intervals computed with `survey` + `mitools` (replicate-weight + imputation-aware):

| own | Label | Mean net worth (USD) | SE | 95% CI |
|---:|---|---:|---:|---|
| 0 | Non-owner (renter/other) | 487,107.06 | 70,743.68 | [348,449.44, 625,764.67] |
| 1 | Owner | 1,150,427.60 | 52,092.55 | [1,048,326.20, 1,252,528.99] |

![Mean net worth by homeownership (owners vs non-owners)](group_by_own_networth.png)

*Figure: Mean household net worth (millions USD) by homeownership status. Error bars show 95% confidence intervals (design + imputation aware).* 

Interpretation: Owner households have a substantially higher mean net worth (~$1.15M) than non-owners (~$0.49M). The 95% CIs do not overlap, indicating a statistically meaningful difference in mean net worth between owners and non-owners in these SCF estimates. Note that homeowners both tend to have higher accumulated wealth and are oversampled in SCF designs; the reported SEs account for design and imputation uncertainty. Smaller subgroup sample sizes and the heavy-tailed wealth distribution mean these CI widths should be considered when comparing subgroups.

## Inequality indicator — mean log(networth + shift) by homeownership

As a simple inequality-related indicator, I computed the mean of log(networth + shift) (shift chosen so the argument of the log is positive) by `own` (homeownership). Mean log networth highlights differences in central tendency on the log scale and is commonly used to summarize skewed wealth distributions.

| own | Mean log(networth + shift) | SE | 95% CI |
|---:|---:|---:|---:|
| 0 | 13.49021 | 0.02718 | [13.43694, 13.54347] |
| 1 | 13.82443 | 0.01772 | [13.78970, 13.85915] |

![Mean log(networth + shift) by homeownership](/group_by_own_log_networth.png)

Interpretation: Owners have a higher mean log(networth + shift) than non-owners, consistent with the earlier level comparisons; the difference is meaningful and the confidence intervals are tight. Mean log differences compress extreme values and therefore give a sense of typical proportional differences rather than absolute-dollar gaps.


