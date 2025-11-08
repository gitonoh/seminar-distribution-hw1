
#-----------------------------------#
#--- SURVEY OF CONSUMER FINANCES ---#
#-----------------------------------#

# The SCF studies net worth across the United States by asking respondents about both active and passive income, mortgages, pensions, credit card debt, even car leases. Administered by the Board of Governors of the Federal Reserve System triennially since 1989, this complex sample survey generalizes to the civilian non-institutional population and comprehensively assesses household wealth.

#This section downloads, imports, and prepares the most current microdata for analysis, then reproduces some statistics and margin of error terms from the Federal Reserve.

#This survey uses a multiply-imputed variance estimation technique described in the 2004 Codebook. Most users do not need to study this function carefully. Define a function specific to only this dataset:

# this is a custom implementation of Rubin's Rule for combining estimates from multiple imputation (MI) datasets
# it combines point estimates and variances from multiple imputed analyses (typically from a with() call on an imputationList) into a single result with proper variance, df and missing information measures

#---------------------------------------------------------------------------------------------------------------------------------------------#
#install.packages("haven")
#install.packages("tidyverse")
#install.packages("survey")
#install.packages("mitools")
#install.packages("dineq")
#install.packages("Hmisc")
#install.packages("convey")

library(haven)
library(tidyverse)
library(survey)
library(mitools)
library(dineq)
library(Hmisc)
library(convey)

#---------------------------------------------------------------------------------------------------------------------------------------------#

# STEP 1 - get the data ####

## 1.A) via zip ####

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
#---------------------------------------------------------------------------------------------------------------------------------------------#
# Download and import the full, summary extract, and replicate weights tables:

# main survey data
scf_2022 <-
  scf_dta_import("https://www.federalreserve.gov/econres/files/scf2022s.zip")
# summary file (summary variables used in FED report)
scf_summary_2022 <-
  scf_dta_import("https://www.federalreserve.gov/econres/files/scfp2022s.zip")
# replicate weights
scf_rw_2022 <-
  scf_dta_import("https://www.federalreserve.gov/econres/files/scf2022rw1s.zip")

# Confirm both the full public data and the summary extract contain five records per family:

stopifnot(nrow(scf_2022) == nrow(scf_rw_2022) * 5)
stopifnot(nrow(scf_2022) == nrow(scf_summary_2022))
#---------------------------------------------------------------------------------------------------------------------------------------------#

## 1.B) Download and save locally ####

# url: https://www.federalreserve.gov/econres/scfindex.htm 

PATH <- "C:/Users/alhuber/OneDrive - WU Wien/Dokumente/Alexander Huber/Teaching/SCF"

# main survey data
scf_2022 <- read_dta(paste0(PATH, "/data/p22i6.dta"))
# summary file (summary variables used in FED report)
scf_summary_2022 <- read_dta(paste0(PATH, "/data/rscfp2022.dta"))
# replicate weights
scf_rw_2022 <- read_dta(paste0(PATH, "/data/p22_rw1.dta"))

# Confirm both the full public data and the summary extract contain five records per family:

stopifnot(nrow(scf_2022) == nrow(scf_rw_2022) * 5)
stopifnot(nrow(scf_2022) == nrow(scf_summary_2022))
#---------------------------------------------------------------------------------------------------------------------------------------------#

# STEP 2 - overview ####

# variable names
names(scf_2022) # full dataset, impossible to read without codebook https://www.federalreserve.gov/econres/files/codebk2022.txt 
names(scf_summary_2022) # variables summarised by FED, easier to work with https://www.federalreserve.gov/econres/files/bulletin.macro.txt 
names(scf_rw_2022) # replicate weights and multiplication factors (how often a PEU was selected for the replicate)

# number of households in the sample
scf_summary_2022$y1 %>% unique() %>% length()  # all HH's (PEU's) * 5 implicates
scf_summary_2022$yy1 %>% unique() %>% length() # total number of HH's in the survey

# total weighted number of households in the US
scf_summary_2022 %>% summarise(sum(wgt)) 

#---------------------------------------------------------------------------------------------------------------------------------------------#

# STEP 3 - working DF ####

scf <- scf_summary_2022

# some examples

# what are we weighting for
scf %>% summarise(mean_networth = mean(networth))                   # sample mean of networth
scf %>% summarise(mean_networth_wtd = sum(networth*wgt) / sum(wgt)) # weighted mean/population mean by hand
scf %>% summarise(mean_networth_wtd = weighted.mean(networth, wgt)) # via function from stats

# quantiles (with Hmisc)
scf %>% summarise(median_income_wtd = Hmisc::wtd.quantile(income, wgt, 0.5))
scf %>% summarise(median_networth_wtd = Hmisc::wtd.quantile(networth, wgt, 0.5))

# from now on we're only interested in weighted statistics

# mean netwealth by sex (of HH reference person)
scf %>% filter(famstruct %in% c(2,3)) %>%         # only for HH not married and without children
  group_by(hhsex) %>% summarise( mean_networth=sum(networth*wgt)/sum(wgt))


## 3.A) merge new variables ####

# NEW VARIABLES: imp.id = Number of the implicate (1:5) and wgt5 = such that HH's in each implicate sum up to total number of HH's in U.S.
scf <- scf %>% mutate(id.imp = as.numeric(str_sub(y1,-1 ,-1)),
                      wgt5 = wgt*5)
scf %>% group_by(id.imp) %>% summarise(n=n(), sum(wgt), sum(wgt5))


# merge additional variables from scf_2022 
# (look them up in codebook https://www.federalreserve.gov/econres/files/codebk2022.txt )

merge_scf <- scf_2022 %>% select(y1, 
                                 x5801, # inheritance received
                                 x5819) # inheritance expected
scf <- left_join(scf, merge_scf, by = "y1")

## 3.B) create quantiles ####

# option 1: CDF (cumulative distribution function)
scf <- scf %>%
  group_by(id.imp) %>%    # separately for each implicate
  arrange(id.imp,networth) %>%   # arrange HH's by networth
  mutate(cdf_networth = cumsum(wgt5)/sum(wgt5),   # computing the cumulative distribution function (CDF)
         cdf_deciles = cut(cdf_networth, breaks = seq(from=0,to=1,by=0.1), labels = 1:10),         # cut into 10 deciles
         cdf_percentiles = cut(cdf_networth, breaks = seq(from=0,to=1,by=0.01), labels = 1:100))   # cut into 100 percentiles

View(scf %>% select(id.imp,y1,networth,cdf_networth,cdf_deciles,cdf_percentiles))

#graphically
ggplot(scf, aes(x=cdf_networth,y=networth, color=factor(id.imp))) + geom_line()
ggplot(scf, aes(x=cdf_networth,y=log(networth), color=factor(id.imp))) + geom_line()
ggplot(scf %>% filter(cdf_networth<0.99), aes(x=cdf_networth,y=networth, color=factor(id.imp))) + geom_line()
ggplot(scf %>% filter(cdf_networth<0.9), aes(x=cdf_networth,y=networth, color=factor(id.imp))) + geom_line()


# assign to regions in the wealth distribution
scf <- scf %>%
  mutate(bottom_50 = ifelse(cdf_networth<=0.5,T,F),
         top_1 = ifelse(cdf_networth>0.99,T,F))

# option 2: weighted quantiles (dineq)
scf <- scf %>%
  group_by(id.imp) %>%
  mutate(cdf_deciles_II = dineq::ntiles.wtd(networth,10,weights=wgt5),
         cdf_percentiles_II = dineq::ntiles.wtd(networth,100,weights=wgt5))
scf %>% group_by(id.imp, cdf_deciles_II) %>% summarise(min = min(networth),
                                                       mean = weighted.mean(networth,wgt5),
                                                       max = max(networth))

stopifnot(scf$cdf_deciles == scf$cdf_deciles_II) # not exactly the same
View(scf %>% select(id.imp,y1,networth,cdf_networth,cdf_deciles,cdf_deciles_II) %>% mutate(same = ifelse(cdf_deciles==cdf_deciles_II,1,0)))
# only 47 observations are not assigned to the same decile --> this is normal due to different methods of computation

#---------------------------------------------------------------------------------------------------------------------------------------------#

# STEP 4 - computations ####

# so far, we only have computed the overall mean or median without taking into account the 5 implicates
# to do it correctly, we have to apply Rubin's rule

## 4.A) Rubin's Rule ####

# mean for each implicate
scf %>% 
  group_by(id.imp) %>%                                     # group by implicate
  summarise(mean_income = sum(income*wgt5)/sum(wgt5),      # weighted mean for income
            mean_networth = sum(networth*wgt5)/sum(wgt5),  # weighted mean for wealth
            number_hhs = sum(wgt5))                        # total number of households in the U.S. population

# Rubin's rule
scf %>% 
  group_by(id.imp) %>%
  summarise(mean_income = sum(income*wgt5)/sum(wgt5),
            mean_networth = sum(networth*wgt5)/sum(wgt5),
            number_hhs = sum(wgt5)) %>%
  summarise_all(mean) # gives us the mean over all 5 implicates

# which weight to use?
scf %>% group_by(id.imp) %>% summarise(mean_income = sum(income*wgt)/sum(wgt),
                                       mean_income_5 = sum(income*wgt5)/sum(wgt5),
                                       total_income = sum(income*wgt),
                                       total_income_5 = sum(income*wgt5),
                                       number_hhs = sum(wgt),
                                       number_hhs_5 = sum(wgt5)) %>%
  summarise_all(mean)
# --> for mean and median: does not matter which weight
# --> for sums and counts: always use the weight that sums up to the total population (wgt5 in this case)

## 4.B) (Weighted) Graphs ####

plot <- scf %>%
  group_by(id.imp, cdf_deciles) %>% 
  summarise(mean_wealth = sum(networth*wgt5)/sum(wgt5)) %>%
  group_by(cdf_deciles) %>%
  summarise_all(mean)

### barplot 1 ####
# with wealth deciles 
plot1 <- ggplot(plot, aes(x=cdf_deciles,y=mean_wealth)) +
  geom_bar(stat="identity",position="dodge",fill="steelblue") +
  labs(
    title = "Mean Wealth by Percentile",
    x = "CDF Wealth Quantiles",
    y = "Mean Wealth"
  ) +
  theme_light()
print(plot1)  

# add mean and median (with Rubin's rule)
mean_wealth <- scf %>% group_by(id.imp) %>% summarise(mean_wealth = sum(networth*wgt5)/sum(wgt5)) %>% summarise_all(mean)
median_wealth <- scf %>% group_by(id.imp) %>% summarise(median_wealth = Hmisc::wtd.quantile(networth, wgt5, 0.5)) %>% summarise_all(mean)

### barplot 2##### 
# with wealth deciles and horizontal lines for mean and median 
plot2 <- ggplot(plot, aes(x=cdf_deciles,y=mean_wealth)) +
  geom_bar(stat="identity",position="dodge",fill="steelblue") +
  geom_hline(yintercept = mean_wealth$mean_wealth, color="red", linetype="dashed", linewidth=1) +
  geom_hline(yintercept = median_wealth$median_wealth, color="darkgreen", linetype="dashed", linewidth=1) +
  annotate("text", x = length(unique(plot$cdf_deciles)) / 4, y = mean_wealth$mean_wealth, label = paste("Mean: ", round(mean_wealth$mean_wealth, 1)), color = "red", vjust = -1) +
  annotate("text", x = length(unique(plot$cdf_deciles)) / 3, y = median_wealth$median_wealth, label = paste("Median: ", round(median_wealth$median_wealth, 1)), color = "darkgreen", vjust = 1.5) +
  labs(
    title = "Mean Wealth by Percentile",
    x = "CDF Wealth Quantiles",
    y = "Mean Wealth"
  ) +
  theme_light()
print(plot2)  

### density plot ##### 
# right-skewed distributions + role of weights
ggplot(scf, aes(x = networth)) +
  geom_density(aes(x = networth, fill = "Unweighted"),  alpha = 0.4,  linetype = "solid") +
  geom_density(aes(x = networth, weight = wgt5, fill = "Weighted"),alpha = 0.6,  linetype = "dashed") +
  theme_classic() +
  scale_x_log10() +
  scale_fill_manual(values = c("Unweighted" = "darkred", "Weighted" = "darkgreen"))

#---------------------------------------------------------------------------------------------------------------------------------------------#

# STEP 5 - Replicate weights ####

# For point estimates (mean, median, etc.) it is fine to work like we did up to now
# For everything that involves standard errors, variances, confidence intervals, etc. we need to take it one step further
# --> now the replicate weights join the equation

names(scf_rw_2022)
head(scf_rw_2022)

# what we see: 999 replicate weights and 999 multipliers/multiplication factors (number of times the case was selected for the replicate) --> by this we end up at the initial number of total households

scf_rw_2022 %>% select(y1,yy1,wt1b1,wt1b2,mm1,mm2)

# first of all: NA's to zeros
scf_rw_2022 <- scf_rw_2022 %>% mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

# sort rw by id and drop id variable
scf_rw_2022_comb <-scf_rw_2022 %>% arrange(y1) %>% select(-y1)

# rep-weights * multiplication factor
scf_rw_2022_comb[, paste0('wgt' , 1:999)] <-
  scf_rw_2022_comb[, paste0('wt1b' , 1:999)] * scf_rw_2022_comb[, paste0('mm' , 1:999)]
scf_rw_2022_comb <- scf_rw_2022_comb[, c('yy1' , paste0('wgt' , 1:999))]

names(scf_rw_2022_comb) # consists of identifier (yy1) and 999 replicate weights (including multiplication factor)


#---------------------------------------------------------------------------------------------------------------------------------------------#

# STEP 6 - Survey Object ####

head(scf) # households are not arranged by their ID's
# arrange households by yy1
scf <- scf %>% arrange(yy1)

# split implicates
scf_list <- split(scf, factor(scf$id.imp))
# merge with large file (give it a thought if you really need all variables at this point --> your device might be substantially slower otherwise)
scf_list <- lapply(scf_list , merge , scf_2022)

# arrange households by yy1
scf_list <-
  lapply(scf_list , function(w)
    w %>% arrange(yy1))

scf_rw_2022_comb <- scf_rw_2022_comb %>% arrange(yy1)

# create a survey object

# option 1: 
# recommended by FED, setting scales of replicate weights (scale/rscales) and computation of variance (mse) manually

scf_design <-
  survey::svrepdesign(
    weights = ~ wgt ,                      # main sampling weight variable
    repweights = scf_rw_2022_comb[,-1] ,   # replicate weights (excluding the first column)
    data = imputationList(scf_list) ,      # multiple imputation object
    scale = 1 ,                            # scaling factor for the point estimates
    rscales = rep(1 / 998 , 999) ,         # scaling factors for the replicate weights (all 999 weights sum up to 131 mio. HH's in the end)
    mse = FALSE ,                          # use classic variance estimation, not MSE (FALSE --> variances computed based on mean of replicates, TRUE --> variances computed based on sum of squares around the point estimate)
    type = "other" ,                       # replication type is not standard (e.g., not BRR or JK)
    combined.weights = TRUE                # indicates replicate weights already include sampling weight
  )

# option 2: 
# to get a feeling on how to work with a survey object, let's reduce the number of replicate weights used, so R runs faster
# also, we use type=bootstrap to handle the replicate weights. why? -> if we want to compute variances and SE by hand, the FED scaling will not give us the same results (although they also use bootstrapping when sampling the replicates)

scf_design_red <- 
  survey::svrepdesign(
    weights = ~wgt,
    repweights = scf_rw_2022_comb[,2:11],  # only the first 10 replicate weights are used
    data = imputationList(scf_list),
    combined.weights = TRUE,
    type = "bootstrap",                    # no custom settings on how to scale the replicates (scale/rscales) and how to compute variances (mse) --> survey applies its own built-in rules for bootstrapping
    mse = FALSE                            # variance computed via replicates, not using MSE (like in the slides)
    )


## convey ####
# important: 
# Run the convey_prep() function on the full design (convey is a package for computing inequality and poverty indicators, before use we need to "prepare" the survey design object):
scf_design$designs <- lapply(scf_design$designs , convey_prep)
scf_design_red$designs <- lapply(scf_design_red$designs , convey_prep)


## 6.A) some examples ####

with(scf_design_red, survey::svymean(~networth)) # gives us mean and SE for all 5 implicates
MIcombine(with(scf_design_red, survey::svymean(~networth))) # overall mean after applying Rubin's rule and SE stemming from total variance (within and between implicates)


#MIcombine(with(scf_design, svymean(~networth)))
#scf_MIcombine(with(scf_design, svymean(~networth))) # se is actually higher when using scf_MIcombine


# that's the inheritance variable from before
scf_design_red$designs[[1]]$variables$x5801

# Add new columns to the data set:

scf_design_red <-
  update(scf_design_red,
         inh_received = factor(x5801, levels = c(1,5), labels = c("yes","no")),
         inh_expected = factor(x5819,  levels = c(1,5), labels = c("yes","no")),
         hhsex = factor(hhsex, levels = 1:2, labels = c("male" , "female")),
         married = as.numeric(married == 1),
         edcl = factor(edcl, levels = 1:4, labels = c("less than high school",
                                                      "high school or GED" ,
                                                      "some college" ,
                                                      "college degree")),
         counts = rep(5)) # this is a helper variable (5 because weights need to be multiplied by 5 within an implicate) 

# Count the unweighted number of records in the survey sample, overall and by groups:

MIcombine(with(scf_design_red , survey::svyby(~ counts , ~ counts , survey::unwtd.count)))
MIcombine(with(scf_design_red , survey::svyby(~ counts , ~ hhsex , survey::unwtd.count)))

# Count the weighted size of the generalizable population, overall and by groups:
MIcombine(with(scf_design_red , survey::svytotal(~ counts)))
MIcombine(with(scf_design_red , 
               survey::svyby(~ counts , ~ inh_received , survey::svytotal)))

# Calculate the mean (average) of a linear variable, overall and by groups:
MIcombine(with(scf_design_red , survey::svymean(~ networth)))
MIcombine(with(scf_design_red , 
               survey::svyby(~ networth , ~ inh_received , survey::svymean)))

# Calculate the distribution of a categorical variable, overall and by groups:
MIcombine(with(scf_design_red , survey::svymean(~ edcl)))
MIcombine(with(scf_design_red ,
               survey::svyby(~ edcl , ~ inh_expected , survey::svymean))) #distribution of edcl within inh_expected

# Calculate the sum of a linear variable, overall and by groups:
MIcombine(with(scf_design_red , survey::svytotal(~ networth)))
MIcombine(with(scf_design_red ,
               survey::svyby(~ networth , ~ hhsex , svytotal)))

# Calculate the weighted sum of a categorical variable, overall and by groups:
MIcombine(with(scf_design_red , survey::svytotal(~ edcl)))
MIcombine(with(scf_design_red ,
               survey::svyby(~ edcl , ~ hhsex , survey::svytotal)))

# Calculate the median (50th percentile) of a linear variable, overall and by groups:

MIcombine(with(
  scf_design_red ,
  survey::svyquantile(~ networth ,
              0.5 , se = TRUE , interval.type = 'quantile')
))

MIcombine(with(
  scf_design_red ,
  survey::svyquantile(~ networth ,
              c(0.25,0.5,0.75) , se = TRUE , interval.type = 'quantile')
))

MIcombine(with(
  scf_design_red ,
  survey::svyby(
    ~ networth ,
    ~ hhsex ,
    survey::svyquantile ,
    0.5 ,
    se = TRUE ,
    interval.type = 'quantile' ,
    ci = TRUE
  )
))

# Estimate a ratio:

MIcombine(with(
  scf_design_red ,
  survey::svyratio(numerator = ~ income , 
           denominator = ~ networth)
))

# Extract the coefficient, standard error, confidence interval, and coefficient of variation from any descriptive statistics function result, overall and by groups:
this_result <-
  MIcombine(with(scf_design_red ,
                     survey::svymean(~ networth)))

coef(this_result)
survey::SE(this_result)
confint(this_result)
survey::cv(this_result) # coefficient of variation

grouped_result <-
  MIcombine(with(scf_design_red ,
                     survey::svyby(~ networth , ~ hhsex , svymean)))

coef(grouped_result)
SE(grouped_result)
confint(grouped_result)
cv(grouped_result)

# Perform a survey-weighted generalized linear model:

glm_result <-
  MIcombine(with(scf_design_red ,
                     survey::svyglm(networth ~ married + edcl + inh_received)))
summary(glm_result)




# Calculate the Gini coefficient with family net worth:

MIcombine(with(scf_design_red , svygini(~ networth)))

# Calculate the Gini coefficient with income:

MIcombine(with(scf_design_red , svygini(~ income)))


#---------------------------------------------------------------------------------------------------------------------------------------------#

# STEP 7 - Intuition of replicate weights ####

## 7.A) Justification ####
# why replicate weights make sense

# create a survey object without the replicate weights
scf_design_naive <- survey::svydesign(~ 1 , data = imputationList(scf_list), weights = ~ wgt) 

# compare the SE

# 1. if we treat it as a large dataset (5x inflated) --> small SE --> we overestimate the reliability of the results
scf %>% ungroup() %>%
  summarise(theta = weighted.mean(networth,wgt),
            se = sqrt(sum((wgt/sum(wgt))^2*(networth-weighted.mean(networth,wgt))^2)))

# 2. if we correct for multiple imputation --> SE more than double
# by hand
scf %>% group_by(id.imp) %>%
  summarise(theta = weighted.mean(networth,wgt),
            se_within = sqrt(sum((wgt/sum(wgt))^2*(networth-weighted.mean(networth,wgt))^2))) %>%
  mutate(within_var = se_within^2,
         W = mean(within_var),                        # the mean within variance (W) is based on the survey weights in this case (not on the replicate weights like in the slides)
         theta_bar = mean(theta),                     # theta_bar is needed to compute the between variances
         B = (1/(5-1))*sum((theta-theta_bar)^2),      # between variance (B) formula you find in the slides
         Tvar = W+(1+1/5)*B,                          # total variance (V) formula you also find in the slides
         se = sqrt(Tvar)) #%>% summarise(theta = mean(theta), se = mean(se))

# with the survey package
with(scf_design_naive, survey::svymean(~networth))
MIcombine(with(scf_design_naive, survey::svymean(~networth)))

# 3. if we correct for multiple imputation AND sampling variance --> correcting for variability helps reducing SE --> but still higher than in 1.
MIcombine(with(scf_design_naive, survey::svymean(~networth)))
MIcombine(with(scf_design_red, survey::svymean(~networth)))
MIcombine(with(scf_design, survey::svymean(~networth)))



## 7.B) Illustration ####
# what happens when we increase number of rep. weights?

# let's make a loop where we compute variance and SE using different amounts of replicate weights

df_var <- data.frame() # create an empty dataframe

for (w in seq(1, 250 , by = 1)) { # 250 separate computations, results will be stored in the dataframe
  
  ## SETUP SURVEY DESIGN
  scf.svy_tmp <- survey::svrepdesign(weights=~wgt, repweights = scf_rw_2022_comb[,2:(w+1)], # remember: first column is yy1 rep-weights start from column 2
                             data=imputationList(scf_list),combined.weights=TRUE,type = "bootstrap")
  
  ## CALCULATE MEAN OF NETWEALTH
  result_tmp <- MIcombine(with(scf.svy_tmp, svymean(~networth)))
  
  row_tmp <- c(w , result_tmp$variance, sqrt(result_tmp$variance))
  df_var <- rbind(df_var, row_tmp)
}

#variance of total net
names(df_var) <- c( "rep_weights" , "Variance", "SE")
head(df_var)
names(df_var) <- c( "rep_weights" , "Variance", "SE")

ggplot(df_var %>% gather(key,val,Variance:SE),
       aes(x=rep_weights,y=val,color=key)) + geom_line() + 
  labs(title = "Increasing the number of replicate weights to get more stable estimates of SE and variance",
       x = "Number of replicate weights used",
       y = "Value") +
  theme_light() + theme(legend.position = "none") +
  facet_wrap(~key, scales = "free_y")

# --> variance and SE do not decrease steadily, but stabilize
# goes from noisy to more stable estimates of the variance and SE

#------------------------------------------------------------------------------#
# if you have time, the same for the total number of replicate weights

# CAUTION: THIS MIGHT TAKE AGES TO RUN

df_var <- data.frame()
for (w in seq(1, 999 , by = 1)) {
  
  
  ## SETUP SURVEY DESIGN
  scf.svy_tmp <- svrepdesign(weights=~wgt, repweights = scf_rw_2022_comb[,2:(w+1)] ,data=imputationList(scf_list),combined.weights=TRUE,type = "bootstrap")
  
  ## CALCULATE MEAN OF NETWEALTH
  result_tmp <- MIcombine(with(scf.svy_tmp, svymean(~networth)))
  
  row_tmp <- c(w , result_tmp$variance, sqrt(result_tmp$variance))
  df_var <- rbind(df_var, row_tmp)
}

#variance of total net
names(df_var) <- c( "rep_weights" , "Variance", "SE")
ggplot(df_var %>% gather(key,val,Variance:SE),
       aes(x=rep_weights,y=val,color=key)) + geom_line() + 
  labs(title = "Increasing the number of replicate weights to get more stable estimates of SE and variance",
       x = "Number of replicate weights used",
       y = "Value") +
  theme_light() + theme(legend.position = "none") +
  facet_wrap(~key, scales = "free_y")

