
#########
# R-setup:

rm(list=ls(all=TRUE)) # Remove all variables from workspace
library(FreqEstimationModel) # Load the haplotype frequency estimation package
library(data.table)
library(plyr)
library(tidyverse)
library(arsenal)
library(kableExtra)



#########
# TO-DO by model user:

# 1. define filename for PDF with results & diagnostics
# 2. Set working directory to source file location
# 3. Load dataset, make sure it has the following sturcture:
# - a data.frame or data.table in wide format, with data per sample stored by row;
# - SNP data is stored in the first columns, with one SNP per column;
# - SNP data is of class numeric, with values 0, 0.5, 1, or 99;
# - row names will be used as sample-ID, do not store sample ID in the first column.
# 4. set MCMC variables
# 5. set MOI variables

#########
# Define file name for model results & set wd:
filename <- "arthur_it50000_maxmoi8"

dir_vis <- "XXX"
dir_dat <- "XXX"


#########
# Load genotyping data:

mydata <- fread(file.path(dir_dat, "ARTHUR_pfmdr1_genotypes_01.csv")) %>%
  filter(!(locus1 == 99 & locus2 == 99 & locus3 == 99)) %>%
  column_to_rownames("ID")

data_summary <- list(Data = mydata)


#########
# MCMC variables:

thinning_interval <- 5 # Change this variable to increase the number of iterations per chain that are not saved in memory
no_traces_preburnin <- 10000 # For more traces but managable pdfs, don't exceed 10k and increase thinning interval instead.

no_mcmc_chains <- 3 # Number of mcmcm chains to run
parallel <- FALSE # Set to true if running code in parallel (this option isn't active yet)
NGS <- FALSE # Set to true if data are in NGS format (this option isn't active yet)
log_like_zero <- FALSE # QC check; set log(p(yi|ai)) to zero (only impacts NGS) and wipes data (see below)
mcmc_variable_list <- list(no_mcmc_chains = no_mcmc_chains,
                           no_traces_preburnin = no_traces_preburnin,
                           thinning_interval = thinning_interval,
                           NGS = NGS,
                           log_like_zero = log_like_zero)


#########
# MOI variables:

if(log_like_zero){ # When log_like_zero == TRUE, code should return the prior, hence use a prior that is easy to eye-ball
  moi_prior <- 'Geometric'
} else {
  moi_prior <- 'Geometric' # Choose between 'Uniform', 'Poisson', 'Geometric' or 'nBinomial'
}

moi_max <- 8 # Maximum MOI regarded as possible by the model (haven't tested beyond 20)
moi_hyperparameter <- 3 # Specify mean MOI (Hyperparameter for the prior on the MOI)
moi_size_hyperparameter <- 0.5 # Only applies if moi_prior == 'nBinomial' (Hyperparameter for the prior on the MOI in the prior is the negative Binomial)
moi_prior_min2 <- NULL # Specify the lower bound for the MOI per individual
moi_initial <- NULL # If null, the initial vector of multiplicties of infection (MOI) is set internally, otherwise set to inputted moi_initial
moi_list <- list(moi_hyperparameter = moi_hyperparameter,
                 moi_size_hyperparameter = moi_size_hyperparameter,
                 moi_prior = moi_prior,
                 moi_max = moi_max,
                 moi_prior_min2 = moi_prior_min2,
                 moi_initial = moi_initial)


#########
# Preprocess data
if(!NGS | log_like_zero){ # Set to true to augment partial observations rather than discard them
  augment_missing_data <- TRUE
} else {
  augment_missing_data <- FALSE
}

processed_data_list <- preprocess_data(data_summary,
                                       log_like_zero,
                                       NGS,
                                       augment_missing_data,
                                       moi_prior_min2)


# frequency variables
frequency_hyperparameter <-rep(1, processed_data_list$no_haplotypes) # The Parameter vector for the dirichlet prior on the frequency estimate
frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to inputted external_frequency_initial
frequency_list <- list(frequency_hyperparameter = frequency_hyperparameter,
                       frequency_initial = frequency_initial)

#########
# Run Markov Chain Monte Carlo (MCMC)

set.seed(1) # Set seed for reproducibility
results <- mcmc_sampling_parallel(processed_data_list,
                                  moi_list,
                                  frequency_list,
                                  mcmc_variable_list,
                                  cores_max = 3)


#########
# Save model results:

arguments <- names(as.list(args(mcmc_sampling)))
arguments <- arguments[-length(arguments)]
if(!log_like_zero){
  dataset<-'Data'
} else {
  dataset<-'QC'
} # Class character: set to dataset name in order to save pdf with data set name if REAL=TRUE and to apend if REAL=FALSE

save(list = c('results', arguments), file = file.path(dir_dat, sprintf('./%s.RData', filename)))


#####
# Load model results:
load(file.path(dir_dat, paste(filename, ".RData", sep = "")))


#####
# Visualize results:
# Generate pdf of results
source(file.path(dir_vis, './visualise_results.R'))
visualise_results(results,
                  data_summary,
                  PDF = sprintf('./%s.pdf', filename),
                  child = FALSE,
                  augment_missing_data,
                  Simulated = F)

# Convergence diagnostics
mcmc_frequency_chains <- mcmc.list(lapply(alply(results$genotype_freq_store_chains, no_mcmc_chains), mcmc))
gelman.diag(mcmc_frequency_chains, transform = TRUE, multivariate = FALSE) # Transformed (using log or logit) to improve the normality of the distribution


#########
# Extract data summary -WvL:

# Extract estimates of haplotype frequency (post burnin)

arr_genfreq <- results$genotype_freq_store_chains
arr_genfreq %>% dim

lst_genfreq <- lapply(alply(arr_genfreq, 3), as.data.table) # transform the 3 dimensional array to a list, with each list element the results of one MCMC chain
lst_genfreq[[1]] %>% dim

dt_genfreq <- rbindlist(lst_genfreq)
dt_genfreq %>% dim

# Have a quick look at the distributions

dt_genfreq %>%
  setNames(paste("haplo_", names(.), sep = "")) %>%
  pivot_longer(cols = starts_with("haplo_"), names_to = "haplotype", values_to = "freq") %>%
  ggplot(.) +
  geom_density(aes(freq, fill = haplotype), alpha = 0.5) +
  theme_bw()

dt_genfreq %>%
  setNames(paste("haplo_", names(.), sep = "")) %>%
  pivot_longer(cols = starts_with("haplo_"), names_to = "haplotype", values_to = "freq") %>%
  ggplot(.) +
  geom_density(aes(freq, fill = haplotype), alpha = 0.5) +
  facet_wrap(~haplotype, scales = "free", ncol = 3) +
  theme_bw()

lst_genfreq[[1]] %>%
  setNames(paste("haplo_", names(.), sep = "")) %>%
  mutate(itnr = 1:nrow(.)) %>%
  pivot_longer(cols = starts_with("haplo_"), names_to = "haplotype", values_to = "freq") %>%
  ggplot(.) +
  geom_line(aes(x = itnr, y = freq, color = haplotype), alpha = 0.75) +
  theme_bw()


# summarize the frequencies by haplotype:

tableby(~., data = dt_genfreq,
        numeric.stats = c("mean", "median", "range", "q1q3")) %>%
  summary(text = T)

dt_genfreq %>%
  setNames(paste("haplo_", names(.), sep = "")) %>%
  pivot_longer(cols = starts_with("haplo_"), names_to = "haplotype", values_to = "freq") %>%
  group_by(haplotype) %>%
  dplyr::summarize(mean = mean(freq, na.rm = T),
                   median = median(freq, na.rm = T),
                   min = min(freq, na.rm = T),
                   max = max(freq, na.rm = T),
                   Q25 = quantile(freq, 0.25, na.rm = T),
                   Q75 = quantile(freq, 0.75, na.rm = T),
                   Q2.5 = quantile(freq, 0.025, na.rm = T),
                   Q97.5 = quantile(freq, 0.975, na.rm = T))


# Final summary table:

dt_genfreq %>%
  setNames(paste("haplo_", names(.), sep = "")) %>%
  pivot_longer(cols = starts_with("haplo_"), names_to = "haplotype", values_to = "freq") %>%
  mutate(haplotype = factor(haplotype, ordered = T,
                            levels = c("haplo_000", "haplo_010", "haplo_100", "haplo_110",
                                       "haplo_101", "haplo_011", "haplo_001", "haplo_111"))) %>%
  group_by(haplotype) %>%
  dplyr::summarize(median = median(freq, na.rm = T),
                   Q2.5 = quantile(freq, 0.025, na.rm = T),
                   Q97.5 = quantile(freq, 0.975, na.rm = T)) %>%
  mutate_if(is.numeric, function(x){ formatC(x*100, digits = 1, format = "f")}) %>%
  mutate(comb = paste(median, " (", Q2.5, ", ", Q97.5, ")", sep = "")) %>%
  dplyr::select(haplotype, comb) %>%
  kable() %>% kable_styling("condensed")


