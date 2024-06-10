#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr) 
library(data.table)

# Avoid nested parallel computation...
Sys.setenv(OPENBLAS_NUM_THREADS=1)

#---- Input/output ----
map_rds <- snakemake@input[["map_rds"]]
gwas_rds <- snakemake@input[["gwas_rds"]]

plot_file <- snakemake@output[["plot_file"]]
plot_file2 <- snakemake@output[["plot_file2"]]


plot_path <- dirname(plot_file)
plot_file <- basename(plot_file)

plot_path2 <- dirname(plot_file2)
plot_file2 <- basename(plot_file2)

#---- Load map_with frequency ----
# map_rds <- "plinkFiles/chrall_map.rds"
map <- readRDS(map_rds)

#---- Load gwas ----
# Remove following lines for local run/debug
# df_beta_rds <- "data/merge_df_beta.rds"
gwas <- readRDS(gwas_rds)

#---- Merge with local genotypes ----
cat("Merging snps...\n")
df_beta <- snp_match(gwas, map)
df_beta <- as.data.table(df_beta)

sd_ldref <- sqrt(2* df_beta$freq * (1 - df_beta$freq))
# sd_ldref <- df_beta[, .(sqrt(2 * freq * (1 - freq)))]

# Quantitative
# TODO: review the current formula (should get an estimation from the original trait used in the GWAS)
# TODO: this is too much dependend by the name of the phenotype, should be standardize
# 
# Estimate trait standard deviation as the first percentile of the GWAS effects
# See: https://privefl.github.io/bigsnpr/articles/LDpred2.html for further details
# sd_trait_est <- quantile(df_beta[, .(sqrt(0.5 * (n_eff * beta_se**2) + beta**2))]$V1, probs=0.99)

sd_ss <- df_beta[, .(1 / sqrt(n_eff * beta_se**2 + beta**2))]
sd_ss <- as.numeric(sd_ss$V1)
sd_ss <- sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)

# Binary case
# sd_ss <- df_beta[, .(2 / sqrt(n_eff * beta_se**2 + beta**2))]

#---- Mark bad snps ----
#sd_ss < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
# sd_ss < 0.1 | sd_af < 0.05)
#
is_bad <- 
  sd_ss < (0.7 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05 #| df_beta$freq < 0.01 | df_beta$freq > 0.99


pp <- ggplot(slice_sample(data.frame(sd_ldref, sd_ss, is_bad), n = 50e3)) +
  geom_point(aes(sd_ldref, sd_ss, color = is_bad), alpha = 0.5) +
  theme_bigstatsr(0.9) +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2) +
  labs(x = "Standard deviations in the reference set",
       y = "Standard deviations derived from the summary statistics",
       color = "To remove?")
ggsave(pp, file=plot_file, path=plot_path)

pp2 <- ggplot(df_beta, aes(x=beta)) + geom_histogram() + theme_bigstatsr() + labs(x="Beta distribution of GWAS for selected SNPs")
ggsave(pp2, file=plot_file2, path=plot_path2)
