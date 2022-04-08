# ==================== healthy_relative_abundances_hmm.R ==================== #
#' description: script for comparing the relative abundances of putative 
#' cysteine degrading and sulfite reducing bacteria in stool samples from cMD.
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
# =========================================================================== #

print("- importing taxa hit information and feature information")
azored_taxa_hits_raw <- read.csv(
  file = "../../results/from-GutFunFind/Azoreductase.taxa_hits_hmm.tsv", 
  sep = "\t", header = TRUE)
azored_feature <- read.csv(
  file = "../../results/from-GutFunFind/Azoreductase.feature_hmm.tsv", 
  sep = "\t", header = TRUE)

uhgg_genid2taxa <- read.csv("../../data/from-xiaofang/spec.txt", sep = '\t', header = FALSE)
uhgg_genid2taxa$V2 <- gsub(";", "|", uhgg_genid2taxa$V2)
colnames(uhgg_genid2taxa) <- c("genome_ID", "taxa_ID")

azored_taxa_hits_present <- filter(azored_taxa_hits_raw, azored_pa == "present")
azored_taxa_hits <- as.vector(t(azored_taxa_hits_present$taxa_name), mode = "character") 
azored_taxa_hits <- gsub(";", "|", azored_taxa_hits)
azored_taxa_hits <- unique(azored_taxa_hits) ## removing some duplicates, (ex. s__Collinsella_aerofaciens_F)

print("- subsetting kraken2 RA controls data for azored bacteria")
azo_k2_RA_controls <- k2_RA_controls[azored_taxa_hits, ]
azo_k2_cs_controls <- colSums(azo_k2_RA_controls)

## prepping data for plotting
df <- data.frame(azored = azo_k2_cs_controls)
dim(df)

print("- plotting abundances of sulfite red bac and cys deg bac")
azored_healthy_RA <- ggplot(melt(df), aes(x = factor(variable), y = value, fill = factor(variable))) +
  geom_violin(scale = "width") +
  # geom_boxplot(outlier.alpha = 0.12) + #geom_jitter(aes(alpha = 0.0001), width = 0.25) +
  scale_shape_manual("", values=c("mean" = "X")) +
  ylab("relative microbial abundances (%)") + labs(fill = "Function") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 15),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 14),
        # plot.title = element_blank(),
        legend.position = "none") +
  ggtitle("Relative Abundances of azo-reducing bacteria")
azored_healthy_RA
# ggsave("", plot = azored_healthy_RA, width = 6, height = 7)
# ggsave("", plot = azored_healthy_RA, width = 6, height = 7)
