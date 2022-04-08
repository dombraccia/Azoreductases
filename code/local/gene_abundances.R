# ============================= gene_abundances.R =========================== #
#' description: script for exploring the abundance of azoreduction genes 
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
#' 
#' Author: Domenick J. Braccia  
#' Last Update: 22 Nov 2021
# =========================================================================== #

## load data
hmp2_gene_counts_raw <- read.csv("../../results/salmon_processed/hmp2/Azoreductase_counts.tsv",
                                 sep = "\t", header = TRUE)
prism_gene_counts_raw <- read.csv("../../results/salmon_processed/prism/Azoreductase_counts.tsv",
                                  sep = "\t", header = TRUE)

## process data
hmp2_rowData <- as.data.frame(do.call(rbind, strsplit(hmp2_gene_counts_raw$Name, ';')))
hmp2_rowData <- cbind(hmp2_rowData, hmp2_gene_counts_raw[, 2])
colnames(hmp2_rowData) <- c("genomeID", "contigID", "sequenceID", "gene_name", "length")
hmp2_gene_counts <- cbind(hmp2_rowData, hmp2_gene_counts_raw[, 4:ncol(hmp2_gene_counts_raw)])

prism_rowData <- as.data.frame(do.call(rbind, strsplit(prism_gene_counts_raw$Name, ';')))
prism_rowData <- cbind(prism_rowData, prism_gene_counts_raw[, 2])
colnames(prism_rowData) <- c("genomeID", "contigID", "sequenceID", "gene_name", "length")
prism_gene_counts <- cbind(prism_rowData, prism_gene_counts_raw[, 4:ncol(prism_gene_counts_raw)])

## subsetting hmp2 metadata by available metagenomic data
hmp2_sample_names <- colnames(hmp2_gene_counts)[-c(1:5)]
hmp2_sample_names <- substr(hmp2_sample_names, 1, nchar(hmp2_sample_names) - 1)
hmp2_sample_names <- as.data.frame(do.call(rbind, strsplit(hmp2_sample_names, '_')))[, 2]
hmp2_pData_subset <- k2_hmp2_pData[k2_hmp2_pData$site_sub_coll %in% hmp2_sample_names, ]

## renaming and reordering columns of hmp2_gene_counts
colnames(hmp2_gene_counts) <- c(colnames(hmp2_rowData), hmp2_sample_names)
hmp2_gene_counts <- hmp2_gene_counts[, c(colnames(hmp2_gene_counts)[1:5], hmp2_pData_subset$site_sub_coll)]

## removing samples with <100k filtered reads
filtered_samples <- filter(hmp2_pData_subset, reads_filtered > 1e5)$site_sub_coll
hmp2_pData_subset <- hmp2_pData_subset[hmp2_pData_subset$site_sub_col %in% filtered_samples, ]
hmp2_gene_counts <- select(hmp2_gene_counts, c(colnames(prism_rowData), filtered_samples))

## normalizing counts by gene length and number of filtered reads 
hmp2_rpk <- sweep(x = hmp2_gene_counts[,-c(1:5)], MARGIN = 1, STATS = (hmp2_gene_counts$length / 1e3), FUN = "/")
# hmp2_rpk_norm_factors <- colSums(hmp2_rpk) ## TPM style scaling
hmp2_nreads <- hmp2_pData_subset$reads_filtered 

prism_rpk <- sweep(x = prism_gene_counts[,-c(1:5)], MARGIN = 1, STATS = (prism_gene_counts$length / 1e3), FUN = "/")
# prism_rpk_norm_factors <- colSums(prism_rpk) ## TPM style scaling
prism_nreads <- k2_prism_pData$spots

### testing out per thousand vs. per million scalling factors
hmp2_counts_normed <- sweep(x = hmp2_rpk, MARGIN = 2, STATS = hmp2_nreads, FUN = "/") * 1e6
prism_counts_normed <- sweep(x = prism_rpk, MARGIN = 2, STATS = prism_nreads, FUN = "/") * 1e6

## aggregating normalized counts
hmp2_counts_normed <- cbind(hmp2_rowData, hmp2_counts_normed)
hmp2_counts_normed %>% 
  select(-c(1:3, 5)) %>% 
  group_by(gene_name) %>% 
  summarise_each(list(sum)) -> hmp2_agg_counts

prism_counts_normed <- cbind(prism_rowData, prism_counts_normed)
prism_counts_normed %>% 
  select(-c(1:3, 5)) %>% 
  group_by(gene_name) %>% 
  summarise_each(list(sum)) -> prism_agg_counts

## prepping data frames for plotting
### hmp2
hmp2_agg_df <- t(hmp2_agg_counts[,-1])
colnames(hmp2_agg_df) <- hmp2_agg_counts$gene_name
hmp2_agg_df <- as.data.frame(hmp2_agg_df)
# hmp2_agg_df <- hmp2_agg_df[match(hmp2_pData_subset$site_sub_coll, rownames(hmp2_agg_df)),]

#### getting desired metadata for plotting
hmp2_agg_pData_df <- cbind(select(hmp2_pData_subset, site_sub_coll, Participant.ID, diagnosis),
                           hmp2_agg_df) 

### prism
prism_agg_df <- t(prism_agg_counts[,-1])
colnames(prism_agg_df) <- prism_agg_counts$gene_name
prism_agg_df <- as.data.frame(prism_agg_df)
prism_agg_df <- prism_agg_df[match(k2_prism_pData$SampleName, rownames(prism_agg_df)),]

#### getting desired metadata for plotting

prism_agg_pData_df <- cbind(select(k2_prism_pData, Diagnosis, SampleName),
                           prism_agg_df)

## summarizing information per individual (only needed for hmp2)
hmp2_agg_pData_df %>% 
  select(-site_sub_coll) %>% 
  group_by(Participant.ID, diagnosis) %>% 
  summarise_each(list(mean)) -> hmp2_agg_pp_mean ## pp = per participant

std_err <- function(x) sd(x) / sqrt(length(x))
hmp2_agg_pData_df %>% 
  select(-site_sub_coll) %>% 
  group_by(Participant.ID, diagnosis) %>% 
  summarise_each(list(std_err)) -> hmp2_agg_pp_se ## pp = per participant

## prepping data frames for plotting

### hmp2
hmp2_agg_pp_mean_long <- melt(hmp2_agg_pp_mean, id.vars = c("Participant.ID", "diagnosis"), value.name = "counts", variable.name = "gene_name")
hmp2_agg_pp_mean_long %>% group_by(Participant.ID) %>% summarise(counts_sum = sum(counts)) %>% arrange(desc(counts_sum)) -> hmp2_counts_pp_sum
hmp2_relevel <- as.character(hmp2_counts_pp_sum$Participant.ID)
hmp2_agg_pp_mean_long$Participant.ID <- as.character(hmp2_agg_pp_mean_long$Participant.ID)
hmp2_agg_pp_mean_long$Participant.ID <- factor(hmp2_agg_pp_mean_long$Participant.ID, levels = c(hmp2_relevel))

hmp2_agg_pData_long <- melt(hmp2_agg_pData_df, id.vars = c("site_sub_coll", "Participant.ID", "diagnosis"), value.name = "counts", variable.name = "gene_name")
hmp2_agg_pData_long %>% group_by(site_sub_coll) %>% summarise(counts_sum = sum(counts)) %>% arrange(desc(counts_sum)) -> hmp2_counts_sum
hmp2_relevel <- as.character(hmp2_counts_sum$site_sub_coll)
hmp2_agg_pData_long$site_sub_coll <- as.character(hmp2_agg_pData_long$site_sub_coll)
hmp2_agg_pData_long$site_sub_coll <- factor(hmp2_agg_pData_long$site_sub_coll, levels = c(hmp2_relevel))


hmp2_agg_pp_se_long <- melt(hmp2_agg_pp_se, id.vars = c("Participant.ID", "diagnosis"), value.name = "stderr", variable.name = "gene_name")
hmp2_agg_pp_se_long$diagnosis <- factor(hmp2_agg_pp_se_long$diagnosis, levels = c("nonIBD", "UC", "CD"))

### prism
prism_agg_long <- melt(prism_agg_pData_df, id.vars = c("Diagnosis", "SampleName"), value.name = "counts", variable.name = "gene_name")
prism_agg_long %>% group_by(SampleName) %>% summarise(counts_sum = sum(counts)) %>% arrange(desc(counts_sum)) -> prism_counts_sum
prism_relevel <- as.character(prism_counts_sum$SampleName)
# prism_agg_long$SampleName <- as.character(hmp2_agg_pp_mean_long$Participant.ID)
prism_agg_long$SampleName <- factor(prism_agg_long$SampleName, levels = c(prism_relevel))

##### ======================= plotting for hmp2 ========================= #####
### per sample
low_count_genes <- names(head(sort(colSums(hmp2_agg_pp_mean[,-c(1:2)])), 5))
hmp2_agg_pData_long %>% 
  filter(!(gene_name %in% low_count_genes)) %>%
  ggplot(aes(x = site_sub_coll, y = counts, fill = gene_name)) +
  geom_bar(position = "stack", stat = "identity") +
  # geom_hline(yintercept = 1e3) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Samples (n = 703)") + ylab("normalized gene abundance") + 
  ggtitle("Normalized abundances of azoreductase genes in HMP2") -> hmp2_gene_abundnaces
hmp2_gene_abundnaces


### per individual
# hmp2_agg_pp_mean_long %>% 
  # mutate(pseudo_counts = counts + 1) %>% 
hmp2_agg_pp_mean_long %>% 
  filter(!(gene_name %in% low_count_genes)) %>%
ggplot(aes(x = Participant.ID, y = counts, fill = gene_name)) +
  geom_bar(position = "stack", stat = "identity") +
  # geom_hline(yintercept = 1e3) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Individuals (n = 104)") + ylab("normalized gene abundance") + 
  ggtitle("Normalized abundances of azoreductase genes in HMP2") -> hmp2_gene_abundnaces_pp
hmp2_gene_abundnaces_pp

top_counts_part <- head(hmp2_relevel)
hmp2_agg_pp_mean_long %>%
  filter(!(Participant.ID %in% top_counts_part)) %>% 
  filter(!(gene_name %in% low_count_genes)) %>% 
  ggplot(aes(x = Participant.ID, y = counts, fill = gene_name)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() + theme(legend.position = "none",
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       text = element_text(size = 16)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("") + ylab("")  -> hmp2_gene_abundnaces_inset
hmp2_gene_abundnaces_inset

hmp2_gene_abundnaces + inset_element(
  hmp2_gene_abundnaces_inset, 
  left = 0.06, 
  bottom = 0.1, 
  right = unit(1, 'npc') - unit(1, 'cm'), 
  top = unit(1, 'npc') - unit(1, 'cm')
)

#### histograms
hmp2_agg_pp_mean_long %>% 
  filter(gene_name == "mdaB", counts < 1000) %>%
ggplot(aes(x = counts)) +
  geom_histogram(binwidth = 10)



### plotting standard errors 
clrs2 <- c("#999999", "#F0E442", "#D55E00")
head(hmp2_agg_pp_se_long)
hmp2_agg_pp_se_long %>% 
  mutate(pseudo_stderr = stderr + 1) %>%
  ggplot() + 
  geom_violin(aes(x = diagnosis, y = pseudo_stderr, fill = diagnosis)) +
  scale_y_continuous(trans = "log2",
                     breaks = c(1,10,100,1000,5000),
                     limits = c(1,5e3)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = clrs2) +
  ylab("log2(SE + 1)")

hmp2_agg_pp_se_long %>% 
  mutate(pseudo_stderr = stderr + 1) %>%
  na.omit() %>% 
  ggplot() + 
  geom_violin(aes(x = diagnosis, y = pseudo_stderr, fill = diagnosis)) +
  scale_y_continuous(trans = "log2",
                     breaks = c(1,10,100,1000,5000),
                     limits = c(1,5e3)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = clrs2)

hmp2_agg_pp_se_long %>% 
  mutate(pseudo_stderr = stderr + 1) %>%
  na.omit() %>% 
  filter(gene_name == "cladeIVa") %>%
  ggplot() + 
  geom_violin(aes(x = diagnosis, y = pseudo_stderr, fill = diagnosis)) +
  scale_y_continuous(trans = "log2",
                     breaks = c(1,10,100,1000,5000),
                     limits = c(1, 5e3)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = clrs2)


##### ======================= plotting for prism ======================== #####
low_count_genes <- names(head(sort(colSums(prism_agg_pData_df[,-c(1:3)])), 5))
prism_agg_long %>% 
  filter(!(gene_name %in% low_count_genes)) %>%
  ggplot(aes(x = SampleName, y = counts, fill = gene_name)) +
  geom_bar(position = "stack", stat = "identity") +
  # geom_hline(yintercept = 1e3) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Individuals (n = 220)") + ylab("normalized gene abundance") + 
  ggtitle("Normalized abundances of azoreductase genes in PRISM") -> prism_gene_abundnaces_pp
prism_gene_abundnaces_pp


##### ====================== calculating AZscore ======================== #####
library(vegan)

## calculating shannon diversity for hmp2 and prism samples
hmp2_sample_diversity <- sapply(hmp2_agg_counts[,-1], diversity)
prism_sample_diversity <- sapply(prism_agg_counts[,-1], diversity) 

hist(hmp2_sample_diversity, breaks = 30)
hist(prism_sample_diversity, breaks = 30)

## calculating pooled azoreductase gene abundances
hmp2_pooled_azored_abundance <- log2(colSums(hmp2_agg_counts[,-1]))
prism_pooled_azored_abundance <- log2(colSums(prism_agg_counts[,-1]))

hist(hmp2_pooled_azored_abundance, breaks = 30)
hist(prism_pooled_azored_abundance, breaks = 30)

## calculating azrScore
hmp2_azr_scores <- hmp2_sample_diversity * hmp2_pooled_azored_abundance
hmp2_azr_scores <- hmp2_azr_scores[!is.na(hmp2_azr_scores)]
prism_azr_scores <- prism_sample_diversity * prism_pooled_azored_abundance
prism_azr_scores <- prism_azr_scores[!is.na(prism_azr_scores)]

hist(hmp2_azr_scores, xlim = range(0:20), breaks = 30)
hist(prism_azr_scores, breaks = 30)

hmp2_azr_scores %>% 
  as.data.frame() %>% 
ggplot(aes(x = hmp2_azr_scores)) + 
  geom_histogram(color = "black", fill = "grey", binwidth = 0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     # axis.text.x = element_blank(),
                     # axis.ticks.x = element_blank(),
                     text = element_text(size = 16)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  xlab("") + ylab("") +
  ggtitle("Distribution of azrScore for HMP2, num. samples = 703")

## incorporating azrScores into dataset info
hmp2_agg_pData_df <- hmp2_agg_pData_df[hmp2_agg_pData_df$site_sub_coll %in% names(hmp2_azr_scores), ]
hmp2_azrScore_pData <- cbind(hmp2_agg_pData_df[,c(1:3)], hmp2_azr_scores)
hmp2_azrScore_pData$diagnosis <- factor(hmp2_azrScore_pData$diagnosis, levels = c("nonIBD", "UC", "CD"))
prism_agg_pData_df <- prism_agg_pData_df[rownames(prism_agg_pData_df) %in% names(prism_azr_scores), ]
prism_azrScore_pData <- cbind(prism_agg_pData_df[,c(1:2)], prism_azr_scores)
prism_azrScore_pData$Diagnosis <- factor(prism_azrScore_pData$Diagnosis, levels = c("HC", "UC", "CD"))

## comparing azrScores across disease types
hmp2_azrScore_pData_long <- melt(hmp2_azrScore_pData, 
                                 id.vars = c("site_sub_coll", "Participant.ID", "diagnosis"), 
                                 value.name = "azrScore")
prism_azrScore_pData_long <- melt(prism_azrScore_pData, 
                                  id.vars = c("SampleName", "Diagnosis"), 
                                  value.name = "azrScore")

clrs2 <- c("#999999", "#F0E442", "#D55E00")
hmp2_azrScore_pData_long %>% 
  ggplot(aes(x = diagnosis, y = azrScore, fill = diagnosis)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = clrs2)

prism_azrScore_pData_long %>% 
  ggplot(aes(x = Diagnosis, y = azrScore, fill = Diagnosis)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = clrs2)

##### =========================== basement ============================== #####

## [may not be necessary?] combining data to summarized experiment for ease of use
hmp2_se <- SummarizedExperiment(assays = list(counts = hmp2_agg_counts[, -1]),
                                rowData = hmp2_agg_counts[, 1],
                                colData = hmp2_pData_subset)

prism_se <- SummarizedExperiment(assays = list(counts = prism_agg_counts[, -1]),
                                 rowData = prism_agg_counts[, 1],
                                 colData = k2_prism_pData)




