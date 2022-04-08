# ========================== azored_heterogeneity.R ========================= #
#' description: script for visualizing and quantifying the heterogeneity of 
#' azoreducing genes and azoreducing species 
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
#' 
#' Author: Domenick J. Braccia  
#' Last Update: 22 Nov 2021
# =========================================================================== #

## data overview

# hmp2_agg_pData_df
# hmp2_agg_pp_mean_long
# hmp2_agg_pData_long
# 
# prism_agg_pData_df
# prism_agg_long
# 
# azored_k2_cMD_df
# azored_k2_hmp2_df
# azored_k2_prism_df
# azored_k2_cib_df


#### gene abundance ####

# HMP2
## per sample 
low_count_genes <- names(head(sort(colSums(hmp2_agg_pp_mean[,-c(1:2)])), 5))
hmp2_agg_pData_long %>% 
  filter(!(gene_name %in% low_count_genes)) %>%
  ggplot(aes(x = site_sub_coll, y = counts, fill = gene_name)) +
  geom_bar(position = "stack", stat = "identity") +
  # geom_hline(yintercept = 1e3) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(legend.position = "bottom",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 12)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Samples (n = 703)") + ylab("normalized gene abundance") + 
  ggtitle("Abundances of azoreductase genes in HMP2 (per sample)") -> hmp2_gene_abundnaces
hmp2_gene_abundnaces
ggsave("figures/figure4_v1/hmp2_gene_abund_per_sample.png", plot = hmp2_gene_abundnaces, width = 8, height = 5)
ggsave("figures/figure4_v1/hmp2_gene_abund_per_sample.svg", plot = hmp2_gene_abundnaces, width = 8, height = 5)

## per individual (samples averaged) 
low_count_genes <- names(head(sort(colSums(hmp2_agg_pp_mean[,-c(1:2)])), 5))
hmp2_agg_pp_mean_long %>% 
  filter(!(gene_name %in% low_count_genes)) %>%
  ggplot(aes(x = Participant.ID, y = counts, fill = gene_name)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  # geom_hline(yintercept = 1e3) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 12)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Subjects (n = 104)") + ylab("normalized gene abundance") -> p1a
p1a
ggsave("figures/figure4_v1/hmp2_gene_abund_per_individual.png", plot = p1a, width = 8, height = 5)
ggsave("figures/figure4_v1/hmp2_gene_abund_per_individual.svg", plot = p1a, width = 8, height = 5)

## relative abundance
azored_k2_hmp2_df %>% 
  ggplot(aes(x = RA, color = diagnosis, fill = diagnosis)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  # geom_density(alpha = 0.2) +
  theme_bw() + theme(legend.position = "bottom",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 12)) +
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("relative abundance") -> hmp2_rel_abund_hist
hmp2_rel_abund_hist
ggsave("figures/figure4_v1/hmp2_rel_abund_hist.png", plot = hmp2_rel_abund_hist, width = 8, height = 5)
ggsave("figures/figure4_v1/hmp2_rel_abund_hist.svg", plot = hmp2_rel_abund_hist, width = 8, height = 5)

azored_k2_hmp2_df %>% 
  ggplot(aes(x = RA, fill = diagnosis)) +
  # geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
  geom_density(alpha = 0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 12)) +
  # scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.06)) +
  scale_x_continuous(limits = c(0,100)) +
  xlab("relative abundance") -> p1b
p1b
ggsave("figures/figure4_v1/hmp2_rel_abund_density.png", plot = p1b, width = 8, height = 5)
ggsave("figures/figure4_v1/hmp2_rel_abund_density.svg", plot = p1b, width = 8, height = 5)

##### PRISM #####
low_count_genes <- names(head(sort(colSums(prism_agg_pData_df[,-c(1:3)])), 5))
prism_agg_long %>% 
  filter(!(gene_name %in% low_count_genes)) %>%
  ggplot(aes(x = SampleName, y = counts, fill = gene_name)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  # geom_hline(yintercept = 1e3) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(legend.position = "none",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     text = element_text(size = 12)) +
  scale_fill_brewer(palette = "Dark2") +
  xlab("Subjects (n = 220)") + ylab("normalized gene abundance") -> p2a
p2a

top_counts_prism <- head(prism_relevel,12)
prism_agg_long %>% 
  filter(!(SampleName %in% top_counts_prism)) %>% 
  filter(!(gene_name %in% low_count_genes)) %>% 
  ggplot(aes(x = SampleName, y = counts, fill = gene_name)) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(legend.position = "none",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title = element_blank(),
                     text = element_text(size = 12)) +
  scale_fill_brewer(palette = "Dark2") -> prism_gene_abundnaces_inset
prism_gene_abundnaces_inset

p2a <- p2a + 
  inset_element(prism_gene_abundnaces_inset, 
  left = 0.1, 
  bottom = 0.12, 
  right =  0.9,
  top = 0.9)
  # right = unit(1, 'npc') - unit(1, 'cm'), 
  # top = unit(1, 'npc') - unit(1, 'cm'))
  
ggsave("figures/figure4_v1/prism_gene_abund_per_individual.png", plot = p2a, width = 8, height = 5)
ggsave("figures/figure4_v1/prism_gene_abund_per_individual.svg", plot = p2a, width = 8, height = 5)

## relative abundance
azored_k2_prism_df %>% 
  ggplot(aes(x = RA, color = diagnosis, fill = diagnosis)) +
  geom_histogram(alpha = 0.5, bins = 25) +
  # geom_density(alpha = 0.2) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 12)) +
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100)) +
  xlab("relative abundance") 

azored_k2_prism_df %>% 
  ggplot(aes(x = RA, fill = diagnosis)) +
  # geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
  geom_density(alpha = 0.5) +
  theme_bw() + theme(legend.position = "bottom",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 12)) +
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.08)) +
  scale_x_continuous(limits = c(0,100)) +
  xlab("relative abundance") -> p2b
p2b
ggsave("figures/figure4_v1/prism_rel_abund_density.png", plot = p2b, width = 8, height = 5)
ggsave("figures/figure4_v1/prism_rel_abund_density.svg", plot = p2b, width = 8, height = 5)

### getting legends
# extract the legend from one of the plots
legend_a <- get_legend(
  # create some space to the left of the legend
  p1a + theme(legend.box.margin = margin(0, 0, 0, 12),
              legend.position = "bottom")
)
legend_b <- get_legend(
  # create some space to the left of the legend
  p1b + theme(legend.box.margin = margin(0, 0, 0, 12),
              legend.position = "bottom")
)

##### cMD #####
azored_k2_cMD_df %>% 
  ggplot(aes(x = RA, fill = diagnosis)) +
  geom_histogram() +
  # geom_density() +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 250)) +
  xlab("relative abundance") + ggtitle("Relative abundances of azoreducing species across cMD")

azored_k2_cMD_df %>% 
  ggplot(aes(x = RA, fill = diagnosis)) +
  # geom_histogram() +
  geom_density(alpha = 0.7) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.065)) +
  xlab("relative abundance") + ggtitle("curatedMetagenomicData") -> cmd_rel_abund_density_pi
cmd_rel_abund_density_pi
ggsave("figures/figure4_v1/cMD_rel_abund_density.png", plot = cmd_rel_abund_density_pi, width = 8, height = 5)
ggsave("figures/figure4_v1/cMD_rel_abund_density.svg", plot = cmd_rel_abund_density_pi, width = 8, height = 5)


##### ggridges plots #####
clrs2_alt <- rev(clrs2)
azored_k2_prism_df$diagnosis <- factor(azored_k2_prism_df$diagnosis, levels = c("CD", "UC", "HC"))
azored_k2_prism_df %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  # geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
  geom_density_ridges(scale = 2) +
  theme_bw() + theme(legend.position = "bottom",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 12)) +
  scale_color_manual(values = clrs2_alt) +
  scale_fill_manual(values = clrs2_alt) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(limits = c(0,100)) +
  xlab("relative abundance") + ylab("") + ggtitle("Relative Abundance of putative azoreducing species - PRISM") -> ridge_plt
ridge_plt

## cMD
azored_k2_cMD_df %>% 
  ggplot(aes(x = RA, fill = diagnosis, y = diagnosis)) +
  geom_density_ridges(scale = 2) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 16)) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("relative abundance") + ggtitle("Relative Abundance of putative azoreducing species (cMD)")

