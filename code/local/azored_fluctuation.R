# ========================== azored_fluctuation.R ========================= #
#' description: script for analyzing/visualizing the flutctuation of 
#' azoreducing genes and species over time (hmp2 data only)
#' 
#' Abbv.: 
#'     1. (ms = multi-sample)
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
#' 
#' Author: Domenick J. Braccia  
#' Last Update: 6 Jan 2022
# =========================================================================== #

## data overview: must have most (or all) of these variables to run the 
## following script
# 
# hmp2_agg_pData_df
# azored_k2_hmp2_df <- as.data.frame(cbind(azored_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(azored_k2_hmp2_cs_risk)]))
# colnames(azored_k2_hmp2_df) <- c("RA", "diagnosis")
# azored_k2_hmp2_df$RA <- as.numeric(azored_k2_hmp2_df$RA)
# azored_k2_hmp2_df$func <- rep("azored", length(azored_k2_hmp2_cs_risk))
# azored_k2_hmp2_df$diagnosis <- factor(azored_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# 
# azored_k2_hmp2_RA
# k2_hmp2_RA <- as.data.frame(k2_hmp2_RA)
# azored_k2_hmp2_RA <- k2_hmp2_RA[azored_taxa_hits, ]
# azored_k2_hmp2_cs <- colSums(azored_k2_hmp2_RA)
# 
# azored_k2_hmp2_RA ## dim: 2076 x 1338
# azored_k2_hmp2_df ## dim: 1317 x    3
# k2_hmp2_pData
# 

##### ============ fluctuation of putative azoreducing spp. ============ #####
azored_k2_hmp2_RA ## relative abundance data (known + putative)
azored_k2_hmp2_df ## diagnosis data (UC, CD, nonIBD) & summed RA for putative
                  ##    azoreducing spp
k2_hmp2_pData     ## all other relevant pData for hmp2 samples


## matching metadata for relative abundance estimates per sample (_m :: matched)
azored_k2_hmp2_RA <- azored_k2_hmp2_RA[, colnames(azored_k2_hmp2_RA) %in% rownames(azored_k2_hmp2_df)]

## reorder plotting dataframe based on External.ID column of hmp2_pData
azored_k2_hmp2_df <- azored_k2_hmp2_df[match(k2_hmp2_pData$External.ID, rownames(azored_k2_hmp2_df)), ]

## add more metadata to the _df
azored_k2_hmp2_df$site_sub_col <- k2_hmp2_pData$site_sub_coll

## extract participant ID and collection number from column `site_sub_col
azored_k2_hmp2_df$Participant.ID <- sapply(strsplit(azored_k2_hmp2_df$site_sub_col, "C"), tail, 2)[1, ]
azored_k2_hmp2_df$coll_num <- sapply(strsplit(azored_k2_hmp2_df$site_sub_col, "C"), tail, 1)

## filter out participants based on number of stool collections obtained (>5)
azored_k2_hmp2_df %>% 
  group_by(Participant.ID) %>% 
  summarise(n_coll = n()) %>% 
  filter(n_coll >= 20) -> hmp2_n_coll

azored_k2_hmp2_df_subset <- azored_k2_hmp2_df[azored_k2_hmp2_df$Participant.ID %in% hmp2_n_coll$Participant.ID, ]
azored_k2_hmp2_df_subset$coll_num <- as.integer(azored_k2_hmp2_df_subset$coll_num)

## pre-pending Subject IDs with 'C' if necessary
azored_k2_hmp2_df_subset$Participant.ID[azored_k2_hmp2_df_subset$Participant.ID == "3015"] <- "C3015"
azored_k2_hmp2_df_subset$Participant.ID[azored_k2_hmp2_df_subset$Participant.ID == "3022"] <- "C3022"
azored_k2_hmp2_df_subset$Participant.ID[azored_k2_hmp2_df_subset$Participant.ID == "3027"] <- "C3027"

#### ============================= plotting v2 ============================== ####

clrs_hmp2 <- c("#999999", "#D6B409", "#D55E00")

azored_k2_hmp2_df_subset %>% 
  group_by(Participant.ID) %>% 
  summarise(mean_RA = mean(RA)) %>% 
  arrange(mean_RA) -> pID_lvls

azored_k2_hmp2_df_subset$Participant.ID <- factor(azored_k2_hmp2_df_subset$Participant.ID, levels = pID_lvls$Participant.ID)

azored_k2_hmp2_df_subset %>% 
  ggplot(aes(x = Participant.ID, y = RA, group = Participant.ID, color = diagnosis)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.5, width = 0.1) +
  # facet_wrap(~diagnosis, ncol = 3) +
  theme_bw() + theme(text = element_text(size = 10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0.03,0.03)) +
  scale_color_manual(values = clrs_hmp2) +
  ylab("relative abundance, %") + xlab("Subjects") -> top
top

azored_k2_hmp2_df_subset %>% 
  ggplot(aes(x = coll_num, y = RA, group = Participant.ID, color = diagnosis)) + 
  geom_point(size = 1) + geom_line() +
  facet_wrap(~diagnosis, ncol = 3) +
  theme_bw() + theme(text = element_text(size = 10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_blank(),
                     strip.text.x = element_blank()) +
  scale_color_manual(values = clrs_hmp2) + 
  scale_x_continuous(breaks = c(1,5,10,15,20,24),
                     limits = c(1,24),
                     expand = c(0.04,0.04)) +
  ylab("relative abundance, %") + xlab("collection number") -> bottom
bottom

azored_k2_hmp2_df_subset %>% 
  filter(diagnosis == "nonIBD") %>% 
  ggplot(aes(x = coll_num, y = RA, group = Participant.ID)) +
  geom_point(aes(shape = Participant.ID), size = 2, color = clrs_hmp2[1]) + geom_line(color = clrs_hmp2[1]) +
  theme_bw() + theme(text = element_text(size = 10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.x = element_blank()) + 
  scale_shape_manual(values = 1:nlevels(azored_k2_hmp2_df_subset$Participant.ID)) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(breaks = c(1,5,10,15,20,24),
                     limits = c(1,24),
                     expand = c(0.04,0.04)) +
  ylab("relative abundance, %") + labs(shape = "nonIBD Subjects") -> bottom_left

azored_k2_hmp2_df_subset %>% 
  filter(diagnosis == "UC") %>% 
  ggplot(aes(x = coll_num, y = RA, group = Participant.ID)) +
  geom_point(aes(shape = Participant.ID), size = 2, color = clrs_hmp2[2]) + geom_line(color = clrs_hmp2[2]) +
  theme_bw() + theme(text = element_text(size = 10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank()) + 
  scale_shape_manual(values = 1:nlevels(azored_k2_hmp2_df_subset$Participant.ID)) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(breaks = c(1,5,10,15,20,24),
                     limits = c(1,24),
                     expand = c(0.04,0.04)) +
  xlab("collection number") + labs(shape = "UC Subjects") -> bottom_mid
# bottom_mid

azored_k2_hmp2_df_subset %>% 
  filter(diagnosis == "CD") %>% 
  ggplot(aes(x = coll_num, y = RA, group = Participant.ID)) +
  geom_point(aes(shape = Participant.ID), size = 2, color = clrs_hmp2[3]) + geom_line(color = clrs_hmp2[3]) +
  theme_bw() + theme(text = element_text(size = 10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank()) + 
  scale_shape_manual(values = 1:nlevels(azored_k2_hmp2_df_subset$Participant.ID)) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(breaks = c(1,5,10,15,20,24),
                     limits = c(1,24),
                     expand = c(0.04,0.04)) + labs(shape = "CD Subjects") -> bottom_right
# bottom_right

top / (bottom_left | bottom_mid | bottom_right) + plot_layout(guides = "collect")
ggsave("../../figures/figure5/figure5_from_R_v2_complex.svg", width = 12, height = 8, units = "in")
ggsave("../../figures/figure5/figure5_from_R_v2_complex.png", width = 12, height = 8, units = "in")
  

top / bottom + plot_layout(guides = "collect")
ggsave("../../figures/figure5/figure5_from_R_v2_simple.svg", width = 15, height = 10, units = "in")
ggsave("../../figures/figure5/figure5_from_R_v2_simple.png", width = 15, height = 10, units = "in")


#### ============================= plotting v1 ============================== ####

clrs_hmp2 <- c("#999999", "#D6B409", "#D55E00")

azored_k2_hmp2_df_subset %>% 
  ggplot(aes(x = coll_num, y = RA, group = Participant.ID, color = diagnosis)) + 
  geom_point(size = 1) + geom_line() +
  facet_wrap(~diagnosis, ncol = 3) +
  theme_bw() + theme(text = element_text(size = 10),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     strip.background = element_blank(),
                     strip.text.x = element_blank()) +
  scale_color_manual(values = clrs_hmp2) + 
  scale_x_continuous(breaks = c(1,5,10,15,20,24),
                     limits = c(1,24),
                     expand = c(0.04,0.04)) +
  ylab("relative abundance (%)") + xlab("collection number") -> top
top

## generating new df for plotting mean and std.err. per individual
azored_k2_hmp2_summary <- azored_k2_hmp2_df_subset %>% 
  group_by(Participant.ID) %>% 
  summarise(diagnosis = diagnosis,
            mean_RA = mean(RA),
            median_RA = median(RA),
            sd_RA = sd(RA),
            N_RA = n(),
            se_RA = sd_RA / sqrt(N_RA)) %>% 
  distinct(diagnosis, Participant.ID, .keep_all = TRUE) %>% 
  arrange(diagnosis) %>% 
  group_by(diagnosis) %>% 
  arrange(mean_RA, .by_group = TRUE)

azored_k2_hmp2_summary %>% 
  filter(diagnosis == "nonIBD") %>% ungroup() -> nonIBD_summary
nonIBD_summary$Participant.ID[5] <- "C3022"
nonIBD_summary$Participant.ID <- factor(nonIBD_summary$Participant.ID, levels = nonIBD_summary$Participant.ID)
nonIBD_summary %>% 
  ggplot(aes(x = Participant.ID, y = mean_RA)) +
  geom_errorbar(aes(ymin = mean_RA - sd_RA, ymax = mean_RA + sd_RA),
                width = 0.2, color = clrs_hmp2[1]) +
  geom_point(shape = "square", size = 2, color = clrs_hmp2[1]) +
  theme_bw() + theme(text = element_text(size = 10),
                     axis.text.x = element_text(angle = 90),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  # scale_color_manual(values = c(clrs_hmp2[1])) +
  ylab("mean relative abundance (%)") + 
  xlab("") -> bottom_left
bottom_left

azored_k2_hmp2_summary %>% 
  filter(diagnosis == "UC") %>% ungroup() -> UC_summary
UC_summary$Participant.ID[4] <- "C3015"
UC_summary$Participant.ID <- factor(UC_summary$Participant.ID, levels = UC_summary$Participant.ID)
UC_summary %>% 
  ggplot(aes(x = Participant.ID, y = mean_RA)) +
  geom_errorbar(aes(ymin = mean_RA - sd_RA, ymax = mean_RA + sd_RA), 
                width = 0.2, color = clrs_hmp2[2]) +
  geom_point(shape = "square", size = 2, color = clrs_hmp2[2])  +
  theme_bw() + theme(text = element_text(size = 10),
                     axis.text.x = element_text(angle = 90),
                     axis.text.y = element_blank(),
                     axis.title.y = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  # scale_color_manual(values = c("red")) +
  xlab("Participant") -> bottom_mid

azored_k2_hmp2_summary %>% 
  filter(diagnosis == "CD") %>% ungroup() -> CD_summary
CD_summary$Participant.ID[3] <- "C3027"
CD_summary$Participant.ID <- factor(CD_summary$Participant.ID, levels = CD_summary$Participant.ID)
CD_summary %>% 
  ggplot(aes(x = Participant.ID, y = mean_RA)) +
  geom_errorbar(aes(ymin = mean_RA - sd_RA, ymax = mean_RA + sd_RA), 
                width = 0.2, color = clrs_hmp2[3]) +
  geom_point(shape = "square", size = 2, color = clrs_hmp2[3]) +
  theme_bw() + theme(text = element_text(size = 10),
                     axis.text.x = element_text(angle = 90),
                     axis.text.y = element_blank(),
                     axis.title.y = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  # scale_color_manual(values = c(clrs_hmp2[3])) +
  xlab("") -> bottom_right

top / (bottom_left + bottom_mid + bottom_right)
ggsave("../../figures/figure5/figure5_from_R_v1.svg", width = 7.5, height = 5, units = "in")
ggsave("../../figures/figure5/figure5_from_R_v1.png", width = 7.5, height = 5, units = "in")


#### ============= statistics performed on figure 5 data ================= ####

## linear mixed effects model ANOVA

# important: (https://www.researchgate.net/post/How-to-compare-two-groups-with-multiple-measurements)

model_nonIBD_CD <- lmer(RA~diagnosis+(1|Participant.ID), data = filter(azored_k2_hmp2_df_subset, diagnosis != "UC"))
summary(model_nonIBD_CD)
anova(model_nonIBD_CD)

model_nonIBD_UC <- lmer(RA~diagnosis+(1|Participant.ID), data = filter(azored_k2_hmp2_df_subset, diagnosis != "CD"))
summary(model_nonIBD_UC)
anova(model_nonIBD_UC)

## getting ranges of median & +/- std. err. for each disease population
azored_k2_hmp2_df_subset %>% 
  filter()

##### ================ basement: gene abundance analysis ================ #####
# 
# hmp2_agg_pData_df %>% 
#   dplyr::count(Participant.ID) %>% 
#   filter(n > 1) %>% 
#   select(Participant.ID) -> multi_sample_subjects
# 
# ## filtereing out subjects for which only one sample is present (ms = multi-sample)
# hmp2_agg_pData_df %>% 
#   filter(Participant.ID %in% multi_sample_subjects$Participant.ID) -> hmp2_agg_pData_ms
# 
# ## adding column containing the collection number
# hmp2_agg_pData_ms$collection <-  as.numeric(substr(hmp2_agg_pData_ms$site_sub_coll, 7, 8))
# hmp2_agg_pData_ms <- hmp2_agg_pData_ms[, c(1:3,17,4:16)]
# 
# ## adding column of gene abundance sums
# hmp2_agg_pData_ms %>% 
#   mutate(sum = rowSums(across(where(is.numeric)))) -> hmp2_agg_pData_ms
# # hmp2_agg_pData_ms <- hmp2_agg_pData_ms[, c(1:4,18,5:17)]
# 
# ## plotting azored gene abundance over time for nonIBD, UC and CD groups
# ## sum of all azored gene abundances
# ggplot(hmp2_agg_pData_ms, aes(x = collection, y = sum, group = Participant.ID)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~diagnosis, ncol = 1) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank()) +
#   scale_y_continuous(limits = c(0,500)) 
# 
# ## clades
# ggplot(hmp2_agg_pData_ms, aes(x = collection, y = cladeIVb, group = Participant.ID)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~diagnosis, ncol = 1) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank()) #+
#   # scale_y_continuous(limits = c(0,100))
# 
# ggplot(hmp2_agg_pData_ms, aes(x = collection, y = cladeIVb, group = Participant.ID)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~diagnosis, ncol = 1) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank()) +
# scale_y_continuous(limits = c(0,100))
# 
# ## other types of visualizations
# 
# hmp2_agg_pData_ms %>% 
#   select(Participant.ID, diagnosis, sum) %>% 
#   group_by(Participant.ID, diagnosis) %>% 
#   summarise_each(funs(mean, sd, se = sd(.)/sqrt(n()))) -> hmp2_sum_stats
# 
# hmp2_sum_stats <- hmp2_sum_stats[order(-hmp2_sum_stats$mean), ]
# hmp2_sum_stats$Participant.ID <- factor(hmp2_sum_stats$Participant.ID, levels = hmp2_sum_stats$Participant.ID)
# 
# hmp2_sum_stats$diagnosis <- factor(hmp2_sum_stats$diagnosis, levels = c("nonIBD", "UC", "CD"))
# 
# ggplot(hmp2_sum_stats, aes(x = Participant.ID, y = mean, fill = diagnosis)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = mean - se, ymax = mean + se)) +
#   theme_bw() + theme(panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      axis.text.x = element_blank(),
#                      axis.ticks.x = element_blank(),
#                      text = element_text(size = 10)) +
#   scale_fill_manual(values = clrs2) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,1210)) +
#   xlab("Subjects with > 1 sample (n = 101)") + ylab("averaged gene abunance") +
#   ggtitle("Averaged abundance of azoreductase genes")
#   # geom_point(hmp2_agg_pData_ms, aes(x = Participant.ID, y = sum))
# 
# ggplot(hmp2_agg_pData_ms, aes(x = Participant.ID, y = sum)) +
#   geom_boxplot()
# 
# ggplot(hmp2_agg_pData_ms, aes(x = Participant.ID, y = sum)) +
#   geom_bar(stat = "identity")

        
