# ======================== risk_relative_abundances.R ======================= #
#' description: script for comparing the relative abundances of azoreducing
#' species in cMD, HMP2, PRISM and CIB
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
#' 
#' Author: Domenick J. Braccia  
#' Last Update: 1 Nov 2021
# =========================================================================== #

##### =================================================================== #####
print("- calculating full column sums")

### cMD
azored_k2_RA <- k2_RA[azored_taxa_hits, ]
azored_k2_cs <- colSums(azored_k2_RA)

### hmp2
azored_k2_hmp2_RA <- k2_hmp2_RA[azored_taxa_hits, ]
azored_k2_hmp2_cs <- colSums(azored_k2_hmp2_RA)

### prism
azored_k2_prism_RA <- k2_prism_RA[azored_taxa_hits, ]
azored_k2_prism_cs <- colSums(azored_k2_prism_RA)

### CIB
azored_k2_cib_RA <- k2_cib_RA[azored_taxa_hits, ]
azored_k2_cib_cs <- colSums(azored_k2_cib_RA)

##### =================================================================== #####
print("- filtering pData to only include populations for plotting")

### cMD
# filter risk samples
k2_pData %>%
  filter(study_condition == c("IBD", "CRC")) -> ibd_crc
unique(ibd_crc$dataset_name) -> ibd_crc_studies
kraken_risk_pData <- k2_pData[k2_pData$dataset_name %in% ibd_crc_studies, ]

#get colSums
azored_k2_cs_risk <- azored_k2_cs[names(azored_k2_cs) %in% rownames(kraken_risk_pData)]

### hmp2
azored_k2_hmp2_cs_risk <- azored_k2_hmp2_cs[names(azored_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]

### prism
azored_k2_prism_cs_risk <- azored_k2_prism_cs[names(azored_k2_prism_cs) %in% k2_prism_pData$Sample]

### cib
## filter risk samples
subset_k2_cib_metadata <- cbind(select(k2_cib_pData, Sample_SRA),
                                select(k2_cib_metadata, disease, Time, Sample.Name))
subset_k2_cib_metadata$Time[is.na(subset_k2_cib_metadata$Time)] <- 1
subset_k2_cib_metadata <- subset_k2_cib_metadata %>%
  filter(Time == c(1)) 
k2_cib_risk_metadata <- subset_k2_cib_metadata$disease
names(k2_cib_risk_metadata) <- subset_k2_cib_metadata$Sample_SRA 

## get colSums
azored_k2_cib_cs_risk <- azored_k2_cib_cs[names(azored_k2_cib_cs) %in% k2_cib_pData$Sample_SRA]


##### =========================== PLOTTING ============================== #####
print("- prepping dataframes for plotting")

### cMD
cMD_study_cond <- as.character(kraken_risk_pData$study_condition)
names(cMD_study_cond) <- rownames(kraken_risk_pData)

azored_k2_df <- as.data.frame(cbind(azored_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(azored_k2_cs_risk)]))
colnames(azored_k2_df) <- c("RA", "population")
azored_k2_df$RA <- as.numeric(azored_k2_df$RA)
azored_k2_df$func <- rep("cysteine", length(azored_k2_cs_risk))
azored_k2_df$population <- factor(azored_k2_df$population, levels = c("control", "IBD", "adenoma", "CRC"))

#### HMP2 ####
subset_k2_pData <- select(k2_hmp2_pData, External.ID, diagnosis)
k2_risk_pData <- subset_k2_pData$diagnosis
names(k2_risk_pData) <- subset_k2_pData$External.ID

azored_k2_hmp2_df <- as.data.frame(cbind(azored_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(azored_k2_hmp2_cs_risk)]))
colnames(azored_k2_hmp2_df) <- c("RA", "diagnosis")
azored_k2_hmp2_df$RA <- as.numeric(azored_k2_hmp2_df$RA)
azored_k2_hmp2_df$func <- rep("cysteine", length(azored_k2_hmp2_cs_risk))
azored_k2_hmp2_df$diagnosis <- factor(azored_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
# colnames(azored_k2_hmp2_df) <- c("RA", "diagnosis", "func")

#### PRISM ####
subset_k2_prism_pData <- select(k2_prism_pData, Sample, Diagnosis)
k2_prism_risk_pData <- subset_k2_prism_pData$Diagnosis
names(k2_prism_risk_pData) <- subset_k2_prism_pData$Sample

azored_k2_prism_df <- as.data.frame(cbind(azored_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(azored_k2_prism_cs_risk)]))
colnames(azored_k2_prism_df) <- c("RA", "diagnosis")
azored_k2_prism_df$RA <- as.numeric(azored_k2_prism_df$RA)
azored_k2_prism_df$func <- rep("cysteine", length(azored_k2_prism_cs_risk))
azored_k2_prism_df$diagnosis <- factor(azored_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

#### CIB ####
azored_k2_cib_df <- as.data.frame(cbind(azored_k2_cib_cs_risk, k2_cib_risk_metadata ))
colnames(azored_k2_cib_df) <- c("RA", "diagnosis")
azored_k2_cib_df$RA <- as.numeric(azored_k2_cib_df$RA)
azored_k2_cib_df$func <- rep("azored", length(azored_k2_cib_cs_risk))
azored_k2_cib_df$diagnosis <- factor(azored_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))


##### ==================== plotting risk cMD sample data ================ #####
print("- plotting risk cMD sample data")

clrs1 <- c("#999999", "#E69F00", "#56B4E9", "#0072B2")

azored_k2_boxplot <- ggplot(azored_k2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = population)) +
  # geom_violin(aes(x = func, y = RA, fill = population), scale = "width") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100),
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs1) +
  ylab("relative abundance, (%)") + xlab("azored") + ggtitle("curatedMetagenomicData")
azored_k2_boxplot

##### ==================== plotting risk hmp2 sample data ================ #####
print("- plotting hmp2 data")

clrs2 <- c("#999999", "#F0E442", "#D55E00")

azored_k2_hmp2_boxplot <- ggplot(azored_k2_hmp2_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs2) +
  ylab("relative abundance (%)") + xlab("azored") + ggtitle("HMP2")
azored_k2_hmp2_boxplot

##### ==================== plotting prism data ========================== #####
print("- plotting prism data")

clrs3 <- c("#999999", "#F0E442", "#D55E00")

azored_k2_prism_boxplot <- ggplot(azored_k2_prism_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs3) +
  ylab("relative abundance, (%)") + xlab("azored") ggtitle("PRISM")
azored_k2_prism_boxplot


##### ==================== plotting cib data ============================ #####
print("- plotting cib data")

clrs4 <- c("#999999", "#D55E00")

azored_k2_cib_boxplot <- ggplot(azored_k2_cib_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = clrs4) +
  ylab("relative abundance, (%)") + xlab("azored") + ggtitle("Lewis et al. 2015")
azored_k2_cib_boxplot
