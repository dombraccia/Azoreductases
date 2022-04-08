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
tic()
##### =================================================================== #####
print("- separating known and putative azoreducing species")

## changing "|" to ";" for taxa names
azored_taxa_hits <- gsub("\\|", ";", azored_taxa_hits)
rownames(k2_RA) <- gsub("\\|", ";", rownames(k2_RA))
rownames(k2_hmp2_RA) <- gsub("\\|", ";", rownames(k2_hmp2_RA))
rownames(k2_prism_RA) <- gsub("\\|", ";", rownames(k2_prism_RA))
rownames(k2_cib_RA) <- gsub("\\|", ";", rownames(k2_cib_RA))
uhgg_genid2taxa$taxa_ID <- gsub("\\|", ";", uhgg_genid2taxa$taxa_ID)

## defining short functions for taxa names search & retreiving RA values
search_azored_spp <- function(name) {str_subset(azored_taxa_hits, name)}
search_all_spp <- function(name) {str_subset(uhgg_genid2taxa$taxa_ID, name)}

#### ===== gathering all known azoreducing species from literature ======= ####

## Suzuki 2019
azored_taxa_known <- character()
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Bacillus_[A-Z,_]*cereus")) ## aka Bacillus sp. B29
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Bacillus_[A-Z,_]*subtilis")) 
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Clostridium_[A-Z,_]*perfringens"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Enterococcus_[A-Z,_]*faecalis"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Enterococcus_[A-Z,_]*faecium")) ## aka Streptococcus faecium 
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Escherichia_[A-Z,_]*coli"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Pseudomonas_[A-Z,_]*putida"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Pseudomonas_[A-Z,_]*aeruginosa"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("Staphylococcus_[A-Z,_]*aureus"))

## Brown 1981
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bacteroides_[A-Z,_]*uniformis"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Erysipelatoclostridium_[A-Z,_]*ramosum")) ## aka Clostridium ramosum (Brown et al. 1981)
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Clostridium_[A-Z,_]*paraputrificum"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Clostridium_[A-Z,_]*sporogenes"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Ruminococcus_[A-Z,_]*fragilis"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bacteroides_[A-Z,_]*bromii"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Veillonella_[A-Z,_]*parvula"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Anaerococcus_[A-Z,_]*prevotii")) ## aka Peptococcus prevotii
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bifidobacterium_[A-Z,_]*adolescentis"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Collinsella_[A-Z,_]*aerofaciens")) ## aka Eubacterium aerofaciens
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Salmonella_[A-Z,_]*enterica")) ## aka Salmonella typhimurium/paratyphi
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Proteus_[A-Z,_]*vulgaris")) 
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Klebsiella_[A-Z,_]*pneumoniae")) 

## Rafii 1990
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Tyzzerella_[A-Z,_]*nexilis")) ## aka Clostridium nexile

## McBain and Macfarlane 1998
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bacteroides_[A-Z,_]*ovatus"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Clostridium_[A-Z,_]*septicum"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Clostridium_[A-Z,_]*butyricum"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Lactobacillus_[A-Z,_]*acidophilus"))

## Xu 2010 (Sudan dyes paper)
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bacteroides_[A-Z,_]*caccae"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bacteroides_[A-Z,_]*thetaiotaomicron"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Parabacteroides_[A-Z,_]*distasonis"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bifidobacterium_[A-Z,_]*adolescentis"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Bifidobacterium_[A-Z,_]*infantis"))

## Saratale 2011
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Lactobacillus_[A-Z,_]*fermentum"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Aeromonas_[A-Z,_]*hydrophila"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Desulfovibrio_[A-Z,_]*desulfuricans"))
azored_taxa_known <- append(azored_taxa_known, search_all_spp("s__Pseudomonas_[A-Z,_]*luteola"))

 
search_azored_spp("_[A-Z,_]*")

search_all_spp("s__Aeromonas_[A-Z,_]*hydrophila")
search_all_spp("")

azored_taxa_putative <- azored_taxa_hits[!grepl(paste(azored_taxa_known, collapse = "|"), azored_taxa_hits)]

print("- calculating full column sums of known and putative azored. spp.")

### cMD
azored_k2_cMD_RA <- k2_RA[azored_taxa_hits, ]
azored_k2_cMD_cs <- colSums(azored_k2_cMD_RA)
azored_known_k2_cMD_RA <- k2_RA[azored_taxa_known, ]
azored_known_k2_cMD_cs <- colSums(azored_known_k2_cMD_RA)
azored_putative_k2_cMD_RA <- k2_RA[azored_taxa_putative, ]
azored_putative_k2_cMD_cs <- colSums(azored_putative_k2_cMD_RA)

### hmp2
k2_hmp2_RA <- as.data.frame(k2_hmp2_RA)
azored_k2_hmp2_RA <- k2_hmp2_RA[azored_taxa_hits, ]
azored_k2_hmp2_cs <- colSums(azored_k2_hmp2_RA)
azored_known_k2_hmp2_RA <- k2_hmp2_RA[azored_taxa_known, ]
azored_known_k2_hmp2_cs <- colSums(azored_known_k2_hmp2_RA)
azored_putative_k2_hmp2_RA <- k2_hmp2_RA[azored_taxa_putative, ]
azored_putative_k2_hmp2_cs <- colSums(azored_putative_k2_hmp2_RA)

### prism
k2_prism_RA <- as.data.frame(k2_prism_RA)
azored_k2_prism_RA <- k2_prism_RA[azored_taxa_hits, ]
azored_k2_prism_cs <- colSums(azored_k2_prism_RA)
azored_known_k2_prism_RA <- k2_prism_RA[azored_taxa_known, ]
azored_known_k2_prism_cs <- colSums(azored_known_k2_prism_RA)
azored_putative_k2_prism_RA <- k2_prism_RA[azored_taxa_putative, ]
azored_putative_k2_prism_cs <- colSums(azored_putative_k2_prism_RA)

### CIB
k2_cib_RA <- as.data.frame(k2_cib_RA)
azored_k2_cib_RA <- k2_cib_RA[azored_taxa_hits, ]
azored_k2_cib_cs <- colSums(azored_k2_cib_RA)
azored_known_k2_cib_RA <- k2_cib_RA[azored_taxa_known, ]
azored_known_k2_cib_cs <- colSums(azored_known_k2_cib_RA)
azored_putative_k2_cib_RA <- k2_RA[azored_taxa_putative, ]
azored_putative_k2_cib_cs <- colSums(azored_putative_k2_cib_RA)

##### =================================================================== #####
print("- filtering pData to only include populations for plotting")

### cMD
# filter risk samples
k2_pData %>%
  filter(study_condition == c("IBD", "CRC")) -> ibd_crc
unique(ibd_crc$dataset_name) -> ibd_crc_studies
kraken_risk_pData <- k2_pData[k2_pData$dataset_name %in% ibd_crc_studies, ]

# get colSums
azored_k2_cs_risk <- azored_k2_cMD_cs[names(azored_k2_cMD_cs) %in% rownames(kraken_risk_pData)]
azored_known_k2_cs_risk <- azored_known_k2_cMD_cs[names(azored_known_k2_cMD_cs) %in% rownames(kraken_risk_pData)]
azored_putative_k2_cs_risk <- azored_putative_k2_cMD_cs[names(azored_putative_k2_cMD_cs) %in% rownames(kraken_risk_pData)]

### hmp2
azored_k2_hmp2_cs_risk <- azored_k2_hmp2_cs[names(azored_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
azored_known_k2_hmp2_cs_risk <- azored_known_k2_hmp2_cs[names(azored_known_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
azored_putative_k2_hmp2_cs_risk <- azored_putative_k2_hmp2_cs[names(azored_putative_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]

### prism
azored_k2_prism_cs_risk <- azored_k2_prism_cs[names(azored_k2_prism_cs) %in% k2_prism_pData$Sample]
azored_known_k2_prism_cs_risk <- azored_known_k2_prism_cs[names(azored_known_k2_prism_cs) %in% k2_prism_pData$Sample]
azored_putative_k2_prism_cs_risk <- azored_putative_k2_prism_cs[names(azored_putative_k2_prism_cs) %in% k2_prism_pData$Sample]

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
azored_known_k2_cib_cs_risk <- azored_known_k2_cib_cs[names(azored_known_k2_cib_cs) %in% k2_cib_pData$Sample_SRA]
azored_putative_k2_cib_cs_risk <- azored_putative_k2_cib_cs[names(azored_putative_k2_cib_cs) %in% k2_cib_pData$Sample_SRA]

##### =========================== PLOTTING ============================== #####
print("- prepping dataframes for plotting")

#### cMD ####
cMD_study_cond <- as.character(kraken_risk_pData$study_condition)
names(cMD_study_cond) <- rownames(kraken_risk_pData)

## all azored spp.
azored_k2_cMD_df <- as.data.frame(cbind(azored_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(azored_k2_cs_risk)]))
colnames(azored_k2_cMD_df) <- c("RA", "diagnosis")
azored_k2_cMD_df$RA <- as.numeric(azored_k2_cMD_df$RA)
azored_k2_cMD_df$func <- rep("azored", length(azored_k2_cs_risk))
azored_k2_cMD_df$diagnosis <- factor(azored_k2_cMD_df$diagnosis, levels = c("control", "IBD", "adenoma", "CRC"))

## known azored spp.
azored_known_k2_cMD_df <- as.data.frame(cbind(azored_known_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(azored_known_k2_cs_risk)]))
colnames(azored_known_k2_cMD_df) <- c("RA", "diagnosis")
azored_known_k2_cMD_df$RA <- as.numeric(azored_known_k2_cMD_df$RA)
azored_known_k2_cMD_df$func <- rep("azored", length(azored_known_k2_cs_risk))
azored_known_k2_cMD_df$diagnosis <- factor(azored_known_k2_cMD_df$diagnosis, levels = c("control", "IBD", "adenoma", "CRC"))

## putative azored spp.
azored_putative_k2_cMD_df <- as.data.frame(cbind(azored_putative_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(azored_putative_k2_cs_risk)]))
colnames(azored_putative_k2_cMD_df) <- c("RA", "diagnosis")
azored_putative_k2_cMD_df$RA <- as.numeric(azored_putative_k2_cMD_df$RA)
azored_putative_k2_cMD_df$func <- rep("azored", length(azored_putative_k2_cs_risk))
azored_putative_k2_cMD_df$diagnosis <- factor(azored_putative_k2_cMD_df$diagnosis, levels = c("control", "IBD", "adenoma", "CRC"))

#### HMP2 ####
subset_k2_pData <- select(k2_hmp2_pData, External.ID, diagnosis)
k2_risk_pData <- subset_k2_pData$diagnosis
names(k2_risk_pData) <- subset_k2_pData$External.ID

## all azored spp (known + putative)
azored_k2_hmp2_df <- as.data.frame(cbind(azored_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(azored_k2_hmp2_cs_risk)]))
colnames(azored_k2_hmp2_df) <- c("RA", "diagnosis")
azored_k2_hmp2_df$RA <- as.numeric(azored_k2_hmp2_df$RA)
azored_k2_hmp2_df$func <- rep("azored", length(azored_k2_hmp2_cs_risk))
azored_k2_hmp2_df$diagnosis <- factor(azored_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))

## known azored spp.
azored_known_k2_hmp2_df <- as.data.frame(cbind(azored_known_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(azored_known_k2_hmp2_cs_risk)]))
colnames(azored_known_k2_hmp2_df) <- c("RA", "diagnosis")
azored_known_k2_hmp2_df$RA <- as.numeric(azored_known_k2_hmp2_df$RA)
azored_known_k2_hmp2_df$func <- rep("azored", length(azored_known_k2_hmp2_cs_risk))
azored_known_k2_hmp2_df$diagnosis <- factor(azored_known_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))

## putative azored spp.
azored_putative_k2_hmp2_df <- as.data.frame(cbind(azored_putative_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(azored_putative_k2_hmp2_cs_risk)]))
colnames(azored_putative_k2_hmp2_df) <- c("RA", "diagnosis")
azored_putative_k2_hmp2_df$RA <- as.numeric(azored_putative_k2_hmp2_df$RA)
azored_putative_k2_hmp2_df$func <- rep("azored", length(azored_putative_k2_hmp2_cs_risk))
azored_putative_k2_hmp2_df$diagnosis <- factor(azored_putative_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))

#### PRISM ####
subset_k2_prism_pData <- select(k2_prism_pData, Sample, Diagnosis)
k2_prism_risk_pData <- subset_k2_prism_pData$Diagnosis
names(k2_prism_risk_pData) <- subset_k2_prism_pData$Sample

## all azored spp (known + putative)
azored_k2_prism_df <- as.data.frame(cbind(azored_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(azored_k2_prism_cs_risk)]))
colnames(azored_k2_prism_df) <- c("RA", "diagnosis")
azored_k2_prism_df$RA <- as.numeric(azored_k2_prism_df$RA)
azored_k2_prism_df$func <- rep("azored", length(azored_k2_prism_cs_risk))
azored_k2_prism_df$diagnosis <- factor(azored_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

## known azored spp.
azored_known_k2_prism_df <- as.data.frame(cbind(azored_known_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(azored_known_k2_prism_cs_risk)]))
colnames(azored_known_k2_prism_df) <- c("RA", "diagnosis")
azored_known_k2_prism_df$RA <- as.numeric(azored_known_k2_prism_df$RA)
azored_known_k2_prism_df$func <- rep("azored", length(azored_known_k2_prism_cs_risk))
azored_known_k2_prism_df %>% 
  mutate(diagnosis = case_when(diagnosis == "HC" ~ "nonIBD",
                               diagnosis == "UC" ~ "UC",
                               diagnosis == "CD" ~ "CD")) -> azored_known_k2_prism_df
azored_known_k2_prism_df$diagnosis <- factor(azored_known_k2_prism_df$diagnosis, levels = c("nonIBD", "UC", "CD"))

## putative azored spp.
azored_putative_k2_prism_df <- as.data.frame(cbind(azored_putative_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(azored_putative_k2_prism_cs_risk)]))
colnames(azored_putative_k2_prism_df) <- c("RA", "diagnosis")
azored_putative_k2_prism_df$RA <- as.numeric(azored_putative_k2_prism_df$RA)
azored_putative_k2_prism_df$func <- rep("azored", length(azored_putative_k2_prism_cs_risk))
azored_putative_k2_prism_df$diagnosis <- factor(azored_putative_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))

# #### CIB ####
# 
# ## all azored spp (known + putative)
# azored_k2_cib_df <- as.data.frame(cbind(azored_k2_cib_cs_risk, k2_cib_risk_metadata ))
# colnames(azored_k2_cib_df) <- c("RA", "diagnosis")
# azored_k2_cib_df$RA <- as.numeric(azored_k2_cib_df$RA)
# azored_k2_cib_df$func <- rep("azored", length(azored_k2_cib_cs_risk))
# azored_k2_cib_df$diagnosis <- factor(azored_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))
# 
# ## known azored spp
# azored_known_k2_cib_df <- as.data.frame(cbind(azored_known_k2_cib_cs_risk, k2_cib_risk_metadata ))
# colnames(azored_known_k2_cib_df) <- c("RA", "diagnosis")
# azored_known_k2_cib_df$RA <- as.numeric(azored_known_k2_cib_df$RA)
# azored_known_k2_cib_df$func <- rep("azored", length(azored_known_k2_cib_cs_risk))
# azored_known_k2_cib_df$diagnosis <- factor(azored_known_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))
# 
# ## putative azored spp
# azored_putative_k2_cib_df <- as.data.frame(cbind(azored_putative_k2_cib_cs_risk, k2_cib_risk_metadata ))
# colnames(azored_putative_k2_cib_df) <- c("RA", "diagnosis")
# azored_putative_k2_cib_df$RA <- as.numeric(azored_putative_k2_cib_df$RA)
# azored_putative_k2_cib_df$func <- rep("azored", length(azored_putative_k2_cib_cs_risk))
# azored_putative_k2_cib_df$diagnosis <- factor(azored_putative_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

## DEFINING PREP FUNCTIONS
prep_cMD <- function(taxa) {
  
  ## subset by taxa (input can be char vector or single char)
  taxa_k2_cMD_cs <- azored_k2_cMD_RA[grep(paste(taxa, collapse = "|"), rownames(azored_k2_cMD_RA), value = TRUE), ] %>% colSums()
  
  ## subset by risk samples
  taxa_k2_cs_risk <- taxa_k2_cMD_cs[names(taxa_k2_cMD_cs) %in% rownames(kraken_risk_pData)]
  
  ## building long dataframe
  taxa_k2_cMD_df <- as.data.frame(cbind(taxa_k2_cs_risk, cMD_study_cond = cMD_study_cond[names(taxa_k2_cs_risk)]))
  colnames(taxa_k2_cMD_df) <- c("RA", "diagnosis")
  taxa_k2_cMD_df$RA <- as.numeric(taxa_k2_cMD_df$RA)
  taxa_k2_cMD_df$func <- rep("azored", length(taxa_k2_cs_risk))
  taxa_k2_cMD_df$diagnosis <- factor(taxa_k2_cMD_df$diagnosis, levels = c("control", "IBD", "adenoma", "CRC"))
  
  return(taxa_k2_cMD_df)
}

prep_hmp2 <- function(taxa) {
  
  ## subset by taxa (input can be char vector or single char)
  taxa_k2_hmp2_cs <- azored_k2_hmp2_RA[grep(paste(taxa, collapse = "|"), rownames(azored_k2_hmp2_RA), value = TRUE), ] %>% colSums()
  
  ## subset by risk samples
  taxa_k2_hmp2_cs_risk <- taxa_k2_hmp2_cs[names(taxa_k2_hmp2_cs) %in% k2_hmp2_pData$External.ID]
  
  ## building long dataframe
  taxa_k2_hmp2_df <- as.data.frame(cbind(taxa_k2_hmp2_cs_risk, k2_risk_pData = k2_risk_pData[names(taxa_k2_hmp2_cs_risk)]))
  colnames(taxa_k2_hmp2_df) <- c("RA", "diagnosis")
  taxa_k2_hmp2_df$RA <- as.numeric(taxa_k2_hmp2_df$RA)
  taxa_k2_hmp2_df$func <- rep("azored", length(taxa_k2_hmp2_cs_risk))
  taxa_k2_hmp2_df$diagnosis <- factor(taxa_k2_hmp2_df$diagnosis, levels = c("nonIBD", "UC", "CD"))
  
  return(taxa_k2_hmp2_df)
}

prep_prism <- function(taxa) {
  
  ## subset by taxa (input can be char vector or single char)
  taxa_k2_prism_cs <- azored_k2_prism_RA[grep(paste(taxa, collapse = "|"), rownames(azored_k2_prism_RA), value = TRUE), ] %>% colSums()
  
  ## subset by risk samples
  taxa_k2_prism_cs_risk <- taxa_k2_prism_cs[names(taxa_k2_prism_cs) %in% k2_prism_pData$Sample]
  
  ## building long dataframe
  taxa_k2_prism_df <- as.data.frame(cbind(taxa_k2_prism_cs_risk, k2_prism_risk_pData = k2_prism_risk_pData[names(taxa_k2_prism_cs_risk)]))
  colnames(taxa_k2_prism_df) <- c("RA", "diagnosis")
  taxa_k2_prism_df$RA <- as.numeric(taxa_k2_prism_df$RA)
  taxa_k2_prism_df$func <- rep("azored", length(taxa_k2_prism_cs_risk))
  taxa_k2_prism_df$diagnosis <- factor(taxa_k2_prism_df$diagnosis, levels = c("HC", "UC", "CD"))
  
  return(taxa_k2_prism_df)
}

prep_cib <- function(taxa) {

  ## subset by taxa (input can be char vector or single char)
  taxa_k2_cib_cs <- azored_k2_cib_RA[grep(paste(taxa, collapse = "|"), rownames(azored_k2_cib_RA), value = TRUE), ] %>% colSums()

  ## subset by risk samples
  taxa_k2_cib_cs_risk <- taxa_k2_cib_cs[names(taxa_k2_cib_cs) %in% k2_cib_pData$Sample_SRA]

  ## building long dataframe
  taxa_k2_cib_df <- as.data.frame(cbind(taxa_k2_cib_cs_risk, k2_cib_risk_metadata ))
  colnames(taxa_k2_cib_df) <- c("RA", "diagnosis")
  taxa_k2_cib_df$RA <- as.numeric(taxa_k2_cib_df$RA)
  taxa_k2_cib_df$func <- rep("azored", length(taxa_k2_cib_cs_risk))
  taxa_k2_cib_df$diagnosis <- factor(taxa_k2_cib_df$diagnosis, levels = c("Control", "Crohn"))

  return(taxa_k2_cib_df)
}

##### ============ STATISTICS FOR ASM ABSTRACT (23 Nov 2021) ============ #####

## cMD + hmp2 + prism, healthy controls
cMD_study_cond <- as.character(k2_pData$study_condition)
names(cMD_study_cond) <- rownames(k2_pData)

azored_k2_cMD_df_tmp <- as.data.frame(cbind(azored_k2_cMD_cs, cMD_study_cond = cMD_study_cond[names(azored_k2_cMD_cs)]))

colnames(azored_k2_cMD_df_tmp) <- c("RA", "diagnosis")
azored_k2_cMD_df_tmp$RA <- as.numeric(azored_k2_cMD_df_tmp$RA)
azored_k2_cMD_df_tmp$func <- rep("azored", length(azored_k2_cMD_cs))

all_k2_df <- rbind(azored_k2_cMD_df_tmp, azored_k2_hmp2_df, azored_k2_prism_df)

all_k2_df %>% 
  filter(diagnosis == "control" | diagnosis == "nonIBD" | diagnosis == "HC") %>% 
  select(RA) %>%  summary()
all_k2_df %>% 
  filter(diagnosis == "control" | diagnosis == "nonIBD" | diagnosis == "HC") %>% 
  nrow()

## cMD, CRC samples
glimpse(azored_k2_cMD_df_tmp)
azored_k2_cMD_df_tmp %>% 
  filter(diagnosis == "CRC") %>% select(RA) %>% summary()
azored_k2_cMD_df_tmp %>% 
  filter(diagnosis == "CRC") %>% nrow()


## hmp2 + prism, IBD samples
azored_k2_ibd_df <- rbind(azored_k2_hmp2_df, azored_k2_prism_df)
azored_k2_ibd_df %>% 
  filter(diagnosis == "UC" | diagnosis == "CD") %>% select(RA) %>% summary()
azored_k2_ibd_df %>% 
  filter(diagnosis == "UC" | diagnosis == "CD") %>% nrow()


##### ===== plotting RA of experimentally validated azored species ====== #####

## defining color vectors
clrs_cMD <- c("#999999", "#E69F00", "#56B4E9", "#0072B2")
clrs_hmp2 <- c("#999999", "#F0E442", "#D55E00")
clrs_prism <- c("#999999", "#F0E442", "#D55E00")


## defining plot functions for each dataset
plot_RA <- function(df, x, fill, clrs) {
  ggplot(df, aes_string(x = x, y = fill, fill = fill)) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
    scale_x_continuous(trans = "log",
                       breaks = c(1e-5,1e-4,1e-3, 1e-2, 0.1, 1, 10, 100),
                       labels = c(-5:2),
                       limits = c(1e-6,100)) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       text = element_text(size = 8),
                       axis.title.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.text.y = element_blank()) +
    scale_y_continuous(expand = c(0,0.01),
                       limits = c(0,0.4)) +
    scale_fill_manual(values = clrs)
}

##### ======================== boxplots ===================== #####
print("- plotting risk cMD sample data")

azored_k2_cMD_boxplot <- ggplot(azored_k2_cMD_df) + 
  geom_boxplot(aes(x = func, y = RA, fill = diagnosis)) +
  # geom_violin(aes(x = func, y = RA, fill = diagnosis), scale = "width") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        # axis.title.y = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 100),
                     minor_breaks = NULL) +
  scale_fill_manual(values = clrs_cMD) +
  ylab("relative abundance, (%)") + xlab("azored") + ggtitle("curatedMetagenomicData")
azored_k2_cMD_boxplot

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
  ylab("relative abundance, %") + xlab("azored") + ggtitle("HMP2")
azored_k2_hmp2_boxplot

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
  ylab("relative abundance, (%)") + xlab("azored") + ggtitle("PRISM")
azored_k2_prism_boxplot 

#### plotting summed relative abundances of azored. spp (density) ####

cMD_putative <- azored_putative_k2_cMD_df %>% 
  # filter(RA > 0) %>% 
  ggplot(aes(x = RA, fill = diagnosis)) +
  # geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
  geom_density(alpha = 0.5) +
  theme_bw() + theme(legend.position = "none",
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size = 8),
                     axis.text = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()
  ) +
  scale_color_manual(values = clrs_cMD) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.17)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1e2)) +
  ggtitle("Putative azoreducing spp.")
cMD_putative

# cMD_known <-  azored_known_k2_cMD_df %>% 
#   ggplot(aes(x = RA, fill = diagnosis)) +
#   geom_density(alpha = 0.5) +
#   theme_bw() + theme(legend.position = "none",
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      text = element_text(size = 8),
#                      axis.text.x = element_blank(),
#                      axis.title = element_blank(),
#                      # axis.ticks.x = element_blank()
#   ) +
#   scale_color_manual(values = clrs_cMD) +
#   scale_fill_manual(values = clrs_cMD) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 0.17)) +
#   scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
#   ggtitle("Known azoreducing spp.")
# cMD_known

# hmp2_putative <- azored_putative_k2_hmp2_df %>% 
#   ggplot(aes(x = RA, fill = diagnosis)) +
#   # geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
#   geom_density(alpha = 0.5) +
#   theme_bw() + theme(legend.position = "none",
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      text = element_text(size = 8),
#                      axis.text = element_blank(),
#                      axis.title = element_blank(),
#                      axis.ticks.y = element_blank()) +
#   scale_color_manual(values = clrs2) +
#   scale_fill_manual(values = clrs2) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 0.17)) +
#   scale_x_continuous(expand = c(0,0), limits = c(0, 100))
# hmp2_putative

# hmp2_known <-  azored_known_k2_hmp2_df %>% 
#   ggplot(aes(x = RA, fill = diagnosis)) +
#   geom_density(alpha = 0.5) +
#   theme_bw() + theme(legend.position = "none",
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      text = element_text(size = 8),
#                      axis.text.x = element_blank(),
#                      axis.title = element_blank()
#   ) +
#   scale_color_manual(values = clrs_hmp2) +
#   scale_fill_manual(values = clrs_hmp2) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 0.17)) +
#   scale_x_continuous(expand = c(0,0), limits = c(0, 100)) #+
#   # ggtitle("Known azoreducing spp.")
# hmp2_known

# prism_putative <- azored_putative_k2_prism_df %>% 
#   ggplot(aes(x = RA, fill = diagnosis)) +
#   # geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
#   geom_density(alpha = 0.5) +
#   theme_bw() + theme(legend.position = "none",
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      text = element_text(size = 8),
#                      axis.title.y = element_blank(),
#                      axis.text.y = element_blank(),
#                      axis.ticks.y = element_blank()
#                      ) +
#   scale_color_manual(values = clrs2) +
#   scale_fill_manual(values = clrs2) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 0.17)) +
#   scale_x_continuous(expand = c(0,0), limits = c(0, 100)) +
#   xlab("relative abundance, %")
# prism_putative

# prism_known <-  azored_known_k2_prism_df %>% 
#   ggplot(aes(x = RA, fill = diagnosis)) +
#   geom_density(alpha = 0.5) +
#   theme_bw() + theme(legend.position = "none",
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      axis.title.y = element_blank(),
#                      # axis.text.y = element_blank(),
#                      # axis.ticks.y = element_blank(),
#                      text = element_text(size = 8)) +
#   scale_color_manual(values = clrs_prism) +
#   scale_fill_manual(values = clrs_prism) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 0.17)) +
#   scale_x_continuous(expand = c(0,0), limits = c(0, 100)) +
# xlab("relative abundance, %")
# prism_known




##### determining highest RA of putative azored spp. by rowMeans #####

cMD_azored_RA_rowmeans <- rowMeans(azored_putative_k2_cMD_RA)
names(cMD_azored_RA_rowmeans) <- rownames(azored_putative_k2_cMD_RA)
tail(sort(cMD_azored_RA_rowmeans), 20)

hmp2_azored_RA_rowmeans <- rowMeans(azored_putative_k2_hmp2_RA)
names(hmp2_azored_RA_rowmeans) <- rownames(azored_putative_k2_hmp2_RA)
tail(sort(hmp2_azored_RA_rowmeans), 20)

prism_azored_RA_rowmeans <- rowMeans(azored_putative_k2_prism_RA)
names(prism_azored_RA_rowmeans) <- rownames(azored_putative_k2_prism_RA)
tail(sort(prism_azored_RA_rowmeans), 20)

##### determining highest RA of putative azored spp. by rowMedians #####

cMD_azored_RA_rowmedians <- rowMedians(as.matrix(azored_putative_k2_cMD_RA))
names(cMD_azored_RA_rowmedians) <- rownames(azored_putative_k2_cMD_RA)
tail(sort(cMD_azored_RA_rowmedians), 20)

hmp2_azored_RA_rowmedians <- rowMedians(as.matrix(azored_putative_k2_hmp2_RA))
names(hmp2_azored_RA_rowmedians) <- rownames(azored_putative_k2_hmp2_RA)
tail(sort(hmp2_azored_RA_rowmedians), 20)

prism_azored_RA_rowmedians <- rowMedians(as.matrix(azored_putative_k2_prism_RA))
names(prism_azored_RA_rowmedians) <- rownames(azored_putative_k2_prism_RA)
tail(sort(prism_azored_RA_rowmedians), 20)

#### compiling known ridge line plots ####

cMD_known_gl <- azored_known_k2_cMD_df %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(from = -5, to = 100,quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     axis.title = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.5))) +
  scale_x_continuous(limits = c(-5,100)) +
  scale_fill_manual(values = clrs_cMD) +
  xlab("relative abundance, %") + ggtitle("Known azoreducing spp.")
cMD_known_gl

hmp2_known_gl <- azored_known_k2_hmp2_df %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(from = -5, to = 100,quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     axis.title = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.2))) +
  scale_x_continuous(limits = c(-5,100)) +
  scale_fill_manual(values = clrs_prism) +
  xlab("relative abundance, %")
hmp2_known_gl

prism_known_gl <- azored_known_k2_prism_df %>%
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(from = -5, to = 100,quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.3))) +
  scale_x_continuous(limits = c(-5,100)) +
  scale_fill_manual(values = clrs_prism) +
  xlab("relative abundance, %")
prism_known_gl

known_gl <- (cMD_known_gl / hmp2_known_gl / prism_known_gl)
known_gl

#### compiling putative ridge line plots ####

cMD_putative_gl <- azored_putative_k2_cMD_df %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.6))) +
  scale_x_continuous(limits = c(-5,100)) +
  scale_fill_manual(values = clrs_cMD) +
  xlab("relative abundance, %") + ggtitle("Putative azoreducing spp.")
cMD_putative_gl

hmp2_putative_gl <- azored_putative_k2_hmp2_df %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.3))) +
  scale_x_continuous(limits = c(-5,100)) +
  scale_fill_manual(values = clrs_hmp2) +
  xlab("relative abundance, %")
hmp2_putative_gl

prism_putative_gl <- azored_putative_k2_prism_df %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.4))) +
  scale_x_continuous(limits = c(-5,100)) +
  scale_fill_manual(values = clrs_prism) +
  xlab("relative abundance, %")
prism_putative_gl

putative_gl <- (cMD_putative_gl / hmp2_putative_gl / prism_putative_gl)
putative_gl

#### compiling A. rectale ridge line plots ####

cMD_subset_tmp <- prep_cMD("s__Agathobacter_rectale")
cMD_Arectale_gl <- cMD_subset_tmp %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     plot.title = element_text(face = "italic"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.4))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  ggtitle("A. rectale")
cMD_Arectale_gl
  
hmp2_subset_tmp <- prep_hmp2("s__Agathobacter_rectale")
hmp2_Arectale_gl <- hmp2_subset_tmp %>% 
  filter(RA > 0) %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_hmp2) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.2))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) 
hmp2_Arectale_gl

prism_subset_tmp <- prep_prism("s__Agathobacter_rectale")
prism_Arectale_gl <- prism_subset_tmp %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_prism) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.4))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  xlab("relative abundance (log(%))")
prism_Arectale_gl

Arectale_gl <- (cMD_Arectale_gl / hmp2_Arectale_gl / prism_Arectale_gl)
Arectale_gl

#### compiling Bacteroides xylanisolvens ridge line plots ####

cMD_subset_tmp <- prep_cMD("s__Bacteroides_xylanisolvens")
cMD_Bxylanisolvens_gl <- cMD_subset_tmp %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     plot.title = element_text(face = "italic"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.6))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  ggtitle("B. xylanisolvens")
cMD_Bxylanisolvens_gl

hmp2_subset_tmp <- prep_hmp2("s__Bacteroides_xylanisolvens")
hmp2_Bxylanisolvens_gl <- hmp2_subset_tmp %>% 
  filter(RA > 0) %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_hmp2) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.6))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) 
hmp2_Bxylanisolvens_gl

prism_subset_tmp <- prep_prism("s__Bacteroides_xylanisolvens")
prism_Bxylanisolvens_gl <- prism_subset_tmp %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_prism) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.3))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  xlab("log10(relative abundance, %)")
prism_Bxylanisolvens_gl

Bxylanisolvens_gl <- (cMD_Bxylanisolvens_gl / hmp2_Bxylanisolvens_gl / prism_Bxylanisolvens_gl)
Bxylanisolvens_gl

#### compiling Fusicatenibacter saccharivorans ridge line plots ####

cMD_subset_tmp <- prep_cMD("s__Fusicatenibacter_saccharivorans")
cMD_Fsaccharivorans_gl <- cMD_subset_tmp %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     plot.title = element_text(face = "italic"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.5))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  ggtitle("F. saccharivorans")
cMD_Fsaccharivorans_gl

hmp2_subset_tmp <- prep_hmp2("s__Fusicatenibacter_saccharivorans")
hmp2_Fsaccharivorans_gl <- hmp2_subset_tmp %>% 
  filter(RA > 0) %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_hmp2) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.2))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) 
hmp2_Fsaccharivorans_gl

prism_subset_tmp <- prep_prism("s__Fusicatenibacter_saccharivorans")
prism_Fsaccharivorans_gl <- prism_subset_tmp %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_prism) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.4))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  xlab("relative abundance (log(%))")
prism_Fsaccharivorans_gl

Fsaccharivorans_gl <- (cMD_Fsaccharivorans_gl / hmp2_Fsaccharivorans_gl / prism_Fsaccharivorans_gl)
Fsaccharivorans_gl

#### compiling Fusobacterium nucleatum ridge line plots ####

cMD_subset_tmp <- prep_cMD("s__Fusobacterium_nucleatum")
cMD_Fnucleatum_gl <- cMD_subset_tmp %>% 
  filter(RA > 0) %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     plot.title = element_text(face = "italic"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_cMD) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.5))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  ggtitle("F. nucleatum")
cMD_Fnucleatum_gl

hmp2_subset_tmp <- prep_hmp2("s__Fusobacterium_nucleatum")
hmp2_Fnucleatum_gl <- hmp2_subset_tmp %>% 
  filter(RA > 0) %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_hmp2) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.2))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) 
hmp2_Fnucleatum_gl

prism_subset_tmp <- prep_prism("s__Fusobacterium_nucleatum")
prism_Fnucleatum_gl <- prism_subset_tmp %>% 
  ggplot(aes(x = RA, y = diagnosis, fill = diagnosis)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 8),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text.y = element_blank(),
                     axis.title = element_blank(),
                     axis.ticks.y = element_blank()) +
  scale_fill_manual(values = clrs_prism) +
  scale_y_discrete(expand = expansion(add = c(0.1, 1.4))) +
  scale_x_continuous(trans = "log",
                     breaks = c(1e-4, 1e-2, 1, 100),
                     labels = c(-4,-2, 0, 2),
                     limits = c(1e-5,100)) +
  xlab("relative abundance (log(%))")
prism_Fnucleatum_gl

Fnucleatum_gl <- (cMD_Fnucleatum_gl / hmp2_Fnucleatum_gl / prism_Fnucleatum_gl)
Fnucleatum_gl

#### finalizing figure 4 ####

figure4_Rout <- (known_gl | putative_gl | Arectale_gl | Bxylanisolvens_gl | Fnucleatum_gl) + 
  plot_layout(widths = c(1.5, 1.5, 1, 1, 1))
figure4_Rout

print("- SAVING TO figures/figure4_Rout/figure4_Rout.[png,svg]")
ggsave("../../figures/figure4_v3/figure4_Rout.png", plot = figure4_Rout, width = 8, height = 5)
ggsave("../../figures/figure4_v3/figure4_Rout.svg", plot = figure4_Rout, width = 8, height = 5)

#######################################################################################

##### stats for manuscript #####

## number of samples and subjects explored in this work
print(paste0("- number of samples in cMD: ", length(unique(k2_pData$sampleID))))
print(paste0("- number of subjects in cMD: ", length(unique(k2_pData$subjectID))))
print(paste0("- number of samples in hmp2: ", length(unique(k2_hmp2_pData$site_sub_coll))))
print(paste0("- number of subjects in hmp2: ", length(unique(k2_hmp2_pData$Participant.ID))))
print(paste0("- number of samples in prism: ", length(unique(k2_prism_pData$Run))))
print(paste0("- number of subjects in prism: ", length(unique(k2_prism_pData$Patient))))

print(paste0("total number of samples: ", (length(unique(k2_pData$sampleID)) + 
                                             length(unique(k2_hmp2_pData$site_sub_coll)) + 
                                             length(unique(k2_prism_pData$Run)))))
print(paste0("total number of subjects: ", (length(unique(k2_pData$subjectID)) + 
                                              length(unique(k2_hmp2_pData$Participant.ID)) +
                                              length(unique(k2_prism_pData$Patient)))))

## comparing RA of known vs putative bacteria within sample individuals
### cMD
azored_known_k2_cMD_df %>% filter(diagnosis == "CRC") %>% select(RA) -> cMD_known_CRC
azored_putative_k2_cMD_df %>% filter(diagnosis == "CRC") %>% select(RA) -> cMD_putative_CRC
wilcox.test(cMD_known_CRC$RA, cMD_putative_CRC$RA)

azored_known_k2_cMD_df %>% filter(diagnosis == "adenoma") %>% select(RA) -> cMD_known_adenoma
azored_putative_k2_cMD_df %>% filter(diagnosis == "adenoma") %>% select(RA) -> cMD_putative_adenoma
wilcox.test(cMD_known_adenoma$RA, cMD_putative_adenoma$RA)

azored_known_k2_cMD_df %>% filter(diagnosis == "IBD") %>% select(RA) -> cMD_known_IBD
azored_putative_k2_cMD_df %>% filter(diagnosis == "IBD") %>% select(RA) -> cMD_putative_IBD
wilcox.test(cMD_known_IBD$RA, cMD_putative_IBD$RA)

azored_known_k2_cMD_df %>% filter(diagnosis == "control") %>% select(RA) -> cMD_known_control
azored_putative_k2_cMD_df %>% filter(diagnosis == "control") %>% select(RA) -> cMD_putative_control
wilcox.test(cMD_known_control$RA, cMD_putative_control$RA)

### HMP2
azored_known_k2_hmp2_df %>% filter(diagnosis == "CD") %>%  select(RA) -> hmp2_known_CD
azored_putative_k2_hmp2_df %>% filter(diagnosis == "CD") %>%  select(RA) -> hmp2_putative_CD
wilcox.test(hmp2_known_CD$RA, hmp2_putative_CD$RA)

azored_known_k2_hmp2_df %>% filter(diagnosis == "UC") %>%  select(RA) -> hmp2_known_UC
azored_putative_k2_hmp2_df %>% filter(diagnosis == "UC") %>%  select(RA) -> hmp2_putative_UC
wilcox.test(hmp2_known_UC$RA, hmp2_putative_UC$RA)

azored_known_k2_hmp2_df %>% filter(diagnosis == "nonIBD") %>%  select(RA) -> hmp2_known_nonIBD
azored_putative_k2_hmp2_df %>% filter(diagnosis == "nonIBD") %>%  select(RA) -> hmp2_putative_nonIBD
wilcox.test(hmp2_known_nonIBD$RA, hmp2_putative_nonIBD$RA)

### prism
azored_known_k2_prism_df %>% filter(diagnosis == "CD") %>%  select(RA) -> prism_known_CD
azored_putative_k2_prism_df %>% filter(diagnosis == "CD") %>%  select(RA) -> prism_putative_CD
wilcox.test(prism_known_CD$RA, prism_putative_CD$RA)

azored_known_k2_prism_df %>% filter(diagnosis == "UC") %>%  select(RA) -> prism_known_UC
azored_putative_k2_prism_df %>% filter(diagnosis == "UC") %>%  select(RA) -> prism_putative_UC
wilcox.test(prism_known_UC$RA, prism_putative_UC$RA)

azored_known_k2_prism_df %>% filter(diagnosis == "nonIBD") %>%  select(RA) -> prism_known_nonIBD
azored_putative_k2_prism_df %>% filter(diagnosis == "HC") %>%  select(RA) -> prism_putative_nonIBD
wilcox.test(prism_known_nonIBD$RA, prism_putative_nonIBD$RA)



azored_k2_cMD_RA[grep("aeruginosa", rownames(azored_k2_cMD_RA)),] %>% 
  colSums() %>% summary()
azored_k2_cMD_RA[grep("Pseudomonas", rownames(azored_k2_cMD_RA)),] %>% 
  colSums() %>% summary()

## number of species from "new" genera with no previous experimental validation
(search_azored_spp("g__Bacillus") %>% length() +
search_azored_spp("g__Brevibacillus") %>% length() +
search_azored_spp("g__Clostridium") %>% length() +
search_azored_spp("g__Enterococcus") %>% length() +
search_azored_spp("g__Escherichia") %>% length() +
search_azored_spp("g__Geobacillus") %>% length() +
search_azored_spp("g__Halomonas") %>% length() +
search_azored_spp("g__Lysinibacillus") %>% length() +
search_azored_spp("g__Paracoccus") %>% length() +
search_azored_spp("g__Parageobacillus") %>% length() +
search_azored_spp("g__Pigmentiphaga") %>% length() +
search_azored_spp("g__Pseudomonas") %>% length() +
search_azored_spp("g__Rhizobium") %>% length() +
search_azored_spp("g__Rhodobacter") %>% length() +
search_azored_spp("g__Rhodococcus") %>% length() +
search_azored_spp("g__Shewanella") %>% length() +
search_azored_spp("g__Shigella") %>% length() +
search_azored_spp("g__Staphylococcus") %>% length() +
search_azored_spp("g__Xenophilus") %>% length() - length(azored_taxa_hits))*-1

toc()


#### BASEMENT ####

# ## cMD
# clrs_cMD <- c("#999999", "#E69F00", "#56B4E9", "#0072B2")

# cMD_Pseudomonas <- prep_cMD("Pseudomonas")
# # cMD_Pseudomonas$diagnosis <- factor(cMD_Pseudomonas$diagnosis, levels = c("IBD", "control", "adenoma", "CRC"))
# a1 <- plot_RA(cMD_Pseudomonas, x = "RA", fill = "diagnosis", clrs = clrs_cMD) +
#   theme(legend.position = "none", 
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.title = element_text(face = "italic")) +
#   # scale_y_continuous(limits = c(0, 0.4))
#   labs(title = "Pseudonomas spp.") 
# a1
# 
# cMD_Bsubtilis <- prep_cMD("subtilis")
# b1 <- plot_RA(cMD_Bsubtilis, x = "RA", fill = "diagnosis", clrs = clrs_cMD) +
#   theme(legend.position = "none",
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         plot.title = element_text(face = "italic")) +
#   labs(title = "Bacillus subtilis") 
# b1
# 
# cMD_Ecoli <- prep_cMD("Escherichia_coli")
# c1 <- plot_RA(cMD_Ecoli, x = "RA", fill = "diagnosis", clrs = clrs_cMD) +
#   theme(legend.position = "right",
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         plot.title = element_text(face = "italic")) +
#   labs(title = "Escherichia coli") 
# # scale_y_continuous(limits = c(0, 0.4))
# c1
# 
# cMD_Pseudomonas <- prep_cMD("Pseudomonas")
# x1 <- plot_RA(cMD_Pseudomonas, x = "RA", fill = "diagnosis", clrs = clrs_cMD)
# 
# 
# ## hmp2
# clrs_hmp2 <- c("#999999", "#F0E442", "#D55E00")
# 
# hmp2_Pseudomonas <- prep_hmp2("Pseudomonas")
# a2 <- plot_RA(hmp2_Pseudomonas, x = "RA", fill = "diagnosis", clrs = clrs_hmp2) +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) 
# a2
# 
# hmp2_Bsubtilis <- prep_hmp2("subtilis")
# b2 <- plot_RA(hmp2_Bsubtilis, x = "RA", fill = "diagnosis", clrs = clrs_hmp2) +
#   theme(legend.position = "none",
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
# b2
# 
# hmp2_Ecoli <- prep_hmp2("Escherichia_coli")
# c2 <- plot_RA(hmp2_Ecoli, x = "RA", fill = "diagnosis", clrs = clrs_hmp2) +
#   theme(legend.position = "right",
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank()) #+
# # scale_y_continuous(expand = c(0,0.005))
# c2
# 
# ## prism
# clrs_prism <- c("#999999", "#F0E442", "#D55E00")
# 
# prism_Pseudomonas <- prep_prism("Pseudomonas")
# a3 <- plot_RA(prism_Pseudomonas, x = "RA", fill = "diagnosis", clrs = clrs_prism) +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) 
# a3
# 
# 
# prism_Bsubtilis <- prep_prism("subtilis")
# b3 <- plot_RA(prism_Bsubtilis, x = "RA", fill = "diagnosis", clrs = clrs_prism) +
#   theme(legend.position = "none",
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank()) 
# b3
# 
# prism_Ecoli <- prep_prism("Escherichia_coli")
# c3 <- plot_RA(prism_Ecoli, x = "RA", fill = "diagnosis", clrs = clrs_prism) +
#   theme(legend.position = "none",
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank()) #+
# # scale_y_continuous(expand = c(0,0.005))
# c3
# 
# 
# ## cib
# clrs_cib <- c("#999999", "#D55E00")
# 
# cib_Pseudomonas <- prep_cib("Pseudomonas")
# a4 <- plot_RA(cib_Pseudomonas, x = "RA", fill = "diagnosis", clrs = clrs_cib) +
#   theme(legend.position = "none",
#         axis.title.x = element_blank())
# 
# cib_Bsubtilis <- prep_cib("subtilis")
# b4 <- plot_RA(cib_Bsubtilis, x = "RA", fill = "diagnosis", clrs = clrs_cib) +
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()) +
#   xlab("relative abundance, log10(%)")
# b4
# 
# cib_Ecoli <- prep_cib("Escherichia_coli")
# c4 <- plot_RA(cib_Ecoli, x = "RA", fill = "diagnosis", clrs = clrs_cib) +
#   theme(legend.position = "none",
#         axis.title = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank()) #+
# # scale_y_continuous(expand = c(0,0.003))
# c4
#
# cMD_RA <- (a1 | b1 | c1)
# hmp2_RA <- (a2 | b2 | c2)
# prism_RA <- (a3 | b3 | c3)
# cib_RA <- (a4 | b4 | c4)