# ============================ gene_hits_hmm_table.R ============================ #
#' description: script for generating a table showing all of the gene hits to 
#' UHGG genomes. going into a table in the manuscript
#' 
#' NOTE: it is assumed that all necessary packages and data are loaded into the 
#' R session where this script is being run.
#' 
#' Author: Domenick J. Braccia
#' Last updated: 29 October 2021
# =========================================================================== #

print("- reading in data files") 

## code below taken from original upset.R script
azored_genes <- read.table("../../results/from-GutFunFind/Azoreductase.genes_hmm.tsv",
                        header = TRUE)
colnames(azored_genes) <- c("genome", "genes")
gene_names <- sort(unique(unlist(unique(strsplit(azored_genes$genes, ',')))))
azored_genes <- mutate(azored_genes, absent = 0, arsH = 0, azo1 = 0, azoR_ropacus = 0, 
                       azoR1_cperf = 0, azoR1_llentus = 0, cladeI = 0, cladeII = 0, 
                       cladeIII = 0, cladeIVa = 0, cladeIVb = 0, ferB = 0, mdaB = 0, yieF = 0)


for (i in 1:dim(azored_genes)[1]) {
  current_genes <- strsplit(azored_genes[i, "genes"], ",")
  gene_table <- table(current_genes)
  for (j in 1:length(gene_table)) {
    if ("absent" %in% names(current_genes)) {break}
    azored_genes[i, names(gene_table)[j]] <- gene_table[j]
  }
}
# select(azored_genes, -c(absent)) # removing "absent" column

