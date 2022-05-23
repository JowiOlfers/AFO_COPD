############ ENRICHMENT ANALYSIS #############

# GSEA online " .... "
library(tidyverse)
library (dplyr)
library (tibble)

###working directories
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir_com <- file.path(".", "Comparison_DEGs")
data.dir <- file.path(data.dir_com, "results")
data.dir_DEG.B <- file.path(".", "Nasal DEG_DE in bronchi")
data.dir_reflect <- file.path(data.dir_DEG.B, "results")
data.dir.Nasal <- file.path(".", "Nasal data")
data.dir.NU <- file.path(data.dir.Nasal, "unique selection")

results.dir <- file.path(".", "Enrichment analysis")



#### load SAMID group annotation
annotation_nasal <- read.table(
  file.path(data.dir, "annotation_nasal.txt"))
annotation_nasal <- tibble::rownames_to_column(annotation_nasal, "SAMID")
annotation_nasal.SC <- annotation_nasal [!(annotation_nasal$group == "1"),]

annotation_bronchial <- read.table(
  file.path(data.dir, "annotation_bronchial.txt"))
annotation_bronchial <- tibble::rownames_to_column(annotation_bronchial, "SAMID")
annotation_bronchial.SC <- annotation_bronchial [!(annotation_bronchial$group == "1"),]
bronchial.SC_IDlist <- annotation_bronchial.SC [,-2]



#### load normalized counts for both nasal and bronchial with all groups
Total.nasal_counts <- read.table (
  file.path(data.dir, "Totalnasal.normalizedcounts.txt"))

Total.bronchial_counts <- read.table(
  file.path(data.dir, "totalbronchial.normalizedcounts.txt"))



######################################################
#Load unique DEG in severe COPD in nasal and bronchi 
######################################################

#Nasal DEG (527 list)
NB_unique.NNLS <- readr::read_csv(
  file.path(data.dir.NU, "nasal_unique-DEGs_severeCOPD-NNLS-FDR0.05FC2.csv"))
NB_unique.NNLS <- NB_unique.NNLS [order(NB_unique.NNLS$logFC.x, decreasing = TRUE), ]


#Nasal DEG in bronchi (104 list)
BB_unique.NNLS <- readr::read_csv(
  file.path(data.dir_reflect, "bronchial_unique-DEGs_severeCOPD-P0.01_NNLS.csv"))
BB_unique.NNLS <- BB_unique.NNLS [order(BB_unique.NNLS$logFC.x, decreasing = TRUE), ]





####### select only genes of the list of interest
bronchial_counts <- Total.bronchial_counts [Total.bronchial_counts$hgnc_symbol %in% BB_unique.NNLS$hgnc_symbol, ]
#readr::write_csv(
#  bronchial_counts, file = file.path (results.dir, "normalizedcounts.bronchial_selection.csv"))
write.table(
  bronchial_counts, file = file.path (results.dir, "normalizedcounts.bronchial_selection.txt"))

#rownames(bronchial_counts) <- bronchial_counts[,1]
#bronchial_counts <- bronchial_counts [,-1]


nasal_counts <- Total.nasal_counts [Total.nasal_counts$hgnc_symbol %in% NB_unique.NNLS$hgnc_symbol, ]
write.table(
  nasal_counts, file = file.path (results.dir, "normalizedcounts.nasalDEG.txt"))
readr::write_csv(
  nasal_counts, file = file.path (results.dir, "normalizedcounts.nasalDEG.csv"))

#rownames(nasal_counts) <- nasal_counts[,1]
#nasal_counts <- nasal_counts [,-1]


