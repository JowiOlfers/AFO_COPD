library(tidyverse)
library (dplyr)
library (tibble)

###working directories
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir.Nasal <- file.path(".", "Nasal data")
data.dir.N <- file.path(data.dir.Nasal, "data")


############### Clinical data ####################
### NASAL
MasterTable.nasal <- readr::read_csv(
  file.path(data.dir.N, "nasal clinical data.csv") 
) %>% 
  dplyr:: select(1:12)

annotation.nasal <- MasterTable.nasal %>%
  dplyr::select(
    SAMID,
    group,
  ) %>% 
  dplyr::mutate(
    group=factor(group))

annotation.nasal <- column_to_rownames(annotation.nasal, var = "SAMID")
annotation.nasal <- annotation.nasal %>%
  dplyr::mutate(
  group = dplyr::case_when(
    group == "0" ~ "Healthy control",
    group == "1" ~ "Mild COPD",
    group == "2" ~ "Severe COPD")
  )


setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis/Comparison_DEGs")
data.dir <- file.path(".", "results")
Nasal <- read.table(
  file.path (data.dir, "Totalnasal.normalizedcounts.txt"))

setwd("~/HVHL/Stages/Afstudeerstage_2022/DE analysis/Nasal data/unique selection")
results.dir <- file.path(".")
results.dir.img <- file.path(".", "img")


Nasal.unique <- readr::read_csv(
  file.path (results.dir, "nasal_unique-DEGs_severeCOPD-NNLS-FDR0.05FC2.csv"))

Nasal.unique_up = Nasal.unique[which(Nasal.unique$logFC.x >= 1), ]
readr::write_csv(Nasal.unique_up, 
                 file = file.path (results.dir, "Nasal.unique_up_FDR0.05FC2.txt"))

Nasal.unique_down = Nasal.unique[which(Nasal.unique$logFC.x <= -1), ]
readr::write_csv(Nasal.unique_down, 
                 file = file.path (results.dir, "Nasal.unique_down_FDR0.05FC2.txt"))


#select top 20
up = arrange (Nasal.unique_up, desc(logFC.x))
up2 = head(up[,1:7], 10)
down = arrange (Nasal.unique_down, desc(-logFC.x))
down2 = head(down[,1:7], 10)
Unique_UD <- rbind (up2, down2)
readr::write_csv(Unique_UD, file = file.path (results.dir, "nasal_top20.up_down_uniqueDEGs.csv"))

Unique <- Nasal [Nasal$hgnc_symbol %in% Unique_UD$hgnc_symbol, ]
rownames(Unique) <- Unique$hgnc_symbol
Unique <- Unique [,-1]

#rownames(Nasal) <- c()
#Nasal <- column_to_rownames(Nasal, var = "hgnc_symbol")



#####################heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

table (annotation.nasal$group)
annotation_col = data.frame(Group = factor(c(rep("Healthy control", 22), rep("Mild COPD", 24), rep("Severe COPD", 113))))
rownames(annotation_col)

data <- Unique
colnames(data)
rownames(annotation_col) <- colnames(data)
head(annotation_col)

data1 <- data [,annotation_col$Group == "Healthy control"]
data2 <- data [,annotation_col$Group == "Mild COPD"]
data3 <- data [,annotation_col$Group == "Severe COPD"]

HC1 <- hclust (dist(t(data1), method = "binary"), method = "complete", members = NULL)
HC2 <- hclust (dist(t(data2), method = "binary"), method = "complete", members = NULL)
HC3 <- hclust (dist(t(data3), method = "binary"), method = "complete", members = NULL)

data1 <- data1[,HC1$order]
data2 <- data2[,HC2$order]
data3 <- data3[,HC3$order]
datax <- cbind(data1, data2, data3)

rownames (annotation_col) <- colnames (datax)

ph <- pheatmap(datax,
               scale="row",
               annotation_col = annotation_col,
               annotation_legend = T,
               annotation_names_col = F,
               number_format="%.2e",
               border="white",  
               color=colorRampPalette(c("#013E80", "white", "orangered"))(100),
               border_color=NA, 
               cellwidth = 9,cellheight = 20, 
               cluster_cols = F, 
               cluster_rows = T,
               show_rownames = T, 
               show_colnames = F,
               legend = T,   
               legend_breaks = -5:5, 
               fontsize = 12,
               fontsize_row = 12, 
               fontsize_col = 10,
               filename = file.path (results.dir.img, "Heatmap_nasal_severe.vs.control_FDR0.05.png") ,dpi=600)
