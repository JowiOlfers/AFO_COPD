library(tidyverse)
library (dplyr)
library (tibble)

###working directories
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir.Bronchial <- file.path(".", "Nasal DEG_DE in bronchi")
data.dir.B <- file.path(data.dir.Bronchial, "data")
results.dir <- file.path(data.dir.Bronchial, "results")
results.dir.img <- file.path(results.dir, "img")


############### Clinical data ####################
### BRONCHIAL
MasterTable.bronchial <- readr::read_csv(
  file.path(data.dir.B, "bronchial clinical data.csv") 
) %>% 
  dplyr:: select(1:12)

annotation.bronchial <- MasterTable.bronchial %>%
  dplyr::select(
    SAMID,
    group,
  ) %>% 
  dplyr::mutate(
    group=factor(group))

annotation.bronchial <- column_to_rownames(annotation.bronchial, var = "SAMID")
annotation.bronchial <- annotation.bronchial %>%
  dplyr::mutate(
  group = dplyr::case_when(
    group == "0" ~ "Healthy control",
    group == "1" ~ "Mild COPD",
    group == "2" ~ "Severe COPD")
  )


setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis/Comparison_DEGs")
data.dir <- file.path(".", "results")
Bronchial <- read.table(
  file.path (data.dir, "totalbronchial.normalizedcounts.txt"))

setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
Bronchial.unique <- readr::read_csv(
  file.path (results.dir, "bronchial_unique-DEGs_severeCOPD-P0.01_NNLS.csv"))


Bronchial.unique_up = Bronchial.unique[which(Bronchial.unique$logFC.x > 0), ]
readr::write_csv(Bronchial.unique, 
                 file = file.path (results.dir, "Bronchial.unique_up_P0.01.txt"))

Bronchial.unique_down = Bronchial.unique[which(Bronchial.unique$logFC.x < 0), ]
readr::write_csv(Bronchial.unique, 
                 file = file.path (results.dir, "Bronchial.unique_down_P0.01.txt"))







#select top 20
up = arrange (Bronchial.unique_up, desc(logFC.x))
up2 = head(up[,1:7], 10)
down = arrange (Bronchial.unique_down, desc(-logFC.x))
down2 = head(down[,1:7], 10)
Unique_UD <- rbind (up2, down2)
readr::write_csv(Unique_UD, file = file.path (results.dir, "bronchial_top20.up_down_uniqueDEGs.csv"))

Unique <- Bronchial [Bronchial$hgnc_symbol %in% Unique_UD$hgnc_symbol, ]
rownames(Unique) <- Unique$hgnc_symbol
Unique <- Unique [,-1]

#rownames(Bronchial) <- c()
#Bronchial <- column_to_rownames(Bronchial, var = "hgnc_symbol")



#####################heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

table (annotation.bronchial$group)
annotation_col = data.frame(Group = factor(c(rep("Healthy control", 23), rep("Mild COPD", 23), rep("Severe COPD", 122))))
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
               filename = file.path (results.dir.img, "Heatmap_bronchial-total.P0.01.png") ,dpi=600)
