library (ggplot2)
library(ggrepel)
library (dplyr)
library(tidyverse)

#Directories 
setwd("~/HVHL/Stages/Afstudeerstage_2022/DE analysis/Nasal DEG_DE in bronchi")
data.dir <- file.path(".", "results")
results.dir <- file.path(".", "results")
results.dir.img <- file.path(results.dir, "img")



### Mild COPD vs. controls 
x <-readr::read_csv(
  file.path(data.dir, "DEG.bronchial_all-DEGs_mild.vs.control-P0.01_NNLS.csv"))
table1 <- data.frame(x)

df1 <- select(table1, hgnc_symbol, logFC, FDR, PValue)
df1_list <- df1 [,-2:-4]

df1_pos = df1[which(df1$logFC > 0),]
df1_pos_list <- df1_pos [,-2:-4]
df1_neg = df1[which(df1$logFC < 0),]
df1_neg_list <- df1_neg [,-2:-4]



### Severe COPD vs. controls 
x <-readr::read_csv(
  file.path(data.dir, "DEG.bronchial_all-DEGs_severe.vs.controls-P0.01_NNLS.csv"))
table2 <- data.frame(x)

df2 <- select(table2, hgnc_symbol, logFC, FDR, PValue)
df2_list <- df2 [,-2:-4]

df2_pos = df2[which(df2$logFC > 0),]
df2_pos_list <- df2_pos [,-2:-4]
df2_neg = df2[which(df2$logFC < 0),]
df2_neg_list <- df2_neg [,-2:-4]


### Severe vs. mild COPD 
x <-readr::read_csv(
  file.path(data.dir, "DEG.bronchial_all-DEGs_severe.vs.mild-P0.01_NNLS.csv"))
table3 <- data.frame(x)

df3 <- select(table3, hgnc_symbol, logFC, FDR, PValue)
df3_list <- df3 [,-2:-4]

df3_pos = df3[which(df3$logFC > 0),]
df3_pos_list <- df3_pos [,-2:-4]
df3_neg = df3[which(df3$logFC < 0),]
df3_neg_list <- df3_neg [,-2:-4]




######################overlap - specific for severe COPD
overlap_mc.sc <- merge (df1, df2, by = "hgnc_symbol" )
#overlap_mc.sc_POS <-  merge (df1_pos, df2_pos, by= "hgnc_symbol")
#overlap_mc.sc_NEG <-  merge (df1_neg, df2_neg, by= "hgnc_symbol")


####### UNIQUE GENES FOR SEVERE COPD ##########
df2.mc.rm <- df2 [ ! df2$hgnc_symbol %in% overlap_mc.sc$hgnc_symbol,]
unique_severe <- merge (df2.mc.rm, df3, by = "hgnc_symbol" )
dim (unique_severe)
readr::write_csv(unique_severe, file = file.path (results.dir, "bronchial_unique-DEGs_severeCOPD-P0.01_NNLS.csv"))

#unique_severe_pos <- unique_severe [which(unique_severe$logFC > 0), ]
#readr::write_csv(unique_severe_pos, file = file.path (results.dir, "bronchial_unique-UP-DEGs_severeCOPD-P0.01_NNLS.csv"))

#unique_severe_neg <- unique_severe [which(unique_severe$logFC < 0), ]
#readr::write_csv(unique_severe_neg, file = file.path (results.dir, "bronchial_unique-NEG-DEGs_severeCOPD-P0.01_NNLS.csv"))


#############venn diagrams
library (VennDiagram)

list <- list (df1_list, df2_list, df3_list)
list.up <- list (df1_pos_list, df2_pos_list, df3_pos_list)
list.down <- list (df1_neg_list, df2_neg_list, df3_neg_list)

# Mild COPD vs. control = C9F5FD
# Severe COPD vs. control = C9FDE5
# Severe vs. mild COPD = C9DEFD

venn.diagram(
  x=list, 
  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Severe vs. mild COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Differentially expressed genes in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#E0E0E0"), alpha('#C9EFFD'),alpha('#EADFFC')),
  cex = .35,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-25, 25, 180),
  cat.dist = c(0.050, 0.050, 0.050),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "SevereCOPD_venndiagram.nnls.png")
)


venn.diagram(
  x=list.up, 
  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Severe vs. mild COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Up-regulated DEGs in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#E0E0E0"), alpha('#C9EFFD'),alpha('#EADFFC')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(200, 180, 180),
  cat.dist = c(0.050, 0.070, 0.050),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Up-regulated_severeCOPD_venndiagram.nnls.png")
)


venn.diagram(
  x=list.down, 
  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Severe vs. mild COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Down-regulated DEGs in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#E0E0E0"), alpha('#C9EFFD'),alpha('#EADFFC')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(180, 0, 180),
  cat.dist = c(0.040, 0.050, 0.050),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Down-regulated_severeCOPD_venndiagram.nnls.png")
)


