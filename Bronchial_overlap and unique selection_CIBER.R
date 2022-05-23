library (ggplot2)
library(ggrepel)
library (dplyr)

#Directories 
setwd("~/HVHL/Stages/Afstudeerstage_2022/DE analysis/Bronchial data")
data.dir <- file.path(".", "results")
results.dir <- file.path(".", "unique selection")
results.dir.img <- file.path(results.dir, "img")


### Mild COPD vs. controls 
#FDR 0.5
x <-readr::read_csv(
  file.path(data.dir , "bronchial_all-DEGs_mild.vs.controls-CIBER-FDR0.05FC2.csv"))
table1 <- data.frame(x)

df1 <- select(table1, hgnc_symbol, logFC, FDR)
df1_list <- df1 [,-2:-3]

df1_pos = df1[which(df1$logFC >= 1),]
df1_pos_list <- df1_pos [,-2:-3]
df1_neg = df1[which(df1$logFC <= -1),]
df1_neg_list <- df1_neg [,-2:-3]


### Severe COPD vs. controls 
x <-readr::read_csv(
  file.path(data.dir, "bronchial_all-DEGs_severe.vs.controls-CIBER-FDR0.05FC2.csv"))
table2 <- data.frame(x)

df2 <- select(table2, hgnc_symbol, logFC, FDR)
df2_list <- df2 [,-2:-3]

df2_pos = df2[which(df2$logFC >= 1),]
df2_pos_list <- df2_pos [,-2:-3]
df2_neg = df2[which(df2$logFC <= -1),]
df2_neg_list <- df2_neg [,-2:-3]



### Severe vs. mild COPD 
x <-readr::read_csv(
  file.path(data.dir, "bronchial_all-DEGs_severevsmild-FDR0.05FC2_ciber.csv"))
table3 <- data.frame(x)

df3 <- select(table3, hgnc_symbol, logFC, FDR)
df3_list <- df3 [,-2:-3]

df3_pos = df3[which(df3$logFC >= 1),]
df3_pos_list <- df3_pos [,-2:-3]
df3_neg = df3[which(df3$logFC <= -1),]
df3_neg_list <- df3_neg [,-2:-3]



### Mild vs. severe COPD
#x <-readr::read_csv(
#  file.path(data.dir, "bronchial_all-DEGs_mildvssevere-FDR0.05FC2_ciber.csv"))
#table4 <- data.frame(x)

#df4 <- select(table4, hgnc_symbol, logFC, FDR)
#df4_list <- df4 [,-2:-3]

#df4_pos = df4[which(df4$logFC > 1),]
#df4_pos_list <- df4_pos [,-2:-3]
#df4_neg = df4[which(df4$logFC < -1),]
#df4_neg_list <- df4_neg [,-2:-3]





###########################################################
## Overlap and selection  
###########################################################

overlap_mc.sc <- merge (df1, df2, by = "hgnc_symbol" )
#overlap_mc.sc_POS <-  merge (df1_pos, df2_pos, by= "hgnc_symbol")
#overlap_mc.sc_NEG <-  merge (df1_neg, df2_neg, by= "hgnc_symbol")


####### UNIQUE GENES FOR SEVERE COPD ##########
df2.mc.rm <- df2 [ ! df2$hgnc_symbol %in% overlap_mc.sc$hgnc_symbol,]
unique_severe <- merge (df2.mc.rm, df3, by = "hgnc_symbol" )
dim (unique_severe)
readr::write_csv(unique_severe, file = file.path (results.dir, "bronchial_unique-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))

#LogFC.x is the fold change from severe vs. control
unique_severe_pos <- unique_severe [which(unique_severe$logFC.x > 1), ]
readr::write_csv(unique_severe_pos, file = file.path (results.dir, "bronchial_unique-UP-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))

unique_severe_neg <- unique_severe [which(unique_severe$logFC.x < -1), ]
readr::write_csv(unique_severe_neg, file = file.path (results.dir, "bronchial_unique-NEG-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))


####### UNIQUE GENES FOR MILD COPD ##########
#df1.mc.rm <- df1 [ ! df1$hgnc_symbol %in% overlap_mc.sc$hgnc_symbol,]
#unique_mild <- merge (df1.mc.rm, df4, by = "hgnc_symbol" )
#dim (unique_mild)
#readr::write_csv(unique_mild, file = file.path (results.dir, "mild_unique-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))

#unique_mild_pos <- unique_mild [which(unique_mild$logFC > 1), ]
#readr::write_csv(unique_mild_pos, file = file.path (results.dir, "mild_unique-UP-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))

#unique_mild_neg <- unique_mild [which(unique_mild$logFC < -1), ]
#readr::write_csv(unique_mild_neg, file = file.path (results.dir, "mild_unique-NEG-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))



############# VennDiagrams
library (VennDiagram)

#######SEVERE COPD
list <- list (df1_list, df2_list, df3_list)
list_up <- list (df1_pos_list, df2_pos_list, df3_pos_list)
list_down <- list (df1_neg_list, df2_neg_list, df3_neg_list)

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
  fill = c(alpha("#C9F5FD"), alpha('#C9FDE5'),alpha('#C9DEFD')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 180),
  cat.dist = c(0.030, 0.040, 0.040),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "severeCOPD_CIBER_venndiagram.png")
)


venn.diagram(
  x=list_up, 
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
  fill = c(alpha("#C9F5FD"), alpha('#C9FDE5'),alpha('#C9DEFD')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 180),
  cat.dist = c(0.030, 0.040, 0.040),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Up-regulated_severeCOPD_CIBER_venndiagram.png")
)


venn.diagram(
  x=list_down, 
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
  fill = c(alpha("#C9F5FD"), alpha('#C9FDE5'),alpha('#C9DEFD')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 180),
  cat.dist = c(0.030, 0.040, 0.040),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Down-regulated_severeCOPD_CIBER_venndiagram.png")
)

#######MILD COPD
#listM <- list (df1_list, df2_list, df4_list)
#listM_up <- list (df1_pos_list, df2_pos_list, df4_pos_list)
#listM_down <- list (df1_neg_list, df2_neg_list, df4_neg_list)

#venn.diagram(
#  x=listM, 
#  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Mild vs. severe COPD"),
#  height = 600,
#  width = 600,
#  resolution = 300,
#  compression = "lzw",
#  lwd = 1.2,
#  main = "Differentially expressed genes in mild COPD", 
#  main.fontface = "bold",
#  main.fontfamily = "sans", 
#  main.cex = .4,
#  fill = c(alpha("#FBE2F1"), alpha('#C9FDE5'),alpha('#C9DEFD')),
#  cex = .4,
#  fontface = "bold", 
#  fontfamily = "sans",
#  cat.cex = 0.35,
#  cat.fontface = "bold",
#  cat.default.pos = "outer",
#  cat.pos = c(0, 0, 180),
#  cat.dist = c(0.030, 0.040, 0.040),
#  cat.fontfamily = "sans",
#  filename =  file.path (results.dir.img, "MildCOPD_CIBER_venndiagram.png")
#)


#venn.diagram(
#  x=listM_up, 
#  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Mild vs. severe COPD"),
#  height = 600,
#  width = 600,
#  resolution = 300,
#  compression = "lzw",
#  lwd = 1.2,
#  main = "Up-regulated DEGs in mild COPD", 
#  main.fontface = "bold",
#  main.fontfamily = "sans", 
#  main.cex = .4,
#  fill = c(alpha("#FBE2F1"), alpha('#C9FDE5'),alpha('#C9DEFD')),
#  cex = .4,
#  fontface = "bold", 
#  fontfamily = "sans",
#  cat.cex = 0.35,
#  cat.fontface = "bold",
#  cat.default.pos = "outer",
#  cat.pos = c(0, 0, 180),
#  cat.dist = c(0.030, 0.040, 0.040),
#  cat.fontfamily = "sans",
#  filename =  file.path (results.dir.img, "Up-regulated_mildCOPD_CIBER_venndiagram.png")
#)


#venn.diagram(
#  x=listM_down, 
#  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Mild vs. severe COPD"),
#  height = 600,
#  width = 600,
#  resolution = 300,
#  compression = "lzw",
#  lwd = 1.2,
#  main = "Down-regulated DEGs in mild COPD", 
#  main.fontface = "bold",
#  main.fontfamily = "sans", 
#  main.cex = .4,
#  fill = c(alpha("#FBE2F1"), alpha('#C9FDE5'),alpha('#C9DEFD')),
#  cex = .4,
#  fontface = "bold", 
#  fontfamily = "sans",
#  cat.cex = 0.35,
#  cat.fontface = "bold",
#  cat.default.pos = "outer",
#  cat.pos = c(0, 0, 180),
#  cat.dist = c(0.030, 0.040, 0.040),
#  cat.fontfamily = "sans",
#  filename =  file.path (results.dir.img, "Down-regulated_mildCOPD_CIBER_venndiagram.png")
#)


