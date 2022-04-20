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
  file.path(data.dir , "bronchial-mildCOPD.vs.control_FDR0.05.txt"))
df1 <- data.frame(x)

tT1 <- select(df1, hgnc_symbol, logFC, FDR)
df1_pos = tT1[which(tT1$logFC > 1),]
df1_pos_list <- df1_pos [,-2:-3]
df1_neg = tT1[which(tT1$logFC < -1),]
df1_neg_list <- df1_neg [,-2:-3]

df1 <- rbind (df1_pos, df1_neg)
df1_list <- df1 [,-2:-3]

##P<0.001
#x <-readr::read_csv(
#  file.path(data.dir, "bronchial-mildCOPD.vs.control_P0.001.txt"))
#df1b <- data.frame(x)
#
#tT1b <- select(df1b, hgnc_symbol, logFC, FDR)
#df1b_pos = tT1b[which(tT1b$logFC > 0),]
#df1b_pos_list <- df1b_pos [,-2:-3]
#df1b_neg = tT1b[which(tT1b$logFC < 0),]
#df1b_neg_list <- df1b_neg [,-2:-3]
#
#df1b <- rbind (df1b_pos, df1b_neg)
#df1b_list <- df1b [,-2:-3]


### Severe COPD vs. controls 
x <-readr::read_csv(
  file.path(data.dir, "bronchial-severeCOPD.vs.control_FDR0.05.txt"))
table2 <- data.frame(x)

tT2 <- select(table2, hgnc_symbol, logFC, FDR)
df2_pos = tT2[which(tT2$logFC > 1),]
df2_pos_list <- df2_pos [,-2:-3]
df2_neg = tT2[which(tT2$logFC < -1),]
df2_neg_list <- df2_neg [,-2:-3]

df2 <- rbind (df2_pos, df2_neg)
df2_list <- df2 [,-2:-3]


### Severe vs. mild COPD 
x <-readr::read_csv(
  file.path(data.dir, "bronchial-severeCOPD.vs.mild_FDR0.05.txt"))
table3 <- data.frame(x)

tT3 <- select(table3, hgnc_symbol, logFC, FDR)
df3_pos = tT3[which(tT3$logFC > 1),]
df3_pos_list <- df3_pos [,-2:-3]
df3_neg = tT3[which(tT3$logFC < -1),]
df3_neg_list <- df3_neg [,-2:-3]

df3 <- rbind (df3_pos, df3_neg)
df3_list <- df3 [,-2:-3]



############ Overlap and selection ############ 

###Mild COPD
overlap_mild_POS <-  merge (df1_pos, df2_pos, by= "hgnc_symbol")
uniqueG_mild_POS <- df1_pos [ ! df1_pos$hgnc_symbol %in% overlap_mild_POS$hgnc_symbol,]
uniqueGenes_mild_POS <-  merge (uniqueG_mild_POS, df3_pos, by= "hgnc_symbol")
readr::write_csv(uniqueG_mild_POS, file = file.path (results.dir, "bronchial_UP-unique-DEGs_mildCOPD-FDR0.05FC2.csv"))

overlap_mild_NEG <-  merge (df1_neg, df2_neg, by= "hgnc_symbol")
uniqueG_mild_NEG <- df1_neg [ ! df1_neg$hgnc_symbol %in% overlap_mild_NEG$hgnc_symbol,]
uniqueGenes_mild_NEG <-  merge (uniqueG_mild_NEG, df3_pos, by= "hgnc_symbol")
readr::write_csv(uniqueG_mild_NEG, file = file.path (results.dir, "bronchial_NEG-unique-DEGs_mildCOPD-FDR0.05FC2.csv"))

uniqueGenes_mild <- rbind (uniqueGenes_mild_POS, uniqueGenes_mild_NEG)
readr::write_csv(uniqueGenes_mild, file = file.path (results.dir, "bronchial_allunique-DEGs_mildCOPD-FDR0.05FC2.csv"))


###Severe COPD
overlap_severe_POS <-  merge (df1_pos, df2_pos, by= "hgnc_symbol")
uniqueG_severe_POS <- df2_pos [ ! df2_pos$hgnc_symbol %in% overlap_severe_POS$hgnc_symbol,]
uniqueGenes_severe_POS <-  merge (uniqueG_severe_POS, df3_pos, by= "hgnc_symbol")
readr::write_csv(uniqueG_severe_POS, file = file.path (results.dir, "bronchial_UP-unique-DEGs_severeCOPD-FDR0.05FC2.csv"))

overlap_severe_NEG <-  merge (df1_neg, df2_neg, by= "hgnc_symbol")
uniqueG_severe_NEG <- df2_neg [ ! df2_neg$hgnc_symbol %in% overlap_severe_NEG$hgnc_symbol,]
uniqueGenes_severe_NEG <-  merge (uniqueG_severe_NEG, df3_pos, by= "hgnc_symbol")
readr::write_csv(uniqueG_severe_NEG, file = file.path (results.dir, "bronchial_NEG-unique-DEGs_severeCOPD-FDR0.05FC2.csv"))

uniqueGenes_severe <- rbind (uniqueGenes_severe_POS, uniqueGenes_severe_NEG)
readr::write_csv(uniqueGenes_severe, file = file.path (results.dir, "bronchial_allunique-DEGs_severeCOPD-FDR0.05FC2.csv"))




#############venn diagrams
library (VennDiagram)

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
  main = "Differentially expressed genes in COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#FFDF00",0.3), alpha('#3090C7',0.3), alpha('#3090C7',0.3)),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(195, -205, 100),
  cat.dist = c(0.030, 0.040, 0.040),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "COPD_venndiagram.png")
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
  fill = c(alpha("#FFDF00",0.3), alpha('#3090C7',0.3), alpha('#3090C7',0.3)),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(195, -205, 100),
  cat.dist = c(0.030, 0.040, 0.040),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Up-regulated_COPD_venndiagram.png")
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
  fill = c(alpha("#FFDF00",0.3), alpha('#3090C7',0.3), alpha('#3090C7',0.3)),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(195, -205, 100),
  cat.dist = c(0.030, 0.040, 0.040),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Down-regulated_COPD_venndiagram.png")
)

