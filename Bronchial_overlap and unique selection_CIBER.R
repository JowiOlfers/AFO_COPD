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
  file.path(data.dir , "bronchial-mildCOPD.vs.control_CIBER_FDR0.05.txt"))
df1 <- data.frame(x)

tT1 <- select(df1, hgnc_symbol, logFC, FDR)
df1_pos = tT1[which(tT1$logFC > 1),]
df1_pos_list <- df1_pos [,-2:-3]
df1_neg = tT1[which(tT1$logFC < -1),]
df1_neg_list <- df1_neg [,-2:-3]

df1 <- rbind (df1_pos, df1_neg)
df1_list <- df1 [,-2:-3]


### Severe COPD vs. controls 
x <-readr::read_csv(
  file.path(data.dir, "bronchial-severeCOPD.vs.control_CIBER_FDR0.05.txt"))
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


### Mild vs. severe COPD
x <-readr::read_csv(
  file.path(data.dir, "bronchial-mildCOPD.vs.severe_FDR0.05.txt"))
table4 <- data.frame(x)

tT4 <- select(table4, hgnc_symbol, logFC, FDR)
df4_pos = tT4[which(tT4$logFC > 1),]
df4_pos_list <- df4_pos [,-2:-3]
df4_neg = tT4[which(tT4$logFC < -1),]
df4_neg_list <- df4_neg [,-2:-3]

df4 <- rbind (df4_pos, df4_neg)
df4_list <- df4 [,-2:-3]




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

unique_severe_pos <- unique_severe [which(unique_severe$logFC > 1), ]
readr::write_csv(unique_severe_pos, file = file.path (results.dir, "bronchial_unique-UP-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))

unique_severe_neg <- unique_severe [which(unique_severe$logFC < -1), ]
readr::write_csv(unique_severe_neg, file = file.path (results.dir, "bronchial_unique-NEG-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))


####### UNIQUE GENES FOR MILD COPD ##########
df1.mc.rm <- df1 [ ! df1$hgnc_symbol %in% overlap_mc.sc$hgnc_symbol,]
unique_mild <- merge (df1.mc.rm, df4, by = "hgnc_symbol" )
dim (unique_mild)
readr::write_csv(unique_mild, file = file.path (results.dir, "mild_unique-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))

unique_mild_pos <- unique_mild [which(unique_mild$logFC > 1), ]
readr::write_csv(unique_mild_pos, file = file.path (results.dir, "mild_unique-UP-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))

unique_mild_neg <- unique_mild [which(unique_mild$logFC < -1), ]
readr::write_csv(unique_mild_neg, file = file.path (results.dir, "mild_unique-NEG-DEGs_severeCOPD-CIBER-FDR0.05FC2.csv"))



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
  filename =  file.path (results.dir.img, "Down-regulated_severeCOPD_CIBER_venndiagram.png")
)

#######MILD COPD
listM <- list (df1_list, df2_list, df4_list)
listM_up <- list (df1_pos_list, df2_pos_list, df4_pos_list)
listM_down <- list (df1_neg_list, df2_neg_list, df4_neg_list)

venn.diagram(
  x=listM, 
  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Mild vs. severe COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Differentially expressed genes in mild COPD", 
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
  filename =  file.path (results.dir.img, "MildCOPD_CIBER_venndiagram.png")
)


venn.diagram(
  x=listM_up, 
  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Mild vs. severe COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Up-regulated DEGs in mild COPD", 
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
  filename =  file.path (results.dir.img, "Up-regulated_mildCOPD_CIBER_venndiagram.png")
)


venn.diagram(
  x=listM_down, 
  category.names = c("Mild COPD vs. control", "Severe COPD vs. control", "Mild vs. severe COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Down-regulated DEGs in mild COPD", 
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
  filename =  file.path (results.dir.img, "Down-regulated_mildCOPD_CIBER_venndiagram.png")
)


