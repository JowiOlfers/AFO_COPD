library (ggplot2)
library(ggrepel)
library (dplyr)

#Directories 
setwd("~/HVHL/Stages/Afstudeerstage_2022/DE analysis/Nasal data")
data.dir <- file.path(".", "results")
results.dir <- file.path(".", "unique selection")
results.dir.img <- file.path(results.dir, "img")

### Mild COPD vs. control
#FDR 0.5
#x <-readr::read_csv(
#  file.path(data.dir, "nasal-mildCOPD.vs.control_FDR0.5.txt"))
#df1a <- data.frame(x)

#tT1a <- select(df1a, hgnc_symbol, logCPM, logFC, FDR)
#df1a_pos = tT1a[which(tT1a$logFC > 0),]
#df1a_pos_list <- df1a_pos [,-2:-4]
#df1a_neg = tT1a[which(tT1a$logFC < 0),]
#df1a_neg_list <- df1a_neg [,-2:-4]

#df1a <- rbind (df1a_pos, df1a_neg)
#df1a_list <- df1a [,-2:-4]



### Severe COPD vs. controls 
x <-readr::read_csv(
  file.path(data.dir, "nasal_all-DEGs_severe.vs.controls-FDR0.05FC2.csv"))
table2 <- data.frame(x)

df2 <- select(table2, hgnc_symbol, logCPM, logFC, FDR)
df2_list <- df2 [,-2:-4]

df2_pos = df2[which(df2$logFC >= 1),]
df2_pos_list <- df2_pos [,-2:-4]
df2_neg = df2[which(df2$logFC <= -1),]
df2_neg_list <- df2_neg [,-2:-4]


### Severe vs. mild COPD 
x <-readr::read_csv(
  file.path(data.dir, "nasal_all-DEGs_severe.vs.mildCOPD-FDR0.05FC2.csv"))
table3 <- data.frame(x)

df3 <- select(table3, hgnc_symbol, logCPM, logFC, FDR)
df3_list <- df3 [,-2:-4]

df3_pos = df3[which(df3$logFC >= 1),]
df3_pos_list <- df3_pos [,-2:-4]
df3_neg = df3[which(df3$logFC <= -1),]
df3_neg_list <- df3_neg [,-2:-4]




######################overlap and excluding
#nothing with mild COPD since nothing significant differential expressed

##SEVERE COPD
overlap <- merge (df2, df3, by= "hgnc_symbol")
overlap_pos <-  merge (df2_pos, df3_pos, by= "hgnc_symbol")
overlap_neg <-  merge (df2_neg, df3_neg, by= "hgnc_symbol")

unique_severe <- overlap
dim (unique_severe)
readr::write_csv(unique_severe, file = file.path (results.dir, "nasal_unique-DEGs_severeCOPD-FDR0.05FC2.csv"))

unique_severe_pos <- overlap_pos
dim (unique_severe_pos)
readr::write_csv(unique_severe_pos, file = file.path (results.dir, "nasal_unique-up-DEGs_severeCOPD-FDR0.05FC2.csv"))

unique_severe_neg <- overlap_neg
dim (unique_severe_neg)
readr::write_csv(unique_severe_neg, file = file.path (results.dir, "nasal_unique-down-DEGs_severeCOPD-FDR0.05FC2.csv"))



#############venn diagrams
library (VennDiagram)

list <- list (df2_list, df3_list)
list.up <- list (df2_pos_list, df3_pos_list)
list.down <- list (df2_neg_list, df3_neg_list)

venn.diagram(
  x=list, 
  category.names = c("Severe COPD vs. control", "Severe vs. mild COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Differentially expressed genes in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#CCEAF7"), alpha('#E3D0FC')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(185, -205),
  cat.dist = c(0.03, 0.045),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "SevereCOPD_venndiagram.png")
)



venn.diagram(
  x=list.up, 
  category.names = c("Severe COPD vs. control", "Severe vs. mild COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Up-regulated DEGs in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#CCEAF7"), alpha('#E3D0FC')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(185, -205),
  cat.dist = c(0.030, 0.045),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Up-regulated_severeCOPD_venndiagram.png")
)


venn.diagram(
  x=list.down, 
  category.names = c("Severe COPD vs. control", "Severe vs. mild COPD"),
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Down-regulated DEGs in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#CCEAF7"), alpha('#E3D0FC')),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.35,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(185, -205),
  cat.dist = c(0.030, 0.045),
  cat.fontfamily = "sans",
  filename =  file.path (results.dir.img, "Down-regulated_severeCOPD_venndiagram.png")
)

