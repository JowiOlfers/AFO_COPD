#First set working directory (go to session)

library (ggplot2)
library(ggrepel)
library (dplyr)

#nasal_mild.DEG = P<0.001 and FC > 0 or < 0 
#nasal_severe.DEG = FDR < 0.01 and FC > 2 or < -2
#bronchial_mild.DEG = FDR < 0.01 and FC > 2 or < -2
#bronchial_severe.DEG = FDR < 0.01 and FC > 2 or < -2


############

#nasal mild vs. control
nasal_mild.DEG <- read.table("nasal_allDEGs_mild.vs.controls-P0.001FC.csv",sep = "",header=T,check.names=F)
dim (nasal_mild.DEG)
nasal_mild.DEG_up <- read.table("nasal_up-DEGs_mild.vs.controls-P0.001FC.csv",sep = "",header=T,check.names=F)
nasal_mild.DEG_down <- read.table("nasal_down-DEGs_mild.vs.controls-P0.001FC.csv",sep = "",header=T,check.names=F)

#nasal severe vs. control
nasal_severe.DEG <- read.table("nasal_allDEGs_severe.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)
dim (nasal_severe.DEG)
nasal_severe.DEG_up <- read.table("nasal_up-DEGs_severe.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)
nasal_severe.DEG_down <- read.table("nasal_down-DEGs_severe.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)

#bronchial mild vs. control
bronchial_mild.DEG <- read.table("bronchial_all-excluded-DEGs_mild.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)
dim (bronchial_mild.DEG)
bronchial_mild.DEG_up <- read.table("bronchial_excluded-up-DEGs_mild.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)
bronchial_mild.DEG_down <- read.table("bronchial_excluded-down-DEGs_mild.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)

#bronchial severe vs. control
bronchial_severe.DEG <- read.table("bronchial_all-excluded-DEGs_severe.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)
dim (bronchial_severe.DEG)
bronchial_severe.DEG_up <- read.table("bronchial_excluded-up-DEGs_severe.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)
bronchial_severe.DEG_down <- read.table("bronchial_excluded-down-DEGs_severe.vs.controls-FDR0.01FC.csv",sep = "",header=T,check.names=F)


########overlap

#all
overlap.mild <- merge (nasal_mild.DEG, bronchial_mild.DEG, by="hgnc_symbol")
overlap.severe <- merge (nasal_severe.DEG, bronchial_severe.DEG, by="hgnc_symbol")

#up-regulated genes
overlap.mild_up <- merge (nasal_mild.DEG_up, bronchial_mild.DEG_up, by="hgnc_symbol")
write.table(overlap.mild_up, "overlap_mild.up.csv")
overlap.severe_up <- merge (nasal_severe.DEG_up, bronchial_severe.DEG_up, by="hgnc_symbol")
write.table(overlap.severe_up, "overlap_severe.up.csv")

#down-regulated genes
overlap.mild_down <- merge (nasal_mild.DEG_down, bronchial_mild.DEG_down, by="hgnc_symbol")
write.table(overlap.mild_down, "overlap_mild.down.csv")
overlap.severe_down <- merge (nasal_severe.DEG_down, bronchial_severe.DEG_down, by="hgnc_symbol")
write.table(overlap.severe_down, "overlap_severe.down.csv")

down <- merge (bronchial_severe.DEG_down, nasal_severe.DEG_down, by="hgnc_symbol")

#############venn diagrams
library (VennDiagram)


#list genes
nasal_mild.up <- nasal_mild.DEG_up [, -1:-6]
nasal_mild.down <- nasal_mild.DEG_down [, -1:-6]
nasal_severe.up <- nasal_severe.DEG_up [, -1:-6]
nasal_severe.down <- nasal_severe.DEG_down [, -1:-6]

bronchial_mild.up <- bronchial_mild.DEG_up [, -1:-6]
bronchial_mild.down <- bronchial_mild.DEG_down [, -1:-6]
bronchial_severe.up <- bronchial_severe.DEG_up [, -1:-6]
bronchial_severe.down <- bronchial_severe.DEG_down [, -1:-6]


list_mild.up <- list (nasal_mild.up, bronchial_mild.up)
list_mild.down <- list (nasal_mild.down, bronchial_mild.down)
list_severe.up <- list (nasal_severe.up, bronchial_severe.up)
list_severe.down <- list(nasal_severe.down, bronchial_severe.down)


#mild COPD - down-regulated genes
venn.diagram(
  x=list_mild.down, 
  category.names = c("Nasal", "Bronchial"),
  filename = 'mildCOPD_down_venndiagram.png',
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.2,
  main = "Overlap down-regulated DEGs in mild COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#FFDF00",0.3), alpha('#3090C7',0.3)),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(180, -180),
  cat.dist = c(0.050, 0.040),
  cat.fontfamily = "sans",
)


#severe COPD - up-regulated genes
venn.diagram(
  x=list_severe.up,
  category.names = c("Nasal", "Bronchial"),
  filename = 'severeCOPD_up_venndiagram.png',
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 0.8,
  main = "Overlap up-regulated DEGs in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#FFDF00",0.3), alpha('#3090C7',0.3)),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0,8),
  cat.dist = c(0.040, 0.040),
  cat.fontfamily = "sans",
)



#severe COPD - down-regulated genes
venn.diagram(
  x=list_severe.down,
  category.names = c("Nasal", "Bronchial"),
  filename = 'severeCOPD_down_venndiagram.png',
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  lwd = 1.3,
  main = "Overlap down-regulated DEGs in severe COPD", 
  main.fontface = "bold",
  main.fontfamily = "sans", 
  main.cex = .4,
  fill = c(alpha("#FFDF00",0.3), alpha('#3090C7',0.3)),
  cex = .4,
  fontface = "bold", 
  fontfamily = "sans",
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 15),
  cat.dist = c(0.040, 0.040),
  cat.fontfamily = "sans",
)
