#First set working directory (go to session)

library (ggplot2)
library(ggrepel)
library (dplyr)

  
table1a <- read.table ("Nasal-mildCOPD.vs.control_FDR0.5.txt", sep ="", header = T)
tT1a <- select(table1a, hgnc_symbol, logFC, FDR)
df_data1a_pos = tT1a[which(tT1a$logFC > 0),]
df_data1a_neg = tT1a[which(tT1a$logFC < 0),]
data1a_neg <- df_data1a_neg [,-2:-3]


table1b <- read.table ("Nasal-mildCOPD.vs.control_P0.001.txt", sep ="", header = T)
tT1b <- select(table1b, hgnc_symbol, logFC, FDR)

df_data1b_pos = tT1b[which(tT1b$logFC > 0),]
data1b_pos <- df_data1b_pos [,-2:-3]
df_data1b_neg = tT1b[which(tT1b$logFC < 0),]
data1b_neg <- df_data1b_neg [,-2:-3]



table2 <- read.table ("nasal-severeCOPD.vs.control_FDR0.05.txt", sep ="", header = T)
tT2 <- select(table2, hgnc_symbol, logFC, FDR)

df_data2_pos = tT2[which(tT2$logFC > 0),]
data2_pos <- df_data2_pos [,-2:-3]
df_data2_neg = tT2[which(tT2$logFC < 0),]
data2_neg <- df_data2_neg [,-2:-3]


table3 <- read.table ("nasal-severe.vs.mild_FDR0.05.txt", sep ="", header = T)
tT3 <- select(table3, hgnc_symbol, logFC, FDR)

df_data3_pos = tT3[which(tT3$logFC > 0),]
data3_pos <- df_data3_pos [,-2:-3]
df_data3_neg = tT3[which(tT3$logFC < 0),]
data3_neg <- df_data3_neg [,-2:-3]





######################overlap and excluding tables
overlap_df.pos <-  merge (df_data2_pos, df_data3_pos, by= "hgnc_symbol")
write.table(overlap_df.pos, "nasal_overlap_df.pos.csv")

overlap_df.neg <-  merge (df_data2_neg, df_data3_neg, by= "hgnc_symbol")
filter(overlap_df.neg, "hgnc_symbol" !=df_data1a_neg$hgnc_symbol)
write.table(overlap_df.neg, "nasal_overlap_df.neg.csv")


#top 20 DEGs
overlap_df.pos = arrange (overlap_df.pos, desc(logFC.x))
top10_overlap.pos = head (overlap_df.pos[,1:3], 10)
overlap_df.neg = arrange (overlap_df.neg, desc(-logFC.x))
top10_overlap.neg = head (overlap_df.neg[,1:3], 10)
top20_overlap = rbind (top10_overlap.neg, top10_overlap.pos)
write.table(top20_overlap, "Nasal_top20.DEGs_uniquesevereCOPD.csv")




####################### Venn diagrams
library ("devtools")
library("ggVennDiagram")
library ("ggnewscale")

list_pos <- list (data2_pos, data3_pos)
list_neg <- list (data1a_neg, data2_neg, data3_neg)


#######Positive FC
ggVennDiagram(list_pos, label_alpha = 0, category.names = c("Severe COPD vs control","Severe vs mild COPD")) +
  guides(fill = guide_legend(title = "DEF genes")) +
  theme(legend.title = element_text(color = "black"), legend.position = "right") + ggplot2::scale_fill_manual(low="white", high = "grey") 
  + png("venndiagram_nasal.pos.png")


#######Negative FC
ggVennDiagram(list_neg, label_alpha = 0, category.names = c("Mild COPD vs control","Severe COPD vs control","Severe vs mild COPD")) +
  guides(fill = guide_legend(title = "DEF genes")) +
  theme(legend.title = element_text(color = "black"), legend.position = "right") +ggplot2::scale_fill_gradient(low="white",high = "grey") 
  + png("venndiagram_nasal.neg.png")
  
  
  

# Jos edit