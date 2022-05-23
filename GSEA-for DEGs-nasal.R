#BiocManager::install("ggupset")
library(clusterProfiler)
library(enrichplot)
library(topGO)
library(GO.db)
library(SparseM)
library(matrixStats)
library(cowplot)
library(grid)
library(ReactomePA)
library(scales)
library(ggplot2)
library(ggrepel)
library(grid)
library(ggupset)
library(ggridges)
library(Rgraphviz)

# set WD and pathways
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir_DEG.B <- file.path(".", "Nasal DEG_DE in bronchi")
data.dir_reflect <- file.path(data.dir_DEG.B, "results")
data.dir.Nasal <- file.path(".", "Nasal data")
data.dir.NU <- file.path(data.dir.Nasal, "unique selection")
results.dir <- file.path(".", "Enrichment analysis")
results.dir_sep <- file.path(results.dir, "results_sep")
results.dir_sep.GO <- file.path(results.dir_sep, "GO")
results.dir_sep.KEGG <- file.path(results.dir_sep, "KEGG")
results.dir_sep.R <- file.path(results.dir_sep, "Reactome")

# DATA
data <- readr::read_csv(
  file.path(data.dir.NU, "nasal_unique-DEGs_severeCOPD-NNLS-FDR0.05FC2.csv"))
names (data)[names(data) == "hgnc_symbol"] <- "SYMBOL"

data1 <- data[which(data$logFC.x >= 1),]
data2 <- data[which(data$logFC.x <= -1),]
gene <- bitr(data$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene1 <- bitr(data1$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene2 <- bitr(data2$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


 
# -------------------------------------------------------------------------------------------
#######GO
ego_up <- enrichGO(gene1$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff = 0.05,qvalueCutoff  = 1, pAdjustMethod = "BH")
ego_up <- setReadable(ego_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(ego_up)
write.table(ego_up,
            file = file.path (results.dir_sep.GO, 
                              "DEGs_GO-UP.nasal.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 

ego_down <- enrichGO(gene2$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",pvalueCutoff = 0.05,qvalueCutoff  = 1, pAdjustMethod = "BH")
ego_down <- setReadable(ego_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(ego_down)
write.table(ego_down,
            file = file.path (results.dir_sep.GO, 
                              "DEGs_GO-DOWN.nasal.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 


### PLOTS ###

## ------------barplot----------------------
plot1_up <-barplot(ego_up, showCategory=10,drop = T)
plot1_up
## ------------dotplot----------------------
plot2_up <-dotplot(ego_up, showCategory=10)
plot2_up
plot3_up <-dotplot(ego_up, split="ONTOLOGY", showCategory = 5) + facet_grid(ONTOLOGY~., scale="free")
plot3_up

ggsave(filename = file.path (results.dir_sep.GO,"barplot.GO-up-Nasal.png"), plot1_up,height=4,width=10,dpi=600)
ggsave(filename = file.path (results.dir_sep.GO,"dotplot-GO-up-Nasal.png"), plot2_up,height=8,width=10,dpi=600)
ggsave(filename = file.path (results.dir_sep.GO,"dotplot2-GO-up-Nasal.pdf"),plot3_up,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.GO,"dotplot2-GO-up-Nasal.png"),plot3_up,height=8,width=10)


## ------------barplot----------------------
plot1_down <-barplot(ego_down, showCategory=10,drop = T)
plot1_down
## ------------dotplot----------------------
plot2_down <-dotplot(ego_down, showCategory=10)
plot2_down
plot3_down <-dotplot(ego_down, split="ONTOLOGY", showCategory = 5) + facet_grid(ONTOLOGY~., scale="free")
plot3_down

ggsave(filename = file.path (results.dir_sep.GO,"barplot.GO-down-Nasal.png"), plot1_down,height=4,width=10,dpi=600)
ggsave(filename = file.path (results.dir_sep.GO,"dotplot-GO-down-Nasal.png"), plot2_down,height=8,width=10,dpi=600)
ggsave(filename = file.path (results.dir_sep.GO,"dotplot2-GO-down-Nasal.pdf"), plot3_down,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.GO,"dotplot2-GO-down-Nasal.png"), plot3_down,height=8,width=10)



######## KEGG
ekegg_up <- enrichKEGG(gene1$ENTREZID, organism='hsa',keyType="ncbi-geneid",pvalueCutoff=0.5,qvalueCutoff =1,pAdjustMethod='BH')
ekegg_up <- setReadable(ekegg_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(ekegg_up)
write.table(ekegg_up,
            file = file.path (results.dir_sep.KEGG, 
                              "DEGs_KEGG-UP.nasal.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 


ekegg_down <- enrichKEGG(gene2$ENTREZID, organism='hsa',keyType="ncbi-geneid",pvalueCutoff=0.5,qvalueCutoff =1,pAdjustMethod='BH')
ekegg_down <- setReadable(ekegg_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(ekegg_down)
write.table(ekegg_down,
            file = file.path (results.dir_sep.KEGG, 
                              "DEGs_KEGG-DOWN.nasal.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 


### PLOTS ###

########UP regulated
## ------------barplot----------------------
plot11_up <-barplot(ekegg_up, showCategory=10)
plot11_up

## ------------dotplot----------------------
plot22_up <-dotplot(ekegg_up, showCategory=10)
plot22_up

KEGG_UP.pw <- pairwise_termsim (ekegg_up)
plot33_up <-emapplot(KEGG_UP.pw, color = "p.adjust", cex_label_category = 0.7)
plot33_up
plot331_up <-emapplot(KEGG_UP.pw, color = "p.adjust", cex_label_category = 0.8, showCategory=10)
plot331_up

plot44_up <-cnetplot(ekegg_up,showCategory=10,readable = T)
plot44_up
plot441_up <-cnetplot(ekegg_up,showCategory=10,circular=TRUE,colorEdge=TRUE,readable = T)
plot441_up
plot55_up <-upsetplot(ekegg_up)
plot55_up
plot6_up <-heatplot(ekegg_up)
plot6_up

ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_barplot-UP-Nasal.png"),plot11_up,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_dotplot-UP-Nasal.png"),plot22_up,height=8,width=10,dpi=600)
ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_emapplot1-UP-Nasal.pdf"),plot33_up,height=10,width=10)
ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_emapplot2-UP-Nasal.pdf"),plot331_up,height=10,width=10)
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_cnetplot1-UP-Nasal.pdf"),plot44_up,height=8,width=10)
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_cnetplot2-UP-Nasal.pdf"),plot441_up,height=8,width=10)      
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_upset-UP-Nasal.pdf"),plot55_up,height=8,width=10)      
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_heatplot-UP-Nasal.pdf"),plot6_up,height=8,width=10)    


######## DOWN regulated
## ------------barplot----------------------
plot11_down <-barplot(ekegg_down, showCategory=10)
plot11_down

## ------------dotplot----------------------
plot22_down <-dotplot(ekegg_down, showCategory=10)
plot22_down

KEGG_DOWN.pw <- pairwise_termsim (ekegg_down)
plot33_down <-emapplot(KEGG_DOWN.pw, color = "p.adjust", cex_label_category = 0.6)
plot33_down
plot331_down <-emapplot(KEGG_DOWN.pw, color = "p.adjust", cex_label_category = 0.8, showCategory=10)
plot331_down

plot44_down <-cnetplot(ekegg_down,showCategory=10,readable = T)
plot44_down
plot441_down <-cnetplot(ekegg_down,showCategory=10,circular=TRUE,colorEdge=TRUE,readable = T)
plot441_down
plot55_down <-upsetplot(ekegg_down)
plot55_down
plot6_down <-heatplot(ekegg_down)
plot6_down

ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_barplot-DOWN-Nasal.png"),plot11_down,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_dotplot-DOWN-Nasal.png"),plot22_down,height=8,width=10,dpi=600)
ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_emapplot1-DOWN-Nasal.pdf"),plot33_down,height=10,width=10)
ggsave(filename = file.path (results.dir_sep.KEGG,"KEGG_emapplot2-DOWN-Nasal.pdf"),plot331_down, height=10,width=10)
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_cnetplot1-DOWN-Nasal.pdf"),plot44_down,height=8,width=10)
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_cnetplot2-DOWN-Nasal.pdf"),plot441_down,height=8,width=10)      
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_upset-DOWN-Nasal.pdf"),plot55_down,height=8,width=10)      
ggsave(filename = file.path (results.dir_sep.KEGG, "KEGG_heatplot-DOWN-Nasal.pdf"),plot6_down,height=8,width=10)    



##ReactomePA
RPA_up <- enrichPathway(gene1$ENTREZID,readable=T, pvalueCutoff = 0.05,qvalueCutoff  = 1,pAdjustMethod = "BH") 
dim(RPA_up)
write.table(RPA_up,
            file = file.path (results.dir_sep.R, 
                              "DEGs_ReactomePA-UP.nasal.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 


RPA_down <- enrichPathway(gene2$ENTREZID,readable=T, pvalueCutoff = 0.05,qvalueCutoff  = 1,pAdjustMethod = "BH") 
dim(RPA_down)
write.table(RPA_down,
            file = file.path (results.dir_sep.R, 
                              "DEGs_ReactomePA-DOWN.nasal.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 


### PLOTS ###

######## UP regulated
## ------------barplot----------------------
plot111_up <- barplot(RPA_up, showCategory=10)
plot111_up
## ------------dotplot----------------------
plot222_up <-dotplot(RPA_up, showCategory=10)
plot222_up

RPA_UP.pw <- pairwise_termsim (RPA_up)
plot333_up <- emapplot(RPA_UP.pw, color = "p.adjust", cex_label_category = 0.8)
plot333_up
plot3331_up <- emapplot(RPA_UP.pw, color = "p.adjust", cex_label_category = 0.8, showCategory=10)
plot3331_up

plot444_up <-cnetplot(RPA_up,categorySize="pvalue")
plot444_up
plot42_up <- cnetplot(RPA_up,showCategory=10,circular=TRUE,colorEdge=TRUE,readable = T)
plot42_up

#plot555_up <- viewPathway("R-HSA-209776", readable=TRUE)
#plot555_up
plot666_up <-heatplot(RPA_up)
plot666_up

ggsave(filename = file.path (results.dir_sep.R,"RPA_barplot-UP-Nasal.png"),plot111_up,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_dotplot-UP-Nasal.png"), plot222_up,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_emapplot-UP-Nasal.png"), plot333_up,height=10,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_emapplot2-UP-Nasal.png"), plot3331_up,height=10,width=10, dpi=600)

#ggsave(filename = file.path (results.dir_sep.R, "RPA_viewPathway-UP-Nasal.png"), plot555_up,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_heatplot-UP-Nasal.png"), plot666_up,height=8,width=10, dpi=600)


######## DOWN regulated
## ------------barplot----------------------
plot111_down <- barplot(RPA_down, showCategory=10)
plot111_down
## ------------dotplot----------------------
plot222_down <-dotplot(RPA_down, showCategory=10)
plot222_down

RPA_DOWN.pw <- pairwise_termsim (RPA_down)
plot333_down <- emapplot(RPA_DOWN.pw, color = "p.adjust", cex_label_category = 0.6)
plot333_down
plot3331_down <- emapplot(RPA_DOWN.pw, color = "p.adjust", cex_label_category = 0.7, showCategory=15)
plot3331_down

plot444_down <-cnetplot(RPA_down,categorySize="pvalue")
plot444_down
plot42_down <- cnetplot(RPA_down,showCategory=10,circular=TRUE,colorEdge=TRUE,readable = T)
plot42_down
#plot555_down <- viewPathway("R-HSA-209776", readable=TRUE)
#plot555_down
plot666_down <-heatplot(RPA_down)
plot666_down

ggsave(filename = file.path (results.dir_sep.R,"RPA_barplot-DOWN-Nasal.png"),plot111_down,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_dotplot-DOWN-Nasal.png"), plot222_down,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_emapplot-DOWN-Nasal.png"), plot333_down,height=10,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_emapplot2-DOWN-Nasal.png"), plot3331_down,height=10,width=10, dpi=600)
#ggsave(filename = file.path (results.dir_sep.R, "RPA_viewPathway-DOWN-Nasal.png"), plot555_down,height=8,width=10, dpi=600)
ggsave(filename = file.path (results.dir_sep.R, "RPA_heatplot-DOWN-Nasal.png"), plot666_down,height=8,width=10, dpi=600)




