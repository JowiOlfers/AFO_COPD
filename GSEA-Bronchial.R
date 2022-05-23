############ ENRICHMENT ANALYSIS #############

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rgraphviz")
#BiocManager::install("topGO")
#BiocManager::install("ReactomePA")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggridges)
library(ggplot2)
library(tidyverse)
library (dplyr)
library (tibble)
#library(cowplot)

# set WD and pathways
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir_DEG.B <- file.path(".", "Nasal DEG_DE in bronchi")
data.dir_reflect <- file.path(data.dir_DEG.B, "results")
data.dir.Nasal <- file.path(".", "Nasal data")
data.dir.NU <- file.path(data.dir.Nasal, "unique selection")
results.dir <- file.path(".", "Enrichment analysis")


#data
data <- readr::read_csv(
  file.path(data.dir_reflect, "bronchial_unique-DEGs_severeCOPD-P0.01_NNLS.csv"))
names (data)[names(data) == "hgnc_symbol"] <- "SYMBOL"
names(data)[names(data) == "logFC.x"] <- "logFC"
data <- select(data, SYMBOL, logFC)

gene <- bitr(data$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
comb_gene <- merge(data, gene, by= "SYMBOL", all=F)
order_gene <- comb_gene[order(comb_gene$logFC, decreasing = T),]

genelist = order_gene$logFC
names(genelist) <- order_gene$ENTREZID
head(genelist)


#########################################################
### GO
#########################################################
ego <- gseGO(geneList  = genelist,
              keyType = "ENTREZID",
              OrgDb  = org.Hs.eg.db,
              ont  = "ALL",
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              verbose  = T)

plot4 <- ridgeplot(ego,10) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "p.adjust", option = "plasma") +
  labs(title = 'GO - pathway enrichment in severe COPD')
plot4
ggsave(filename = file.path (results.dir,"GO-all_bronchial.png"),plot4,height=10,width=12,dpi=600)

ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ego.res <- ego@result
dim(ego)

write.table(ego,
            file = file.path (results.dir, "GSEA_GO-all.bronchial_FDR0.05.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 

plot1 <- dotplot(ego, showCategory = 2)
plot1
plot2 <-dotplot(ego, showCategory = 5, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
plot2
plot21 <- dotplot(ego, showCategory = 10, split=".sign") + facet_grid(~.sign) 
plot21
plot3 <- gseaplot2(ego,1:5,color="red",pvalue_table = F) 
plot3

ggsave(filename = file.path (results.dir,"go-all_dotplot1_bronchial.pdf"),plot1,height=8,width=10)
ggsave(filename = file.path (results.dir,"go-all-dotplot2_bronchial.png"),plot2,height=8,width=10,dpi=600)
ggsave(filename = file.path (results.dir,"go-all_gseaplot_bronchial.pdf"),plot3,height=8,width=10)



###GO - ont = biological process or "BP"
egoA <- gseGO(geneList  = genelist,
             keyType = "ENTREZID",
             OrgDb  = org.Hs.eg.db,
             ont  = "BP",
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             verbose  = T)

plot4A <- ridgeplot(egoA,10) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "p.adjust", option = "plasma") +
  labs(title = 'GO - Biological Process - pathway enrichment in severe COPD')
plot4A
ggsave(filename = file.path (results.dir,"GO-BP_bronchial.png"),plot4A,height=10,width=12,dpi=600)

egoA <- setReadable(egoA, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ego.resA <- egoA@result
dim(egoA)

write.table(egoA,
            file = file.path (results.dir, "GSEA_GO-BP.bronchial_FDR0.05.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 



###GO - ont = Molecular Function or "MF"
egoB <- gseGO(geneList  = genelist,
             keyType = "ENTREZID",
             OrgDb  = org.Hs.eg.db,
             ont  = "MF",
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             verbose  = T)

plot4B <- ridgeplot(egoB ,10) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "p.adjust", option = "plasma") +
  labs(title = 'GO - Molecular Function - pathway enrichment in severe COPD')
plot4B
ggsave(filename = file.path (results.dir,"GO-MF.bronchial.png"),plot4B,height=10,width=12,dpi=600)


egoB <- setReadable(egoB, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ego.resB <- egoB@result
dim(egoB)

write.table(egoB,
            file = file.path (
              results.dir, "GSEA_GO-MF.bronchial_FDR0.05.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 



###GO - ont = CC for Cellular Component.
egoC <- gseGO(geneList  = genelist,
              keyType = "ENTREZID",
              OrgDb  = org.Hs.eg.db,
              ont  = "CC",
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              verbose  = T)

plot4C <- ridgeplot(egoC, 10) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "p.adjust", option = "plasma") +
  labs(title = 'GO - Cellular Component - pathway enrichment in severe COPD')
plot4C
ggsave(filename = file.path (results.dir,"GO-MF.bronchial.png"),plot4C,height=10,width=12,dpi=600)


egoC <- setReadable(egoC, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ego.resC <- egoC@result
dim(egoC)
write.table(egoC,
            file = file.path (
              results.dir, "GSEA_GO-CC.bronchial_FDR0.05.txt") ,sep = "\t",quote = F,col.names = T,row.names = F) 




#########################################################
### KEGG
#########################################################
kk <- gseKEGG(geneList = genelist,
              keyType  = 'kegg',
              organism = 'hsa',
              pvalueCutoff = 1,
              pAdjustMethod = "BH",
              verbose  = T)

plot23 <- ridgeplot(kk,10) + 
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "p.adjust", option = "plasma") +
  labs(title = 'KEGG pathway enrichment in severe COPD')
plot23
ggsave(filename = file.path (results.dir, "kegg.bronchial.png"),plot23,height=8,width=10,dpi=600)

kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(kk)

write.table(kk, 
            file = file.path (results.dir,"GSEA_kegg-bronchial_FDR0.05.txt"),sep = "\t",quote = F,col.names = T,row.names = F) 


plot21 <- dotplot(kk)
plot21
plot211 <- dotplot(kk, showCategory = 5, split=".sign") + facet_grid(~.sign) +
  labs(title = 'KEGG - pathway enrichment in severe COPD')
plot211
plot22 <- gseaplot2(kk,1:5,color="red",pvalue_table = F)
plot22
plot221 <- gseaplot2(kk,"hsa05164",color="green",pvalue_table = T)
plot221

ggsave(filename = file.path (results.dir,"kegg_dotplot1.bronchial.pdf"),plot21,height=8,width=10)
ggsave(filename = file.path (results.dir, "kegg_dotplot2.bronchial.png"),plot211,height=8,width=10,dpi=600)
ggsave(filename = file.path (results.dir,"kegg_gseaplot.bronchial.pdf"),plot22,height=8,width=10)
ggsave(filename = file.path (results.dir,"kegg_gseaplot1.bronchial.png"),plot221,height=8,width=10,dpi=600)




#########################################################
### reactomePA
#########################################################
Reactome <- gsePathway(geneList = genelist,
                       organism = 'human',
                       pvalueCutoff = 0.5,
                       pAdjustMethod = "BH",
                       verbose = T)
plot33 <- ridgeplot(Reactome,10) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "p.adjust", option = "plasma") +
  labs(title = 'Reactome pathway enrichment in severe COPD')
plot33
ggsave(filename = file.path (results.dir,"reactomePA.bronchial.png"), plot33,height=8,width=10,dpi=600)


Reactome <- setReadable(Reactome, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
dim(Reactome)

write.table(Reactome,
            file = file.path (results.dir,"GSEA_reactomePA-bronchial_FDR0.5.txt"),sep = "\t",quote = F,col.names = T,row.names = F) 


plot31 <- dotplot(Reactome, showCategory = 10)
plot31
plot311 <- dotplot(Reactome,showCategory = 9, split=".sign")+facet_grid(~.sign) 
plot311
plot32 <- gseaplot2(Reactome,1:4,color="red",pvalue_table = F)
plot32
#plot321 <- gseaplot2(Reactome,"R-HSA-5683826",color="green",pvalue_table = T) 
#plot322 <- gseaplot2(Reactome,"R-HSA-1474244",color="green",pvalue_table = T)
#plot323 <- gseaplot2(Reactome,"R-HSA-6805567",color="green",pvalue_table = T)


ggsave(filename = file.path (results.dir,"reactomepathway_dotplot1.bronchial.pdf"),plot31,height=8,width=20)
ggsave(filename = file.path (results.dir,"reactomepathway_dotplot2.bronchial.png"),plot311,height=8,width=20,dpi=600)
ggsave(filename = file.path (results.dir,"reactomepathway_gseaplot.bronchial.pdf"),plot32,height=8,width=10)
#ggsave(filename = file.path (results.dir,"RPA_gseaplot1.bronchial.png"),plot321,height=8,width=10,dpi=600)
#ggsave(filename = file.path (results.dir,"RPA_gseaplot2.bronchial.png"),plot322,height=8,width=10,dpi=600)
#ggsave(filename = file.path (results.dir,"RPA_gseaplot3.bronchial.png"),plot323,height=8,width=10,dpi=600)

