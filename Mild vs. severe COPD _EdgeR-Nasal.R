library (limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(heatmap3)
library(gplots)
library(ggfortify)
library(tidyverse)
library (dplyr)

#Directories 
setwd("~/HVHL/Stages/Afstudeerstage_2022/DE analysis/Nasal data")
library(here)

data.dir <- file.path(".", "data")
results.dir <- file.path(here::here(), "results")
dir.create(results.dir, recursive=TRUE)
results.dir.img <- file.path(results.dir, "img")


x <- readr::read_csv(
  file.path(data.dir, "nasal rawcounts.csv"))
array <- data.frame (x)

y <- readr::read_csv(
  file.path(data.dir, "nasal clinical data_CD_mild.csv"))
s <- data.frame (y)

row.names(array)=array[,1]
array=array[,-1]
array=array[,s$SAMID]

array1=array[,which(s$severe.vs.mild==1)]
s1=s[which(s$severe.vs.mild==1),]



# EdgeR: disease

###############################array3 - Mild vs. severe COPD
expressionDataToUse= array1
samplesToUse=s1

##########create factors
group=as.factor(samplesToUse$group)
smoking=as.factor(samplesToUse$smoking)
age=as.numeric(samplesToUse$age)
gender=as.factor(samplesToUse$gender)
packyears=as.numeric(samplesToUse$packyears)
years.cessation=as.numeric(samplesToUse$years.of.cessation)

#CD method: CIBERSORT
#proportion_dendritic=as.numeric(samplesToUse$proportion_ciber_Dendritic)
#proportion_goblet1N=as.numeric(samplesToUse$proportion_ciber_Goblet_1N)

#CD method: NNLS
#proportion_basal=as.numeric(samplesToUse$proportion_nnls_Basal)
#proportion_fibro=as.numeric(samplesToUse$proportion_nnls_Fibroblasts)
#proportion_goblet1N=as.numeric(samplesToUse$proportion_nnls_Goblet_1N)


###replace years.cessation NA with average value of corresponding group
df = data.frame(years.cessation, group)
df1 <- subset(df, group == "1")
df2 <- subset(df, group == "2")

#calculate average of group 1 and 2
dfx <- df[!(is.na(df$years.cessation)), ]

dfx1 <- subset(dfx, group == "1")
mean.cessation_1 <- mean(dfx1$years.cessation)
print (mean.cessation_1)

dfx2 <- subset(dfx, group == "2")
mean.cessation_2 <- mean(dfx2$years.cessation)
print (mean.cessation_2)

#replace NA with average cessation of group
df1[is.na(df1)]<- mean.cessation_1
df2[is.na(df2)]<- mean.cessation_2
df_na.rm <- rbind(df1,df2)
years.of.cessation <- as.numeric(df_na.rm$years.cessation)



######List and filter on expression
total_nasaldata_DGEL <- DGEList(expressionDataToUse)
keep <- filterByExpr(total_nasaldata_DGEL,group=group)
nasaldata_DGEL <- total_nasaldata_DGEL[keep,, keep.lib.sizes=FALSE]
dim(nasaldata_DGEL)
nasaldata_DGEL$samples

##############Normalization - TMM
nasaldata_DGEL <- calcNormFactors(nasaldata_DGEL, method="TMM")
nasaldata_DGEL$samples
#plotMDS(nasaldata_DGEL)

#count per million (cpm) read
normalized_counts <- cpm(nasaldata_DGEL, log = TRUE)
normalized_counts <- as.data.frame (normalized_counts)
readr::write_csv(normalized_counts, 
                 file = file.path (results.dir, "nasal-log2CPM.mildvssevere.csv"))
readr::write_csv(normalized_counts, 
                 file = file.path (results.dir, "nasal-log2CPM.mildvssevere.txt"))


##############Differential expression analysis
design <- model.matrix(~group + age + gender + packyears + years.of.cessation)
nasaldata_DGEL <- estimateDisp(nasaldata_DGEL, design)
#plotBCV(nasaldata_DGEL)

fit <- glmFit(nasaldata_DGEL, design)
lrt <- glmLRT(fit, coef = 2)
results <- topTags(lrt,n=nrow(nasaldata_DGEL))
#plotMeanVar(nasaldata_DGEL, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)


#Ensembl - hgnc_symbols
ensembl.dir <- file.path("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/CD")
source(file.path(ensembl.dir, "_helpers.r"), verbose=TRUE, chdir=TRUE)

tT1=cbind(results$table,row.names(results$table))

geneData <- getGenedataByEnsemblId(
  ensemblIds = tT1$row.names(results$table) %>% unique(),
  file.location = data.dir
) %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::filter(
    dplyr::row_number() == 1,
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    hgnc_symbol,
    ensembl_gene_id,
  )

#mart = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
#                  dataset="hsapiens_gene_ensembl")

#genesID = as.character(rownames(tT1))
#G_list = getBM(
#  filters = "ensembl_gene_id",
#  attributes = c(
#    "ensembl_gene_id",
#    "hgnc_symbol"
#  ),
#  values = genesID,
#  mart = mart
#)
#G_list=as.matrix(G_list)


tT1=merge(
  x = tT1,
  y = geneData,
  all.x = T,
  all.y = F,
  by.x = "row.names(results$table)",
  by.y ="ensembl_gene_id")


rownames(tT1)=tT1[,1]
#standard an error due to DUPLICATES:'ENSG00000187510', 'ENSG00000276085'  and 'ENSG00000255374'
names(tT1)[1] <- "ENSGid"
tT1 <- tT1 [!(is.na(tT1$hgnc_symbol) | tT1$hgnc_symbol == ""), ]
readr::write_csv(tT1, file = file.path (results.dir, "nasal-mild.vs.severeCOPD.txt"))

tT2=tT1[which(tT1$FDR<0.05), ]
names(tT2)[1] <- "ENSGid"
readr::write_csv(tT2, file = file.path (results.dir, "nasal-mild.vs.severeCOPD_FDR0.05.txt"))

tT3=tT1[which(tT1$FDR<0.01), ]
names(tT3)[1] <- "ENSGid"
readr::write_csv(tT3, file = file.path (results.dir, "nasal-mild.vs.severeCOPD_FDR0.01.txt"))

#tT4=tT1[which(tT1$PValue<0.05),]
#names(tT4)[1] <- "ENSGid"
#readr::write_csv(tT4, file = file.path (results.dir, "nasal-mild.vs.severeCOPD_P0.05.txt"))

#tT5=tT1[which(tT1$PValue<0.01),]
#names(tT5)[1] <- "ENSGid"
#readr::write_csv(tT5, file = file.path (results.dir, "nasal-mild.vs.severeCOPD_P0.01.txt"))

#tT6=tT1[which(tT1$PValue<0.001),]
#names(tT6)[1] <- "ENSGid"
#readr::write_csv(tT6, file = file.path (results.dir, "nasal-mild.vs.severeCOPD_P0.001.txt"))




#Selection and presentation significant up and down regulated genes
######volcano plot
library(ggpubr)
library(ggthemes)

#FDR 0.05 selected genes
tT1$logFDR <- -log10(tT1$FDR)
tT1$Group <- "No diff"
tT1$Group[which((tT1$FDR < 0.05)&(tT1$logFC > 1))] = "Up"
tT1$Group[which((tT1$FDR < 0.05)&(tT1$logFC < -1))] = "Down"

#amount of DEG with FDR < 0.05 and LogFC > 1 of < -1
table(tT1$Group)

tT1$label = ""
tT1 <- tT1[order(tT1$FDR),]
#select that the top 10 DEG show names in the plot
up_gene <- head(tT1$hgnc_symbol[which(tT1$Group == "Up")],10)
down_gene <- head(tT1$hgnc_symbol[which(tT1$Group == "Down")],10)

#present top 10 up and down regulated genes 
up_gene
down_gene

tT1_gene <- c(as.character(up_gene), as.character(down_gene))
tT1$label[match(tT1_gene, tT1$hgnc_symbol)] <- tT1_gene

p <- ggscatter(tT1, x = "logFC", y = "logFDR", 
               color = "Group", 
               palette = c("#2f5688","#BBBBBB","#CC0000"), 
               size = 1,
               label = tT1$label, 
               font.label = 8, 
               repel = T,
               xlab = "log2FoldChange", ylab = "-log10(FDR)") + theme_base() +
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = c(-1, 1), linetype="dashed")
p
ggsave (filename = file.path (results.dir.img, "Volcano-mildvssevere-FDR0.05.png"), width=25,height=25,units="cm",dpi=600 )
ggsave (filename = file.path (results.dir.img, "gene_diff-mildvssevere-FDR0.05.pdf"), width=25,height=25,units="cm")



##### Pvalue selected genes - P<0.001
#tT1$logPValue <- -log10(tT1$PValue)
#tT1$GroupP <- "No diff"
#tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC > 1))] = "Up"
#tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC < -1))] = "Down"

#amount of DEG with Pvalue < 0.001 and LogFC > 0 of < 0
#table(tT1$Group)

#tT1$labelP = ""
#tT1 <- tT1[order(tT1$PValue),]

#select that the top 10 DEG show names in the plot
#up_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Up")],10)
#down_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Down")],10)

#present top 10 up and down regulated genes 
#up_geneP
#down_geneP

#tT1_geneP <- c(as.character(up_geneP), as.character(down_geneP))
#tT1$labelP[match(tT1_geneP, tT1$hgnc_symbol)] <- tT1_geneP

#pp <- ggscatter(tT1, x = "logFC", y = "logPValue", 
#                color = "GroupP", 
#                palette = c("#2f5688","#BBBBBB","#CC0000"), 
#                size = 1,
#                label = tT1$labelP, 
#                font.label = 8, 
#                repel = T,
#                xlab = "log2FoldChange", ylab = "-log10(PValue)") + theme_base() +
#  geom_hline(yintercept = -log10(0.001),linetype="dashed")+
#  geom_vline(xintercept = c(-1, 1),linetype="dashed")
#pp
#ggsave (filename = file.path (results.dir.img, "Volcano-WC-mildvssevere-PValue.png"), width=30,height=25,units="cm",dpi=600 )
#ggsave (filename = file.path (results.dir.img, "gene_diff-mildvssevere-PValue.pdf"), width=25,height=25,units="cm")




######Create data.frame with normalized expression for genes with FDR sig. DEGs
#select up and down regulated genes with FDR < 0.05 and FC > 2 or < -2
tT21=tT2[which(tT2$logFC >= 1),]
dim (tT21)
readr::write_csv(tT21, file = file.path (results.dir, "nasal_up-DEGs_mild.vs.severeCOPD-FDR0.05FC2.csv"))

tT22=tT2[which(tT2$logFC <= -1), ]
dim (tT22)
readr::write_csv(tT22, file = file.path (results.dir, "nasal_down-DEGs_mild.vs.severeCOPD-FDR0.05FC2.csv"))

tT23 <- rbind (tT21, tT22)
dim (tT23)
readr::write_csv(tT23, file = file.path (results.dir, "nasal_allDEGs_mild.vs.severeCOPD-FDR0.05FC2.csv"))


#select top 20 - FDR < 0.05
tT21 = arrange (tT21, desc(logFC))
tT21a = head(tT21[,1:7], 10)
tT22 = arrange (tT22, desc(-logFC))
tT22a = head(tT22[,1:7], 10)
tT23a <- rbind (tT21a, tT22a)
readr::write_csv(tT23a, file = file.path (results.dir, "nasal_top20.DEGs_mild.vs.severeCOPD-FDR0.05FC2.csv"))

#select top 10 - FDR < 0.05
tT21 = arrange (tT21, desc(logFC))
tT221a = head(tT21[,1:7], 5)
tT22 = arrange (tT22, desc(-logFC))
tT222a = head(tT22[,1:7], 5)
tT223a <- rbind (tT221a, tT222a)
readr::write_csv(tT223a, file = file.path (results.dir, "nasal_top10.DEGs_mild.vs.severeCOPD-FDR0.05FC2.csv"))


#select up and down regulated genes with FDR < 0.01 and FC > 2 or < -2
#tT31=tT3[which(tT3$logFC >= 1),] 
#readr::write_csv(tT31, file = file.path (results.dir, "nasal_up-DEGs_mild.vs.severeCOPD-FDR0.01FC2.csv"))

#tT32=tT3[which(tT3$logFC <= -1), ]
#readr::write_csv(tT31, file = file.path (results.dir, "nasal_down-DEGs_mild.vs.severeCOPD-FDR0.01FC2.csv"))

#tT33 <- rbind (tT31, tT32)
#write.table(tT33, "nasal_allDEGs_severe.vs.mildCOPD-FDR0.01FC2.csv")

#select top 20 FDR < 0.01 and FC > 2 or < -2
#tT31 = arrange (tT31, desc(logFC))
#tT31a = head(tT31[,1:7], 10)
#tT32 = arrange (tT32, desc(-logFC))
#tT32a = head(tT32[,1:7], 10)
#tT33a <- rbind (tT31a, tT32a)
#readr::write_csv(tT33a, file = file.path (results.dir, "nasal_top20.DEGs_mild.vs.severeCOPD-FDR0.01FC2.csv"))

#select top 10 FDR < 0.01 and FC > 2 or < -2
#tT31 = arrange (tT31, desc(logFC))
#tT321a = head(tT31[,1:7], 5)
#tT32 = arrange (tT32, desc(-logFC))
#tT322a = head(tT32[,1:7], 5)
#tT323a <- rbind (tT321a, tT322a)
#readr::write_csv(tT323a, file = file.path (results.dir, "nasal_top20.DEGs_mild.vs.severeCOPD-FDR0.01FC2.csv"))


##############Combine normalized raw counts with expression data
d <- normalized_counts
ENSGid <- rownames(d)
rownames(d) <- NULL
normalized_counts2 <- cbind (ENSGid, d)

normalized_counts.FDR <- merge(tT23, normalized_counts2, by="ENSGid" )
normalized_counts.FDRa <- merge(tT23a, normalized_counts2, by="ENSGid" )

row.names(normalized_counts.FDR) <- normalized_counts.FDR$hgnc_symbol
normalized_counts.FDR$LR <- normalized_counts.FDR$ENSGid <- normalized_counts.FDR$logFC <- normalized_counts.FDR$logCPM <- normalized_counts.FDR$FDR <- normalized_counts.FDR$PValue <- normalized_counts.FDR$hgnc_symbol <- NULL
readr::write_csv(normalized_counts.FDR, file = file.path (results.dir, "nasal-normalized_counts.FDR0.05.mild.vs.severeCOPD.csv"))

row.names(normalized_counts.FDRa) <- normalized_counts.FDRa$hgnc_symbol
normalized_counts.FDRa$LR <- normalized_counts.FDRa$ENSGid <- normalized_counts.FDRa$logFC <- normalized_counts.FDRa$logCPM <- normalized_counts.FDRa$FDR <- normalized_counts.FDRa$PValue <- normalized_counts.FDRa$hgnc_symbol <- NULL
readr::write_csv(normalized_counts.FDRa, file = file.path (results.dir, "nasal-normalized_counts.top20-FDR0.05.mild.vs.severeCOPD.csv"))



#######################heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

table (s1$group)
annotation_col = data.frame(Group = factor(c(rep("Severe COPD", 113),rep("Mild COPD", 24))))
rownames(annotation_col)

data <- normalized_counts.FDRa
colnames(data)
rownames(annotation_col) <- colnames(data)
head(annotation_col)

data1 <- data [,annotation_col$Group == "Mild COPD"]
data2 <- data [,annotation_col$Group == "Severe COPD"]

HC1 <- hclust (dist(t(data1), method = "manhattan"), method = "complete", members = NULL)
HC2 <- hclust (dist(t(data2), method = "manhattan"), method = "complete", members = NULL)

data1 <- data1[,HC1$order]
data2 <- data2[,HC2$order]
datax <- cbind(data2, data1)

rownames (annotation_col) <- colnames (datax)

ph <- pheatmap(datax,
              scale="row",
              annotation_col = annotation_col,
              annotation_legend = T,
              annotation_names_col = F,
              number_format="%.2e",
              border="white",  
              color=colorRampPalette(c("navy", "white", "red"))(100),
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
              filename = file.path (results.dir.img, "Heatmap_nasal_mild.vs.severeCOPD_FDR0.05.png") ,dpi=600)
