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
  file.path(data.dir, "nasal clinical data_CD.csv"))
s <- data.frame (y)

row.names(array)=array[,1]
array=array[,-1]
array=array[,s$SAMID]


array1=array[,which(s$mild.vs.control==1)]
s1=s[which(s$mild.vs.control==1),]

array2=array[,which(s$severe.vs.control==1)]
s2=s[which(s$severe.vs.control==1),]

array3=array[,which(s$severe.vs.mild==1)]
s3=s[which(s$severe.vs.mild==1),]

#array1= contains groups control and mild COPD (n=46)
#array2= contains groups control and severe COPD (n=135)
#array3= contains groups mild and severe COPD (n=137)



#EdgeR: disease

################################## array1 - Mild vs. control
expressionDataToUse= array1
samplesToUse=s1

#create factors
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

#calculate average of group 0
dfx <- df[!(is.na(df$years.cessation)), ]
dfx0 <- subset(dfx, group == "0")
mean.cessation_0 <- mean(dfx0$years.cessation)
print (mean.cessation_0)

#replace NA with average cessation of group
df[is.na(df)]<- mean.cessation_0
years.of.cessation <- as.numeric(df$years.cessation)
print(years.of.cessation)



#List and filter on expression
total_nasaldata_DGEL <- DGEList(expressionDataToUse)
keep <- filterByExpr(total_nasaldata_DGEL,group=group)
nasaldata_DGEL <- total_nasaldata_DGEL[keep,, keep.lib.sizes=FALSE]
dim(nasaldata_DGEL)
nasaldata_DGEL$samples

################Normalization - TMM
nasaldata_DGEL <- calcNormFactors(nasaldata_DGEL, method="TMM")
nasaldata_DGEL$samples
#plotMDS(nasaldata_DGEL)

#count per million (cpm) read
normalized_counts <- cpm(nasaldata_DGEL, log = TRUE)
normalized_counts <- as.data.frame (normalized_counts)
readr::write_csv(normalized_counts, 
  file = file.path (results.dir, "nasal-log2CPM.mildvscontrol.csv"))
readr::write_csv(normalized_counts, 
  file = file.path (results.dir, "nasal-log2CPM.mildvscontrol.txt"))




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
#standard an error due to DUPLICATES:'ENSG00000187510', 'ENSG00000255374' and ENSG00000276085
names(tT1)[1] <- "ENSGid"
tT1 <- tT1 [!(is.na(tT1$hgnc_symbol) | tT1$hgnc_symbol == ""), ]
readr::write_csv(tT1, file = file.path (results.dir, "nasal-mildCOPD.vs.control.txt"))

tT2=tT1[which(tT1$FDR<0.05),]
names(tT2)[1] <- "ENSGid"
readr::write_csv(tT2, file = file.path (results.dir, "nasal-mildCOPD.vs.control_FDR0.05.txt"))

#FDR 0.25 again 0 DEG

tT3=tT1[which(tT1$FDR<0.5),]
names(tT3)[1] <- "ENSGid"
readr::write_csv(tT3, file = file.path (results.dir, "nasal-mildCOPD.vs.control_FDR0.5.txt"))

tT4=tT1[which(tT1$PValue<0.05),]
names(tT4)[1] <- "ENSGid"
#readr::write_csv(tT4, file = file.path (results.dir, "nasal-mildCOPD.vs.control_P0.05.txt"))

tT5=tT1[which(tT1$PValue<0.01),]
names(tT5)[1] <- "ENSGid"
#readr::write_csv(tT5, file = file.path (results.dir, "nasal-mildCOPD.vs.control_P0.01.txt"))

tT6=tT1[which(tT1$PValue<0.001),]
names(tT6)[1] <- "ENSGid"
readr::write_csv(tT6, file = file.path (results.dir, "nasal-mildCOPD.vs.control_P0.001.txt"))



#Selection and presentation significant up and down regulated genes
###volcano plot
library(ggpubr)
library(ggthemes)

############ FDR 0.5 selected genes
tT1$logFDR <- -log10(tT1$FDR)
tT1$Group <- "No diff"
tT1$Group[which((tT1$FDR < 0.5)&(tT1$logFC > 0))] = "Up"
tT1$Group[which((tT1$FDR < 0.5)&(tT1$logFC < 0))] = "Down"

#amount of DEG with FDR < 0.5 and LogFC > 0 of < 0
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
  geom_hline(yintercept = -log10(0.5),linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")
p
ggsave (filename = file.path (results.dir.img, "Volcano-WC-mildvscontrol-FDR-0.5.png"), width=25,height=25,units="cm",dpi=600 )
ggsave (filename = file.path (results.dir.img, "gene_diff-mildvscontrol-FDR-0.5.pdf"), width=25,height=25,units="cm")


###Pvalue selected genes - P<0.001
tT1$logPValue <- -log10(tT1$PValue)
tT1$GroupP <- "No diff"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC > 0))] = "Up"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC < 0))] = "Down"

#amount of DEG with Pvalue < 0.001 and LogFC > 0 of < 0
table(tT1$GroupP)

tT1$labelP = ""
tT1 <- tT1[order(tT1$PValue),]

#select that the top DEG show names in the plot
up_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Up")],10)
down_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Down")],15)

#present top 10 up and down regulated genes 
up_geneP
down_geneP

tT1_geneP <- c(as.character(up_geneP), as.character(down_geneP))
tT1$labelP[match(tT1_geneP, tT1$hgnc_symbol)] <- tT1_geneP

p1 <- ggscatter(tT1, x = "logFC", y = "logPValue", 
               color = "GroupP", 
               palette = c("#2f5688","#BBBBBB","#CC0000"), 
               size = 1,
               label = tT1$labelP, 
               font.label = 8, 
               repel = T,
               xlab = "log2FoldChange", ylab = "-log10(PValue)") + theme_base() +
  geom_hline(yintercept = -log10(0.001),linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")
p1
ggsave (filename = file.path (results.dir.img, "Volcano-WC-mildvscontrol-PValue.png"), width=30,height=30,units="cm",dpi=600 )
ggsave (filename = file.path (results.dir.img, "gene_diff-mildvscontrol-PValue.pdf"), width=25,height=30,units="cm")



#select up and down regulated genes with P < 0.001 and FC > 2 or < -2
tT61=tT6[which(tT6$logFC > 0),] 
dim (tT61)
readr::write_csv(tT61, file = file.path (results.dir, "nasal_up-DEGs_mild.vs.controls-P0.001FC2.csv"))

tT62=tT6[which(tT6$logFC < 0), ]
dim (tT62)
readr::write_csv(tT62, file = file.path (results.dir, "nasal_down-DEGs_mild.vs.controls-P0.001FC2.csv"))

tT63 <- rbind (tT61, tT62)
dim (tT63)
readr::write_csv(tT63, file = file.path (results.dir, "nasal_allDEGs_mild.vs.controls-P0.001FC2.csv"))

#select top 20 - top 10 up and top 10 down regulated genes
library(dplyr)
tT61 = arrange (tT61, desc(logFC))
tT61a = head(tT61[,1:7], 10)
tT62 = arrange (tT62, desc(-logFC))
tT62a = head(tT62[,1:7], 10)
tT63a <- rbind (tT61a, tT62a)
readr::write_csv(tT63a, file = file.path (results.dir, "nasal_top20.DEGs_mild.vs.controls-P0.001.csv"))


########Create data.frame with normalized expression for genes with significant FC
d <- normalized_counts
ENSGid <- rownames(d)
rownames(d) <- NULL
normalized_counts2 <- cbind (ENSGid, d)

normalized_counts.P0.001 <- merge(tT6, normalized_counts2, by="ENSGid" )
row.names(normalized_counts.P0.001) <- normalized_counts.P0.001$hgnc_symbol
normalized_counts.P0.001$LR <- normalized_counts.P0.001$ENSGid <- normalized_counts.P0.001$logFC <- normalized_counts.P0.001$logCPM <- normalized_counts.P0.001$FDR <- normalized_counts.P0.001$PValue <- normalized_counts.P0.001$hgnc_symbol <- NULL 
readr::write_csv(normalized_counts.P0.001, file = file.path (results.dir, "nasal-normalized_counts.P0.001.mildvscontrol.csv"))


#####################heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

table (s1$group)
annotation_col = data.frame(Group = factor(c(rep("Non-COPD", 22), rep("Mild COPD", 24))))
rownames(annotation_col)

data <- normalized_counts.P0.001 
colnames(data)
rownames(annotation_col) <- colnames(data)
head(annotation_col)

data1 <- data [,annotation_col$Group == "Non-COPD"]
data2 <- data [,annotation_col$Group == "Mild COPD"]

HC1 <- hclust (dist(t(data1), method = "manhattan"), method = "complete", members = NULL)
HC2 <- hclust (dist(t(data2), method = "manhattan"), method = "complete", members = NULL)

data1 <- data1[,HC1$order]
data2 <- data2[,HC2$order]
datax <- cbind(data1, data2)

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
              cellwidth = 15,cellheight = 30, 
              cluster_cols = F, 
              cluster_rows = T,
              show_rownames = T, 
              show_colnames = F,
              legend = T,   
              legend_breaks = -4:4, 
              fontsize = 11,
              fontsize_row = 12, 
              fontsize_col = 10,
              filename = file.path (results.dir.img, "Heatmap_nasal_mild.vs.control_P0.001.png") ,dpi=600)

