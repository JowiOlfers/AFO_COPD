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
library (tibble)

#Directories 
setwd("~/HVHL/Stages/Afstudeerstage_2022/DE analysis/Nasal DEG_DE in bronchi")
library(here)

data.dir <- file.path(".", "data")
results.dir <- file.path(here::here(), "results")
dir.create(results.dir, recursive=TRUE)
results.dir.img <- file.path(results.dir, "img")

x <- readr::read_csv(
  file.path(data.dir, "bronchial rawcounts.csv"))
array <- data.frame (x)

y <- readr::read_csv(
  file.path(data.dir, "bronchial clinical data.csv"))
s <- data.frame (y)

#s1= contains groups mild COPD and controls (n=46)
#s2= contains groups severe COPD and controls (n=145)
#s3= contains groups severe and mild COPD (n=145)



row.names(array)=array[,1]
array=array[,-1]
array=array[,s$SAMID]

array <- rownames_to_column(array, var = "ensembl_gene_id")

#Ensembl - hgnc_symbols
ensembl.dir <- file.path("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/CD")
source(file.path(ensembl.dir, "_helpers.r"), verbose=TRUE, chdir=TRUE)

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

DEG <- left_join(array, geneData, by = "ensembl_gene_id")

unique.dir <- file.path("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis/Nasal data/unique selection")
DEG.nasal <- readr::read_csv(
  file.path(unique.dir, "nasal_unique-DEGs_severeCOPD-FDR0.05FC2.csv"))

DEG <- DEG [DEG$hgnc_symbol %in% DEG.nasal$hgnc_symbol, ]
rownames(DEG) <- DEG$hgnc_symbol
DEG <- DEG [,-1]

#check if these are actually overlapping
DEG_list <- DEG$hgnc_symbol
DEG_nasal_list <- DEG.nasal$hgnc_symbol


DEG1=DEG[,which(s$mild.vs.control==1)]
s1=s[which(s$mild.vs.control==1),]

DEG2=DEG[,which(s$severe.vs.control==1)]
s2=s[which(s$severe.vs.control==1),]

DEG3=DEG[,which(s$severe.vs.mild==1)]
s3=s[which(s$severe.vs.mild==1),]





# EdgeR: disease

###########################DEG1 mild COPD vs control
expressionDataToUse= DEG1
samplesToUse=s1

#create factors
group=as.factor(samplesToUse$group)
smoking=as.factor(samplesToUse$smoking)
age=as.numeric(samplesToUse$age)
gender=as.factor(samplesToUse$gender)
packyears=as.numeric(samplesToUse$packyears)
years.cessation=as.numeric(samplesToUse$years.of.cessation)

###replace years.cessation NA with average value of corresponding group
df = data.frame(years.cessation, group)
df0 <- subset(df, group == "0")
df1 <- subset(df, group == "1")


#calculate average of group 0 and 1
dfx <- df[!(is.na(df$years.cessation)), ]

dfx0 <- subset(dfx, group == "0")
mean.cessation_0 <- mean(dfx0$years.cessation)
print (mean.cessation_0)

dfx1 <- subset(dfx, group == "1")
mean.cessation_1 <- mean(dfx1$years.cessation)
print (mean.cessation_1)

#replace NA with average cessation of group
df0[is.na(df0)]<- mean.cessation_0
df1[is.na(df1)]<- mean.cessation_1
df_na.rm <- rbind(df0,df1)
years.of.cessation <- as.numeric(df_na.rm$years.cessation)



#List and filter on expression
total_bronchialdata_DGEL <- DGEList(expressionDataToUse)
keep <- filterByExpr(total_bronchialdata_DGEL,group=group)
bronchialdata_DGEL <- total_bronchialdata_DGEL[keep,, keep.lib.sizes=FALSE]
dim(bronchialdata_DGEL)
bronchialdata_DGEL$samples

##############Normalization - TMM
bronchialdata_DGEL <- calcNormFactors(bronchialdata_DGEL, method="TMM")
bronchialdata_DGEL$samples
#plotMDS(bronchialdata_DGEL)

#count per million (cpm) read
normalized_counts <- cpm(bronchialdata_DGEL, log = TRUE)
normalized_counts <- as.data.frame (normalized_counts)
write.csv(normalized_counts, 
                 file = file.path (results.dir, "bronchial-log2CPM.mildvscontrol.csv"))
write.csv(normalized_counts, 
                 file = file.path (results.dir, "bronchial-log2CPM.mildvscontrol.txt"))



##############Differential expression analysis
design <- model.matrix(~group + age + gender + packyears + years.of.cessation)
bronchialdata_DGEL <- estimateDisp(bronchialdata_DGEL, design)
#plotBCV(bronchialdata_DGEL)

fit <- glmFit(bronchialdata_DGEL, design)
lrt <- glmLRT(fit, coef = 2)
results <- topTags(lrt,n=nrow(bronchialdata_DGEL))
#plotMeanVar(bronchialdata_DGEL, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)


tT1 <- results$table
tT1$hgnc_symbol <- row.names(tT1)
readr::write_csv(tT1, file = file.path (results.dir, "DEG.bronchial-mildCOPD.vs.control.txt"))

tT2=tT1[which(tT1$PValue<0.01),]
readr::write_csv(tT2, file = file.path (results.dir, "DEG.bronchial-mildCOPD.vs.control-P0.01.txt"))

tT3=tT1[which(tT1$PValue<0.001),]
#readr::write_csv(tT3, file = file.path (results.dir, "DEG.bronchial-mildCOPD.vs.control-P0.001.txt"))


###### select up and down regulated genes with P < 0.01 and logFC < or > 0 #######
tT21=tT2[which(tT2$logFC > 0),]
dim (tT21)
readr::write_csv(tT21, file = file.path (results.dir, "DEG.bronchial_up-DEGs_mild.vs.control-P0.01.csv"))

tT22=tT2[which(tT2$logFC < 0), ]
dim (tT22)
readr::write_csv(tT22, file = file.path (results.dir, "DEG.bronchial_down-DEGs_mild.vs.control-P0.01.csv"))

tT23 <- rbind (tT21, tT22)
dim (tT23)
readr::write_csv(tT23, file = file.path (results.dir, "DEG.bronchial_all-DEGs_mild.vs.control-P0.01.csv"))

#select top 20 - P < 0.01
tT21 = arrange (tT21, desc(logFC))
tT21a = head(tT21[,1:6], 10)
tT22 = arrange (tT22, desc(-logFC))
tT22a = head(tT22[,1:6], 10)
tT23a <- rbind (tT21a, tT22a)
readr::write_csv(tT23a, file = file.path (results.dir, "DEG.bronchial_top20.DEGs_mild.vs.control-P0.01.csv"))

#select top 10 - P < 0.01
tT21 = arrange (tT21, desc(logFC))
tT221a = head(tT21[,1:6], 5)
tT22 = arrange (tT22, desc(-logFC))
tT222a = head(tT22[,1:6], 5)
tT223a <- rbind (tT221a, tT222a)
readr::write_csv(tT223a, file = file.path (results.dir, "DEG.bronchial_top10.DEGs_mild.vs.control-P0.01.csv"))



###### select up and down regulated genes with P < 0.001 and logFC < or > 0 ######
#tT31=tT3[which(tT3$logFC > 0),]
#dim (tT31)
#readr::write_csv(tT31, file = file.path (results.dir, "DEG.bronchial_up-DEGs_mild.vs.control-P0.001.csv"))
#
#tT32=tT3[which(tT3$logFC < 0), ]
#dim (tT32)
#readr::write_csv(tT32, file = file.path (results.dir, "DEG.bronchial_down-DEGs_mild.vs.control-P0.001.csv"))
#
#tT33 <- rbind (tT31, tT32)
#dim (tT33)
#readr::write_csv(tT33, file = file.path (results.dir, "DEG.bronchial_allDEGs_mild.vs.control-P0.001.csv"))
#
##select top 20 - P < 0.001
#tT31 = arrange (tT31, desc(logFC))
#tT31a = head(tT31[,1:6], 10)
#tT32 = arrange (tT32, desc(-logFC))
#tT32a = head(tT32[,1:6], 10)
#tT33a <- rbind (tT31a, tT32a)
#readr::write_csv(tT33a, file = file.path (results.dir, "DEG.bronchial_top20.DEGs_mild.vs.control-P0.001.csv"))
#
##select top 10 - P < 0.001
#tT31 = arrange (tT31, desc(logFC))
#tT321a = head(tT31[,1:6], 5)
#tT32 = arrange (tT32, desc(-logFC))
#tT322a = head(tT32[,1:6], 5)
#tT323a <- rbind (tT321a, tT322a)
#readr::write_csv(tT323a, file = file.path (results.dir, "DEG.bronchial_top10.DEGs_mild.vs.control-P0.001.csv"))



########## Selection and presentation significant up and down regulated genes
#volcano plot
library(ggpubr)
library(ggthemes)

################### Pvalue < 0.01 selected genes with logFC < or > 0
tT1$logPValue <- -log10(tT1$PValue)
tT1$Group <- "No diff"
tT1$Group[which((tT1$PValue < 0.01)&(tT1$logFC > 0))] = "Up"
tT1$Group[which((tT1$PValue < 0.01)&(tT1$logFC < 0))] = "Down"

#amount of DEG with P < 0.01 and LogFC > 0 < 
table(tT1$Group)

tT1$label = ""
tT1 <- tT1[order(tT1$PValue),]
#select that the top 10 DEG show names in the plot
up_gene <- head(tT1$hgnc_symbol[which(tT1$Group == "Up")],10)
down_gene <- head(tT1$hgnc_symbol[which(tT1$Group == "Down")],10)

#present top 10 up and down regulated genes 
up_gene
down_gene

tT1_gene <- c(as.character(up_gene), as.character(down_gene))
tT1$label[match(tT1_gene, tT1$hgnc_symbol)] <- tT1_gene

p <- ggscatter(tT1, x = "logFC", y = "logPValue", 
               color = "Group", 
               palette = c("#013E80","#BBBBBB","#CC4D00"),
               size = 1,
               label = tT1$label, 
               font.label = 8, 
               repel = T,
               xlab = "log2FoldChange", ylab = "-log10(PValue)") + theme_base() +
  geom_hline(yintercept = -log10(0.01),linetype="dashed")+
  geom_vline(xintercept = c(0,0),linetype="dashed")
p
ggsave (filename = file.path (results.dir.img, "Volcano-bronchial-mildvscontrol-P0.01.png"), width=25,height=25,units="cm",dpi=600 )
ggsave (filename = file.path (results.dir.img, "gene_diff-bronchial-mildvscontrol-P0.01.pdf"), width=25,height=25,units="cm")




########Create data.frame with normalized expression for genes with significant FC
d <- normalized_counts
hgnc_symbol <- rownames(d)
rownames(d) <- NULL
normalized_counts2 <- cbind (hgnc_symbol, d)

normalized_counts.P <- merge(tT23, normalized_counts2, by="hgnc_symbol" )
normalized_counts.Pa <- merge(tT23a, normalized_counts2, by="hgnc_symbol" )

row.names(normalized_counts.P) <- normalized_counts.P$hgnc_symbol
normalized_counts.P$LR <- normalized_counts.P$ENSGid <- normalized_counts.P$logFC <- normalized_counts.P$logCPM <- normalized_counts.P$FDR <- normalized_counts.P$PValue <- normalized_counts.P$hgnc_symbol <- NULL
write.csv(normalized_counts.P, file = file.path (results.dir, "DEG_bronchial-normalized_counts.P0.01.mildvscontrol.csv"))

row.names(normalized_counts.Pa) <- normalized_counts.Pa$hgnc_symbol
normalized_counts.Pa$LR <- normalized_counts.Pa$ENSGid <- normalized_counts.Pa$logFC <- normalized_counts.Pa$logCPM <- normalized_counts.Pa$FDR <- normalized_counts.Pa$PValue <- normalized_counts.Pa$hgnc_symbol <- NULL
write.csv(normalized_counts.Pa, file = file.path (results.dir, "DEG_bronchial-normalized_counts.top20-P0.01.mildvscontrol.csv"))



#####################heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

table (s1$group)
annotation_col = data.frame(Group = factor(c(rep("Non-COPD", 23),rep("Mild COPD", 23))))
rownames(annotation_col)

#normalized_counts.Pa = top 20
data <- normalized_counts.Pa
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
               color=colorRampPalette(c("#013E80", "white", "orangered"))(100),
               border_color=NA, 
               cellwidth =10,cellheight = 22, 
               cluster_cols = F, 
               cluster_rows = T,
               show_rownames = T, 
               show_colnames = F,
               legend = T,   
               legend_breaks = -4:4, 
               fontsize = 10,
               fontsize_row = 12, 
               fontsize_col = 10,
               filename = file.path (results.dir.img, "Heatmap_bronchial_mild.vs.control_P0.01.png") ,dpi=600)
