library("DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
####

library("edgeR")
library("ggplot2")
library("ggrepel")
library("biomaRt")
library("heatmap3")
library("gplots")
library("ggfortify")

array=read.table("bronchial_rawcounts.csv",sep = ";",header=T,check.names=F)
s=read.csv("bronchial clinical data.csv", sep=";",header=T,check.names=F)

row.names(array)=array[,1]
array=array[,-1]
array=array[,s$SAMID]

array1=array[,which(s$mild.vs.control==1)]
s1=s[which(s$mild.vs.control==1),]

array2=array[,which(s$severe.vs.control==1)]
s2=s[which(s$severe.vs.control==1),]

array3=array[,which(s$severe.vs.mild==1)]
s3=s[which(s$severe.vs.mild==1),]

#array1= contains groups mild COPD and controls (n=46)
#array2= contains groups severe COPD and controls (n=146)
#array3= contains groups severe and mild COPD (n=146)


# EdgeR: disease


################################## array1 - Mild COPD vs. non-COPD
expressionDataToUse= array1
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

#calculate average of group 0
dfx <- df[!(is.na(df$years.cessation)), ]
dfx0 <- subset(dfx, group == "0")
mean.cessation_0 <- mean(dfx0$years.cessation)
print (mean.cessation_0)

#replace NA with average cessation of group
df[is.na(df)]<- mean.cessation_0

#create years.cessation factor without NA
years.of.cessation <- as.numeric(df$years.cessation)
print(years.of.cessation)


#######
total_bronchialdata_DGEL <- DGEList(expressionDataToUse)
keep <- filterByExpr(total_bronchialdata_DGEL,group=group)
bronchialdata_DGEL <- total_bronchialdata_DGEL[keep,, keep.lib.sizes=FALSE]
dim(bronchialdata_DGEL)
bronchialdata_DGEL$samples


################Normalization
bronchialdata_DGEL <- calcNormFactors(bronchialdata_DGEL, method="TMM")
bronchialdata_DGEL$samples
plotMDS(bronchialdata_DGEL)


#count per million (cpm) read, also called normalized count
normalized_counts <- cpm(bronchialdata_DGEL, log = TRUE)
write.table(normalized_counts,file="bronchial-log2CPM.mildvscontrol.txt",sep="\t",quote=F)
write.table(normalized_counts,file="bronchial-log2CPM.mildvscontrol.csv",sep=",")



##############Differential expression analysis
design <- model.matrix(~group + age + gender + packyears + years.of.cessation)
bronchialdata_DGEL <- estimateDisp(bronchialdata_DGEL, design)
plotBCV(bronchialdata_DGEL)

fit <- glmFit(bronchialdata_DGEL, design)
lrt <- glmLRT(fit, coef = 2)
results <- topTags(lrt,n=nrow(bronchialdata_DGEL))

#plotMeanVar(nasaldata_DGEL, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)



tT1=cbind(results$table,row.names(results$table))

mart = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl")

genesID = as.character(rownames(tT1))
G_list = getBM(
  filters = "ensembl_gene_id",
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol"
  ),
  values = genesID,
  mart = mart
)
G_list=as.matrix(G_list)


tT1=merge(
  x = tT1,
  y = G_list,
  all.x = T,
  all.y = F,
  by.x = "row.names(results$table)",
  by.y ="ensembl_gene_id")


rownames(tT1)=tT1[,1]
#standard an error due to DUPLICATES:'ENSG00000254876', 'ENSG00000276085
names(tT1)[1] <- "ENSGid"
tT1 <- tT1 [!(is.na(tT1$hgnc_symbol) | tT1$hgnc_symbol == ""), ]
tT1 <- tT1 [!(tT1$hgnc_symbol == "CCL3L3"), ]
write.table(tT1, "bronchial-mildCOPD.vs.control.txt")

tT2=tT1[which(tT1$FDR<0.05),]
names(tT2)[1] <- "ENSGid"
write.table(tT2, "bronchial-mildCOPD.vs.control_FDR0.05.txt")

tT3=tT1[which(tT1$FDR<0.01),]
names(tT3)[1] <- "ENSGid"
write.table(tT3, "bronchial-mildCOPD.vs.control_FDR0.01.txt")

tT4=tT1[which(tT1$PValue<0.01),]
names(tT4)[1] <- "ENSGid"
write.table(tT4, "bronchial-mildCOPD.vs.control_P0.05.txt")

tT5=tT1[which(tT1$PValue<0.001),]
names(tT5)[1] <- "ENSGid"
write.table(tT5, "bronchial-mildCOPD.vs.control_P0.001.txt")





########## Selection and presentation significant up and down regulated genes

#volcano plot
library(ggpubr)
library(ggthemes)



################### FDR 0.01 selected genes
tT1$logFDR <- -log10(tT1$FDR)
tT1$Group <- "No diff"
tT1$Group[which((tT1$FDR < 0.01)&(tT1$logFC > 1))] = "Up"
tT1$Group[which((tT1$FDR < 0.01)&(tT1$logFC < -1))] = "Down"

#amount of DEG with FDR < 0.01 and LogFC > 1 or < -1
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
  geom_hline(yintercept = -log10(0.01),linetype="dashed")+
  geom_vline(xintercept = c(-1,1),linetype="dashed")
p

ggsave(filename = "Volcano-bronchial-mildvscontrol-FDR0.01.png", width=21,height=21,units="cm",dpi=600)
ggsave(filename = "gene_diff-bronchial-mildvscontrol-FDR0.01.pdf", width=21,height=21,units="cm")




###################### Pvalue selected genes - P<0.001

tT1$logPValue <- -log10(tT1$PValue)
tT1$GroupP <- "No diff"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC > 1))] = "Up"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC < -1))] = "Down"

#amount of DEG with Pvalue < 0.001 and LogFC > 0 of < 0
table(tT1$GroupP)

tT1$labelP = ""
tT1 <- tT1[order(tT1$PValue),]

#select that the top 10 DEG show names in the plot
up_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Up")],10)
down_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Down")],10)

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
  geom_vline(xintercept = c(-1,1),linetype="dashed")

p1
ggsave(filename = "Volcano-bronchial-mildvscontrol-PValue.png", width=25,height=25,units="cm",dpi=600)
ggsave(filename = "gene_diff-bronchial-mildvscontrol-PValue.pdf", width=25,height=25,units="cm")




######################

#select up and down regulated genes with FDR < 0.01 and FC > 2 or < -2

tT31_all=tT3[which(tT3$logFC >= 1),] 
#dim (tT31_all)
tT32_all=tT3[which(tT3$logFC <= -1), ]
#dim (tT32_all)
tT33_all <- rbind (tT31_all, tT32_all)
#dim (tT33_all)
write.table(tT33_all, "bronchial_allDEGs_mild.vs.controls-FDR0.01FC.csv")



#exclude bronchial severe vs control DEGs (FDR 0.01 FC > 2 or < -2)
severevscontrol <- read.table("bronchial_allDEGs_severe.vs.controls-FDR0.01FC.csv")
severe_DEGs <- severevscontrol$hgnc_symbol
overlap <- intersect(tT33_all$hgnc_symbol, severe_DEGs)
tT33 <- tT33_all[ ! tT33_all$hgnc_symbol %in% overlap,]
dim (tT33)
write.table(tT33, "bronchial_all-excluded-DEGs_mild.vs.controls-FDR0.01FC.csv")


tT31 =tT33 [which(tT33$logFC >= 1),] 
dim (tT31)
write.table(tT31, "bronchial_excluded-up-DEGs_mild.vs.controls-FDR0.01FC.csv")


tT32 =tT33[which(tT33$logFC <= -1), ]
dim (tT32)
write.table(tT32, "bronchial_excluded-down-DEGs_mild.vs.controls-FDR0.01FC.csv")




#select top 20 up and down regulated genes with FDR < 0.01 and FC > 2 or < -2
library(dplyr)

tT31 = arrange (tT31, desc(logFC))
tT31a = head(tT31[,1:7], 10)
tT32 = arrange (tT32, desc(-logFC))
tT32a = head(tT32[,1:7], 10)
tT33a <- rbind (tT31a, tT32a)
write.table(tT33a, "bronchial_top20.DEGs_mild.vs.controls-FDR0.01.csv")


#select top 10 up and down regulated genes with FDR < 0.01 and 
tT31 = arrange (tT31, desc(logFC))
tT321a = head(tT31[,1:7], 5)
tT32 = arrange (tT32, desc(-logFC))
tT322a = head(tT32[,1:7], 5)
tT323a <- rbind (tT321a, tT322a)
write.table(tT323a, "bronchial_top10.DEGs_mild.vs.controls-FDR0.01.csv")




########Create data.frame with normalized expression for genes with significant FC
d <- normalized_counts
ENSGid <- rownames(d)
rownames(d) <- NULL
normalized_counts2 <- cbind (ENSGid, d)

normalized_counts.FDR <- merge(tT33, normalized_counts2, by="ENSGid" )
normalized_counts.FDRa <- merge(tT33a, normalized_counts2, by="ENSGid" )

row.names(normalized_counts.FDR) <- normalized_counts.FDR$hgnc_symbol
normalized_counts.FDR$LR <- normalized_counts.FDR$ENSGid <- normalized_counts.FDR$logFC <- normalized_counts.FDR$logCPM <- normalized_counts.FDR$FDR <- normalized_counts.FDR$PValue <- normalized_counts.FDR$hgnc_symbol <- NULL
write.table(normalized_counts.FDR,file="bronchial-normalized_counts.FDR0.01.mildvscontrol.csv",sep=",")

row.names(normalized_counts.FDRa) <- normalized_counts.FDRa$hgnc_symbol
normalized_counts.FDRa$LR <- normalized_counts.FDRa$ENSGid <- normalized_counts.FDRa$logFC <- normalized_counts.FDRa$logCPM <- normalized_counts.FDRa$FDR <- normalized_counts.FDRa$PValue <- normalized_counts.FDRa$hgnc_symbol <- NULL
write.table(normalized_counts.FDRa,file="bronchial-normalized_counts.FDRa0.01.mildvscontrol.csv",sep=",")





#####################heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

data <- read.table("bronchial-normalized_counts.FDRa0.01.mildvscontrol.csv",sep=",")

table (s1$group)
annotation_col = data.frame(Group = factor(c(rep("Non-COPD", 23), rep("Mild COPD", 23))))
rownames(annotation_col)
colnames(data)
rownames(annotation_col) <- colnames(data)
head(annotation_col)


ph <- pheatmap(data,
              scale="row",
              annotation_col = annotation_col,
              annotation_legend = T,
              annotation_names_col = F,
              number_format="%.2e",
              border="white",  
              color=colorRampPalette(c("navy", "white", "red"))(100),
              border_color=NA, 
              cellwidth = 18,cellheight = 11, 
              cluster_cols = F, 
              cluster_rows = T,
              show_rownames = T, 
              show_colnames = F,
              legend = T,   
              legend_breaks = -4:4, 
              fontsize = 10,
              fontsize_row = 10, 
              fontsize_col = 10,
              filename = "Heatmap_bronchial_mild.vs.control_FDR0.01.png",dpi=600)

