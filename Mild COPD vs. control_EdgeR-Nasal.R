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



array=read.table("nasal rawcounts.csv",sep = ";",header=T,check.names=F)
s=read.csv("nasal clinical data.csv", sep=";",header=T,check.names=F)

row.names(array)=array[,1]
array=array[,-1]
array=array[,s$SAMID]

array1=array[,which(s$control.vs.mild==1)]
s1=s[which(s$control.vs.mild==1),]

array2=array[,which(s$control.vs.severe==1)]
s2=s[which(s$control.vs.severe==1),]

array3=array[,which(s$mild.vs.severe==1)]
s3=s[which(s$mild.vs.severe==1),]

#array1= contains groups control and mild COPD (n=46)
#array2= contains groups control and severe COPD (n=136)
#array3= contains groups mild and severe COPD (n=138)




# EdgeR: disease

################################## array1 - Mild vs. non-COPD
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
total_nasaldata_DGEL <- DGEList(expressionDataToUse)
keep <- filterByExpr(total_nasaldata_DGEL,group=group)
nasaldata_DGEL <- total_nasaldata_DGEL[keep,, keep.lib.sizes=FALSE]
dim(nasaldata_DGEL)
nasaldata_DGEL$samples



################Normalization
nasaldata_DGEL <- calcNormFactors(nasaldata_DGEL, method="TMM")
nasaldata_DGEL$samples

#count per million (cpm) read, also called normalized count
normalized_counts <- cpm(nasaldata_DGEL, log = TRUE)
write.table(normalized_counts,file="nasal-log2CPM.mildvscontrol.txt",sep="\t",quote=F)
write.table(normalized_counts,file="nasal-log2CPM.mildvscontrol.csv",sep=",")

plotMDS(nasaldata_DGEL)





##############Differential expression analysis
design <- model.matrix(~group + age + gender + packyears + years.of.cessation)
nasaldata_DGEL <- estimateDisp(nasaldata_DGEL, design)
plotBCV(nasaldata_DGEL)

fit <- glmFit(nasaldata_DGEL, design)
lrt <- glmLRT(fit, coef = 2)
results <- topTags(lrt,n=nrow(nasaldata_DGEL))

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
write.table(tT1, "Nasal-mildCOPD.vs.control.txt")

tT2=tT1[which(tT1$FDR<0.05),]
names(tT2)[1] <- "ENSGid"
write.table(tT2, "Nasal-mildCOPD.vs.control_FDR.txt")

#FDR 0.25 again 0 DEG

tT3=tT1[which(tT1$FDR<0.5),]
names(tT3)[1] <- "ENSGid"
write.table(tT3, "Nasal-mildCOPD.vs.control_FDR0.5.txt")

tT4=tT1[which(tT1$PValue<0.05),]
names(tT4)[1] <- "ENSGid"
write.table(tT4, "Nasal-mildCOPD.vs.control_P0.05.txt")

tT5=tT1[which(tT1$PValue<0.01),]
names(tT5)[1] <- "ENSGid"
write.table(tT5, "Nasal-mildCOPD.vs.control_P0.01.txt")

tT6=tT1[which(tT1$PValue<0.001),]
names(tT6)[1] <- "ENSGid"
write.table(tT6, "Nasal-mildCOPD.vs.control_P0.001.txt")





########## Selection and presentation significant up and down regulated genes

#volcano plot
library(ggpubr)
library(ggthemes)



################### FDR 0.5 selected genes
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
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")
p
ggsave(filename = "Volcano-WC-mildvscontrol-FDR-0.5.png", width=21,height=21,units="cm",dpi=600)
ggsave(filename = "gene_diff-mildvscontrol-FDR-0.5.pdf", width=21,height=21,units="cm")




###################### Pvalue selected genes - P<0.001

tT1$logPValue <- -log10(tT1$PValue)
tT1$GroupP <- "No diff"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC > 0))] = "Up"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC < 0))] = "Down"

#amount of DEG with Pvalue < 0.001 and LogFC > 0 of < 0
table(tT1$GroupP)

tT1$labelP = ""
tT1 <- tT1[order(tT1$PValue),]

#select that the top 10 DEG show names in the plot
up_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Up")],10)
down_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Down")],12)

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
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")
p1
ggsave(filename = "Volcano-WC-mildvscontrol-PValue.png", width=25,height=25,units="cm",dpi=600)
ggsave(filename = "gene_diff-mildvscontrol-PValue.pdf", width=25,height=25,units="cm")



#select up and down regulated genes with P < 0.001 and FC > 2 or < -2
tT61=tT6[which(tT6$logFC > 0),] 
dim (tT61)
write.table(tT61, "nasal_up-DEGs_mild.vs.controls-P0.001FC.csv")

tT62=tT6[which(tT6$logFC < 0), ]
dim (tT62)
write.table(tT63, "nasal_down-DEGs_mild.vs.controls-P0.001FC.csv")

tT63 <- rbind (tT61, tT62)
dim (tT63)
write.table(tT63, "nasal_allDEGs_mild.vs.controls-P0.001FC.csv")

#select top 20
library(dplyr)

tT61 = arrange (tT61, desc(logFC))
tT61a = head(tT61[,1:7], 10)
tT62 = arrange (tT62, desc(-logFC))
tT62a = head(tT62[,1:7], 10)
tT63a <- rbind (tT61a, tT62a)
write.table(tT63a, "Nasal_top20.DEGs_mild.vs.controls-P0.001.csv")



########Create data.frame with normalized expression for genes with significant FC
d <- normalized_counts
ENSGid <- rownames(d)
rownames(d) <- NULL
normalized_counts2 <- cbind (ENSGid, d)

normalized_counts.P0.001 <- merge(tT6, normalized_counts2, by="ENSGid" )
row.names(normalized_counts.P0.001) <- normalized_counts.P0.001$hgnc_symbol
normalized_counts.P0.001$LR <- normalized_counts.P0.001$ENSGid <- normalized_counts.P0.001$logFC <- normalized_counts.P0.001$logCPM <- normalized_counts.P0.001$FDR <- normalized_counts.P0.001$PValue <- normalized_counts.P0.001$hgnc_symbol <- NULL 

write.table(normalized_counts.P0.001,file="nasal-normalized_counts.P0.001.mildvscontrol.csv",sep=",")



#####################heatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

data <- read.table("nasal-normalized_counts.P0.001.mildvscontrol.csv",sep=",")

table (s1$group)
annotation_col = data.frame(Group = factor(c(rep("Non-COPD", 22), rep("Mild COPD", 24))))
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
              filename = "Heatmap_nasal_mild.vs.control_P0.001.png",dpi=600)

