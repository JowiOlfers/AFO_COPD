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

#array1= contains the groups mild COPD and control (n=46)
#array2= contains the groups severe COPD and control (n=135)
#array3= contains the groups severe COPD and control (n=137)




# EdgeR: disease

################################### array2 - Severe COPD vs. controls
expressionDataToUse= array2
samplesToUse=s2


##########create factors
group=as.factor(samplesToUse$group)
smoking=as.factor(samplesToUse$smoking)
age=as.numeric(samplesToUse$age)
gender=as.factor(samplesToUse$gender)
packyears=as.numeric(samplesToUse$packyears)
years.cessation=as.numeric(samplesToUse$years.of.cessation)


###replace years.cessation NA with average value of corresponding group
df = data.frame(years.cessation, group)
df0 <- subset(df, group == "0")
df2 <- subset(df, group == "2")


#calculate average of group 0 and 2
dfx <- df[!(is.na(df$years.cessation)), ]

dfx0 <- subset(dfx, group == "0")
mean.cessation_0 <- mean(dfx0$years.cessation)
print (mean.cessation_0)

dfx2 <- subset(dfx, group == "2")
mean.cessation_2 <- mean(dfx2$years.cessation)
print (mean.cessation_2)


#replace NA with average cessation of group
df0[is.na(df0)]<- mean.cessation_0
df2[is.na(df2)]<- mean.cessation_2
df_na.rm <- rbind(df0,df2)

#create years.cessation factor without NA
years.of.cessation <- as.numeric(df_na.rm$years.cessation)



######List and filter on expression
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
write.table(normalized_counts,file="nasal-log2CPM.Severevscontrols.txt",sep="\t",quote=F)
write.table(normalized_counts,file="nasal-log2CPM.Severevscontrols.csv",sep=",")

plotMDS(nasaldata_DGEL)





##############Differential expression analysis
design <- model.matrix(~group + age + gender + packyears + years.of.cessation)
nasaldata_DGEL <- estimateDisp(nasaldata_DGEL, design)
plotBCV(nasaldata_DGEL)


fit <- glmFit(nasaldata_DGEL, design)
lrt <- glmLRT(fit, coef = 2)
results <- topTags(lrt,n=nrow(nasaldata_DGEL))


#plotMeanVar(nasaldata_DGEL, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE)


#add gene names
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
write.table(tT1, "nasal-severeCOPD.vs.control.txt")

tT2=tT1[which(tT1$FDR<0.05), ]
names(tT2)[1] <- "ENSGid"
write.table(tT2, "nasal-severeCOPD.vs.control_FDR0.05.txt")

tT3=tT1[which(tT1$FDR<0.01), ]
names(tT3)[1] <- "ENSGid"
write.table(tT3, "nasal-severeCOPD.vs.control_FDR0.01.txt")

tT4=tT1[which(tT1$PValue<0.05),]
names(tT4)[1] <- "ENSGid"
write.table(tT4, "nasal-severeCOPD.vs.control_P0.05.txt")

tT5=tT1[which(tT1$PValue<0.01),]
names(tT5)[1] <- "ENSGid"
write.table(tT5, "nasal-severeCOPD.vs.control_P0.01.txt")

tT6=tT1[which(tT1$PValue<0.001),]
names(tT6)[1] <- "ENSGid"
write.table(tT6, "nasal-severeCOPD.vs.control_P0.001.txt")




###########volcano plot
library(ggpubr)
library(ggthemes)


################### FDR 0.01 selected genes
tT1$logFDR <- -log10(tT1$FDR)
tT1$Group <- "No diff"
tT1$Group[which((tT1$FDR < 0.01)&(tT1$logFC > 1))] = "Up"
tT1$Group[which((tT1$FDR < 0.01)&(tT1$logFC < -1))] = "Down"

#amount of DEG with FDR < 0.01 and LogFC > 1 of < -1
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
  geom_vline(xintercept = c(-1, 1), linetype="dashed")
p

ggsave(filename = "Volcano-WC-severevscontrol-FDR0.01.png", width=21,height=21,units="cm",dpi=600)
ggsave(filename = "gene_diff.severevcontrol-FDR0.01.pdf", width=21,height=21,units="cm")



###################### Pvalue selected genes - P<0.001

tT1$logPValue <- -log10(tT1$PValue)
tT1$GroupP <- "No diff"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC > 1))] = "Up"
tT1$GroupP [which((tT1$PValue < 0.001)&(tT1$logFC < -1))] = "Down"

#amount of DEG with Pvalue < 0.001 and LogFC > 0 of < 0
table(tT1$Group)

tT1$labelP = ""
tT1 <- tT1[order(tT1$PValue),]

#select that the top 10 DEG show names in the plot
up_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Up")],10)
down_geneP <- head(tT1$hgnc_symbol[which(tT1$GroupP == "Down")],10)

#present top 10 up and down regulated genes 
#up_geneP
#down_geneP

tT1_geneP <- c(as.character(up_geneP), as.character(down_geneP))
tT1$labelP[match(tT1_geneP, tT1$hgnc_symbol)] <- tT1_geneP

pp <- ggscatter(tT1, x = "logFC", y = "logPValue", 
                color = "GroupP", 
                palette = c("#2f5688","#BBBBBB","#CC0000"), 
                size = 1,
                label = tT1$labelP, 
                font.label = 8, 
                repel = T,
                xlab = "log2FoldChange", ylab = "-log10(PValue)") + theme_base() +
  geom_hline(yintercept = -log10(0.001),linetype="dashed")+
  geom_vline(xintercept = c(-1, 1),linetype="dashed")
pp
ggsave(filename = "Volcano-WC-severevscontrol-PValue.png", width=25,height=25,units="cm",dpi=600)
ggsave(filename = "gene_diff-severevscontrol-PValue.pdf", width=25,height=25,units="cm")





######Create data.frame with normalized expression for genes with FDR sig. DEGs

library(dplyr)
#select up and down regulated genes with FDR < 0.05 and FC > 2 or < -2
tT21=tT2[which(tT2$logFC >= 1),]
tT22=tT2[which(tT2$logFC <= -1), ]
tT23 <- rbind (tT21, tT22)


#select top 20 - FDR < 0.05
tT21 = arrange (tT21, desc(logFC))
tT21a = head(tT21[,1:7], 10)
tT22 = arrange (tT22, desc(-logFC))
tT22a = head(tT22[,1:7], 10)
tT23a <- rbind (tT21a, tT22a)
write.table(tT23a, "Nasal_top20.DEGs_severe.vs.controls-FDR0.05.csv")

#select top 10 - FDR < 0.05
tT21 = arrange (tT21, desc(logFC))
tT221a = head(tT21[,1:7], 5)
tT22 = arrange (tT22, desc(-logFC))
tT222a = head(tT22[,1:7], 5)
tT223a <- rbind (tT221a, tT222a)
write.table(tT223a, "Nasal_top10.DEGs_severe.vs.controls-FDR0.05.csv")



#select up and down regulated genes with FDR < 0.01 and FC > 2 or < -2
tT31=tT3[which(tT3$logFC >= 1),] 
dim (tT31)
write.table(tT31, "nasal_up-DEGs_severe.vs.controls-FDR0.01FC.csv")

tT32=tT3[which(tT3$logFC <= -1), ]
dim (tT32)
write.table(tT32, "nasal_down-DEGs_severe.vs.controls-FDR0.01FC.csv")

tT33 <- rbind (tT31, tT32)
dim (tT33)
write.table(tT33, "nasal_allDEGs_severe.vs.controls-FDR0.01FC.csv")


#select top 20 FDR < 0.01 and FC > 2 or < -2
tT31 = arrange (tT31, desc(logFC))
tT31a = head(tT31[,1:7], 10)
tT32 = arrange (tT32, desc(-logFC))
tT32a = head(tT32[,1:7], 10)
tT33a <- rbind (tT31a, tT32a)
write.table(tT33a, "Nasal_top20.DEGs_severe.vs.controls-FDR0.01.csv")

#select top 10 FDR < 0.01 and FC > 2 or < -2
tT31 = arrange (tT31, desc(logFC))
tT321a = head(tT31[,1:7], 5)
tT32 = arrange (tT32, desc(-logFC))
tT322a = head(tT32[,1:7], 5)
tT323a <- rbind (tT321a, tT322a)
write.table(tT323a, "Nasal_top10.DEGs_severe.vs.controls-FDR0.01.csv")






##############Combine normalized raw counts with expression data
d <- normalized_counts
ENSGid <- rownames(d)
rownames(d) <- NULL
normalized_counts2 <- cbind (ENSGid, d)

normalized_counts.FDR <- merge(tT33, normalized_counts2, by="ENSGid" )
normalized_counts.FDRa <- merge(tT33a, normalized_counts2, by="ENSGid" )


######normalized_counts.FDR
#error: row.names duplicate 'POLR2J4'and. 'TBCE'.
normalized_counts.FDR [normalized_counts.FDR$hgnc_symbol == "TBCE",]
normalized_counts.FDR$hgnc_symbol [normalized_counts.FDR$ENSGid == "ENSG00000284770"] <- "TBCE2"
normalized_counts.FDR [normalized_counts.FDR$hgnc_symbol == "POLR2J4",]
normalized_counts.FDR$hgnc_symbol [normalized_counts.FDR$ENSGid == "ENSG00000272655"] <- "POLR2J42"
#gene name TBCE2 is equal to TBCE and gene name POLR2J42 = POLR2J4

row.names(normalized_counts.FDR) <- normalized_counts.FDR$hgnc_symbol
normalized_counts.FDR$LR <- normalized_counts.FDR$ENSGid <- normalized_counts.FDR$logFC <- normalized_counts.FDR$logCPM <- normalized_counts.FDR$FDR <- normalized_counts.FDR$PValue <- normalized_counts.FDR$hgnc_symbol <- NULL
write.table(normalized_counts.FDR,file="nasal-normalized_counts.FDR-0.01.severevscontrol.csv",sep=",")

######normalized_counts.FDRa - tT23a
#error: row.names duplicate 'POLR2J4'and. 'TBCE'.
normalized_counts.FDRa [normalized_counts.FDRa$hgnc_symbol == "TBCE",]
normalized_counts.FDRa$hgnc_symbol [normalized_counts.FDRa$ENSGid == "ENSG00000284770"] <- "TBCE2"
normalized_counts.FDRa [normalized_counts.FDRa$hgnc_symbol == "POLR2J4",]
normalized_counts.FDRa$hgnc_symbol [normalized_counts.FDRa$ENSGid == "ENSG00000272655"] <- "POLR2J42"
#gene name TBCE2 is equal to TBCE and gene name POLR2J42 = POLR2J4

row.names(normalized_counts.FDRa) <- normalized_counts.FDRa$hgnc_symbol
normalized_counts.FDRa$LR <- normalized_counts.FDRa$ENSGid <- normalized_counts.FDRa$logFC <- normalized_counts.FDRa$logCPM <- normalized_counts.FDRa$FDR <- normalized_counts.FDRa$PValue <- normalized_counts.FDRa$hgnc_symbol <- NULL
write.table(normalized_counts.FDRa,file="nasal-normalized_counts.FDRa-0.01.severevscontrol.csv",sep=",")





############################ heatmap

library(pheatmap)
library(RColorBrewer)
library(ggplot2)

data <- read.table("nasal-normalized_counts.FDRa-0.01.severevscontrol.csv",sep=",")

table (s2$group)
annotation_col = data.frame(Group = factor(c(rep("Non-COPD", 22),rep("Severe COPD", 113))))
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
              cellwidth = 9,cellheight = 20, 
              cluster_cols = F, 
              cluster_rows = T,
              show_rownames = T, 
              show_colnames = F,
              legend = T,   
              legend_breaks = -4:4, 
              fontsize = 12,
              fontsize_row = 12, 
              fontsize_col = 10,
              filename = "Heatmap_nasal_severevscontrol_FDR-0.01.png",dpi=600)

