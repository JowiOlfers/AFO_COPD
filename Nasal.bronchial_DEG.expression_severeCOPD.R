library (limma)
library(edgeR)
library(DESeq2)
library(tidyverse)
library (dplyr)
library (tibble)

###working directories
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir.Nasal <- file.path(".", "Nasal data")
data.dir.N <- file.path(data.dir.Nasal, "data")
data.dir.Bronchial <- file.path(".", "Bronchial data")
data.dir.B <- file.path(data.dir.Bronchial, "data")



############### Clinical data ####################
### NASAL
MasterTable.nasal <- readr::read_csv(
  file.path(data.dir.N, "nasal clinical data.csv") 
) %>% 
  dplyr:: select(1:12)

annotation.nasal <- MasterTable.nasal %>%
  dplyr::select(
    SAMID,
    group,
  ) %>% 
  dplyr::mutate(
    group=factor(group))

annotation.nasal <- column_to_rownames(annotation.nasal, var = "SAMID")

### BRONCHIAL
MasterTable.bronchial <- readr::read_csv(
  file.path(data.dir.B, "bronchial clinical data.csv")
) %>% 
  dplyr:: select(1:12)

annotation.bronchial <- MasterTable.bronchial %>%
  dplyr::select(
    SAMID,
    group,
  ) %>% 
  dplyr::mutate(
    group=factor(group))

annotation.bronchial <- column_to_rownames(annotation.bronchial, var = "SAMID")

#create group factor
group.N=as.factor(MasterTable.nasal$group)
group.B=as.factor(MasterTable.bronchial$group)


############## Load normalized expression data #################
############## NASAL
x <- readr::read_csv(
  file.path(data.dir.N, "nasal rawcounts.csv"))
array.N <- data.frame (x)
row.names(array.N)=array.N[,1]
array.N=array.N[,-1]


##List and filter on expression
total_nasaldata <- DGEList(array.N)
keep <- filterByExpr(total_nasaldata,group=group.N)
nasaldata_DGEL <- total_nasaldata[keep,, keep.lib.sizes=FALSE]
dim(nasaldata_DGEL)
nasaldata_DGEL$samples

##Normalization - TMM
nasaldata_DGEL <- calcNormFactors(nasaldata_DGEL, method="TMM")
nasaldata_DGEL$samples
#plotMDS(nasaldata_DGEL)

#count per million (cpm) read
normalized_counts.N <- cpm(nasaldata_DGEL, log = TRUE)
normalized_counts.N <- as.data.frame (normalized_counts.N)


############## BRONCHIAL
x <- readr::read_csv(
  file.path(data.dir.B, "bronchial rawcounts.csv"))
array.B <- data.frame (x)
row.names(array.B)=array.B[,1]
array.B=array.B[,-1]

##List and filter on expression
total_bronchialdata <- DGEList(array.B)
keep <- filterByExpr(total_bronchialdata,group=group.B)
bronchialdata_DGEL <- total_bronchialdata[keep,, keep.lib.sizes=FALSE]
dim(bronchialdata_DGEL)
bronchialdata_DGEL$samples

##Normalization - TMM
bronchialdata_DGEL <- calcNormFactors(bronchialdata_DGEL, method="TMM")
bronchialdata_DGEL$samples
#plotMDS(bronchialdata_DGEL)

#count per million (cpm) read
normalized_counts.B <- cpm(bronchialdata_DGEL, log = TRUE)
normalized_counts.B <- as.data.frame (normalized_counts.B)


###save normalized counts nasal and bronchial
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir.com <- file.path(".", "Comparison_DEGs")
data.dir <- file.path(data.dir.com, "data")
results.dir <- file.path(data.dir.com, "results")
dir.create(results.dir, recursive=TRUE)

write.csv(normalized_counts.N, 
          file = file.path (results.dir, "nasal-log2CPM.csv"))
write.csv(normalized_counts.N, 
          file = file.path (results.dir, "nasal-log2CPM.txt"))

write.csv(normalized_counts.B, 
          file = file.path (results.dir, "bronchial-log2CPM.csv"))
write.csv(normalized_counts.B, 
          file = file.path (results.dir, "bronchial-log2CPM.txt"))

write.table (annotation.nasal, file = file.path (results.dir, "annotation_nasal.txt"))
write.table (annotation.bronchial, file = file.path (results.dir, "annotation_bronchial.txt"))



############ Ensembl - hgnc_symbols ###########
ensembl.dir <- file.path("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/CD")
gene.dir <- file.path(ensembl.dir, "data.nosync")
source(file.path(ensembl.dir, "_helpers.r"), verbose=TRUE, chdir=TRUE)

geneData.N <- getGenedataByEnsemblId(
  ensemblIds = row.names (normalized_counts.N) %>% unique(),
  file.location = gene.dir
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


geneData.B <- getGenedataByEnsemblId(
  ensemblIds = row.names (normalized_counts.B) %>% unique(),
  file.location = gene.dir
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




########## Data selection and organisation ##########
### NASAL
normalized_counts.N <- rownames_to_column(normalized_counts.N, var = "ensembl_gene_id" )
Nasal_exp <- left_join(geneData.N, normalized_counts.N, by = "ensembl_gene_id")
Nasal_exp <- na.omit(Nasal_exp)
Nasal_exp <- Nasal_exp [,-(2), drop=FALSE]
write.table (Nasal_exp, file = file.path  (results.dir, "Totalnasal.normalizedcounts.txt"))

### BRONCHIAL
normalized_counts.B <- rownames_to_column(normalized_counts.B, var = "ensembl_gene_id" )
Bronchial_exp <- left_join(geneData.B, normalized_counts.B, by = "ensembl_gene_id")
Bronchial_exp <- na.omit(Bronchial_exp)
Bronchial_exp <- Bronchial_exp [,-(2), drop=FALSE]
write.table (Bronchial_exp, file = file.path  (results.dir, "totalbronchial.normalizedcounts.txt"))




#### ------- selection ---------
#### Nasal
Nasal <- Nasal_exp  %>% filter (Nasal_exp$hgnc_symbol == "ALOX15B" | 
                                  Nasal_exp$hgnc_symbol == "KCNA3" | 
                                  Nasal_exp$hgnc_symbol == "SPN" | 
                                  Nasal_exp$hgnc_symbol == "P2RY10" | 
                                  Nasal_exp$hgnc_symbol == "ZNF683" | 
                                  Nasal_exp$hgnc_symbol == "TRBC1" )

Nasal <- column_to_rownames(Nasal, var = "hgnc_symbol")
Nasal <- as.data.frame(t(Nasal))
write.csv(Nasal, 
          file = file.path (results.dir, "Nasal.exp_selection.csv"))

N.SG <- merge(annotation.nasal, Nasal, by="row.names")
N.SG  <- N.SG %>%
  rename(
    nasalSAMID = Row.names
  )

N <- N.SG [,-1]
N <- as.data.frame(N)


#### Bronchial
Bronchial <- Bronchial_exp  %>% filter (Bronchial_exp$hgnc_symbol == "ALOX15B" | 
                                        Bronchial_exp$hgnc_symbol == "KCNA3" |
                                        Bronchial_exp$hgnc_symbol == "SPN" |
                                        Bronchial_exp$hgnc_symbol == "P2RY10" |
                                        Bronchial_exp$hgnc_symbol == "ZNF683" |
                                        Bronchial_exp$hgnc_symbol == "TRBC1" )

Bronchial <- column_to_rownames(Bronchial, var = "hgnc_symbol")
Bronchial <- as.data.frame(t(Bronchial))
write.csv(Bronchial, 
          file = file.path (results.dir, "Bronchial.exp_selection.csv"))


B.SG <- merge(annotation.bronchial, Bronchial, by="row.names")
B.SG  <- B.SG %>%
  rename(
    bronchialSAMID = Row.names
  )

B <- B.SG [,-1]
B <- as.data.frame(B)






################# PLOTS #####################
library(ggplot2)
library(ggpubr)
library(ggprism)
library(hrbrthemes)
library(viridis)

########create datasets
### NASAL
sample.N = N %>%
  group_by(group) %>%
  summarize(num=n())

#ALOX15B = 1
N1 <- dplyr::select(N, group, ALOX15B)

#KCNA3 = 2 
N2 <- dplyr::select(N, group, KCNA3)

#SPN = 3
N3 <- dplyr::select(N, group, SPN)

#P2RY10 = 4
N4 <- dplyr::select(N, group, P2RY10)

#ZNF683 = 5
N5 <- dplyr::select(N, group, ZNF683)

#TRBC1 = 6
N6 <- dplyr::select(N, group, TRBC1)


### BRONCHIAL ###
sample.B = B %>%
  group_by(group) %>%
  summarize(num=n())

#ALOX15B = 1
B1 <- dplyr::select(B, group, ALOX15B)

#KCNA3 = 2 
B2 <- dplyr::select(B, group, KCNA3)

#SPN = 3
B3 <- dplyr::select(B, group, SPN)

#P2RY10 = 4
B4 <- dplyr::select(B, group, P2RY10)

#ZNF683 = 5
B5 <- dplyr::select(B, group, ZNF683)

#TRBC1 = 6
B6 <- dplyr::select(B, group, TRBC1)






###################################################
# EXPRESSION PLOTS 
###################################################

#set colours
blue <- "#49BFD5"
yellow <- "#E7E134"
green   <- "#9CD94E"
red <- "#EA5B2D"
magenta <- "#E531CC"


### NASAL ###
#N1 = ALOX15B
N1 <- N1 %>%
  dplyr::mutate(
  group = dplyr::case_when(
    group == "0" ~ "Healthy",
    group == "1" ~ "Mild COPD",
    group == "2" ~ "Severe COPD"
  )
)


my_comparisons = list( c("Healthy", "Mild COPD"), 
                       c("Healthy", "Severe COPD"), 
                       c ("Mild COPD", "Severe COPD"))

## load source code
devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")
## sourced from github "dgrtwo/geom_flat_violin.R



N1_boxplot <- ggplot(N1, aes(x = group, 
              y = ALOX15B, 
              fill = group)
       )+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
              size = 0.5, color="black"
              ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
  #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
             #color = "white",
            #position = position_identity(),
            #position = position_jitterdodge()
            ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Nasal gene expression of ALOX15B") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(7.65, 8.3, 9.1))  #default = wilcox.test

N1_boxplot
ggsave (filename = file.path (results.dir, "Nasal.boxplot_ALOX15B-expression.png"), width=25,height=25,units="cm",dpi=600)


####N2 = KCNA3 
N2 <- N2 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    )
  )

N2_boxplot <- ggplot(N2, aes(x = group, 
                             y = KCNA3, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Nasal gene expression of KCNA3") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(4.8, 5.7, 6.6))  #default = wilcox.test

N2_boxplot
ggsave (filename = file.path (results.dir, "Nasal.boxplot_KCNA3-expression.png"), width=25,height=25,units="cm",dpi=600)


####N3 = SPN
N3 <- N3 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    )
  )

N3_boxplot <- ggplot(N3, aes(x = group, 
                             y = SPN, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Nasal gene expression of SPN") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(6.5, 7.3, 8.1))  #default = wilcox.test

N3_boxplot
ggsave (filename = file.path (results.dir, "Nasal.boxplot_SPN-expression.png"), width=25,height=25,units="cm",dpi=600)


####N4 = P2RY10
N4 <- N4 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    )
  )

N4_boxplot <- ggplot(N4, aes(x = group, 
                             y = P2RY10, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Nasal gene expression of P2RY10") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.3, 6, 6.8))  #default = wilcox.test

N4_boxplot
ggsave (filename = file.path (results.dir, "Nasal.boxplot_P2RY10-expression.png"), width=25,height=25,units="cm",dpi=600)


####N5 = ZNF683
N5 <- N5 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    )
  )

N5_boxplot <- ggplot(N5, aes(x = group, 
                             y = ZNF683, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Nasal gene expression of ZNF683") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(5, 5.9, 6.8))  #default = wilcox.test

N5_boxplot
ggsave (filename = file.path (results.dir, "Nasal.boxplot_ZNF683-expression.png"), width=25,height=25,units="cm",dpi=600)


####N6 = TRBC1
N6 <- N6 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    )
  )

N6_boxplot <- ggplot(N6, aes(x = group, 
                             y = TRBC1, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Nasal gene expression of TRBC1") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(5, 5.8, 6.7))  #default = wilcox.test

N6_boxplot
ggsave (filename = file.path (results.dir, "Nasal.boxplot_TRBC1-expression.png"), width=25,height=25,units="cm",dpi=600)




### ------ BRONCHIAL ------ ###

#########B1 = ALOX15B
B1 <- B1 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    ))


B1_boxplot <- ggplot(B1, aes(x = group, 
                             y = ALOX15B, 
                             fill = group)
  )+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Bronchial gene expression of ALOX15B") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(6, 6.8, 7.6))  #default = wilcox.test

B1_boxplot
ggsave (filename = file.path (results.dir, "Bronchial.boxplot_ALOX15B-expression.png"), width=25,height=25,units="cm",dpi=600)



####### B2 = KCNA3
B2 <- B2 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    ))


B2_boxplot <- ggplot(B2, aes(x = group, 
                             y = KCNA3, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Bronchial gene expression of KCNA3") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(5, 5.8, 6.6))  #default = wilcox.test

B2_boxplot
ggsave (filename = file.path (results.dir, "Bronchial.boxplot_KCNA3-expression.png"), width=25,height=25,units="cm",dpi=600)


####### B3 = SPN
B3 <- B3 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    ))


B3_boxplot <- ggplot(B3, aes(x = group, 
                             y = SPN, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Bronchial gene expression of SPN") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(6.8, 7.6, 8.4))  #default = wilcox.test

B3_boxplot
ggsave (filename = file.path (results.dir, "Bronchial.boxplot_SPN-expression.png"), width=25,height=25,units="cm",dpi=600)

#####B4 = P2RY10
B4 <- B4 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    ))


B4_boxplot <- ggplot(B4, aes(x = group, 
                             y = P2RY10, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Bronchial gene expression of P2RY10") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.3, 6, 6.8))  #default = wilcox.test

B4_boxplot
ggsave (filename = file.path (results.dir, "Bronchial.boxplot_P2RY10-expression.png"), width=25,height=25,units="cm",dpi=600)


##### B5 = ZNF683
B5 <- B5 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    ))


B5_boxplot <- ggplot(B5, aes(x = group, 
                             y = ZNF683, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Bronchial gene expression of ZNF683") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.1, 5.9, 6.6))  #default = wilcox.test

B5_boxplot
ggsave (filename = file.path (results.dir, "Bronchial.boxplot_ZNF683-expression.png"), width=25,height=25,units="cm",dpi=600)

####B6 = TRBC1
B6 <- B6 %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    ))


B6_boxplot <- ggplot(B6, aes(x = group, 
                             y = TRBC1, 
                             fill = group)
)+ 
  scale_fill_manual(values = c(blue, magenta, red)) + 
  geom_flat_violin (alpha = 0.3, position = position_dodge(width = .7), 
                    size = 0.5, color="black"
  ) +
  geom_dotplot(binaxis = "y", dotsize = 0.6, stackdir = "down", binwidth = 0.2
               #geom_point(shape = 23, size=1.5, alpha = 0.3, fill = "black"
               #color = "white",
               #position = position_identity(),
               #position = position_jitterdodge()
  ) +
  geom_boxplot(alpha = 1, notch = FALSE, width = 0.1, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7, 
               #fill = "white"
  ) +
  theme_classic() + 
  ylab("Expression level") +
  xlab(" ") +
  ggtitle("Bronchial gene expression of TRBC1") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.1, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 15, face = "bold")) +
  
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.1, 5.9, 6.6))  #default = wilcox.test

B6_boxplot
ggsave (filename = file.path (results.dir, "Bronchial.boxplot_TRBC1-expression.png"), width=25,height=25,units="cm",dpi=600)



# ------------ Combine all boxplots ------------ 
#ALOX15B = 1
#KCNA3 = 2 
#SPN = 3
#P2RY10 = 4
#ZNF683 = 5
#TRBC1 = 6

box.ALOX15B <- ggarrange (N1_boxplot, B1_boxplot, labels = c("A", "B"))
box.ALOX15B
ggsave (
  filename = file.path (results.dir, "Boxplot_N-B_ALOX15B.png"), 
  plot = box.ALOX15B,
  width=58,
  height=35,
  units="cm",
  dpi=600)

box.KCNA3 <- ggarrange (N2_boxplot, B2_boxplot, labels = c("C", "D"))
box.KCNA3
ggsave (
  filename = file.path (results.dir, "Boxplot_N-B_KCNA3.png"), 
  plot = box.KCNA3,
  width=58,
  height=35,
  units="cm",
  dpi=600)

box.SPN <- ggarrange (N3_boxplot, B3_boxplot, labels = c("E", "F"))
box.SPN
ggsave (
  filename = file.path (results.dir, "Boxplot_N-B_SPN.png"), 
  plot = box.SPN,
  width=58,
  height=35,
  units="cm",
  dpi=600)

box.P2RY10 <- ggarrange (N4_boxplot, B4_boxplot, labels = c("G", "H"))
box.P2RY10
ggsave (
  filename = file.path (results.dir, "Boxplot_N-B_P2RY10.png"), 
  plot = box.P2RY10,
  width=58,
  height=35,
  units="cm",
  dpi=600)

box.ZNF683 <- ggarrange (N5_boxplot, B5_boxplot, labels = c("I", "J"))
box.ZNF683
ggsave (
  filename = file.path (results.dir, "Boxplot_N-B_ZNF683.png"), 
  plot = box.ZNF683,
  width=58,
  height=35,
  units="cm",
  dpi=600)

box.TRBC1 <- ggarrange (N6_boxplot, B6_boxplot, labels = c("K", "L"))
box.TRBC1
ggsave (
  filename = file.path (results.dir, "Boxplot_N-B_TRBC1.png"), 
  plot = box.TRBC1,
  width=58,
  height=35,
  units="cm",
  dpi=600)


##### combine all
box.all <- ggarrange (box.ALOX15B, box.KCNA3, box.SPN, box.P2RY10, box.ZNF683,box.TRBC1,
                       ncol = 2, nrow = 3)
box.all
ggsave (
  filename = file.path (results.dir, "Boxplot_total6.png"), 
  plot = box.all,
  width=85,
  height=70,
  units="cm",
  dpi=600)

ggsave (
  filename = file.path (results.dir, "Boxplot_total6.pdf"), 
  plot = box.all,
  width=85,
  height=70,
  units="cm",
  dpi=600)



box.all2 <- ggarrange (box.ALOX15B, box.KCNA3, box.SPN, box.P2RY10, box.ZNF683,box.TRBC1,
                      ncol = 1, nrow = 6)
box.all2
ggsave (
  filename = file.path (results.dir, "Boxplot2_total6.png"), 
  plot = box.all2,
  width=30,
  height=93,
  units="cm",
  dpi=600)

ggsave (
  filename = file.path (results.dir, "Boxplot2_total6.pdf"), 
  plot = box.all2,
  width=30,
  height=93,
  units="cm",
  dpi=600)


################################################
# Bronchial versus nasal
################################################
###paired samples
setwd("C:/Users/Jowi/Documents/HVHL/Stages/Afstudeerstage_2022/DE analysis")
data.dir.NB_DEG <- file.path(".", "Nasal DEG_DE in bronchi")
data.dir.NB <- file.path(data.dir.NB_DEG, "results")

Match <- readr::read_csv( 
  file.path(data.dir, "Matched samples_SAMID.csv"))
Match <- Match [,-1]
Match <- Match [,-2]

Unique.NB <- readr::read_csv( 
  file.path(data.dir.NB, "bronchial_unique-DEGs_severeCOPD-P0.01_NNLS.csv"))
Unique.NB <- Unique.NB [,-(2:7)]


# --- correlation between all nasal and bronchial samples, which genes to select? -------
exp.N <- inner_join(Nasal_exp, Unique.NB, by = "hgnc_symbol" )
exp.N <- exp.N %>% 
  arrange(hgnc_symbol)
readr::write_csv(exp.N,
                 file = file.path (results.dir, "Nasal_selection.normalizedexpression.csv"))


exp.N <- exp.N[!(exp.N$hgnc_symbol== "AMMECR1-IT1"),]
exp.N <- column_to_rownames(exp.N, var = "hgnc_symbol")
exp.N <- t(as.matrix(exp.N))
exp.N <- as.data.frame(exp.N)
exp.N <- rownames_to_column(exp.N, var = "nasalSAMID")

exp.N <- inner_join(Match, exp.N, by = "nasalSAMID")
exp.N <- exp.N [,-(1:2)]

#exp.N <- column_to_rownames(exp.N, var = "ID")


exp.B <- inner_join(Bronchial_exp, Unique.NB, by = "hgnc_symbol" )
exp.B <- exp.B %>% 
  arrange(hgnc_symbol)
readr::write_csv(exp.B,
                 file = file.path (results.dir, "Bronchial-selection_normalizedexpression.csv"))

exp.B <- column_to_rownames(exp.B, var = "hgnc_symbol")
exp.B <- t(as.matrix(exp.B))
exp.B <- as.data.frame(exp.B)
exp.B <- rownames_to_column(exp.B, var = "bronchialSAMID")

exp.B <- inner_join(Match, exp.B, by = "bronchialSAMID")
exp.B <- exp.B [,-(1:2)]
#exp.B <- column_to_rownames(exp.B, var = "ID")

exp.total <- full_join(exp.N, exp.B, by = "ID")
readr::write_csv(exp.total,
                 file = file.path (results.dir, "Nasal_bronchial-selection-normalizedexpression.csv"))
exp.total <- column_to_rownames(exp.total, var = "ID") 


library(readxl)
#an_dan <- as.data.frame(read_exel("data.xlsx",sheet=2))

an_dan <- exp.total
cor_an <- matrix(NA, nrow = 103, ncol = 103)
colnames(cor_an) <- colnames (an_dan[1:103])

p_an <- matrix(NA, nrow = 103, ncol = 103)
colnames(p_an) <- colnames (an_dan[1:103])
  
k <- 1
for (t in 1:103) {
  i <- 1
  for (j in 104:206) {
    cor_an1 <- cor.test(an_dan[,t],an_dan[,j],alternative = "two.sided",method = "spearman",conf.level = 0.95)
    cor_an[i,k] <- cor_an1$estimate
    p_an[i,k] <- cor_an1$p.value
    
    i <- i+1
    j <- j+1
  }
  t <- t+1
  k <- k+1
}

cor_an <- as.data.frame(cor_an)
readr::write_csv(cor_an,
                 file = file.path (results.dir, "cor_an.csv"))

p_an <- as.data.frame(p_an)
readr::write_csv(p_an,
                 file = file.path (results.dir, "p_an.csv"))




############----------- selection - select and correlation plot --------------##########
#ALOX15B = 1
#KCNA3 = 2 
#SPN = 3
#P2RY10 = 4
#ZNF683 = 5
#TRBC1 = 6

### Nasal
Nasal <- rownames_to_column(Nasal, var = "nasalSAMID")
Total.N <- merge(Nasal, Match, by = "nasalSAMID")

Total.N <- dplyr::select(
  Total.N, ALOX15B, KCNA3, SPN, P2RY10, ZNF683, TRBC1, ID)

Total.N <- Total.N %>%
  rename(
    ALOX15B.nasal = ALOX15B,
    KCNA3.nasal = KCNA3,
    SPN.nasal = SPN,
    P2RY10.nasal = P2RY10,
    ZNF683.nasal = ZNF683,
    TRBC1.nasal = TRBC1)

Total.N_ALOX15B <- dplyr::select(
  Total.N, ALOX15B.nasal, ID)

Total.N_KCNA3 <- dplyr::select(
  Total.N, KCNA3.nasal, ID)

Total.N_SPN <- dplyr::select(
  Total.N, SPN.nasal, ID)

Total.N_P2RY10 <- dplyr::select(
  Total.N, P2RY10.nasal, ID)

Total.N_ZNF683 <- dplyr::select(
  Total.N, ZNF683.nasal, ID)

Total.N_TRBC1 <- dplyr::select(
  Total.N, TRBC1.nasal, ID)


### bronchial 
Bronchial <- rownames_to_column(Bronchial, var = "bronchialSAMID") 
Total.B <- merge(Bronchial, Match, by = "bronchialSAMID")
Total.B <- dplyr::select(
  Total.B, ALOX15B, KCNA3, SPN, P2RY10, ZNF683, TRBC1, ID)

Total.B <- Total.B %>%
  rename(
    ALOX15B.bronchial = ALOX15B,
    KCNA3.bronchial = KCNA3,
    SPN.bronchial = SPN,
    P2RY10.bronchial = P2RY10,
    ZNF683.bronchial = ZNF683,
    TRBC1.bronchial = TRBC1)

Total.B_ALOX15B <- dplyr::select(
  Total.B, ALOX15B.bronchial, ID)

Total.B_KCNA3 <- dplyr::select(
  Total.B, KCNA3.bronchial, ID)

Total.B_SPN <- dplyr::select(
  Total.B, SPN.bronchial, ID)

Total.B_P2RY10 <- dplyr::select(
  Total.B, P2RY10.bronchial, ID)

Total.B_ZNF683 <- dplyr::select(
  Total.B, ZNF683.bronchial, ID)

Total.B_TRBC1 <- dplyr::select(
  Total.B, TRBC1.bronchial, ID)


#combine paired nasal and bronchial expression in 1 table
Total_ALOX15B <- merge(Total.N_ALOX15B, Total.B_ALOX15B, by = "ID")
Total_KCNA3 <- merge(Total.N_KCNA3, Total.B_KCNA3, by = "ID")
Total_SPN <- merge(Total.N_SPN, Total.B_SPN, by = "ID")
Total_P2RY10 <- merge(Total.N_P2RY10, Total.B_P2RY10, by = "ID")
Total_ZNF683 <- merge(Total.N_ZNF683, Total.B_ZNF683, by = "ID")
Total_TRBC1 <- merge(Total.N_TRBC1, Total.B_TRBC1, by = "ID")

#group categatorization
N.SG <- N.SG[,-(3:4)]
B.SG <- B.SG[,-(3:4)]

MatchG <- merge(N.SG, Match, by = "nasalSAMID")
MatchG <- dplyr::select(MatchG,
  group, ID 
  ) %>%
  dplyr::mutate(
  group = dplyr::case_when(
    group == "0" ~ "Healthy",
    group == "1" ~ "Mild COPD",
    group == "2" ~ "Severe COPD"
  ))


Total_ALOX15B <- merge (Total_ALOX15B, MatchG)
Total_ALOX15B$group <- factor( Total_ALOX15B$group, levels = c("Healthy", "Mild COPD", "Severe COPD"))
readr::write_csv(Total_ALOX15B,
                 file = file.path (results.dir, "Total.expression_ALOX15B.csv"))

Total_KCNA3 <- merge (Total_KCNA3, MatchG)
Total_KCNA3$group <- factor(Total_KCNA3$group, levels = c("Healthy", "Mild COPD", "Severe COPD"))
readr::write_csv(Total_KCNA3,
                 file = file.path (results.dir, "Total.expression_KCNA3.csv"))

Total_SPN <- merge (Total_SPN, MatchG)
Total_SPN$group <- factor(Total_SPN$group, levels = c("Healthy", "Mild COPD", "Severe COPD"))
readr::write_csv(Total_SPN,
                 file = file.path (results.dir, "Total.expression_SPN.csv"))

Total_P2RY10 <- merge (Total_P2RY10, MatchG)
Total_P2RY10$group <- factor(Total_P2RY10$group, levels = c("Healthy", "Mild COPD", "Severe COPD"))
readr::write_csv(Total_P2RY10,
                 file = file.path (results.dir, "Total.expression_P2RY10.csv"))

Total_ZNF683 <- merge (Total_ZNF683, MatchG)
Total_ZNF683$group <- factor(Total_ZNF683$group, levels = c("Healthy", "Mild COPD", "Severe COPD"))
readr::write_csv(Total_ZNF683,
                 file = file.path (results.dir, "Total.expression_ZNF683.csv"))

Total_TRBC1 <- merge (Total_TRBC1, MatchG)
Total_TRBC1$group <- factor(Total_TRBC1$group, levels = c("Healthy", "Mild COPD", "Severe COPD"))
readr::write_csv(Total_TRBC1,
                 file = file.path (results.dir, "Total.expression_TRBC1.csv"))


##########################################
# SCATTERPLOTS/correlation
#########################################
library("ggpubr")
#ALOX15B = 1
#KCNA3 = 2 
#SPN = 3
#P2RY10 = 4
#ZNF683 = 5
#TRBC1 = 6

#blue <- "#49BFD5", magenta <- "#E531CC" and red <- "#EA5B2D"

##### ALOX15B expression nasal vs bronchial #####
ALOX15B_reg <- ggscatter(Total_ALOX15B, 
                   x = "ALOX15B.nasal", 
                   y = "ALOX15B.bronchial", 
                   fill = "group",
                   add = "reg.line", 
                   conf.int = FALSE,
                   cor.coef = TRUE, 
                   cor.method = "spearman",
                   title = "The correlation between nasal and bronchial ALOX15B gene expression", 
                   xlab = "Nasal expression", 
                   ylab = "Bronchial expression"
                   ) +
              ggplot2::theme(
                axis.title = element_text(size = 12.5),
                axis.text = element_text(size = 12),
                plot.title = element_text(size = 14.5, face = "bold")) 

ALOX15B_reg
ggsave(
  filename = file.path (results.dir, "ALOX15B-expression_correlation_regressiononly.png"), 
  plot = ALOX15B_reg,
  width=25,
  height=25,
  units="cm",
  dpi=600)

#49BFD5", magenta <- "#E531CC" and red <- "#EA5B2D"

###with color per group and no reg.line
ALOX15B <- ggscatter(Total_ALOX15B, 
                     x = "ALOX15B.nasal", 
                     y = "ALOX15B.bronchial", 
                     color = "group",
                     palette = c("#49BFD5", "#E531CC", "#EA5B2D"), 
                     #add = "reg.line", 
                     conf.int = FALSE,
                     cor.coef = TRUE, 
                     cor.method = "spearman",
                     cor.coef.size = 5.2,
                     cor.coef.coord = c(0.32, 4.9),
                     title = "A  Nasal and bronchial ALOX15B gene expression", 
                     xlab = "Nasal expression", 
                     ylab = "Bronchial expression"
) +
  ggplot2::theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(size = 16.5, vjust =4.5, face = "bold"),
    plot.margin = unit(c(1,0.7,0.7,0.7), "cm"))

ALOX15B
ggsave(
  filename = file.path (results.dir, "ALOX15B-expression_correlation.png"), 
  plot = ALOX15B,
  width=25,
  height=25,
  units="cm",
  dpi=600)



###with color per group and no reg.line
KCNA3 <- ggscatter(Total_KCNA3, 
                   x = "KCNA3.nasal", 
                   y = "KCNA3.bronchial", 
                   color = "group",
                   palette = c("#49BFD5", "#E531CC", "#EA5B2D"), 
                   #add = "reg.line", 
                   conf.int = FALSE,
                   cor.coef = TRUE, 
                   cor.method = "spearman",
                   cor.coef.size = 5.2,
                   cor.coef.coord = c(-3.5, 2.8),
                   title = "B  Nasal and bronchial KCNA3 gene expression", 
                   xlab = "Nasal expression", 
                   ylab = "Bronchial expression"
                   ) +
  ggplot2::theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(size = 16.5, vjust =4.5, face = "bold"),
    plot.margin = unit(c(1,0.7,0.7,0.7), "cm"))


KCNA3
ggsave(
  filename = file.path (results.dir, "KCNA3-expression_correlation.png"), 
  plot = KCNA3,
  width=25,
  height=25,
  units="cm",
  dpi=600)


SPN <- ggscatter(Total_SPN, 
                   x = "SPN.nasal", 
                   y = "SPN.bronchial", 
                   color = "group",
                   palette = c("#49BFD5", "#E531CC", "#EA5B2D"), 
                   #add = "reg.line", 
                   conf.int = FALSE,
                   cor.coef = TRUE, 
                   cor.method = "spearman",
                   cor.coef.size = 5.2,
                   cor.coef.coord = c(-2, 5.1),
                   title = "C  Nasal and bronchial SPN gene expression", 
                   xlab = "Nasal expression", 
                   ylab = "Bronchial expression"
                 ) +
  ggplot2::theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(size = 16.5, vjust =4.5, face = "bold"),
    plot.margin = unit(c(1,0.7,0.7,0.7), "cm"))

SPN
ggsave(
  filename = file.path (results.dir, "SPN-expression_correlation.png"), 
  plot = SPN,
  width=25,
  height=25,
  units="cm",
  dpi=600)



P2RY10 <- ggscatter(Total_P2RY10, 
                 x = "P2RY10.nasal", 
                 y = "P2RY10.bronchial", 
                 color = "group",
                 palette = c("#49BFD5", "#E531CC", "#EA5B2D"), 
                 #add = "reg.line", 
                 conf.int = FALSE,
                 cor.coef = TRUE, 
                 cor.method = "spearman",
                 cor.coef.size = 5.2,
                 cor.coef.coord = c(-3.5, 2.5),
                 title = "D  Nasal and bronchial P2RY10 gene expression", 
                 xlab = "Nasal expression", 
                 ylab = "Bronchial expression"
) +
  ggplot2::theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(size = 16.5, vjust =4.5, face = "bold"),
    plot.margin = unit(c(1,0.7,0.7,0.7), "cm"))


P2RY10
ggsave(
  filename = file.path (results.dir, "P2RY10-expression_correlation.png"), 
  plot = P2RY10,
  width=25,
  height=25,
  units="cm",
  dpi=600)



ZNF683 <- ggscatter(Total_ZNF683, 
                    x = "ZNF683.nasal", 
                    y = "ZNF683.bronchial", 
                    color = "group",
                    palette = c("#49BFD5", "#E531CC", "#EA5B2D"), 
                    #add = "reg.line", 
                    conf.int = FALSE,
                    cor.coef = TRUE, 
                    cor.method = "spearman",
                    cor.coef.size = 5.2,
                    cor.coef.coord = c(-3.5, 2.1),
                    title = "E  Nasal and bronchial ZNF683 gene expression", 
                    xlab = "Nasal expression", 
                    ylab = "Bronchial expression"
) +
  ggplot2::theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(size = 16.5, vjust =4.5, face = "bold"),
    plot.margin = unit(c(1,0.7,0.7,0.7), "cm"))


ZNF683
ggsave(
  filename = file.path (results.dir, "ZNF683-expression_correlation.png"), 
  plot = ZNF683,
  width=25,
  height=25,
  units="cm",
  dpi=600)


TRBC1 <- ggscatter(Total_TRBC1, 
                    x = "TRBC1.nasal", 
                    y = "TRBC1.bronchial", 
                    color = "group",
                    palette = c("#49BFD5", "#E531CC", "#EA5B2D"), 
                    #add = "reg.line", 
                    conf.int = FALSE,
                    cor.coef = TRUE, 
                    cor.method = "spearman",
                    cor.coef.size = 5.2,
                    cor.coef.coord = c(-2.8, 3),
                    title = "F  Nasal and bronchial TRBC1 gene expression", 
                    xlab = "Nasal expression", 
                    ylab = "Bronchial expression"
                   ) +
  ggplot2::theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(size = 16.5, vjust =4.5, face = "bold"),
    plot.margin = unit(c(1,0.7,0.7,0.7), "cm"))


TRBC1
ggsave(
  filename = file.path (results.dir, "TRBC1-expression_correlation.png"), 
  plot = TRBC1,
  width=25,
  height=25,
  units="cm",
  dpi=600)

#legend.text = element_text(size = 12),
#legend.direction = "vertical", 
#legend.position = "right", 


### arrange into one figure
corrplot_all <- ggarrange (ALOX15B, KCNA3, SPN, P2RY10, ZNF683, TRBC1,
                       ncol = 2, nrow = 3)
corrplot_all
ggsave (
  filename = file.path (results.dir, "Corplot_total6.png"), 
  plot = corrplot_all,
  width=50,
  height=70,
  units="cm",
  dpi=600)


corrplot_all2 <- ggarrange (ALOX15B, KCNA3, SPN, P2RY10, ZNF683, TRBC1,
                           ncol = 3, nrow = 2)
corrplot_all2
ggsave (
  filename = file.path (results.dir, "Corplot2_total6.png"), 
  plot = corrplot_all2,
  width=80,
  height=50,
  units="cm",
  dpi=600)

