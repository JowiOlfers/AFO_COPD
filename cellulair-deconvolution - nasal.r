library(tidyverse)
library (dplyr)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")
library(limma)
library(edgeR)


setwd("~/HVHL/Stages/Afstudeerstage_2022/CD")
source(file.path("_helpers.r"), verbose=TRUE, chdir=TRUE)

 # Global Vars, basically
data.dir <- file.path(".", "data.nosync")
data.dir.BB.mRNA <- file.path(data.dir, "RNA-Seq-BB-mRNA-v38")
data.dir.NB.mRNA <- file.path(data.dir, "RNA-Seq-NB-mRNA-v38")
data.dir.patients <- file.path(data.dir, "patients")
data.dir.samples <- file.path(data.dir, "samples")

results.dir <- file.path(here::here(), "results.nosync")
dir.create(results.dir, recursive=TRUE)


################## Load data #####################
#### Expression data + gene description
expressionData <- readxl::read_xlsx(
  file.path(data.dir.NB.mRNA, "nasal rawcounts.xlsx")
  #col_types = readr::cols()
)

geneData <- getGenedataByEnsemblId(
  ensemblIds = expressionData$Gene %>% unique(),
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

expression <- expressionData %>%
  tibble::column_to_rownames("Gene") %>%
  cpm() %>%                     # normalized
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::left_join(
    y = geneData %>%
      dplyr::select(hgnc_symbol,ensembl_gene_id),
    by = c("Gene" = "ensembl_gene_id")
  ) %>%
  filter(!is.na(hgnc_symbol)) %>%
  dplyr::select(-Gene) %>%
  tibble::column_to_rownames("hgnc_symbol")

genes <- row.names(expression)


#### Patient/clinical data
masterTable <- readxl::read_xlsx(
  file.path(data.dir.patients, "nasal clinical data.xlsx")
) %>%
  dplyr::select(
    1:68
  )



#############################################################
#cellulair-deconvolution
#############################################################
library(nlme)
library(nnls)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(hrbrthemes)

#Directory nasal
out.dir <- file.path(results.dir, "nasal_400")
out.data.dir <- file.path(out.dir, "data")
out.img.dir <- file.path(out.dir, "img")
out.bplot.img.dir <- file.path(out.img.dir, "boxplots")

dir.create(out.data.dir, recursive=TRUE)
dir.create(out.bplot.img.dir, recursive=TRUE)

#### Load scRNAseq data ####
data.dir.scRNAseq_nasal <- file.path (".", "Nasal Epithelial")
RNAseq <- read.csv(file.path(data.dir.scRNAseq_nasal, "nasal_epithelial_means_woRare_400_0.csv"))
rownames(RNAseq)<-RNAseq$X
RNAseq <- RNAseq[,-1]


###########################################
# Deconvolution nasal                                               
###########################################
## source cibersort
source("CIBERSORT.R", verbose=TRUE)

theme_set(theme_prism())

keep <- intersect(rownames(RNAseq), rownames(expression))
Bulk <- expression[keep,]
Ref <- RNAseq[keep,] %>%
  dplyr::rename_all(
    function(col.name) {
      col.name <- gsub(
        x = col.name,
        pattern = "X3_",
        replacement = ""
      )

      col.name <- gsub(
        x = col.name,
        pattern = "\\.",
        replacement = " "
      )
    }
  )


######  1. CIBERSORT (nu-svr) ######
xf <- file.path(out.data.dir, 'reference.tsv')
Ref %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rownames")%>%
  readr::write_tsv(
    file = xf
  )

yf <- file.path(out.data.dir, 'mixture.tsv')
Bulk %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rownames") %>%
  readr::write_tsv(
    file = yf
  )

results <- CIBERSORT(sig_matrix = xf, mixture_file = yf, QN = FALSE,perm=100)
res.cibersort <- t(results[,1:(ncol(results)-3)])


df.ciber <- res.cibersort %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_type") %>%
  tidyr::pivot_longer(
    -cell_type,
    names_to = "sample",
    values_to = "proportion"
  ) %>%
  dplyr::mutate(
    method = "CIBERSORT"
  ) %>%
  readr::write_csv(
    file = file.path(out.dir, "ciber_result.csv")
  )

# If I need to make a heatmap with z-scores instead:
#cal_z_score <- function(x){
#  (x - mean(x)) / sd(x)}
#data_subset_norm <- t(apply(data_subset, 1, cal_z_score))


annotation_col <- masterTable %>%
  dplyr::select(
    SAMID,
    group
  ) %>% 
  dplyr::mutate(
  group=factor(group)) %>%
  tibble::column_to_rownames("SAMID")

hm.ciber <- pheatmap(
  res.cibersort [,rownames(annotation_col)],
  main = "CIBERSORT",
  show_colnames = FALSE,
  color = colorRampPalette(c("grey92", "grey40", "grey0"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  #gaps_col = ,
  fontsize = 20,
  annotation_col = annotation_col
)

saveRDS(
  hm.ciber,
  file = file.path(out.img.dir, "hm_cibersort.rds")
)
png(
  filename = file.path(out.img.dir, "hm_cibersort.png"),
  width = 1000,
  height = 800,
  units = "px"
)
print(hm.ciber)
dev.off()


#add certain groups together?


###### 2.nnls ######
Results_nnls = do.call(
    cbind.data.frame,
    lapply(
      apply(
        Bulk,
        2,
        function(x) nnls::nnls(as.matrix(Ref),x)
      ),
      function(y) y$x
    )
  )
res.nnls = apply(Results_nnls,2,function(x) x/sum(x)) #explicit STO constraint????
rownames(res.nnls) <- stringr::str_trim(colnames(Ref))


df.nnls <- res.nnls %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_type") %>%
  tidyr::pivot_longer(
    -cell_type,
    names_to = "sample",
    values_to = "proportion"
  ) %>%
  dplyr::mutate(
    method = "NNLS"
  ) %>%
  readr::write_csv(
    file = file.path(out.dir, "nnls_result.csv")
  )


hm.nnls <- pheatmap(
  res.nnls [,rownames(annotation_col)],
  main = "NNLS",
  show_colnames = FALSE,
  color = colorRampPalette(c("grey92", "grey40", "grey0"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 20,
  annotation_col = annotation_col
) 

saveRDS(
  hm.nnls,
  file = file.path(out.img.dir, "hm_nnls.rds")
)
png(
  filename = file.path(out.img.dir, "hm_nnls.png"),
  width = 1000,
  height = 800,
  units = "px"
)
print(hm.nnls)
dev.off()



###########
summry_table <- df.nnls %>%
  dplyr::add_row( # weg?
    df.ciber
  ) %>%
  dplyr::left_join(
    y = masterTable %>%
      dplyr::select(SAMID, group),
    by = c("sample" = "SAMID")
  ) %>%
  dplyr::select(-sample) %>%
  dplyr::nest_by(method) %>%
  dplyr::mutate(
    summary = data %>%
      dplyr::group_by(cell_type, group) %>%
      dplyr::summarize(
        .groups = "keep",
        summary = paste0(
          base::mean(proportion) %>%
            round(digits = 2) %>%
            format(nsmall = 2),
          " (",
          stats::sd(proportion) %>%
            round(digits = 2) %>%
            format(nsmall = 2),
          ")"
        )
      ) %>%
      list()
  ) %>%
  dplyr::select(-data) %>%
  tidyr::unnest(summary) %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = summary
  ) %>%
  readr::write_csv(
    file = file.path(out.dir, "result_summary.nasal.csv")
  )


#summary table2 - no groups, 1 combined group
summry_table2 <- df.nnls %>%
    dplyr::add_row( # weg?
      df.ciber
    ) %>%
    dplyr::select(-sample) %>%
    dplyr::nest_by(method) %>%
    dplyr::mutate(
      summary = data %>%
        dplyr::group_by(cell_type) %>%
        dplyr::summarize(
          .groups = "keep",
          summary = paste0(
            base::mean(proportion) %>%
              round(digits = 2) %>%
              format(nsmall = 2),
            " (",
            stats::sd(proportion) %>%
              round(digits = 2) %>%
              format(nsmall = 2),
            ")"
          )
        ) %>%
        list()
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(summary) %>%
    readr::write_csv(
      file = file.path(out.dir, "result_summary_2.nasal.csv")
    )



######### Compare results of methods and groups #############


#############################################
# Boxplots of two methods per cell type                       
#############################################
df <- rbind(df.ciber, df.nnls)
box.plt <- df %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = cell_type,
      y = proportion,
      fill = method
    )
  ) +
  ggplot2::geom_boxplot(
    outlier.alpha = 0
  ) +
  ggplot2::geom_jitter(
    color = "black",
    size = 0.4,
    alpha = 0.35
  ) +
  ggplot2::theme(
    legend.position = "bottom"
  ) +
  ggplot2::coord_flip() +
  ggplot2::xlab(label = "")

df %>%
  readr::write_csv(
    file = file.path(out.bplot.img.dir, "boxplots.nasal.csv")
  )

saveRDS(
  box.plt,
  file = file.path(out.bplot.img.dir, "boxplots.nasal.rds")
)

ggsave(
  filename = file.path(out.bplot.img.dir, "boxplots.nasal.png"),
  plot = box.plt,
  width = 20,
  height = 12.5,
  units = "cm"
)


################################################################
# Correlation plots of two methods per cell types            
################################################################
if (!identical(colnames(res.cibersort),colnames(res.nnls))) {
  stop("Not all cell types in both data frames.")
}

df.cor <- data.frame(
    cell_type = NA_character_,
    coef = NA_real_,
    pvalue = NA_real_
  )

#i = iterator
for (i in rownames(res.nnls)) {
  t <- cor.test(
      res.cibersort[stringr::str_trim(i),],
      res.nnls[i,],
      method = "spearman"
    )

  df.cor <- df.cor %>%
    dplyr::add_row(
      data.frame(
        cell_type = i,
        coef = t$estimate,
        pvalue = t$p.value
      )
    ) %>%
    dplyr::filter(
      !is.na(cell_type)
    ) 
  }

df.cor %>%
  readr::write_csv(
    file = file.path(out.dir, "spearman_correlation_two_methods.csv")
  )

sc_plots <- list()
for (i in 1:min(nrow(res.cibersort), 16)) {
  df.i <- data.frame(
      CIBERSORT = res.cibersort[i,],
      NNLS = res.nnls[i,]
    )
  sc_plots[[i]] <- df.i %>%
    ggplot(
      mapping = aes(
        x = CIBERSORT,
        y = NNLS
      )
    ) +
    geom_point() +
    ggtitle(NULL, subtitle = rownames(res.cibersort)[i]) +
    rremove("ylab") +
    rremove("xlab") +
    theme(
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 7),
      plot.title = element_text(size = 9)
    )
}

sc_plots[['ncol']] <- 4
sc_plots[['nrow']] <- 4

plt <- do.call(
  ggarrange,
  sc_plots
)

plt <- annotate_figure(plt,
    bottom = text_grob("CIBERSORT"),
    left = text_grob("NNLS", rot = 90)
  )

saveRDS(
  plt,
  file = file.path(out.img.dir, "scatters.nasal.rds")
)

ggsave(
  filename = file.path(out.img.dir, "scatters.nasal.png"),
  plot = plt
)


###################################################
# Compare Cell Types over Disease Groups             
###################################################
big.df <- df.nnls %>%
  dplyr::add_row(df.ciber) %>%
  dplyr::left_join(
    y = masterTable %>%
      dplyr::select(
        SAMID,
        group
      ),
    by = c("sample" = "SAMID")
  )

big.df2 <- big.df %>%
  dplyr::group_by(
    cell_type, method
  ) %>%
  dplyr::filter(
    mean(proportion) >= 0.01
  ) %>%
  dplyr::ungroup()


kruskal.wallis.tests <- big.df2 %>%
  dplyr::nest_by(
    cell_type, method
  ) %>%
  dplyr::mutate(
    kruskal.wallis.test = kruskal.test(
      proportion ~ group,
      data = data
    ) %>%
    broom::tidy() %>%
    dplyr::rename(test.method = method) %>%
    list()
  ) %>%
  tidyr::unnest(cols = kruskal.wallis.test) %>%
  dplyr::select(-data)


celltype_bp <- list()
for (c.method in unique(big.df2$method)) {
  for (cell.type in unique(big.df2$cell_type)) {
    kruskal.wallis.pval <- kruskal.wallis.tests %>%
      dplyr::filter(
        cell_type == cell.type &
        method == c.method
      ) %>%
      dplyr::pull(p.value)

    big.figure.data <- big.df2 %>%
      dplyr::filter(
        cell_type == cell.type &
        method == c.method
      ) %>%
      dplyr::group_by(
        cell_type
      ) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group == "0" ~ "Healthy",
      group == "1" ~ "Mild COPD",
      group == "2" ~ "Severe COPD"
    )
  )
    

    big.figure <- big.figure.data %>%
      ggplot2::ggplot(
        mapping =  ggplot2::aes(
          x = group,
          y = proportion,
          fill = group
        )
      ) +
      ggplot2::geom_boxplot(
        outlier.alpha = 0
      ) +
      ggplot2::geom_jitter(
        color = "black",
        #size = 0.4,
        alpha = 0.35
      ) +
      ggpubr::stat_compare_means(
        comparisons = list(
          c("Healthy", "Mild COPD"),
          c("Healthy", "Severe COPD"),
          c("Mild COPD", "Severe COPD")
         
        ),
        #label = "p.format",
        label.y.npc = "top",
        size = 2.8
      ) +
      ggprism::theme_prism() +
      ggplot2::theme(
        axis.text.x = element_text(
          angle = 45
        ),
        legend.position = "none"
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("Proportion") +
      ggplot2::ggtitle(
        cell.type,
        subtitle = paste0("Kruskal-Wallis-test: ", format.pval(kruskal.wallis.pval))
      ) +
      ggprism::scale_colour_prism() +
      ggplot2::scale_fill_manual(
        values = c(
          "Healthy" = "#ff3333",
          "Mild COPD" = "#33ff33",
          "Severe COPD" = "#3333ff"
        )
      )

    saveRDS(
      big.figure,
      file = file.path(out.bplot.img.dir, paste0("boxplots-cell-types-",make.names(cell.type),"-",c.method,".rds")))
    

    ggsave(
      filename = file.path(out.bplot.img.dir, paste0("boxplots-cell-types-",make.names(cell.type),"-",c.method,".png")),
      plot = big.figure,
      width = 15,
      height = 15,
      units = "cm"
    )

    big.figure.data %>%
      readr::write_csv(
        file = file.path(out.bplot.img.dir, paste0("boxplots-cell-types-",make.names(cell.type),"-",c.method,".csv"))
      )

    celltype_bp[[make.names(cell.type)]] <- big.figure
  }
}


########################################################################
# Compare Cell Types over Disease Groups: Test groups
#########################################################################
t.test.res <- data.frame(
  cell.type = character(),
  method = character(),
  group.1 = character(),
  group.2 = character(),
  p.val = numeric(),
  statistic = numeric(),
  stderr = numeric()
)

big.df <- big.df %>%
  dplyr::mutate(
  group=factor(group))
  
combinations <- combn(levels(big.df$group), 2)

for (cell.type_ in unique(big.df$cell_type)) {
  for (method_ in unique(big.df$method)) {
    t.test.groups <- big.df %>%
      dplyr::filter(
        cell_type == cell.type_,
        method == method_
      ) %>%
      dplyr::nest_by(group) %>%
      dplyr::mutate(
        data = list(data$proportion)
      )

    for (combination in 1:ncol(combinations)) {
      group.1.name <- combinations[1, combination]
      group.2.name <- combinations[2, combination]

      group.1 <- t.test.groups %>%
        dplyr::filter(group == group.1.name) %>%
        dplyr::pull(data) %>%
        dplyr::first()
      group.2 <- t.test.groups %>%
        dplyr::filter(group == group.2.name) %>%
        dplyr::pull(data) %>%
        dplyr::first()

      c.t.test.res <- t.test(
        group.1,
        group.2
      )

      t.test.res <- rbind(
        t.test.res,
        data.frame(
          cell.type = cell.type_,
          method = method_,
          group.1 = group.1.name,
          group.2 = group.2.name,
          p.val = c.t.test.res$p.value,
          statistic = c.t.test.res$statistic,
          stderr = c.t.test.res$stderr
        )
      )
      rm(group.1.name, group.2.name, c.t.test.res)
    }
    rm(method_)
  }
  rm(cell.type_)
}

t.test.res %>%
  readr::write_csv(file.path(out.dir, "cell-types-vs-disease-t-test.csv"))

big.df %>%
  readr::write_csv(file.path(out.dir, "cell-types-vs-disease-values.csv"))

style <- ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    axis.title.x = ggplot2::element_text(size = 8),
    axis.title.y = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(size = 15),
    plot.subtitle = ggplot2::element_text(size = 8)
  )

figure_2_row_1 <- box.plt + ggtitle("Cellular deconvolution using NNLS (?)") + theme(legend.position = "none")
figure_2_row_2 <- ggpubr::ggarrange(
    celltype_bp[[1]] + style,
    celltype_bp[[2]] + style,
    celltype_bp[[3]] + style,
    labels = c("B", "C", "D"),
    ncol = 3,
    nrow = 1
  )

figure_2 <- ggpubr::ggarrange(
    figure_2_row_1,
    figure_2_row_2,
    labels = c("A", ""),
    ncol = 1,
    nrow = 2
  ) %>%
  ggpubr::annotate_figure(
    top = text_grob(
      "",
      face = "bold",
      size = 14
    )
  )

ggsave(
  filename = file.path(out.img.dir, "Figure_2.png"),
  plot = figure_2,
  width = 20,
  height = 20,
  units = "cm"
)
