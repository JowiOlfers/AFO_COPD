#@ _helpers.r
library(biomaRt)

####
# Helper functions
####
getGenedataByEnsemblId <- function(ensemblIds = NULL, file.location = ".", hg_version = 38) {
 
  file.name <- file.path(file.location, paste0("genes_info_hg", hg_version, ".csv"))
  print(file.name)
  if (!file.exists(file.name)) {
    if (!("mart" %in% ls())) {
      args <- list(
        biomaRt = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl"
      )

      if (hg_version != 38) {
        args <- append(args, list(GRCh = hg_version))
      }

      assign("mart", do.call(biomaRt::useEnsembl, args))
    }

    args <- list(
      attributes = c(
        "hgnc_symbol",
        "ensembl_gene_id"  ),
      mart = mart
    )

    if (!is.null(ensemblIds)) {
      args <- append(args,
        list(
          filters = "ensembl_gene_id",
          values = as.character(ensemblIds)
        )
      )
    }
    gene.list <- do.call(biomaRt::getBM, args)
    readr::write_csv(gene.list, file = file.name)
  }
  return(
    readr::read_csv(
      file.name,
      col_types = readr::cols()
    )
  )
}

select.columns.in.order <- function(dataframe, columns) {
  dataframe[, columns]
}

drop.columns.if.all.same.value <- function(dataframe) {
  for (name in colnames(dataframe)) {
    is.all.same <- (dataframe[, name] %>% unique() %>% length()) <= 1
    if (is.all.same) {
      dataframe <- dataframe %>%
        dplyr::select(
          -tidyselect::one_of(name)
        )
    }
  }
  dataframe
}

remove.empty.lists <- function(src.list) {
  dest.list <- list()
  for (name in names(src.list)) {
    if (length(src.list[[name]]) > 0) {
      dest.list[[name]] <- src.list[[name]]
    }
  }
  dest.list
}

intersect.all <- function(lists) {
  intersecting <- NA
  for (name in names(lists)) {
    if ((!is.na(intersecting))[1]) {
      intersecting <- intersect(intersecting, lists[[name]])
    } else {
      intersecting <- lists[[name]]
    }
  }
  intersecting
}

unique.all <- function(lists) {
  unique_ <- NA
  for (name in names(lists)) {
    if ((!is.na(unique_))[1]) {
      unique_ <- unique_[!(unique_ %in% intersect(unique_, lists[[name]]))]
    } else {
      unique_ <- lists[[name]]
    }
  }
  unique_
}

meta_zscore <- function(fc1, fc2, pval1, pval2, n1, n2) {
  nsum <- n1 + n2
  w1 <- sqrt(n1) / sqrt(nsum);
  w2 <- sqrt(n2) / sqrt(nsum);

  convert.pvalue.to.zscore <- function(pval, beta) {
    if ( beta > 0 ) {
      z <- qnorm( pval / 2 );
    } else {
      z <- -(qnorm( pval / 2 ));
    }

    return(z);
  }

  p.to.z <- function(pval1, pval2, fc1, fc2, w1, w2) {
    z1 <- convert.pvalue.to.zscore(pval1, fc1);
    z2 <- convert.pvalue.to.zscore(pval2, fc2);
    zsum <- ( w1*z1 ) + ( w2*z2 )

    return(zsum);
  }

  zmeta <- p.to.z(pval1, pval2, fc1, fc2, w1, w2);
  pmeta <- pnorm(-(abs(zmeta))) * 2;
  return(pmeta)
}


meta_zscore_df <- function(df, l.p.val, r.p.val, l.fc, r.fc, l.n, r.n) {
  for (i in 1:nrow(df)) {
    row <- df[i, ]

    if (!is.na(row[l.p.val]) && !is.na(row[r.p.val])) {
      df[i, "meta.p"] <- meta_zscore(
        as.numeric(row[l.fc]), as.numeric(row[r.fc]),
        as.numeric(row[l.p.val]), as.numeric(row[r.p.val]),
        as.numeric(l.n), as.numeric(r.n)
      )
    }
  }

  df[, "meta.p.adj.bh"] <- p.adjust(
    p = df$meta.p,
    method = "BH",
    n = nrow(df)
  )

  df
}

