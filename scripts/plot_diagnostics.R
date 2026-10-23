#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# GWASFlow 诊断图（参数手册 §9.7 契约）
#
# 按需消费真实产物：
#   --pcs + --kinship        → 项目级 pca.png + kinship_heatmap.png
#   --phenotype + --sample-inclusion → 每性状 phenotype_distribution.png + sample_inclusion.png
# 上游产物缺失（模块关闭）→ 由规则层跳过该图并在 plot_manifest 标注原因。

suppressPackageStartupMessages({
  library(ggplot2)
})

parse_args <- function(argv) {
  args <- list(trait = NULL, out_dir = NULL, pcs = NULL, kinship = NULL,
               phenotype = NULL, sample_inclusion = NULL, unit = "unspecified")
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (i + 1 <= length(argv)) {
      args[[gsub("-", "_", key)]] <- argv[i + 1]
      i <- i + 2
    } else {
      stop(paste0("参数缺少值: ", argv[i]))
    }
  }
  args
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$out_dir)) stop("必须提供 --out-dir")
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  made <- c()

  # ---- PCA 图（项目级）----
  if (!is.null(args$pcs) && file.exists(args$pcs)) {
    pca <- read.table(args$pcs, header = TRUE, comment.char = "#",
                      stringsAsFactors = FALSE)
    colnames(pca)[1:2] <- c("FID", "IID")
    pc_cols <- grep("^PC", colnames(pca), value = TRUE)
    if (length(pc_cols) >= 2) {
      pca$group <- as.factor(pca$FID)
      p <- ggplot(pca, aes_string(x = pc_cols[1], y = pc_cols[2], color = "group")) +
        geom_point(size = 1.5, alpha = 0.8) +
        labs(title = "PCA (LD-pruned marker set)",
             x = pc_cols[1], y = pc_cols[2], color = "FID") +
        theme_minimal(base_size = 10)
      ggsave(file.path(args$out_dir, "pca.png"), p, width = 6, height = 5, dpi = 150)
      made <- c(made, "pca.png")
    }
  }

  # ---- K 热图（项目级） ----
  if (!is.null(args$kinship) && file.exists(args$kinship)) {
    k <- as.matrix(read.table(args$kinship, header = FALSE, comment.char = "#"))
    n_show <- min(nrow(k), 200)
    idx <- seq_len(n_show)
    k_sub <- k[idx, idx]
    k_df <- data.frame(
      row = rep(idx, times = n_show),
      col = rep(idx, each = n_show),
      value = as.vector(k_sub)
    )
    p <- ggplot(k_df, aes(x = col, y = row, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "#F5F5F5", high = "#2C6E9B") +
      labs(title = paste0("Kinship heatmap (first ", n_show, " samples)"),
           x = "Sample index", y = "Sample index", fill = "K") +
      theme_minimal(base_size = 10) +
      theme(axis.text = element_blank())
    ggsave(file.path(args$out_dir, "kinship_heatmap.png"), p,
           width = 6, height = 5, dpi = 150)
    made <- c(made, "kinship_heatmap.png")
  }

  # ---- 表型分布（每性状）----
  if (!is.null(args$phenotype) && file.exists(args$phenotype)) {
    pheno <- scan(args$phenotype, what = character(), quiet = TRUE)
    values <- suppressWarnings(as.numeric(pheno))
    valid <- values[!is.na(values)]
    if (length(valid) > 0) {
      df <- data.frame(phenotype = valid)
      p <- ggplot(df, aes(x = phenotype)) +
        geom_histogram(bins = 30, fill = "#2C6E9B", alpha = 0.8, color = "white") +
        labs(
          title = paste0("Phenotype distribution — ", args$trait),
          subtitle = paste0("N valid = ", length(valid), "；unit: ", args$unit),
          x = paste0(args$trait, " [", args$unit, "]"), y = "Count"
        ) +
        theme_minimal(base_size = 10)
      ggsave(file.path(args$out_dir, "phenotype_distribution.png"), p,
             width = 6, height = 4, dpi = 150)
      made <- c(made, "phenotype_distribution.png")
    }
  }

  # ---- 样本纳入概况（每性状）----
  if (!is.null(args$sample_inclusion) && file.exists(args$sample_inclusion)) {
    incl <- read.delim(args$sample_inclusion, header = TRUE, stringsAsFactors = FALSE)
    incl$included <- tolower(as.character(incl$included)) == "true"
    reason_col <- if ("exclusion_reason" %in% colnames(incl)) {
      ifelse(incl$included, "included",
             paste0("excluded: ", incl$exclusion_reason))
    } else {
      ifelse(incl$included, "included", "excluded")
    }
    df <- data.frame(category = reason_col, stringsAsFactors = FALSE)
    df$category <- factor(df$category)
    p <- ggplot(df, aes(x = category, fill = category)) +
      geom_bar() +
      labs(title = paste0("Sample inclusion — ", args$trait),
           x = "Status", y = "Count") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      theme(legend.position = "none")
    ggsave(file.path(args$out_dir, "sample_inclusion.png"), p,
           width = 7, height = 4, dpi = 150)
    made <- c(made, "sample_inclusion.png")
  }

  cat(sprintf("[plot_diagnostics] %s: 生成 %s\n", args$trait,
              if (length(made)) paste(made, collapse = ", ") else "0 张图（输入未提供）"))
}

main()
