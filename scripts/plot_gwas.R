#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# GWASFlow 曼哈顿/QQ 图（参数手册 §9.7 契约）
#
# 契约：
# - 图中阈值线数值与该性状 metadata 一致（唯一来源，不另抄）
# - 绘图抽样保留全部显著位点（point-budget 仅控数据量，记录抽样审计）
# - p_underflow 显示上限显式标注（不修改统计表原始值）
# - 染色体顺序来自 chromosome_table 配置
# - 输出 svg/pdf/png 三格式 + 绘图数据

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

parse_args <- function(argv) {
  args <- list(
    sumstats = NULL, metadata = NULL, chromosome_table = NULL,
    trait = NULL, model = NULL, test = NULL,
    point_budget = 200000, preserve_significant = TRUE,
    p_zero_display_cap = "1e-300", formats = "svg,pdf,png",
    seed = 20261023, out_dir = NULL
  )
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (key %in% c("preserve-significant")) {
      args[[gsub("-", "_", key)]] <- tolower(argv[i + 1]) != "false"
      i <- i + 2
    } else {
      if (i + 1 <= length(argv)) {
        args[[gsub("-", "_", key)]] <- argv[i + 1]
        i <- i + 2
      } else {
        stop(paste0("参数缺少值: ", argv[i]))
      }
    }
  }
  args$point_budget <- as.numeric(args$point_budget)
  args$seed <- as.numeric(args$seed)
  args
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$sumstats) || is.null(args$out_dir)) {
    stop("必须提供 --sumstats 与 --out-dir")
  }

  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  gwas <- read.delim(args$sumstats, header = TRUE, check.names = FALSE,
                     na.strings = c("NA"), stringsAsFactors = FALSE)
  meta <- fromJSON(args$metadata)
  chrom_table <- read.delim(args$chromosome_table, header = TRUE,
                            stringsAsFactors = FALSE)

  threshold <- meta$multiple_testing$threshold
  threshold_num <- suppressWarnings(as.numeric(threshold))
  m_tests <- meta$counts$n_valid_p

  # 显示用 -log10p（不改统计表原始值）
  p_display_cap <- as.numeric(args$p_zero_display_cap)
  cap_neglog <- -log10(p_display_cap)
  gwas$p_value_num <- suppressWarnings(as.numeric(gwas$p_value))
  gwas$neglog <- ifelse(is.finite(gwas$p_value_num) & gwas$p_value_num > 0,
                        -log10(gwas$p_value_num), NA)
  is_underflow <- !is.na(gwas$p_value_text) & gwas$p_value_text == "0"
  gwas$neglog_display <- ifelse(is_underflow, cap_neglog, gwas$neglog)
  gwas$status <- gwas$status

  significant <- !is.na(gwas$p_value_num) & !is.na(threshold_num) &
    ((gwas$p_value_num > 0 & gwas$p_value_num <= threshold_num) |
       (is_underflow & !is.na(threshold_num) & threshold_num > 0))
  gwas$significant <- significant

  # 抽样：非显著点超过预算时抽样；显著点全保留（可审计）
  set.seed(args$seed)
  non_sig_idx <- which(!gwas$significant & !is.na(gwas$neglog_display))
  keep_idx <- which(gwas$significant | is.na(gwas$neglog_display))
  if (length(non_sig_idx) > args$point_budget) {
    sampled <- sample(non_sig_idx, args$point_budget)
  } else {
    sampled <- non_sig_idx
  }
  plot_idx <- sort(unique(c(keep_idx, sampled)))
  plot_data <- gwas[plot_idx, ]

  # 染色体顺序（chromosome_table 配置驱动）
  chrom_order <- as.character(chrom_table[[1]])
  plot_data$chrom <- as.character(plot_data$chrom)
  plot_data <- plot_data[plot_data$chrom %in% chrom_order, ]
  plot_data$chrom_factor <- factor(plot_data$chrom, levels = chrom_order)
  plot_data$pos_num <- suppressWarnings(as.numeric(plot_data$pos))

  # 基因组坐标（按染色体表长度累计偏移）
  length_map <- setNames(as.numeric(chrom_table[[2]]), chrom_order)
  offsets <- c(0, cumsum(as.numeric(length_map[chrom_order])))
  names(offsets) <- c(chrom_order, "END")
  plot_data$gpos <- sapply(seq_len(nrow(plot_data)), function(i) {
    chr <- as.character(plot_data$chrom_factor[i])
    if (is.na(chr)) return(NA)
    offsets[chr] + plot_data$pos_num[i]
  })
  chrom_midpoints <- sapply(chrom_order, function(chr) {
    (offsets[chr] + offsets[chr + 1]) / 2
  })

  # ---- 曼哈顿图 ----
  manhattan <- ggplot(
    plot_data[!is.na(plot_data$neglog_display), ],
    aes(x = gpos, y = neglog_display, color = chrom_factor)
  ) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_color_manual(
      values = rep(c("#2C6E9B", "#D1825E"), length.out = length(chrom_order))
    ) +
    scale_x_continuous(breaks = chrom_midpoints, labels = chrom_order, expand = c(0.01, 0)) +
    labs(
      title = paste0("Manhattan plot — trait: ", args$trait,
                     " (model: ", args$model, ")"),
      subtitle = paste0("Bonferroni threshold = ", threshold,
                        " (alpha/(family_size x M), M = ", m_tests, ")"),
      x = "Chromosome", y = expression(-log[10](P)), color = "Chromosome"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none")

  if (!is.na(threshold_num) && threshold_num > 0) {
    manhattan <- manhattan +
      geom_hline(yintercept = -log10(threshold_num), linetype = "dashed",
                 color = "#C0392B", linewidth = 0.4)
  }
  if (any(is_underflow[plot_idx], na.rm = TRUE)) {
    manhattan <- manhattan +
      geom_hline(yintercept = cap_neglog, linetype = "dotted",
                 color = "#7F8C8D", linewidth = 0.3) +
      annotate("text", x = 0, y = cap_neglog, hjust = 0, vjust = -0.5,
               label = paste0("P=0 display cap (", args$p_zero_display_cap,
                              ")；统计表保留原值"),
               size = 2.5, color = "#7F8C8D")
  }

  # ---- QQ 图 ----
  qq_data <- plot_data[!is.na(plot_data$neglog), ]
  qq_data <- qq_data[order(qq_data$p_value_num), ]
  n_qq <- nrow(qq_data)
  expected <- -log10((seq_len(n_qq)) / (n_qq + 1))
  qq_df <- data.frame(
    expected = expected,
    observed = qq_data$neglog_display[order(qq_data$p_value_num)],
    significant = qq_data$significant
  )
  qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(aes(color = significant), size = 0.5, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "#C0392B", linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "#2C6E9B", "TRUE" = "#C0392B")) +
    labs(
      title = paste0("QQ plot — trait: ", args$trait),
      x = expression(Expected ~ -log[10](P)),
      y = expression(Observed ~ -log[10](P)), color = "Significant"
    ) +
    theme_minimal(base_size = 10)

  # ---- 输出三格式 ----
  formats <- strsplit(args$formats, ",")[[1]]
  for (fmt in formats) {
    ggsave(file.path(args$out_dir, paste0("manhattan.", fmt)), manhattan,
           width = 10, height = 5, dpi = 150, limitsize = FALSE)
    ggsave(file.path(args$out_dir, paste0("qq.", fmt)), qq_plot,
           width = 6, height = 6, dpi = 150, limitsize = FALSE)
  }

  # ---- 绘图数据（抽样审计） ----
  plot_export <- data.frame(
    chrom = plot_data$chrom, pos = plot_data$pos_num,
    variant_id = plot_data$variant_id,
    neg_log10_p_display = plot_data$neglog_display,
    p_underflow = is_underflow[plot_idx],
    significant = plot_data$significant,
    stringsAsFactors = FALSE
  )
  write.table(plot_export, file.path(args$out_dir, "plot_data.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  cat(sprintf(
    "[plot_gwas] %s: %d 点绘制（显著 %d 全保留；非显著抽样 %d/%d，预算 %g）\n",
    args$trait, nrow(plot_export), sum(plot_data$significant, na.rm = TRUE),
    length(sampled), length(non_sig_idx), args$point_budget
  ))
}

main()
