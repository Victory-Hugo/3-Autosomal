#!/usr/bin/env Rscript

# 加载常用包，保持与模板脚本一致的写法
library(optparse)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrastr)
library(stringr)


# 颜色配置，与 Jupyter 可视化保持一致
SIG_COLOR <- "#FF5858"
BASE_COLOR <- "#8A8A8A"

# 参数设置
option_list <- list(
  make_option(c("-i", "--input"),
              default="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/output/bugwas_biallelic_lmmout_allSNPs.txt",
              help="GWAS 结果文件路径"),
  make_option(c("-g", "--gff"),
              default="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/conf/NC_000915.gff",
              help="GFF 注释文件路径"),
  make_option(c("-o", "--outdir"),
              default="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/output/plots",
              help="输出目录"),
  make_option(c("-t", "--threshold"),
              default=5, type="double",
              help="-log10(p) 阈值"),
  make_option(c("-l", "--label_threshold"),
              default=7, type="double",
              help="显著 SNP 标注所需的 -log10(p) 阈值"),
  make_option(c("-f", "--format"),
              default="png",
              help="输出格式：png/pdf/svg"),
  make_option(c("--title"),
              default="GWAS Manhattan Plot",
              help="图标题")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (!file.exists(opt$input)) stop("找不到输入文件: ", opt$input)
if (!file.exists(opt$gff)) stop("找不到 GFF 文件: ", opt$gff)

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# 读取并校验数据
message("读取 GWAS 结果: ", opt$input)
gwas <- read_delim(opt$input, delim="\t", col_types=cols())
required_cols <- c("chr", "rs", "ps", "p_lrt")
missing_cols <- setdiff(required_cols, colnames(gwas))
if (length(missing_cols) > 0) stop("缺少必要列: ", paste(missing_cols, collapse=", "))

# 数据预处理
gwas <- gwas |>
  filter(!is.na(p_lrt), p_lrt > 0, p_lrt <= 1) |>
  mutate(chr = as.character(chr))
if (!"negLog10" %in% names(gwas)) {
  gwas$negLog10 <- -log10(gwas$p_lrt)
}
gwas$negLog10[is.infinite(gwas$negLog10)] <- NA
gwas <- drop_na(gwas, negLog10, ps)
if (nrow(gwas) == 0) stop("清洗后没有可用 SNP。")

palette_colors <- c("#D55E00", "#E39400", "#E0C318", "#14A16A", "#0072B2")
max_pos <- max(gwas$ps, na.rm=TRUE)
physical_limits <- c(0, max(1600000, max_pos))
x_breaks <- seq(0, physical_limits[2], by=200000)

# 绘制曼哈顿图
manhattan_path <- file.path(opt$outdir, paste0("GWAS_Manhattan.", opt$format))
manhattan_data <- gwas
y_breaks <- seq(0, ceiling(max(manhattan_data$negLog10, na.rm=TRUE) / 2) * 2, by=2)
if (tail(y_breaks, 1) < max(manhattan_data$negLog10, na.rm=TRUE)) {
  y_breaks <- c(y_breaks, tail(y_breaks, 1) + 2)
}
threshold_line <- opt$threshold
below_threshold <- manhattan_data |> filter(negLog10 < threshold_line)
above_threshold <- manhattan_data |> filter(negLog10 >= threshold_line)

theme_manhattan <- theme_classic(base_size=10) +
  theme(
    text = element_text(family="Arial"),
    axis.line = element_line(color="black", linewidth=0.4),
    axis.ticks = element_line(color="black", linewidth=0.3),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(hjust=0.5, face="bold"),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8),
    panel.border = element_rect(color="black", fill=NA, linewidth=0.4)
  )

manhattan_core <- ggplot() +
  geom_point_rast(
    data=below_threshold,
    aes(x=ps, y=negLog10),
    color="#B0B0B0",
    size=0.8,
    alpha=0.8,
    show.legend=FALSE
  ) +
  geom_point_rast(
    data=above_threshold,
    aes(x=ps, y=negLog10, color=negLog10),
    size=1.1,
    alpha=0.95,
    show.legend=FALSE
  ) +
  scale_color_gradientn(colors=palette_colors, limits=c(threshold_line, max(manhattan_data$negLog10, na.rm=TRUE))) +
  geom_hline(yintercept=threshold_line, color="black", linetype="dashed", linewidth=0.4) +
  scale_x_continuous(
    breaks=x_breaks,
    limits=physical_limits,
    expand=c(0, 0)
  ) +
  scale_y_continuous(
    breaks=y_breaks,
    limits=c(0, max(y_breaks)),
    expand=c(0, 0)
  ) +
  labs(x="Genomic position (bp)", y="-log10(p-value)", title=opt$title) +
  theme_manhattan

# 读取 GFF 构建注释范围供标注使用
gff_cols <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff_data <- read_delim(
  opt$gff,
  delim="\t",
  comment="#",
  col_names=gff_cols,
  col_types=cols(.default = "c")
)
gene_ranges <- gff_data |>
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    gene_name = str_remove(str_extract(attributes, "gene=[^;]+"), "gene="),
    product = str_remove(str_extract(attributes, "product=[^;]+"), "product=")
  ) |>
  filter(!is.na(start), !is.na(end), type %in% c("gene", "CDS")) |>
  transmute(
    start,
    end,
    product_label = coalesce(product, gene_name)
  )

# 显著 SNP 注释
annotation_df <- NULL
if (nrow(above_threshold) > 0) {
  annotation_df <- above_threshold |>
    rowwise() |>
    mutate(
      product_label = {
        hits <- gene_ranges$product_label[gene_ranges$start <= ps & gene_ranges$end >= ps]
        if (length(hits) == 0) NA_character_ else hits[1]
      }
    ) |>
    ungroup() |>
    mutate(product_label = ifelse(is.na(product_label), paste0("Pos:", ps), product_label))
}

if (!is.null(annotation_df) && nrow(annotation_df) > 0) {
  annotation_df <- annotation_df |> filter(negLog10 >= opt$label_threshold)
  annotation_df <- annotation_df |>
    filter(!is.na(product_label), product_label != "") |>
    arrange(desc(negLog10)) |>
    group_by(product_label) |>
    slice_head(n = 1) |>
    ungroup()
}

if (!is.null(annotation_df) && nrow(annotation_df) > 0) {
  manhattan_core <- manhattan_core +
    geom_label(
      data=annotation_df,
      aes(x=ps, y=negLog10, label=product_label),
      size=2.5,
      color="black",
      fill=NA,
      label.size=NA,
      label.padding=unit(0.05, "lines")
    )
}

message("绘制曼哈顿图 -> ", manhattan_path)
ggsave(manhattan_path, manhattan_core, width=7, height=4, dpi=300, bg="white")

# 绘制 QQ 图
pvals <- sort(gwas$p_lrt)
n <- length(pvals)
expected <- -log10(ppoints(n))
observed <- -log10(pvals)
conf_int <- sapply(seq_len(n), function(i) qbeta(c(0.025, 0.975), i, n - i + 1))
qq_df <- tibble(
  expected = expected,
  observed = observed,
  lower = -log10(conf_int[1, ]),
  upper = -log10(conf_int[2, ])
)

theme_qq <- theme_classic(base_size=10) +
  theme(
    text = element_text(family="Arial"),
    axis.line = element_line(color="#000000", linewidth=0.4),
    axis.ticks = element_line(color="#262626", linewidth=0.3),
    plot.title = element_text(hjust=0.5, face="bold")
  )

qq_plot <- ggplot(qq_df, 
                 aes(x=expected, y=observed)
                 ) +
  geom_ribbon(
              aes(ymin=lower, ymax=upper), 
              fill="#E0E0E0", 
              alpha=0.8
              ) +
  geom_abline(
              slope=1, 
              intercept=0, 
              color="#99dee8", 
              linewidth=0.4) +
  geom_point(
              color="#7AC6DB", 
              size=0.5) +
  labs(
              x="Expected -log10(p)", 
              y="Observed -log10(p)", 
              title="GWAS QQ Plot"
      ) +
  theme_qq

qq_path <- file.path(opt$outdir, paste0("GWAS_QQ.", opt$format))
message("绘制 QQ 图 -> ", qq_path)
ggsave(qq_path, qq_plot, width=4.5, height=4.5, dpi=300, bg="white")

message("全部完成！")
