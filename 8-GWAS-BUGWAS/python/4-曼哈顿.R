#!/usr/bin/env Rscript

# 加载常用包，保持与模板脚本一致的写法
library(optparse)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrastr)
library(stringr)
library(ggrepel)
library(patchwork)


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
  make_option(c("-s", "--fst"),
              default="",
              help="FST 结果文件路径，留空则跳过 FST 绘图"),
  make_option(c("-o", "--outdir"),
              default="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/output/plots",
              help="输出目录"),
  make_option(c("-t", "--threshold"),
              default=5, type="double",
              help="GWAS -log10(p) 阈值"),
  make_option(c("-l", "--label_threshold"),
              default=7, type="double",
              help="GWAS 标注阈值 (-log10(p))"),
  make_option(c("--fst_lower"),
              default=0.4, type="double",
              help="Fst 低阈值"),
  make_option(c("--fst_upper"),
              default=0.6, type="double",
              help="Fst 高阈值"),
  make_option(c("--fst_label"),
              default=0.7, type="double",
              help="Fst 标注阈值"),
  make_option(c("-f", "--format"),
              default="pdf",
              help="输出格式：png/pdf/svg"),
  make_option(c("--title"),
              default="GWAS Manhattan Plot",
              help="图标题")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (!file.exists(opt$input)) stop("找不到输入文件: ", opt$input)
if (!file.exists(opt$gff)) stop("找不到 GFF 文件: ", opt$gff)
has_fst <- nzchar(opt$fst)
if (has_fst && !file.exists(opt$fst)) stop("找不到 FST 文件: ", opt$fst)

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

ggsave_wrapper <- function(filename, plot_obj, width, height) {
  device_fun <- if (tolower(opt$format) == "pdf") cairo_pdf else NULL
  ggsave(
    filename,
    plot_obj,
    width=width,
    height=height,
    dpi=300,
    bg="white",
    device=device_fun
  )
}

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

# 读取并处理 FST 结果（可选）
fst_data <- NULL
if (has_fst) {
  fst_raw <- read_delim(opt$fst, delim=",", col_types=cols())
  fst_data <- fst_raw |>
    rename(Location = 1, Fst = 2) |>
    mutate(
      Location = as.numeric(Location),
      Fst = as.numeric(Fst)
    ) |>
    drop_na(Location, Fst)
  if (nrow(fst_data) == 0) {
    warning("FST 文件没有有效数据，跳过 FST 绘图。")
    fst_data <- NULL
  }
}

palette_colors <- c("#0072B2", "#14A16A", "#E0C318", "#E39400", "#D55E00")
if (!is.null(fst_data)) {
  max_pos <- max(gwas$ps, fst_data$Location, na.rm=TRUE)
} else {
  max_pos <- max(gwas$ps, na.rm=TRUE)
}
physical_limits <- c(0, max(1600000, max_pos))
x_breaks <- seq(0, physical_limits[2], by=200000)
if (!is.null(fst_data)) {
  fst_data <- fst_data |> filter(Location >= physical_limits[1], Location <= physical_limits[2])
}
fst_lower <- opt$fst_lower
fst_upper <- opt$fst_upper
fst_label_threshold <- opt$fst_label

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

fst_theme <- theme_classic(base_size=10) +
  theme(
    text = element_text(family="Arial"),
    axis.line = element_line(color="black", linewidth=0.4),
    axis.ticks = element_line(color="black", linewidth=0.3),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(hjust=0.5, face="bold"),
    panel.border = element_rect(color="black", fill=NA, linewidth=0.4),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  )

fst_plot <- NULL
fst_annotation <- NULL
if (!is.null(fst_data)) {
  fst_ymax <- max(max(fst_data$Fst, na.rm=TRUE), fst_upper) * 1.05
  fst_below <- fst_data |> filter(Fst < fst_lower)
  fst_above <- fst_data |> filter(Fst >= fst_lower)
  fst_color_upper <- if (nrow(fst_above) > 0) max(fst_above$Fst, na.rm=TRUE) else fst_upper
} else {
  fst_ymax <- NULL
  fst_below <- NULL
  fst_above <- NULL
  fst_color_upper <- fst_upper
}

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

if (!is.null(fst_data)) {
  fst_high <- fst_data |> filter(Fst >= fst_label_threshold)
  if (nrow(fst_high) > 0) {
    fst_annotation <- fst_high |>
      rowwise() |>
      mutate(
        product_label = {
          hits <- gene_ranges$product_label[gene_ranges$start <= Location & gene_ranges$end >= Location]
          if (length(hits) == 0) NA_character_ else hits[1]
        }
      ) |>
      ungroup() |>
      filter(!is.na(product_label), product_label != "") |>
      arrange(desc(Fst)) |>
      group_by(product_label) |>
      slice_head(n = 1) |>
      ungroup()
  }
}

if (!is.null(fst_data)) {
  fst_plot <- ggplot() +
    geom_point_rast(
      data=fst_below,
      aes(x=Location, y=Fst),
      color="#B0B0B0",
      size=0.8,
      alpha=0.8
    ) +
    geom_point_rast(
      data=fst_above,
      aes(x=Location, y=Fst, color=Fst),
      size=1.0,
      alpha=0.95,
      show.legend=FALSE
    ) +
    scale_x_continuous(
      breaks=x_breaks,
      limits=physical_limits,
      expand=c(0, 0)
    ) +
    scale_y_continuous(
      limits=c(0, 1),
      breaks=seq(0, 1, 0.2),
      expand=c(0, 0)
    ) +
    scale_color_gradientn(
      colors=palette_colors,
      limits=c(fst_lower, fst_color_upper)
    ) +
    geom_hline(yintercept=fst_lower, color="black", linetype="dashed", linewidth=0.4) +
    geom_hline(yintercept=fst_upper, color="black", linetype="dotted", linewidth=0.4) +
    labs(x="Genomic position (bp)", y="Fst", title="Fst across genome") +
    fst_theme

  if (!is.null(fst_annotation) && nrow(fst_annotation) > 0) {
    fst_plot <- fst_plot +
      geom_label_repel(
        data=fst_annotation,
        aes(x=Location, y=Fst, label=product_label),
        size=2.5,
        color="black",
        fill=NA,
        label.size=0,
        box.padding=0.3,
        point.padding=0.2,
        min.segment.length=0,
        segment.color="black",
        segment.size=0.3
      )
  }
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
    geom_label_repel(
      data=annotation_df,
      aes(x=ps, y=negLog10, label=product_label),
      size=2.5,
      color="black",
      fill=NA,
      label.size=0,
      box.padding=0.3,
      point.padding=0.2,
      min.segment.length=0,
      segment.color="black",
      segment.size=0.3
    )
}

if (!is.null(fst_plot)) {
  combined_plot <- manhattan_core / fst_plot + plot_layout(heights = c(1, 1))
  message("绘制曼哈顿+Fst 图 -> ", manhattan_path)
  ggsave_wrapper(manhattan_path, combined_plot, width=7, height=7)
} else {
  message("绘制曼哈顿图 -> ", manhattan_path)
  ggsave_wrapper(manhattan_path, manhattan_core, width=7, height=4.5)
}

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

qq_plot <- ggplot(qq_df, aes(x=expected, y=observed)) +
  geom_ribbon(
    aes(ymin=lower, ymax=upper),
    fill="#E0E0E0",
    alpha=0.8
  ) +
  geom_abline(
    slope=1,
    intercept=0,
    color="#99dee8",
    linewidth=0.4
  ) +
  geom_point_rast(
    color="#7AC6DB",
    size=0.8
  ) +
  labs(
    x="Expected -log10(p)",
    y="Observed -log10(p)",
    title="GWAS QQ Plot"
  ) +
  theme_qq

qq_path <- file.path(opt$outdir, paste0("GWAS_QQ.", opt$format))
message("绘制 QQ 图 -> ", qq_path)
ggsave_wrapper(qq_path, qq_plot, width=4.5, height=4.5)

message("全部完成！")
