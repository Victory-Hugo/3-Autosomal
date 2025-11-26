#!/usr/bin/env Rscript

# 加载常用包，保持与模板脚本一致的写法
library(optparse)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggrastr)


# 颜色配置，与 Jupyter 可视化保持一致
SIG_COLOR <- "#FF5858"
BASE_COLOR <- "#8A8A8A"

# 参数设置
option_list <- list(
  make_option(c("-i", "--input"),
              default="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/output/bugwas_biallelic_lmmout_allSNPs.txt",
              help="GWAS 结果文件路径"),
  make_option(c("-o", "--outdir"),
              default="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/output/plots",
              help="输出目录"),
  make_option(c("-t", "--threshold"),
              default=5, type="double",
              help="-log10(p) 阈值"),
  make_option(c("-f", "--format"),
              default="png",
              help="输出格式：png/pdf/svg"),
  make_option(c("--title"),
              default="GWAS Manhattan Plot",
              help="图标题")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (!file.exists(opt$input)) stop("找不到输入文件: ", opt$input)

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

# 绘制曼哈顿图
manhattan_path <- file.path(opt$outdir, paste0("GWAS_Manhattan.", opt$format))
manhattan_data <- gwas |>
  mutate(point_col = ifelse(negLog10 >= opt$threshold, SIG_COLOR, BASE_COLOR))
y_breaks <- seq(0, ceiling(max(manhattan_data$negLog10, na.rm=TRUE) / 2) * 2, by=2)
if (tail(y_breaks, 1) < max(manhattan_data$negLog10, na.rm=TRUE)) {
  y_breaks <- c(y_breaks, tail(y_breaks, 1) + 2)
}
x_breaks <- seq(0, 1600000, by=200000)
threshold_line <- 5

theme_manhattan <- theme_classic(base_size=10) +
  theme(
    text = element_text(family="Arial"),
    axis.line = element_line(color="black", linewidth=0.4),
    axis.ticks = element_line(color="black", linewidth=0.3),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(hjust=0.5, face="bold"),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  )

manhattan_plot <- ggplot(manhattan_data, aes(x=ps, y=negLog10)) +
  geom_point_rast(aes(color=point_col), size=0.9, alpha=0.9, show.legend=FALSE) +
  scale_color_identity() +
  geom_hline(yintercept=threshold_line, color="black", linetype="dashed", linewidth=0.4) +
  scale_x_continuous(
    breaks=x_breaks,
    expand=c(0, 0)
  ) +
  scale_y_continuous(
    breaks=y_breaks,
    limits=c(0, max(y_breaks)),
    expand=c(0, 0)
  ) +
  labs(x="Genomic position (bp)", y="-log10(p-value)", title=opt$title) +
  theme_manhattan


message("绘制曼哈顿图 -> ", manhattan_path)
ggsave(manhattan_path, manhattan_plot, width=7, height=4, dpi=300, bg="white")

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
    axis.line = element_line(color="#1f6d8f", linewidth=0.4),
    axis.ticks = element_line(color="#4cd2de", linewidth=0.3),
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
