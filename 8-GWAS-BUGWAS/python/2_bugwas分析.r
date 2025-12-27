library(bugwas)

# BUGWAS 分析脚本
# 安装 bugwas 包（如果尚未安装）：
#   install.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos = NULL, type = "source")
#   install.packages("https://s1-software.oss-cn-chengdu.aliyuncs.com/bugwas_1.0.tar.gz", repos = NULL, type = "source")
# 安装 gemma 软件（如果尚未安装）：
#   conda install -c bioconda gemma

# 命令行参数：
#   1) GFF 文件路径
#   2) 输入目录（包含 geno_biallelic_SNP.txt、pheno.txt、treeFull.nwk）
#   3) 输出目录
#   4) GEMMA 可执行文件路径
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("用法: Rscript 2_bugwas分析.r <GFF_FILE_PATH> <INPUT_DIR> <OUTPUT_DIR> <GEMMA_PATH>")
}

GFF_FILE_PATH <- args[[1]]
INPUT_DIR     <- args[[2]]
OUTPUT_DIR    <- args[[3]]
gem.path      <- args[[4]]

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
setwd(OUTPUT_DIR)

# 输入文件路径
gen    <- file.path(INPUT_DIR, "geno_biallelic_SNP.txt")
pheno  <- file.path(INPUT_DIR, "pheno.txt")
phylo  <- file.path(INPUT_DIR, "treeFull.nwk")
prefix <- "bugwas"  # 分析结果文件前缀

# 简单工具函数：将子目录内容搬到输出目录后删除子目录
move_dir_contents <- function(src_dir, dest_dir) {
  if (!dir.exists(src_dir)) return(invisible(NULL))
  src_files <- list.files(src_dir, full.names = TRUE)
  if (length(src_files) > 0) {
    file.copy(src_files, dest_dir, overwrite = TRUE)
  }
  unlink(src_dir, recursive = TRUE, force = TRUE)
}

# 执行 GWAS 分析（不生成默认图）
data <- lin_loc(
  gen = gen,
  pheno = pheno,
  phylo = phylo,
  prefix = prefix,
  gem.path = gem.path,
  creatingAllPlots = FALSE
)

# BUGWAS 会在输出目录生成额外子目录，统一搬回主输出目录
move_dir_contents(file.path(OUTPUT_DIR, "output"), OUTPUT_DIR)
move_dir_contents(file.path(OUTPUT_DIR, "bugwas_PCloadings"), OUTPUT_DIR)

# 如需保留中间结果，可启用：
# save(data, file = "data_GWAS.RData")

# 生成图像
all_plots(
  biallelic = data$biallelic,
  triallelic = data$triallelic,
  genVars = data$genVars,
  treeInfo = data$treeInfo,
  config = data$config
)
# 再次清理可能在绘图阶段生成的子目录
move_dir_contents(file.path(OUTPUT_DIR, "output"), OUTPUT_DIR)
move_dir_contents(file.path(OUTPUT_DIR, "bugwas_PCloadings"), OUTPUT_DIR)

###################################
## 数据后处理
###################################
setwd(OUTPUT_DIR)
library(ggplot2)
library(readr)

##############################
# 1. 读取并处理 GFF 文件
##############################
if (!file.exists(GFF_FILE_PATH)) {
  stop("GFF 文件不存在，请检查路径！")
}

gff_temp <- read_delim(
  GFF_FILE_PATH,
  delim = "\t",
  escape_double = FALSE,
  col_names = FALSE,
  trim_ws = TRUE,
  skip = 2
)

colnames(gff_temp) <- c(
  "seq_id", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

gff_df <- as.data.frame(gff_temp)
gff_df$start <- as.numeric(gff_df$start)
gff_df$end   <- as.numeric(gff_df$end)

gff_CDS    <- gff_df[gff_df$type == "CDS", ]
gff_nonCDS <- gff_df[gff_df$type != "CDS", ]

##############################
# 2. 读取 GWAS 分析结果
##############################
gwas_file <- file.path(OUTPUT_DIR, "bugwas_biallelic_lmmout_allSNPs.txt")
data <- read.table(
  gwas_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

data$negLog10 <- -log10(data$p_lrt)

write.table(
  data,
  file = gwas_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 如果分隔符不是空格，则需要加参数 sep="\t" 或其他
gwas_data <- read.table(
  gwas_file,
  header = TRUE,
  stringsAsFactors = FALSE
)
