#* BUGWAS分析脚本
#* 安装bugwas包（如果尚未安装）
# 使用备份：install.packages("https://1835545223.v.123pan.cn/1835545223/24078748", repos=NULL, type="source")
# install.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")
#* 安装gemma软件（如果尚未安装）
# conda install -c bioconda gemma


library(bugwas)
# 定义输入输出目录
GFF_FILE_PATH <- "/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/conf/NC_000915.gff"
BASE_DIR <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/8-临床结局相关"
INPUT_DIR <- file.path(BASE_DIR, "data/GWAS/Atrophy-萎缩")
OUTPUT_DIR <- file.path(BASE_DIR, "output/GWAS/Atrophy-萎缩")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
setwd(OUTPUT_DIR)

# 定义输入文件路径
gen <- file.path(INPUT_DIR, "geno_biallelic_SNP.txt")
pheno <- file.path(INPUT_DIR, "pheno.txt")
phylo <- file.path(INPUT_DIR, "treeFull.nwk")
prefix <- "bugwas"  # 分析结果文件前缀
gem.path <- "/home/luolintao/miniconda3/envs/BigLin/bin/gemma"  # GEMMA软件路径

data <- lin_loc(
  gen = gen, 
  pheno = pheno, 
  phylo = phylo, 
  prefix = prefix, 
  gem.path = gem.path, 
  creatingAllPlots = FALSE
)  # 执行GWAS分析

# save(data, file = "data_GWAS.RData")  # 保存GWAS分析结果

all_plots(
  biallelic = data$biallelic, 
  triallelic = data$triallelic, 
  genVars = data$genVars, 
  treeInfo = data$treeInfo, 
  config = data$config
)  # 生成相关图像

###################################
## 数据后处理
###################################
setwd(OUTPUT_DIR)
library(ggplot2)
library(readr)

##############################
# 1. 读取并处理 GFF 文件
##############################
# 检查 GFF 文件是否存在
if (!file.exists(GFF_FILE_PATH)) {
  stop("GFF 文件不存在，请检查路径！")
}

# 读取 GFF 文件，跳过前两行注释
gff_temp <- read_delim(
  GFF_FILE_PATH,
  delim = "\t",
  escape_double = FALSE,
  col_names = FALSE,
  trim_ws = TRUE,
  skip = 2
)


# 指定列名
colnames(gff_temp) <- c("seq_id", "source", "type", "start", "end", 
                        "score", "strand", "phase", "attributes")

# 转换为数据框并确保坐标列为数值型
gff_df <- as.data.frame(gff_temp)
gff_df$start <- as.numeric(gff_df$start)
gff_df$end   <- as.numeric(gff_df$end)

# 拆分出 CDS 与非 CDS
gff_CDS   <- gff_df[gff_df$type == "CDS", ]
gff_nonCDS <- gff_df[gff_df$type != "CDS", ]

##############################
# 2. 读取 GWAS 分析结果
##############################
# 确保你的 "bugwas_biallelic_lmmout_allSNPs.txt" 中存在列 ps、negLog10
# 读取数据（注意根据实际文件分隔符调整 sep 参数，这里假定为制表符）
data <- read.table(file.path(OUTPUT_DIR, "bugwas_biallelic_lmmout_allSNPs.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 针对每一行，计算 -log10(p_lrt) 并替换原本的 negLog10 列
data$negLog10 <- -log10(data$p_lrt)

# 将处理后的数据写入新的文件中（使用制表符分隔，不包含行名）
write.table(data, file = file.path(OUTPUT_DIR, "bugwas_biallelic_lmmout_allSNPs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
# 如果分隔符不是空格，则需要加参数 sep="\t" 或其他
gwas_data <- read.table(file.path(OUTPUT_DIR, "bugwas_biallelic_lmmout_allSNPs.txt"), 
                        header = TRUE, stringsAsFactors = FALSE)
