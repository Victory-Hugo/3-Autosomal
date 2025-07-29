#!/usr/bin/env Rscript
# 通过命令行参数传入变量
# 使用方法：
# Rscript script.R <csv_file_path> <src_directory> <dest_alignment_file> <popgen_working_directory> <popgen_output_path>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript script.R <csv_file_path> <src_directory> <dest_alignment_file> <popgen_working_directory> <popgen_output_path>\n")
  quit(status = 1)
}

csv_file_path        <- args[1]
src_directory        <- args[2]
dest_alignment_file  <- args[3]
popgen_working_dir   <- args[4]
popgen_output_path   <- args[5]

# 合并 FASTA 文件的函数
mergeFastaFiles <- function(csv_file, src_dir, dest_file, id_column = "strain") {
  # 读取 CSV 文件，要求包含菌株ID的列（默认列名为 "strain"）
  strains_info <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
  if (!(id_column %in% names(strains_info))) {
    stop(paste("CSV 文件中不存在列：", id_column))
  }
  ids <- strains_info[[id_column]]
  
  # 构建每个 ID 对应的 FASTA 文件路径
  fasta_files <- file.path(src_dir, paste0(ids, ".fasta"))
  
  # 筛选出存在的 FASTA 文件
  existing_files <- fasta_files[file.exists(fasta_files)]
  if (length(existing_files) == 0) {
    stop("没有找到任何存在的 FASTA 文件!")
  }
  
  # 如果目标文件已存在，先删除以免重复写入
  if (file.exists(dest_file)) {
    file.remove(dest_file)
  }
  
  # 构建合并命令，使用 shQuote 确保路径中的特殊字符正确处理
  cmd <- paste("cat", paste(shQuote(existing_files), collapse = " "), ">", shQuote(dest_file))
  cat("执行命令:\n", cmd, "\n")
  
  # 执行合并命令
  system(cmd)
}

# PopGenome 分析函数
performFSTAnalysis <- function(csv_file, popgen_input_dir, popgen_working_dir, popgen_output_file) {
  # 设置工作目录
  setwd(popgen_working_dir)
  
  # 加载 PopGenome 包
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("需要安装 PopGenome 包，请先安装后再运行！")
  }
  library(PopGenome)
  
  # 读取样本信息数据文件
  meta <- read.table(csv_file, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  if (!("fs_pop" %in% names(meta))) {
    stop("CSV 文件中不存在 'fs_pop' 列!")
  }
  # 筛选出属于 hsp1 和 hsp2 的菌株
  meta <- meta[meta$fs_pop %in% c("hsp1", "hsp2"), ]
  
  # 读取核心区域的 SNP 数据（FASTA 格式），指定目录下应仅存放一个 .aln 文件
  data <- readData(popgen_input_dir,
                   populations = FALSE,
                   outgroup = FALSE,
                   include.unknown = FALSE,
                   gffpath = FALSE,
                   format = "fasta",
                   parallized = FALSE,
                   progress_bar_switch = FALSE,
                   FAST = TRUE,
                   big.data = TRUE,
                   SNP.DATA = FALSE)
  
  # 展示基本数据
  cat("个体信息：\n")
  print(get.individuals(data))
  
  cat("数据槽结构：\n")
  show.slots(data)
  
  cat("数据汇总：\n")
  print(get.sum.data(data))
  
  # 设置种群分组（要求 CSV 中包含 'strain' 列）
  if (!("strain" %in% names(meta))) {
    stop("CSV 文件中不存在 'strain' 列!")
  }
  pop1 <- meta$strain[meta$fs_pop == "hsp1"]
  pop2 <- meta$strain[meta$fs_pop == "hsp2"]
  data <- set.populations(data, list(pop1, pop2))
  
  # 滑动窗口转换，用于局部 FST 分析
  slide.data <- sliding.window.transform(data,
                                         width = 1,
                                         jump = 1,
                                         type = 2,
                                         start.pos = FALSE,
                                         end.pos = FALSE,
                                         whole.data = TRUE)
  
  # 计算每个位点的核苷酸水平 FST 值
  slide.data <- F_ST.stats(slide.data, mode = "nucleotide")
  cat("FST 结果结构：\n")
  str(slide.data@nuc.F_ST.pairwise)
  fst_values <- slide.data@nuc.F_ST.pairwise
  
  # 先创建空文件，然后写入结果
  if (!file.exists(popgen_output_file)) {
    if (!file.create(popgen_output_file)) {
      stop("无法创建文件: ", popgen_output_file)
    }
  }
  
  write.table(fst_values, file = popgen_output_file, sep = "\t", row.names = FALSE)
  cat("FST 计算结果已保存到：", popgen_output_file, "\n")
}

# 主函数：整合合并 FASTA 与 PopGenome 分析
runFSTPipeline <- function(csv_file, src_dir, dest_file, popgen_working_dir, popgen_output_file) {
  mergeFastaFiles(csv_file, src_dir, dest_file)
  
  # 假设 PopGenome 的输入目录为目标文件所在目录
  popgen_input_dir <- dirname(dest_file)
  performFSTAnalysis(csv_file, popgen_input_dir, popgen_working_dir, popgen_output_file)
}

# 调用主函数
runFSTPipeline(csv_file_path, src_directory, dest_alignment_file, popgen_working_dir, popgen_output_path)
