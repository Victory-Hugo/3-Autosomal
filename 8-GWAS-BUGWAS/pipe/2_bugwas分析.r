dir.create("/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/output/临床指标/临床海拔")
setwd("/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/output/临床指标/临床海拔")  # 设置工作目录

library(bugwas)
# install.packages("https://raw.githubusercontent.com/sgearle/bugwas/master/build/bugwas_1.0.tar.gz", repos=NULL, type="source")
gen <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/input/临床指标/临床海拔/geno_biallelic_SNP.txt"  # 基因型文件,通过/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/script/s4_input_bugwas.py生成
pheno <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/input/临床指标/临床海拔/pheno.txt"  # 表型文件，类似如下
phylo <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/input/临床指标/临床海拔/treeFull.nwk"  # 系统发育树文件
prefix <- "bugwas"  # 分析结果文件前缀
gem.path <- "/home/luolintao/miniconda3/bin/gemma"  # GEMMA软件路径

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
setwd("/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/output/临床指标/临床海拔")
library(ggplot2)
library(readr)

##############################
# 1. 读取并处理 GFF 文件
##############################
# 检查 GFF 文件是否存在
gff_file <- "/mnt/d/幽门螺旋杆菌/Annotation/参考序列/NC_000915.gff"
if (!file.exists(gff_file)) {
  stop("GFF 文件不存在，请检查路径！")
}

# 读取 GFF 文件，跳过前两行注释
gff_temp <- read_delim(
  gff_file,
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
data <- read.table('/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/output/临床指标/临床海拔/bugwas_biallelic_lmmout_allSNPs.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 针对每一行，计算 -log10(p_lrt) 并替换原本的 negLog10 列
data$negLog10 <- -log10(data$p_lrt)

# 将处理后的数据写入新的文件中（使用制表符分隔，不包含行名）
write.table(data, file = '/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/output/临床指标/临床海拔/bugwas_biallelic_lmmout_allSNPs.txt', sep = "\t", row.names = FALSE, quote = FALSE)
# 如果分隔符不是空格，则需要加参数 sep="\t" 或其他
gwas_data <- read.table("/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/output/临床指标/临床海拔/bugwas_biallelic_lmmout_allSNPs.txt", 
                        header = TRUE, stringsAsFactors = FALSE)

# # 绘制曼哈顿图：横轴 = ps，纵轴 = negLog10
# plot(gwas_data$ps, gwas_data$negLog10, 
#      xlab = "Position", ylab = "-log10(p-value)", main = "Manhattan Plot")
# abline(h = 10, col = "green")

#########################################
# 3. 按基因组区间划分并绘图
#########################################
# 根据你所需的区间遍历，示例为 0~160*10000，每次步长 5 => 5万bp为窗口
# for (i in seq(from = 0, to = 160, by = 5)) {
#     xmin <- i * 10000
#     xmax <- xmin + 50000
    
#     tmp <- gwas_data[gwas_data$negLog10 != Inf & 
#                     gwas_data$ps > xmin & 
#                     gwas_data$ps <= xmax, ]
    
#     if ( any(tmp$negLog10 > 5) ) {
#       cat("Generating plot for region:", xmin, "-", xmax, "\n")
    
#     # 从 CDS 注释中取到该区间的基因
#     # 注意，这里的条件可能要改成 >= 或 <=；视实际需求而定
#     region_CDS <- gff_CDS[gff_CDS$start > xmin & gff_CDS$end < xmax, ]
    
#     # 画图
#     p <- ggplot(tmp, aes(x = ps / 1000, y = negLog10)) +
#       geom_point() +
#       geom_hline(yintercept = 5, col = "green") +
#       geom_segment(data = region_CDS, 
#                    aes(x = start / 1000, y = -2, 
#                        xend = end / 1000, yend = -2),
#                    color = "blue") +
#       labs(x = "位置 (kb)", y = "-log10(p-value)") +
#       ggtitle(paste0("Region: ", xmin, " - ", xmax))
    
#     # 保存图像
#     out_png <- paste0("xxx", xmin, "_", xmax, ".png")
#     ggsave(out_png, plot = p, width = 19.6, height = 7.64)
#   }
# }


#########################################
#! 5. 筛选显著 SNP 并与基因注释匹配，使用/mnt/d/幽门螺旋杆菌/Script/分析结果/GWAS/script/3_基因注释.ipynb
#########################################
# # 筛选 -log(p-value) 大于 5 的显著 SNP（排除 logp 为 Inf 的记录）
# sig <- gwas_data[gwas_data$logp > 5 & gwas_data$logp != Inf, ]

# # 计算相邻 SNP 之间的距离（用于进一步过滤）
# delta <- c(0, abs(sig$position[-nrow(sig)] - sig$position[-1]))
# sig <- sig[delta < 1000, ]  # 保留距离小于 1kb 的 SNP

# # 初始化用于存储匹配到 CDS 的结果
# sigCDS <- NULL
# for (x in sig$position) {
#   # 从 CDS 注释中筛选出 SNP 所在区间（start < x 且 end > x）
#   matched <- gff_CDS[gff_CDS$start < x & gff_CDS$end > x, ]
#   if (nrow(matched) > 0) {
#     sigCDS <- rbind(sigCDS, matched)
#   }
# }
# # 去重
# sigCDS <- sigCDS[!duplicated(sigCDS), ]

# # 从 attributes 字段中提取 product 信息
# sigCDS$product <- unlist(lapply(sigCDS$attributes, function(x) {
#   parts <- strsplit(x, split = "product=")[[1]]
#   if (length(parts) > 1) {
#     return(parts[2])
#   } else {
#     return(NA)
#   }
# }))

# # 按照 CDS 的起始位置排序并保存结果
# write.table(sigCDS[order(sigCDS$start), c("start", "end", "product")],
#             file = "significantCDS.csv", quote = FALSE, row.names = FALSE, 
#             col.names = FALSE, sep = ",")
# write.table(sig, file = "significantSNPs.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)


