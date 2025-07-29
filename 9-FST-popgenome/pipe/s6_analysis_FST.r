# 设置文件路径
csv_file <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/FST/output/临床指标/AML-阿莫西林/META.csv"
src_dir <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/merged_fasta/Whole_Genome_NoInDel"
dest_file <- "/mnt/d/幽门螺旋杆菌/Script/分析结果/FST/input/temp.aln"

# 读取CSV文件（文件中包含标题行）
strains_info <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
# 提取第一列的ID（假设列名为"strain"）
ids <- strains_info$strain
# 构建每个ID对应的fasta文件的完整路径
fasta_files <- file.path(src_dir, paste0(ids, ".fasta"))


# 筛选出存在的fasta文件
existing_files <- fasta_files[file.exists(fasta_files)]
if (length(existing_files) == 0) {
  stop("没有找到任何存在的fasta文件!")
}

# 如果目标文件已存在，先删除以免重复写入
if (file.exists(dest_file)) {
  file.remove(dest_file)
}

# 构建cat命令，将所有存在的文件一次性合并到目标文件中
# 使用shQuote确保路径中的特殊字符被正确处理
cmd <- paste("cat", paste(shQuote(existing_files), collapse = " "), ">", shQuote(dest_file))
cat("执行命令:\n", cmd, "\n")

# 执行命令
system(cmd)


## R脚本，用于计算以下比较的FST值（种群分化指数）：
##################################################
## 仅使用双等位SNP进行计算
##################################################
# #todo 安装程序
# install.packages("devtools")
# library(devtools)
# devtools::install_github("pievos101/PopGenome")

# install.packages("路径/PopGenome_x.x.x.tar.gz", repos = NULL, type = "source")
# # 或者使用
# library(devtools)
# devtools::install_local("/mnt/c/Users/Administrator/Desktop/PopGenome_2.7.7.tar.gz")
# # 设置工作目录
setwd("/mnt/d/幽门螺旋杆菌/Script/分析结果/FST/output/临床指标/AML-阿莫西林")  # 设置工作目录
# 加载PopGenome包
library(PopGenome)
# 读取样本信息数据文件
meta <- read.table(csv_file, sep = ",", header = TRUE)  # 读取菌株信息
meta <- meta[meta$fs_pop %in% c("hsp1", "hsp2"),]  # 筛选出属于hsp1和hsp2的菌株

## 使用双等位核心SNP（包括CDS和基因间区域）进行FST计算
## SNP集合与GWAS分析中使用的相同
#* 读取核心区域的SNP数据（FASTA格式）
#! 给出路径即可！该路径下存放1个aln文件，其中包括多条序列！不要存放任何其他文件在该目录下！
data <- readData('/mnt/d/幽门螺旋杆菌/Script/分析结果/FST/input/',
        populations=FALSE,outgroup=FALSE,include.unknown=FALSE,
        gffpath=FALSE,format="fasta",parallized=FALSE,
        progress_bar_switch=FALSE, FAST=TRUE,big.data=TRUE,
        SNP.DATA=FALSE) # 读取核心区域的SNP数据（FASTA格式）

#* 统计数据
get.individuals(data)
#* 展示数据分析的各种方法
show.slots(data) 
#* 统计数据
get.sum.data(data)
# write.csv(get.sum.data(data), file = "/mnt/c/Users/victo/Desktop/1.csv")

#todo 先计算一下Pi和TajimD值
basic_data <- F_ST.stats(data,
                        FAST=FALSE,
                        mode="haplotype",  # mode="haplotype" or mode="nucleotide"
                        ) # If FAST is switched on, this module only calculates nuc.diversity.within, hap.diversity.within, haplotype.F_ST, nucleotide.F_ST and pi.
basic_data@n.sites
basic_data@Pi
basic_data@haplotype.F_ST
basic_data <- neutrality.stats(data)
basic_data@Tajima.D 


## 设置种群，基于菌株信息中的fs_pop分组
data <- set.populations(data, list(
  meta$strain[meta$fs_pop == "hsp1"],  # hsp1菌株
  meta$strain[meta$fs_pop == "hsp2"]   # hsp2菌株
))




## 如果需要每个位点的FST值
## 需要创建窗口大小为1的滑动窗口
## （在每个窗口内计算相关统计量）
slide.data <- sliding.window.transform(data, 
                                      width=1, 
                                      jump=1, 
                                      type=2, # 1 scan only biallelic positions (SNPs), 2 scan the genome. default:1
                                      start.pos=FALSE,
                                      end.pos=FALSE, 
                                      whole.data=TRUE) 

# slide.data <- diversity.stats(slide.data)  # 计算种群内多样性
# str(slide.data@nuc.diversity.within)  # 查看种群内多样性结果的结构
# x <- slide.data@nuc.diversity.within[,2]  # 提取第二列多样性数据

## 计算每个位点的FST值
slide.data <- F_ST.stats(slide.data, mode = "nucleotide")  # 计算核苷酸水平的FST统计量
 str(slide.data@nuc.F_ST.pairwise)  # 查看FST结果的结构
x <- slide.data@nuc.F_ST.pairwise  # 提取FST值

## 保存滑动窗口分析的结果
# save(slide.data, file = "FST_highLow.RData")  # 将结果保存为RData文件
write.table(x, file = "/mnt/d/幽门螺旋杆菌/Script/分析结果/FST/output/hspLAEAsia1_hspLAEAsia3/FST_highLow.txt", sep = "\t", row.names = FALSE)
