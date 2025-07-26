#===============================================================================
# R 脚本：GWAS 可视化分析
#===============================================================================

# 1. 环境设置
#===============================================================================
# 加载必要的包
library(gap)     # 用于QQ图和基因组控制
library(qqman)   # 用于曼哈顿图和QQ图
library(dplyr)   # 用于数据处理

# 设定工作目录
setwd('/mnt/c/Users/Administrator/Desktop/output/')
INPUT_FILE_PLINK1 <- "Han_SNP4_BETA.assoc.logistic" #TODO 替换为使用PLINK1输出的结果文件
INPUT_FILE_PLINK2 <- "Han_SNP4_plink2.PHENO1.glm.logistic.hybrid" #TODO 替换为使用PLINK2输出的结果文件

#TODO 颜色根据自己的喜好进行调整
manhattan_colors <- c("#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#66C2A5","#FC8D62")

# 2. QQ 图分析（使用 gap 包）
# 读取第一个数据集
gwas_data1 <- read.table(INPUT_FILE_PLINK1, header=TRUE)
# 计算基因组控制参数 lambda
gc_result <- gcontrol2(gwas_data1$P)
cat("QQ lambda 值:", gc_result$lambda, "\n")
dev.off()
# 生成 QQ 图
tiff("QQ_V1.tif",
     width=1200,
     height=1200, 
     units="px",
     res=120)
gap::qqunif(gwas_data1$P, main="QQ Plot")
dev.off()

# 3. 曼哈顿图分析
#===============================================================================
# 3.1 第一个曼哈顿图（使用相同数据）
tiff("Manhattan_V1.tif", width=3600, height=900, 
     units="px", res=120)
manhattan(gwas_data1, col=manhattan_colors, 
          main="Manhattan Plot - V1", cex=0.6, cex.axis=0.8)
dev.off()


# 3.2 第二个曼哈顿图（使用新格式数据）
# 读取 PLINK2 格式的数据
gwas_data2 <- read.table(INPUT_FILE_PLINK2, 
                         header=TRUE, 
                         check.names=FALSE,
                         comment.char="")

cat("原始数据列名:", paste(names(gwas_data2), collapse=", "), "\n")

# 标准化列名以符合 qqman 包要求
gwas_data2 <- gwas_data2 |>
  rename(CHR = `#CHROM`,
         BP = POS,
         SNP = ID)

cat("重命名后列名:", paste(names(gwas_data2), collapse=", "), "\n")

# 生成第二个曼哈顿图
tiff("Manhattan_V2.tif", width=3600, height=900, 
     units="px", res=120)
manhattan(gwas_data2, col=manhattan_colors,
          main="Manhattan Plot - PLINK2 Results", cex=0.6, cex.axis=0.8)
dev.off()

# 生成对应的 QQ 图
tiff("QQ_V2.tif", width=1200, height=1200, 
     units="px", res=120)
qqman::qq(gwas_data2$P, main="QQ Plot - PLINK2 Results")
dev.off()


