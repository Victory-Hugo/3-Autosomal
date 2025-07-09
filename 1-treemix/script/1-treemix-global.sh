#!/bin/bash

# Admixture graphs between the different populations were obtained using the software Treemix v. 1.1212. 
# Treemix was run with a number of migration edges between 0 to 15, with 10 replicates for each number of edge and hpAfrica2 was set as the outgroup. 
# The final number of migration edges was chosen as being the smallest number that allowed 99.8% of the variance to be explained.

# 脚本用途: 使用TreeMix分析全球群体之间的迁移和混合关系(不循环版本)
# 使用方法: bash 1-treemix-global.sh [输入文件] [群体文件] [迁移边数量]

# 确保输入文件存在
if [ $# -lt 2 ]; then
    echo "用法: bash $0 [输入文件] [群体文件] [迁移边数量(可选,默认5)]"
    echo "例如: bash $0 input.frq.gz HX19pops 5"
    exit 1
fi

INPUT_FILE=$1
POP_FILE=$2
MIGRATION_EDGES=${3:-5}  # 默认迁移边数量为5，可通过第三个参数指定

# 检查输入文件
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 $INPUT_FILE 不存在"
    exit 1
fi

if [ ! -f "$POP_FILE" ]; then
    echo "错误: 群体文件 $POP_FILE 不存在"
    exit 1
fi

echo "开始运行全球群体TreeMix分析..."

##======= 计算种群间的FST值 =======
echo "计算种群间的FST值..."
/usr/bin/perl /home/biosoftware/ppgv1/pairwiseFst/pairwise.perl $POP_FILE 0 0

##======= 全球群体分析 =======
echo "进行全球群体分析(不循环)..."
# 全球群体分析, 以hpAfrica2为根, 不循环因为全球群体太多循环会非常慢
echo "运行单次TreeMix分析，迁移边数量m=${MIGRATION_EDGES}..."

OUTPUT_DIR="../output/treemix_results/global"
mkdir -p $OUTPUT_DIR

# 运行TreeMix分析
nohup treemix -i $INPUT_FILE -root hpAfrica2 -o ${OUTPUT_DIR}/GlobalTreemix_m${MIGRATION_EDGES} -m ${MIGRATION_EDGES} -se -bootstrap -k 500 -global -noss &

# 等待TreeMix分析完成
echo "TreeMix分析已启动，运行在后台..."
echo "可以使用以下命令查看进度:"
echo "tail -f ${OUTPUT_DIR}/GlobalTreemix_m${MIGRATION_EDGES}.log"

##======= 结果可视化 =======
# 创建R脚本用于可视化结果
cat << 'EOF' > ${OUTPUT_DIR}/plot_treemix_results.R
# R脚本用于可视化TreeMix全球群体分析结果
library(RColorBrewer)
library(R.utils)

# 检查是否已安装所需的包
required_packages <- c("RColorBrewer", "R.utils")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 加载TreeMix绘图函数
source("https://raw.githubusercontent.com/jyanglab/treemix/master/plotting_funcs.R")

# 设置工作目录(根据需要修改)
# setwd("../output/treemix_results/global")

# 设置输出PDF文件
pdf("treemix_results_plot.pdf", width=12, height=8)

# 读取群体顺序文件(如果存在)
pop_order <- NULL
if(file.exists("pop.order.file")) {
  pop_order <- "pop.order.file"
}

# 参数设置
prefix <- "GlobalTreemix_m5"  # 根据实际输出文件名调整
title <- "Global Population TreeMix Analysis (m=5)"

# 绘制迁移图
plot_tree(prefix, title=title, pop_order=pop_order)

# 绘制残差热图
plot_resid(prefix, pop_order=pop_order)

# 关闭PDF设备
dev.off()

# 保存残差矩阵为CSV文件
cov <- as.matrix(read.table(paste0(prefix, ".cov.gz")))
rownames(cov) <- colnames(cov)
se <- as.matrix(read.table(paste0(prefix, ".covse.gz")))
rownames(se) <- colnames(se)

# 计算标准化残差
resid <- (cov - se)/se
write.csv(resid, file="treemix_residuals.csv")

cat("TreeMix结果可视化已完成，输出文件:\n")
cat("1. treemix_results_plot.pdf - 包含迁移图和残差热图\n")
cat("2. treemix_residuals.csv - 残差矩阵数据\n")
EOF

echo "已创建R脚本用于可视化结果: ${OUTPUT_DIR}/plot_treemix_results.R"
echo "TreeMix分析完成后，可以使用以下命令绘制图形:"
echo "cd ${OUTPUT_DIR} && Rscript plot_treemix_results.R"
