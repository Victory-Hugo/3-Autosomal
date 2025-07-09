#!/bin/bash

# Admixture graphs between the different populations were obtained using the software Treemix v. 1.1212. 
# Treemix was run with a number of migration edges between 0 to 15, with 10 replicates for each number of edge and hpAfrica2 was set as the outgroup. 
# The final number of migration edges was chosen as being the smallest number that allowed 99.8% of the variance to be explained.

# 脚本用途: 使用TreeMix分析中国及周边群体之间的迁移和混合关系(循环版本)
# 使用方法: bash 1-treemix-china.sh [输入文件] [群体文件]

# 确保输入文件存在
if [ $# -lt 2 ]; then
    echo "用法: bash $0 [输入文件] [群体文件]"
    echo "例如: bash $0 input.frq.gz HX19pops"
    exit 1
fi

INPUT_FILE=$1
POP_FILE=$2

# 检查输入文件
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 $INPUT_FILE 不存在"
    exit 1
fi

if [ ! -f "$POP_FILE" ]; then
    echo "错误: 群体文件 $POP_FILE 不存在"
    exit 1
fi

echo "开始运行TreeMix分析..."

##======= 计算种群间的FST值 =======
echo "计算种群间的FST值..."
/usr/bin/perl /home/biosoftware/ppgv1/pairwiseFst/pairwise.perl $POP_FILE 0 0

##======= 中国及周边群体分析 =======
echo "进行中国及周边群体分析(循环寻找最优模型)..."
# 中国及周边群体分析, 以非洲为根, 循环测试不同m值以寻找最优模型
##======= 多次分析以评估最佳m值=======
## 比如m取1-10(常用1-5,1-10)，每个m值重复5次(至少两次)

# 可以通过参数指定m值的范围
MIN_M=${3:-0}
MAX_M=${4:-10}
REPEATS=${5:-5}

echo "迁移边数量范围: $MIN_M 到 $MAX_M, 每个m值重复 $REPEATS 次"

for m in $(seq $MIN_M $MAX_M)
do
    echo "测试迁移边数量 m=$m, $REPEATS 次重复..."
    for i in $(seq 1 $REPEATS)
    do
        {
         echo "运行 TreeMix m=${m}, 重复 ${i}..."
         nohup treemix -i $INPUT_FILE -root hpAfrica -o ChinaTreemix${m}.${i} -m ${m} -se -bootstrap -k 500 -global -noss &
        }
    done
    # 等待当前m值的所有重复完成
    wait
    echo "迁移边数量 m=$m 的所有重复已完成"
done

##======= 结果分析 =======
echo "TreeMix分析完成，请使用R脚本分析方差解释率以确定最佳的迁移边数量"
echo "推荐使用以下R命令评估结果:"
echo "R代码示例:"
echo "source(\"https://raw.githubusercontent.com/jyanglab/treemix/master/plotting_funcs.R\")"
echo "plot_resid(\"ChinaTreemix.5.1\", pop_order=\"pop.order.file\")"

# 创建结果汇总目录
mkdir -p ../output/treemix_results
echo "所有结果文件已保存到对应输出目录"

# 提供绘制最优模型选择图的R代码
cat << 'EOF' > ../output/treemix_results/plot_optimal_m.R
# R代码用于分析和选择最优的迁移边数量
# 使用方法: Rscript plot_optimal_m.R

# 加载所需的包
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

# 设置工作目录(根据需要修改)
# setwd("/path/to/treemix/results")

# 读取每个m值的方差解释率
read_variance <- function(prefix, m, rep=1) {
  filename <- paste0(prefix, m, ".", rep, ".llik")
  if(file.exists(filename)) {
    data <- read.table(filename)
    return(data[nrow(data), 2]) # 返回最后一行的方差解释率
  } else {
    return(NA)
  }
}

# 收集所有m值的方差解释率
collect_variance <- function(prefix, min_m=0, max_m=10, reps=5) {
  results <- data.frame(m=integer(), rep=integer(), variance=numeric())
  
  for(m in min_m:max_m) {
    for(rep in 1:reps) {
      variance <- read_variance(prefix, m, rep)
      if(!is.na(variance)) {
        results <- rbind(results, data.frame(m=m, rep=rep, variance=variance))
      }
    }
  }
  
  return(results)
}

# 收集结果
results <- collect_variance("ChinaTreemix", 0, 10, 5)

# 计算每个m值的平均方差解释率和标准差
summary_stats <- aggregate(variance ~ m, results, 
                          function(x) c(mean=mean(x), sd=sd(x), n=length(x)))
summary_df <- data.frame(m=summary_stats$m, 
                         mean=summary_stats$variance[,1],
                         sd=summary_stats$variance[,2],
                         n=summary_stats$variance[,3])

# 绘制方差解释率随m值的变化
p <- ggplot(summary_df, aes(x=m, y=mean)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  labs(title="TreeMix: 方差解释率随迁移边数量的变化",
       x="迁移边数量 (m)",
       y="方差解释率") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        text = element_text(size=14))

# 保存图表
ggsave("optimal_m_selection.pdf", p, width=10, height=7)
ggsave("optimal_m_selection.png", p, width=10, height=7, dpi=300)

# 打印表格形式的结果
print(summary_df)

# 识别"拐点"
# 计算方差解释率的增长率
summary_df$increase <- c(NA, diff(summary_df$mean))
summary_df$rel_increase <- summary_df$increase / summary_df$mean * 100

# 输出结果
cat("\n最优m值选择指南:\n")
cat("1. 查看方差解释率是否达到高值(通常>0.99或99%)\n")
cat("2. 观察增长率是否开始显著下降(通常<1%)\n\n")

cat("方差解释率及其增长率:\n")
print(summary_df[, c("m", "mean", "increase", "rel_increase")])

# 建议最优m值
threshold <- 0.01  # 1%的增长率阈值
candidates <- which(summary_df$rel_increase < threshold & !is.na(summary_df$rel_increase))
if(length(candidates) > 0) {
  suggested_m <- summary_df$m[min(candidates)]
  cat("\n建议的最优m值:", suggested_m, 
      "\n(此时方差解释率为:", round(summary_df$mean[summary_df$m == suggested_m]*100, 2), 
      "%, 相对增长率低于1%)\n")
} else {
  max_var_m <- summary_df$m[which.max(summary_df$mean)]
  cat("\n未发现明显的拐点, 建议使用最大方差解释率对应的m值:", max_var_m, 
      "\n(最大方差解释率为:", round(max(summary_df$mean, na.rm=TRUE)*100, 2), "%)\n")
}
EOF

echo "已创建R脚本用于分析最优m值选择: ../output/treemix_results/plot_optimal_m.R"
echo "运行完所有TreeMix分析后，可以使用以下命令分析结果:"
echo "cd ../output/treemix_results && Rscript plot_optimal_m.R"
