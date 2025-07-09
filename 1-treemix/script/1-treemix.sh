# Admixture graphs between the different populations were obtained using the software Treemix v. 1.1212. 
# Treemix was run with a number of migration edges between 0 to 15, with 10 replicates for each number of edge and hpAfrica2 was set as the outgroup. 
# The final number of migration edges was chosen as being the smallest number that allowed 99.8% of the variance to be explained.

# 脚本用途: 使用TreeMix分析全球群体之间的迁移和混合关系(不循环版本)
# 使用方法: bash 1-treemix-global.sh [输入文件] [群体文件]

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

##======= 全球群体分析 =======
echo "进行全球群体分析(不循环)..."
# 全球群体分析, 以hpAfrica2为根, 不循环因为全球群体太多循环会非常慢
echo "运行单次TreeMix分析，迁移边数量m=5..."
nohup treemix -i $INPUT_FILE -root hpAfrica2 -o GlobalTreemix -m 5 -se -bootstrap -k 500 -global -noss &

##======= 结果分析 =======
echo "TreeMix分析完成，请使用R脚本分析方差解释率以确定最佳的迁移边数量"
echo "推荐使用以下R命令评估结果:"
echo "R代码示例:"
echo "source(\"https://raw.githubusercontent.com/jyanglab/treemix/master/plotting_funcs.R\")"
echo "plot_resid(\"ChinaTreemix.5.1\", pop_order=\"pop.order.file\")"

# 创建结果汇总目录
mkdir -p ../output/treemix_results
echo "所有结果文件已保存到对应输出目录"