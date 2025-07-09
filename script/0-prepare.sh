#!/bin/bash
set -e

# 设置工作目录
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/"
SCRIPT_PY="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/src/0-plink→treemix.py"

# 创建输出目录
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "Step 1: 首先从VCF文件创建PLINK格式文件...这一步已经在vcf→geno,ind,snp完成"

echo "Step 2: 准备群体信息文件（关键步骤）..."
# 如果已有pop.cov文件，可跳过此步骤
# 格式应为：
#   sample1  sample1  POP_A
#   sample2  sample2  POP_B

# echo "Step 3: 生成按群体分层的频率文件..."
# plink \
#   --bfile /mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data \
#   --freq \
#   --within "/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/conf/pop.cov" \
#   --allow-extra-chr \
#   --out /mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/group_freq

# echo "Step 4: 压缩频率文件..."
# gzip -f /mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/group_freq.frq.strat

echo "Step 5: 执行转换为TreeMix格式..."
# 检查输入文件是否存在
INPUT_FILE="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/group_freq.frq.strat.gz"
OUTPUT_FILE="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/treemix_input"

if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 $INPUT_FILE 不存在!"
    echo "请先运行步骤3和4生成此文件，或确认文件路径正确。"
    exit 1
fi

echo "正在运行转换脚本..."
# 使用已有的0-plink→treemix.py脚本
/home/luolintao/miniconda3/envs/pyg/bin/python3 \
   "$SCRIPT_PY" \
    "$INPUT_FILE" \
    > "$OUTPUT_FILE" 2> "${OUTPUT_FILE}.log"

# 检查结果
if [ $? -ne 0 ] || [ ! -s "$OUTPUT_FILE" ]; then
    echo "错误: 转换失败或输出文件为空!"
    echo "查看日志文件了解详情: ${OUTPUT_FILE}.log"
    exit 1
else
    echo "转换成功！输出文件: $OUTPUT_FILE"
    echo "日志文件: ${OUTPUT_FILE}.log"
    echo "输出文件前5行预览:"
    head -5 "$OUTPUT_FILE"
fi
# 压缩输出文件
echo "压缩输出文件..."
gzip -f "$OUTPUT_FILE"
echo "转换完成：${OUTPUT_FILE}.gz"

echo "正确的TreeMix输入格式应如下："
echo "POP_A POP_B POP_C  # 首行为群体名称"
echo "12,5 20,3 ?,?     # 每行为SNP数据，格式为"REF等位计数,ALT等位计数""
echo ""
echo "执行zcat命令查看前几行输出："
zcat "${OUTPUT_FILE}.gz" | head -5

echo ""
echo "如果你需要运行TreeMix分析，请使用以下命令："
echo "bash 1-treemix-china.sh ${OUTPUT_FILE}.gz <群体文件>"