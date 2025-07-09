#!/bin/bash
set -e

# 设置工作目录
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/"
SCRIPT_PY="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/src/0-plink→treemix.py"
POP_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/conf/pop.cov"
# 创建输出目录
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "Step 1: 从现有的PLINK文件创建带有正确SNP ID的版本..."
PLINK_PREFIX="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/geno-ind-snp/7544_filtered_pruned_data"
TEMP_PREFIX="${OUTPUT_DIR}/7544_with_snp_ids"

# # 使用染色体、位置和等位基因信息创建唯一的SNP ID
# echo "为SNP设置唯一ID（格式：chr:pos:A1:A2）..."
# plink \
#   --bfile "$PLINK_PREFIX" \
#   --set-missing-var-ids '@:#:\$1:\$2' \
#   --allow-extra-chr \
#   --make-bed \
#   --out "$TEMP_PREFIX"

# # 检查结果
# if [ $? -ne 0 ]; then
#     echo "错误：PLINK命令执行失败"
#     cat "${TEMP_PREFIX}.log"
    
#     echo "尝试备选方案：使用位置和随机数..."
#     echo "创建新的SNP ID文件..."
    
#     # 使用 awk 创建新的 SNP ID
#     awk -v seed=$RANDOM 'BEGIN{srand(seed)} 
#          NR==FNR{pos[$4]++; cnt[$4]=pos[$4]; next} 
#          {if(cnt[$4]>1) 
#              print $1"\t"$1":"$4"_"cnt[$4]"\t"$3"\t"$4"\t"$5"\t"$6; 
#           else 
#              print $1"\t"$1":"$4"\t"$3"\t"$4"\t"$5"\t"$6}' \
#         "$PLINK_PREFIX.bim" "$PLINK_PREFIX.bim" > "${TEMP_PREFIX}.bim.new"
    
#     if [ -f "${TEMP_PREFIX}.bim.new" ]; then
#         mv "${TEMP_PREFIX}.bim.new" "${TEMP_PREFIX}.bim"
#         # 复制其他必要的文件
#         cp "${PLINK_PREFIX}.bed" "${TEMP_PREFIX}.bed"
#         cp "${PLINK_PREFIX}.fam" "${TEMP_PREFIX}.fam"
#         echo "成功创建新的SNP ID"
#     else
#         echo "错误：无法创建新的SNP ID"
#         exit 1
#     fi
# fi

echo "Step 2: 准备群体信息文件（关键步骤）..."
# 检查群体文件是否存在

if [ ! -f "$POP_FILE" ]; then
    echo "错误: 群体文件 $POP_FILE 不存在!"
    echo "请创建群体文件，格式如下:"
    echo "sample1  sample1  POP_A"
    echo "sample2  sample2  POP_B"
    exit 1
fi

# echo "Step 3: 生成按群体分层的频率文件..."
# plink \
#   --bfile "$TEMP_PREFIX" \
#   --freq \
#   --within "$POP_FILE" \
#   --allow-extra-chr \
#   --out "${OUTPUT_DIR}/group_freq"

# if [ $? -ne 0 ]; then
#     echo "错误：生成频率文件失败"
#     cat "${OUTPUT_DIR}/group_freq.log"
#     exit 1
# fi

# echo "Step 4: 压缩频率文件..."
# gzip -f "${OUTPUT_DIR}/group_freq.frq.strat"

# # 检查生成的频率文件
# echo "检查生成的频率文件前几行..."
# zcat "${OUTPUT_DIR}/group_freq.frq.strat.gz" | head -10

# # 检查是否所有SNP ID都是唯一的
# echo "检查SNP ID的唯一性..."
# TOTAL_SNPS=$(zcat "${OUTPUT_DIR}/group_freq.frq.strat.gz" | wc -l)
# UNIQUE_SNPS=$(zcat "${OUTPUT_DIR}/group_freq.frq.strat.gz" | cut -f2 | sort | uniq | wc -l)
# echo "总行数: $TOTAL_SNPS"
# echo "唯一SNP数: $UNIQUE_SNPS"

echo "Step 5: 执行转换为TreeMix格式..."
# 检查输入文件是否存在
INPUT_FILE="${OUTPUT_DIR}/group_freq.frq.strat.gz"
OUTPUT_FILE="${OUTPUT_DIR}/treemix_input"

if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 $INPUT_FILE 不存在!"
    echo "请检查前面的步骤是否成功完成。"
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
    cat "${OUTPUT_FILE}.log"
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