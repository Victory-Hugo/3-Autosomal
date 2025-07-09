#!/bin/bash

# 处理VCF文件中重复ID问题的预处理脚本
# 用法: ./preprocess_vcf.sh <input_vcf> <output_vcf>

set -e

if [ $# -ne 2 ]; then
    echo "用法: $0 <input_vcf> <output_vcf>"
    echo "说明: 预处理VCF文件，解决重复ID问题"
    exit 1
fi

INPUT_VCF="$1"
OUTPUT_VCF="$2"

if [ ! -f "$INPUT_VCF" ]; then
    echo "错误: 输入VCF文件不存在: $INPUT_VCF"
    exit 1
fi

echo "开始预处理VCF文件..."
echo "输入文件: $INPUT_VCF"
echo "输出文件: $OUTPUT_VCF"

# 创建临时文件
TEMP_VCF="${OUTPUT_VCF}.tmp"

# 方法1: 使用bcftools重新生成唯一ID
if command -v bcftools &> /dev/null; then
    echo "使用bcftools重新生成SNP ID..."
    bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' "$INPUT_VCF" -O z -o "$TEMP_VCF"
    
    if [ $? -eq 0 ]; then
        mv "$TEMP_VCF" "$OUTPUT_VCF"
        echo "使用bcftools成功处理完成"
        exit 0
    else
        echo "bcftools处理失败，尝试其他方法..."
        rm -f "$TEMP_VCF"
    fi
fi

# 方法2: 使用zcat和awk手动处理
echo "使用awk手动处理重复ID..."
{
    # 保留头部信息
    zcat "$INPUT_VCF" | grep '^#'
    
    # 处理数据行，为每行生成唯一ID
    zcat "$INPUT_VCF" | grep -v '^#' | awk -F'\t' '
    {
        # 生成唯一ID: 染色体_位置_参考等位基因_替代等位基因_行号
        unique_id = $1 "_" $2 "_" $4 "_" $5 "_" NR
        $3 = unique_id
        print $0
    }' OFS='\t'
} | gzip > "$TEMP_VCF"

if [ $? -eq 0 ]; then
    mv "$TEMP_VCF" "$OUTPUT_VCF"
    echo "手动处理完成"
else
    echo "错误: 手动处理失败"
    rm -f "$TEMP_VCF"
    exit 1
fi

echo "VCF预处理完成！"
echo "输出文件: $OUTPUT_VCF"

# 检查输出文件
echo "检查输出文件完整性..."
ORIGINAL_LINES=$(zcat "$INPUT_VCF" | grep -v '^#' | wc -l)
PROCESSED_LINES=$(zcat "$OUTPUT_VCF" | grep -v '^#' | wc -l)

echo "原始数据行数: $ORIGINAL_LINES"
echo "处理后数据行数: $PROCESSED_LINES"

if [ "$ORIGINAL_LINES" -eq "$PROCESSED_LINES" ]; then
    echo "✓ 行数检查通过"
else
    echo "✗ 警告: 行数不一致，请检查处理结果"
fi
