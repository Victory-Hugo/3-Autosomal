#!/bin/bash

# 定义输入目录和待处理的主VCF文件
INPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/fineSTRUCTURE/INPUT/New_Chunk"
VCF_FILE="/mnt/d/幽门螺旋杆菌/Script/分析结果/fineSTRUCTURE/INPUT/二等位_核心SNP_99genomes_修正_bgzip_noIndel_CDS_noN.vcf.gz"

# 定义处理单个样本文件的函数
process_file() {
    SAMPLE_FILE="$1"
    BASENAME=$(basename "$SAMPLE_FILE" .txt)
    OUTPUT_VCF="${INPUT_DIR}/${BASENAME}.vcf.gz"

    echo "正在处理样本文件: $SAMPLE_FILE"

    # 使用 bcftools 提取 VCF 子集
    bcftools view \
        --force-samples \
        --threads 16 \
        --samples-file "$SAMPLE_FILE" \
        -Oz \
        "$VCF_FILE" \
        -o "$OUTPUT_VCF"

    # 为生成的 VCF 文件建立索引
    bcftools index "$OUTPUT_VCF"
}

# 导出函数和变量，以便 GNU parallel 在子进程中可以访问
export -f process_file
export INPUT_DIR
export VCF_FILE

# 使用 GNU parallel 并行处理所有 .txt 文件，-j 9 表示同时最多运行 9 个任务
parallel -j 11 process_file {} ::: "$INPUT_DIR"/*.txt

echo "全部文件处理完毕！"
