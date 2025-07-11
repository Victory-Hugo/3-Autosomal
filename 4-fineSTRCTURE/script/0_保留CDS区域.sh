#!/bin/bash

# ====================================================================
# 脚本功能：根据 CDS_region.csv 中的基因区域过滤 VCF 文件，仅保留这些区域内的变异
# 使用工具：awk, bcftools, bgzip, tabix
# ====================================================================

# 设置错误检测
set -euo pipefail

# -----------------------------------
# 定义输入输出文件路径
# -----------------------------------
VCF_FILE="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.NoN.maf99.WGS.recode.SNP.noN.vcf.gz"
REGION_CSV="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/4-fineSTRCTURE/conf/CDS_region.csv"
OUTPUT_VCF="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/merged_biallelic_7544.NoN.maf99.WGS.recode.SNP.noN.CDS.vcf.gz"

# 定义染色体名称（根据实际情况修改）
# 请确保此名称与 VCF 文件中的染色体名称一致
CHROM="NC_000915.1"  # 示例为1，请根据 VCF 文件实际染色体名称调整

# -----------------------------------
# 步骤1：将 CSV 转换为 BED 格式
# -----------------------------------
echo "转换 CSV 文件为 BED 格式..."

awk -F, -v chrom="$CHROM" 'NR>1 {print chrom "\t" ($1 - 1) "\t" $2}' "$REGION_CSV" > regions.bed

echo "BED 文件已生成：regions.bed"

# -----------------------------------
# 步骤2：使用 bcftools 过滤 VCF 文件
# -----------------------------------
echo "开始使用 bcftools 过滤 VCF 文件..."

bcftools view -R regions.bed -Oz -o "$OUTPUT_VCF" "$VCF_FILE" && echo "VCF 文件过滤成功：$OUTPUT_VCF" || { echo "VCF 文件过滤失败"; exit 1; }

# -----------------------------------
# 步骤3：索引输出的 VCF 文件
# -----------------------------------
echo "为输出的 VCF 文件建立索引..."

bcftools index "$OUTPUT_VCF" && echo "索引建立成功：${OUTPUT_VCF}.csi" || { echo "索引建立失败"; exit 1; }

# -----------------------------------
# 步骤4：清理临时文件
# -----------------------------------
echo "删除临时 BED 文件..."

rm regions.bed && echo "临时文件已删除。"

echo "所有步骤已成功完成！最终输出文件为：$OUTPUT_VCF"
