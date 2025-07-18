#!/usr/bin/env bash

# 输入、输出文件路径
INPUT="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/WGS.aln.fasta.CDS_地理东亚及周边_3378.vcf.gz"
OUTPUT="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/WGS.aln.fasta.CDS_地理东亚及周边_3378_filter.vcf.gz"

# 使用 bcftools 去除 ALT 为 "." 的记录
# -e 'ALT="."' 表示剔除 ALT 字段精确为 "." 的行
bcftools view -e 'ALT="."' "$INPUT" -Oz -o "$OUTPUT"

# 为输出文件建立索引（可选，但推荐用于后续快速查询）
bcftools index "$OUTPUT"
