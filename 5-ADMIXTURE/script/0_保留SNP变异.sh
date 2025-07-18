#!/usr/bin/env bash
# 路径定义
INPUT="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/WGS.aln.fasta.CDS_hpEAsia_2792.vcf.gz"
OUTPUT="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge/WGS.aln.fasta.CDS_hpEAsia_2792_SNP.vcf.gz"

# 方法一：使用 bcftools
# 仅保留 SNP（-v snps），并用 -Oz 指定输出为 bgzip 压缩格式
bcftools view -v snps "$INPUT" -Oz -o "$OUTPUT"
# 生成索引，方便后续快速查询
bcftools index "$OUTPUT"
