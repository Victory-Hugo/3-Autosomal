#!/usr/bin/env bash
set -euo pipefail

# 输入 VCF 文件
vcf="/home/luolintao/Helicopter/Script/分析结果/merged_fasta/WGS.aln.fasta.CDS.vcf.gz"
# 输出前缀（可根据需要修改）
out_prefix="/home/luolintao/Helicopter/Script/分析结果/ADMIXTURE/data/7544_CDS/ADMIXTURE_7544_CDS"

# （1）用 PLINK 生成二进制文件
plink --vcf "$vcf" \
      --double-id \
      --make-bed \
      --allow-extra-chr \
      --out "$out_prefix"

# 备份原始 .bim
cp "${out_prefix}.bim" "${out_prefix}.bim.bak"

# （2）修正染色体名：把 NC_000915.1 全部替换成 1
awk 'BEGIN{OFS="\t"} 
     $1=="NC_000915.1"{$1="1"} 
     {print}' \
    "${out_prefix}.bim" > "${out_prefix}.tmp.bim" \
&& mv "${out_prefix}.tmp.bim" "${out_prefix}.bim"

echo "✅ ${out_prefix}.bim 已替换，原文件保存在 ${out_prefix}.bim.bak"

# （3）可选：对缺失率较高的位点/样本做过滤
# 如果你用的是 PLINK1.9，可以去掉 --alleleACGT；PLINK2.0 可保留
plink --bfile "$out_prefix" \
      --geno 0.10 \
      --mind 0.10 \
      --indep-pairwise 50 10 0.1 \
      --max-maf 0.01 \
      --make-bed \
      --alleleACGT \
      --out "${out_prefix}_filtered"

echo "✅ 过滤完成，结果在 ${out_prefix}_filtered.{bed,bim,fam}"

