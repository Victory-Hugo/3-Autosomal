#!/usr/bin/env bash
set -euo pipefail

# 输入 VCF 文件
vcf="/home/luolintao/Helicopter/Script/分析结果/ADMIXTURE/data/TXQ/merged_dedup_biallelic.CORE.SNP.recode_CDS.vcf.gz"
# 输出前缀（可根据需要修改）
out_prefix="/home/luolintao/Helicopter/Script/分析结果/ADMIXTURE/data/TXQ/WGS.aln.fasta.CDS_地理东亚及周边_3378"

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
#! --maf 0.01仅保留MAF大于0.01的频率的变异，这些变异在群体中更常见，减少了低频变异的噪音。
#! --geno 0.10 表示过滤缺失率超过10%的位点，这些位点由于缺失数据过多，可能会影响后续分析的准确性。
#! --mind 0.10 表示过滤缺失率超过10%的样本，这样可以确保分析中只包含数据完整的样本。
#! --indep-pairwise 500 10 0.2 表示进行连锁不平衡(LD)过滤，保留每500个SNP中相关性小于0.2的SNP对。
plink --bfile "$out_prefix" \
      --indep-pairwise 100 10 0.2 \
      --maf 0.01   \
      --make-bed \
      --alleleACGT \
      --out "${out_prefix}_filtered"

echo "✅ 过滤完成，结果在 ${out_prefix}_filtered.{bed,bim,fam}"

