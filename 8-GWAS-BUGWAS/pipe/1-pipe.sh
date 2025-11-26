#!/bin/bash
set -euo pipefail

# 路径配置
# ├── input
# │   ├── list.txt
# │   ├── pheno.txt
# └── output
#* ======软件配置区======
PYTHON3_PATH="/home/luolintao/miniconda3/envs/BigLin/bin/python3"
CODE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON_FASTA_SCRIPT="${CODE_DIR}/python/0_从aln提取需要的样本.py"
PYTHON_BUGWAS_SCRIPT="${CODE_DIR}/python/0_准备bugwas输入文件.py"
PYTHON3_ANNO="${CODE_DIR}/python/3_基因注释.py"
R_SCRIPT="${CODE_DIR}/python/2_bugwas分析.r"
GFF_FILE="${CODE_DIR}/conf/NC_000915.gff"
GEMMA_PATH="${CODE_DIR}/bin/gemma"
VERYFASTTREE_PATH="VeryFastTree"

#* ======文件配置区======
INPUT_DIR="${CODE_DIR}/tmp/海拔/input"  
OUTPUT_DIR="${CODE_DIR}/tmp/海拔/output"
FASTA_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/1-序列处理流/output/merge_fasta/不考虑InDel/7544_Sample/"

# 第一步：准备序列与基因型输入
"${PYTHON3_PATH}" "${PYTHON_FASTA_SCRIPT}" \
    "${INPUT_DIR}/list.txt" \
    "${FASTA_DIR}/" \
    "${INPUT_DIR}/alignmentSub.aln"

snp-sites -c -o \
     "${INPUT_DIR}/alignmentSub.temp" \
     "${INPUT_DIR}/alignmentSub.aln"

nohup "${VERYFASTTREE_PATH}" -nt \
    -threads 24 \
    "${INPUT_DIR}/alignmentSub.temp" \
    > "${INPUT_DIR}/treeFull.nwk" \
    2> "${INPUT_DIR}/treeFull.log" &

snp-sites -c -v -o \
     "${INPUT_DIR}/tmp1.vcf" \
     "${INPUT_DIR}/alignmentSub.aln"

sed "s/#CHROM/CHROM/" \
    "${INPUT_DIR}/tmp1.vcf" \
    | grep -v "#" | cut -f 2,4-5,10- | sed "1s/POS/ps/" \
    > "${INPUT_DIR}/tmp2.txt"

"${PYTHON3_PATH}" "${PYTHON_BUGWAS_SCRIPT}" \
     "${INPUT_DIR}/tmp2.txt" \
     --output_file "${INPUT_DIR}/geno_biallelic_SNP.txt"

# 第二步：运行BUGWAS分析
Rscript "${R_SCRIPT}" \
  "${GFF_FILE}" \
  "${INPUT_DIR}" \
  "${OUTPUT_DIR}" \
  "${GEMMA_PATH}"

# 第三步：整理注释显著位点
#? --dist_threshold   越大，保留的显著位点越多
#? --logp_threshold 2 意味着只保留 p < 0.01 的 SNP 进入后续注释
"${PYTHON3_PATH}" "${PYTHON3_ANNO}" \
    --gwas_file "${OUTPUT_DIR}/bugwas_biallelic_lmmout_allSNPs.txt" \
    --gff_file "${GFF_FILE}" \
    --output_dir "${OUTPUT_DIR}/func" \
    --logp_threshold 2 \
    --dist_threshold 1000           
