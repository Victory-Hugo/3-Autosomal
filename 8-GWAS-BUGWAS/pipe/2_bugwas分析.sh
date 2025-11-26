#!/bin/bash
set -euo pipefail

# 配置路径参数（根据需要自行修改）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GFF_FILE="${SCRIPT_DIR}/../conf/NC_000915.gff"
INPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/input"
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/tmp/海拔/output"
GEMMA_PATH="/home/luolintao/miniconda3/envs/BigLin/bin/gemma"
R_SCRIPT="${SCRIPT_DIR}/2_bugwas分析.r"


Rscript "${R_SCRIPT}" \
  "${GFF_FILE}" \
  "${INPUT_DIR}" \
  "${OUTPUT_DIR}" \
  "${GEMMA_PATH}" 
