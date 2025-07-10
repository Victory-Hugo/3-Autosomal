#!/bin/bash

# 使用TreeMix分析全球群体间的迁移和混合关系 (0-15条迁移边)
# 用法: bash 1-treemix-global.sh

set -euo pipefail

INPUT_FILE="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data/TreeMix.treeout.gz"
OUTPUT_DIR="/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/output/global"
OUTGROUP="hpAfrica2"


echo "开始运行全球群体TreeMix分析..."

mkdir -p "$OUTPUT_DIR"
echo "迁移边数量,方差解释量(%)" > "${OUTPUT_DIR}/variance_explained.txt"

for m in {0..15}; do
  echo "分析迁移边数量 m=${m}..."
  treemix -i "$INPUT_FILE" \
    -root "$OUTGROUP" \
    -o "${OUTPUT_DIR}/GlobalTreemix_m${m}" \
    -m "$m" \
    -se -bootstrap \
    -global -noss
done

echo "TreeMix分析完成。结果摘要:"
