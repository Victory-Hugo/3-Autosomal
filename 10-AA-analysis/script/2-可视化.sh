#!/bin/bash
# 激活archetypal_analysis环境
# conde activate archetypal_analysis

INPUT_Q="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/10-AA-analysis/output/ADMIXTURE_800_top5_AA.5.Q"
OUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/10-AA-analysis/output"
mkdir -p ${OUT_DIR}

cd ${OUT_DIR}

archetypal-plot \
    -i ${INPUT_Q} \
    -p bar_simple

archetypal-plot \
    -i ${INPUT_Q} \
    -p plot_simplex
