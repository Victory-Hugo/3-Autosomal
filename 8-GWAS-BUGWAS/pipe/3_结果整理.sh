#! /bin/bash
PYTHON3_PATH="/home/luolintao/miniconda3/bin/python3" #todo python3的路径
PYTHON3_ANNO='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/python/3_基因注释.py'
OUTPUT_DIR='/mnt/d/幽门螺旋杆菌/Script/分析结果/8-临床结局相关/output/GWAS/Atrophy-萎缩/'
GFF_FILE='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/8-GWAS-BUGWAS/conf/NC_000915.gff'


# 整理结果
 ${PYTHON3_PATH} ${PYTHON3_ANNO} \
    --gwas_file "${OUTPUT_DIR}/bugwas_biallelic_lmmout_allSNPs.txt" \
    --gff_file ${GFF_FILE} \
    --output_dir "${OUTPUT_DIR}/func" \
    --logp_threshold 2 \
    --dist_threshold 1
