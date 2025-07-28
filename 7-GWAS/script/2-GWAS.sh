#!/usr/bin/env bash
#===============================================================================
# 脚本：优化后的 GWAS 流水线片段（仅变量提取，功能与注释不变）
#===============================================================================

set -euo pipefail
IFS=$'\n\t'

#-------------------------------------------------------------------------------
# 定义路径和文件前缀变量
#-------------------------------------------------------------------------------
DATA_DIR="/home/luolintao/09_GWAS_Autosomal/data"
CONF_DIR="/home/luolintao/09_GWAS_Autosomal/conf/"
OUTPUT_DIR="/home/luolintao/09_GWAS_Autosomal/output"
PREFIX="Han_SNP4"


#*========================== 进行 logistic 回归分析 ==========================
#*========================== 推荐使用第三种方式计算 ==========================

# 如下 3 种计算方式。区别详情见：`7-GWAS/markdown/3-逻辑斯蒂回归.md`

# #*======== 第一种方式 ========
# plink --bfile "${DATA_DIR}/${PREFIX}" \
#       --logistic \
#       --covar "${CONF_DIR}/${PREFIX}.cov" \
#       --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
#       --hide-covar \
#       --sex \
#       --out "${OUTPUT_DIR}/${PREFIX}"
# # 输出结果类似：
# #  CHR                  SNP         BP   A1       TEST    NMISS         OR         STAT            P 
# #    1           rs12127425     794332    A        ADD     3135     0.9973     -0.03465       0.9724
# #    1            rs9697294     832359    T        ADD     3090       0.67       -2.081      0.03743

# #*======== 第二种方式 ========
# #! 如果在 `--logistic` 命令中添加 `--beta` 参数，则会返回回归系数（BETA）而不是比值比。
# plink --bfile "${DATA_DIR}/${PREFIX}" \
#     --logistic \
#     --beta \
#     --covar "${CONF_DIR}/${PREFIX}.cov" \
#     --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
#     --hide-covar \
#     --sex \
#     --out "${OUTPUT_DIR}/${PREFIX}_BETA"
# # 输出结果类似：
# #  CHR                  SNP         BP   A1       TEST    NMISS       BETA         STAT            P 
# #    1           rs12127425     794332    A        ADD     3135  -0.002684     -0.03465       0.9724
# #    1            rs9697294     832359    T        ADD     3090    -0.4005       -2.081      0.03743

#*======== 第三种方式 ========
# 这条命令是用最新的 PLINK 2 来做全基因组的广义线性模型（GLM）关联检验
# 整体流程和 PLINK 1.9 的 --logistic 类似，但更灵活、性能更好。
plink2 \
  --bfile "${DATA_DIR}/${PREFIX}" \
  --covar "${CONF_DIR}/${PREFIX}.cov" \
  --covar-variance-standardize \
  --glm hide-covar sex \
  --threads 32 \
  --out "${OUTPUT_DIR}/${PREFIX}_plink2"
# 输出结果类似：
# #CHROM  POS     ID      REF     ALT     PROVISIONAL_REF?        A1      OMITTED A1_FREQ FIRTH?  TEST    OBS_CT  OR      LOG(OR)_SE      Z_STAT  P       ERRCODE
# 1       794332  rs12127425      G       A       Y       A       G       0.137959        N       ADD     3135    0.997319        0.0774623       -0.0346522      0.972357        .
# 1       832359  rs9697294       C       T       Y       T       C       0.0224919       N       ADD     3090    0.669987        0.192449        -2.08106        0.0374286       .
# 1       834056  rs28482280      A       C       Y       C       A       0.0188864       N       ADD     3071    0.688674        0.210703        -1.7702 0.0766931       .
