#!/usr/bin/env bash
set -euo pipefail

# 输入文件前缀（不带 .bed/.bim/.fam 后缀）
bed_prefix="/mnt/d/幽门螺旋杆菌/Script/分析结果/ADMIXTURE/output/ADMIXTURE_800_filtered"

# 参数设置
K_MIN=2       # 最小 K
K_MAX=12       # 最大 K
THREADS=1     # 线程数
CV_FOLDS=10    # 交叉验证折数
BOOTSTRAPS=100 # bootstrap 次数
SEED=12345     # 随机种子（可选，保证可重复性）

# 逐 K 循环运行
for K in $(seq $K_MIN $K_MAX); do
  echo "=== Running ADMIXTURE for K=${K} ==="
  admixture \
    --cv=$CV_FOLDS \
    -j$THREADS   \
    -B$BOOTSTRAPS \
    -s$SEED      \
    "${bed_prefix}.bed" $K \
    | tee "${bed_prefix}_K${K}.log"
done

echo "All done. 结果保存在："
echo "  - ${bed_prefix}_K{10..12}.Q  （群体构成矩阵）"
echo "  - ${bed_prefix}_K{10..12}.P  （等位基因频率矩阵）"
echo "  - ${bed_prefix}_K{10..12}.log（运行日志，包括 CV error）"
