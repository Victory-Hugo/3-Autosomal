#!/usr/bin/env bash
set -euo pipefail

# ==== 配置区 ====

# 输入文件前缀（不带 .bed/.bim/.fam 后缀）
input_prefix="/home/luolintao/Helicopter/Script/分析结果/ADMIXTURE/data/TXQ/WGS.aln.fasta.CDS_地理东亚及周边_3378_filtered"

# 指定 ADMIXTURE 输出结果目录
output_dir="/home/luolintao/Helicopter/Script/分析结果/ADMIXTURE/output"
mkdir -p "$output_dir"

# 存放各个子脚本的目录
script_dir="/home/luolintao/Helicopter/Script/分析结果/ADMIXTURE/script"
mkdir -p "$script_dir"

# ADMIXTURE 参数设置
K_MIN=6         # 最小 K
K_MAX=15        # 最大 K
THREADS=8       # 线程数
CV_FOLDS=10     # 交叉验证折数
BOOTSTRAPS=100  # bootstrap 次数
SEED=12345      # 随机种子

# ==== 脚本生成 ====

for K in $(seq $K_MIN $K_MAX); do
  script_path="${script_dir}/run_ADMIX_K${K}.sh"
  cat > "$script_path" <<EOF
#!/usr/bin/env bash
set -euo pipefail

# 自动生成脚本：K=${K}

# 输入/输出前缀
input_prefix="${input_prefix}"
output_dir="${output_dir}"
mkdir -p "\$output_dir"

# ADMIXTURE 参数
THREADS=${THREADS}
CV_FOLDS=${CV_FOLDS}
BOOTSTRAPS=${BOOTSTRAPS}
SEED=${SEED}

echo "=== Running ADMIXTURE for K=${K} ==="
# 调用 admixture 执行
/home/bin/admixture \\
  --cv=\$CV_FOLDS \\
  -j\$THREADS \\
  -B\$BOOTSTRAPS \\
  -s\$SEED \\
  "\${input_prefix}.bed" ${K} \\
  | tee "\${output_dir}/7544_CDS_K${K}.log"

echo "Finished K=${K}. 结果："
echo "  - \${output_dir}/7544_CDS_K${K}.Q"
echo "  - \${output_dir}/7544_CDS_K${K}.P"
echo "  - \${output_dir}/7544_CDS_K${K}.log"
EOF

  # 赋予执行权限
  chmod +x "$script_path"
  echo "Generated script: $script_path"
done

echo "所有子脚本已生成在：${script_dir}，可以逐个执行或批量提交。"
