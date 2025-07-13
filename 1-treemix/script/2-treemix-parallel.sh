#!/bin/bash
# 生成16个独立的TreeMix运行脚本，对应迁移边数量m从0到15，每种m做10次重复
#*===============================
# todo 设定输入文件路径，输入文件通过./python/plink2treemix.py生成
INPUT_FILE="/home/luolintao/Helicopter/Script/分析结果/treemix/data/TreeMix_China3622.gz"
# todo 设定输出目录
OUTPUT_DIR="/home/luolintao/Helicopter/Script/分析结果/treemix/output/China/"
# todo 设定脚本输出目录
SCRIPT_DIR="/home/luolintao/Helicopter/Script/分析结果/treemix/script/China/"
# todo 设定outgroup
OUTGROUP="hpAfrica2"
#*===============================
echo "生成TreeMix分析脚本..."

# 确保脚本输出目录存在
mkdir -p "$OUTPUT_DIR"
mkdir -p "$SCRIPT_DIR"



# 为每个迁移边数量生成脚本
for m in {0..15}; do
  script_name="${SCRIPT_DIR}/treemix_m${m}.sh"
  
  cat > "$script_name" << EOF
#!/bin/bash
#
# TreeMix分析脚本 - 迁移边数量 m=${m}, 重复 rep=1..10
# 用法: bash treemix_m${m}.sh

set -euo pipefail

INPUT_FILE="$INPUT_FILE"
OUTPUT_DIR="$OUTPUT_DIR"
OUTGROUP="$OUTGROUP"

echo "开始分析迁移边数量 m=${m}..."

mkdir -p "\$OUTPUT_DIR"

# 对每个 m 做 10 次重复
for rep in {1..10}; do
  prefix="\${OUTPUT_DIR}/Treemix${m}.\${rep}"
  echo "  [m=${m} rep=\${rep}] 运行 TreeMix..."
  treemix -i "\$INPUT_FILE" \\
    -root "\$OUTGROUP" \\
    -o "\$prefix" \\
    -m ${m} \\
    -se -bootstrap \\
    -global \\
    -k 500 \\
    -threads 8

  echo "  [m=${m} rep=\${rep}] 分析完成"
done

echo "全部重复实验完成：m=${m}"
EOF

  chmod +x "$script_name"
  echo "已生成脚本: $script_name"
done

echo "所有脚本已生成完毕，存放在 $SCRIPT_DIR"
