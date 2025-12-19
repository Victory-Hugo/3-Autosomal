#!/bin/bash
# 使用自定义archetypal-plot脚本的示例

# 设置路径
SCRIPT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/10-AA-analysis/script"
INPUT_Q="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/10-AA-analysis/output/ADMIXTURE_800_top5_AA.5.Q"
FAM_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/10-AA-analysis/input/ADMIXTURE_800_filtered_top5.fam"
OUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/10-AA-analysis/output"

# 进入输出目录
cd ${OUT_DIR}

echo "使用自定义archetypal-plot脚本生成可视化图形..."

# 生成所有类型的图形
python3 ${SCRIPT_DIR}/archetypal_plot_custom.py \
    -i ${INPUT_Q} \
    -f ${FAM_FILE} \
    -p all \
    -title "Population_Colored" \
    -c tab20 \
    -o aa_population_colored \
    -dpi 300

echo "生成完成！输出文件位于: ${OUT_DIR}"

# 也可以单独生成特定类型的图形
echo "生成单纯形图..."
python3 ${SCRIPT_DIR}/archetypal_plot_custom.py \
    -i ${INPUT_Q} \
    -f ${FAM_FILE} \
    -p plot_simplex \
    -title "Simplex_PopColored" \
    -c Set3 \
    -o aa_simplex_only \
    -dpi 500

echo "生成条形图..."
python3 ${SCRIPT_DIR}/archetypal_plot_custom.py \
    -i ${INPUT_Q} \
    -f ${FAM_FILE} \
    -p bar_simple \
    -title "Bar_PopColored" \
    -c viridis \
    -o aa_bar_only \
    -dpi 300

echo "所有可视化已完成！"