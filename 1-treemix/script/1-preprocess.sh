#!/bin/bash
set -e  # 遇到错误立即退出

WORK_DIR='/mnt/d/幽门螺旋杆菌/Script/分析结果/4-treemix/data'
PAIRE_WISE_PY='/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/1-treemix/src/pairwise_fst.py'
GROUP_FILE="${WORK_DIR}/group.csv"

# bash "/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/1-treemix/src/makecluster.sh" \
#  "${WORK_DIR}/7544_filtered_pruned_data" \
#     0 \
#     0

echo "======== FST 计算前准备 ========"

# 确保 group.csv 文件可访问
echo "检查群体文件是否存在..."
if [ -f "${GROUP_FILE}" ]; then
    echo "找到群体文件: ${GROUP_FILE}"
    
    # 检查文件内容
    echo "文件内容预览:"
    head -n 5 "${GROUP_FILE}"
    
    # 检查文件格式，确保没有空格和特殊字符
    echo "检查文件格式..."
    if grep -q "[[:space:]]" "${GROUP_FILE}"; then
        echo "警告: 群体文件包含空格或制表符，这可能会导致问题"
        echo "建议修复文件格式"
    fi
    
    GROUPS_COUNT=$(wc -l < "${GROUP_FILE}")
    echo "群体文件包含 ${GROUPS_COUNT} 个群体"
else
    echo "错误: ${GROUP_FILE} 文件不存在"
    echo "请确保群体文件存在并命名为 group.csv"
    exit 1
fi

# 检查输入文件
echo "检查输入文件..."
if [ -f "${WORK_DIR}/7544_filtered_pruned_data.bed" ] && \
   [ -f "${WORK_DIR}/7544_filtered_pruned_data.bim" ] && \
   [ -f "${WORK_DIR}/7544_filtered_pruned_data.fam" ]; then
    echo "找到PLINK输入文件"
else
    echo "错误: PLINK输入文件不完整"
    exit 1
fi

# 检查分组文件
CSV_FILE="${WORK_DIR}/7544_filtered_pruned_data.csv"
if [ -f "${CSV_FILE}" ]; then
    echo "找到CSV分组文件: ${CSV_FILE}"
else
    echo "警告: 未找到CSV分组文件 ${CSV_FILE}"
    echo "将尝试运行makecluster.sh创建它"
    
    # 运行makecluster.sh创建CSV文件
    bash "/mnt/f/OneDrive/文档（科研）/脚本/Download/3-Autosomal/1-treemix/src/makecluster.sh" \
      "${WORK_DIR}/7544_filtered_pruned_data" \
      0 \
      0
fi

# 切换到工作目录
cd $WORK_DIR
echo "当前工作目录: $(pwd)"

echo "======== 开始计算FST ========"
echo "将执行: python $PAIRE_WISE_PY ${WORK_DIR}/7544_filtered_pruned_data --threads 8 --group-file ${GROUP_FILE}"

# 运行FST计算
python $PAIRE_WISE_PY \
    "${WORK_DIR}/7544_filtered_pruned_data" \
    --threads 8 \
    --group-file "${GROUP_FILE}"

echo "======== FST 计算完成 ========"
if [ -f "Fst_mat.csv" ]; then
    echo "成功生成结果文件: Fst_mat.csv"
    echo "结果预览:"
    head -n 10 Fst_mat.csv
else
    echo "警告: 未找到结果文件 Fst_mat.csv"
    echo "FST计算可能失败"
fi